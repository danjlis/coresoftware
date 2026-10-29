#include "JetSmear.h"

#include "Jet.h"
#include "JetContainerv1.h"

#include <cdbobjects/CDBTF.h>  // for CDBTF1 (legacy method)

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TSystem.h>

#include <algorithm>  // for std::max, std::min
#include <cmath>
#include <iostream>  // for operator<<, basic_ostream
#include <string>

JetSmear::JetSmear(const std::string &name)
  : SubsysReco(name)
{
}

JetSmear::~JetSmear()
{
  delete m_JerSmearBand;
}

int JetSmear::InitRun(PHCompositeNode *topNode)
{
  // Node names and jet radius have no defaults and must be set by the user.
  if (m_inputRecoNode.empty() || m_inputCalibNode.empty() || m_inputTruthNode.empty() || m_outputRecoNode.empty() || jet_radius <= 0)
  {
    std::cout << PHWHERE << " JetSmear::InitRun() : set_InputRecoNode, set_InputCalibNode, set_InputTruthNode, set_OutputRecoNode and set_JetRadius must all be set! Exiting." << std::endl;
    gSystem->Exit(1);
  }

  // Create smeared jet node.
  int ret = CreateNodeTree(topNode);

  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  // Load JER smear band histogram.
  ret = initSmearing();

  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int JetSmear::initSmearing()
{
  std::string calibFile = m_calibFileOverride.empty() ? fetchCalibDir("JERSmear") : m_calibFileOverride;
  if (calibFile.empty())
  {
    std::cout << PHWHERE << " JetSmear::initSmearing() : No JERSmear calibration file found! Exiting." << std::endl;
    gSystem->Exit(1);
  }

  TFile *smearFile = TFile::Open(calibFile.c_str(), "READ");
  if (!smearFile || smearFile->IsZombie())
  {
    std::cout << PHWHERE << " JetSmear::initSmearing() : Could not open calibration file " << calibFile << "! Exiting." << std::endl;
    gSystem->Exit(1);
  }

  int r_index = (int) (jet_radius * 10 + 0.1);

  const std::string r_suffix = "_r0" + std::to_string(r_index);
  const std::string jes_suffix = m_useEMfracCalib ? "_EMfracJES":"_legacyJES";
  const std::string pileup_suffix = m_usePileup ? "_pileup":"_nop";
  std::string systematic_suffix = "_nominal";

  if (m_sys < 0)
    {
      systematic_suffix = "_sysdown";
    }
  else if (m_sys > 0)
    {
      systematic_suffix = "_sysup";
    }
  const std::string histo_name = "h_jer_smear" + r_suffix + pileup_suffix + jes_suffix + systematic_suffix;

  // Histogram is cloned and detached from the file so it can be closed below.
  delete m_JerSmearBand;
  m_JerSmearBand = nullptr;

  // JER smear band (required).
  TH1 *jerSmearBand = dynamic_cast<TH1 *>(smearFile->Get(histo_name.c_str()));
  if (!jerSmearBand)
  {
    std::cout << PHWHERE << " JetSmear::initSmearing() : Histogram " << histo_name << " not found in " << calibFile << "! Exiting." << std::endl;
    smearFile->Close();
    delete smearFile;
    gSystem->Exit(1);
  }
  m_JerSmearBand = dynamic_cast<TH1 *>(jerSmearBand->Clone());
  m_JerSmearBand->SetDirectory(nullptr);

  smearFile->Close();
  delete smearFile;

  return Fun4AllReturnCodes::EVENT_OK;
}


int JetSmear::process_event(PHCompositeNode *topNode)
{
  m_calib_jets.clear();
  m_raw_jets.clear();
  m_truth_jets.clear();
  m_matched_jets.clear();
  
  // Get raw, calib, truth and output jet nodes.
  // Run and event number for the smearing seed.
  EventHeader *evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!evtHeader)
  {
    std::cout << PHWHERE << " WARNING::Missing node EventHeader! Exiting." << std::endl;
    gSystem->Exit(1);
  }
  m_run = evtHeader->get_RunNumber();
  m_event = evtHeader->get_EvtSequence();

  m_reco_jetcontainer = findNode::getClass<JetContainer>(topNode, m_inputRecoNode);
  if (!m_reco_jetcontainer)
  {
    std::cout << PHWHERE << " WARNING::Missing node " << m_inputRecoNode << "! Exiting." << std::endl;
    gSystem->Exit(1);
  }

  m_calib_jetcontainer = findNode::getClass<JetContainer>(topNode, m_inputCalibNode);
  if (!m_calib_jetcontainer)
  {
    std::cout << PHWHERE << " WARNING::Missing node " << m_inputCalibNode << "! Exiting." << std::endl;
    gSystem->Exit(1);
  }

  m_truth_jetcontainer = findNode::getClass<JetContainer>(topNode, m_inputTruthNode);
  if (!m_truth_jetcontainer)
  {
    std::cout << PHWHERE << " WARNING::Missing node " << m_inputTruthNode << "! Exiting." << std::endl;
    gSystem->Exit(1);
  }

  m_output_jetcontainer = findNode::getClass<JetContainer>(topNode, m_outputRecoNode);
  if (!m_output_jetcontainer)
  {
    std::cout << PHWHERE << " WARNING::Missing node " << m_outputRecoNode << "! Exiting." << std::endl;
    gSystem->Exit(1);
  }

  m_hasVertex = m_reco_jetcontainer->has_zvertex() && std::isfinite(m_reco_jetcontainer->get_vertex_z());

  // Collect calib jets with E > 0 (as in ppg08), each paired with its raw jet.
  // JetCalib keeps calib jets 1:1 in order with the raw jets; calib E is not set
  // by JetCalib, so the raw E is used for the cut when the calib E is NaN.
  int ijet = 0;
  m_calib_jets.reserve(m_calib_jetcontainer->size());
  m_raw_jets.reserve(m_calib_jetcontainer->size());
  for (auto jet : *m_calib_jetcontainer)
  {
    Jet *rawjet = m_reco_jetcontainer->get_jet(ijet);
    ijet++;

    float e = jet->get_e();
    if (!std::isfinite(e) && rawjet)
    {
      e = rawjet->get_e();
    }
    if (!(e > 0))
    {
      continue;
    }

    m_calib_jets.push_back(jet);
    m_raw_jets.push_back(rawjet);
  }

  m_truth_jets.reserve(m_truth_jetcontainer->size());
  for (auto jet : *m_truth_jetcontainer)
  {
    m_truth_jets.push_back(jet);
  }

  match_reco_truth();

  smear_reco_jets();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetSmear::CreateNodeTree(PHCompositeNode *topNode)
{
  // Check nodes.
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "JetSmear::CreateNodeTree : DST Node missing, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // Create the ANTIKT and TOWER nodes if required.
  PHCompositeNode *antiktNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "ANTIKT"));
  if (!antiktNode)
  {
    antiktNode = new PHCompositeNode("ANTIKT");
    dstNode->addNode(antiktNode);
  }
  PHCompositeNode *towerNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "TOWER"));
  if (!towerNode)
  {
    towerNode = new PHCompositeNode("TOWER");
    antiktNode->addNode(towerNode);
  }

  JetContainer *test_jets = findNode::getClass<JetContainer>(topNode, m_outputRecoNode);
  if (!test_jets)
  {
    JetContainer *smear_jets = new JetContainerv1();
    PHIODataNode<PHObject> *smearjetNode;
    if (Verbosity() > 0)
    {
      std::cout << "JetSmear::CreateNode : creating " << m_outputRecoNode << std::endl;
    }
    smearjetNode = new PHIODataNode<PHObject>(smear_jets, m_outputRecoNode, "PHObject");
    towerNode->addNode(smearjetNode);
  }
  else
  {
    std::cout << "JetSmear::CreateNode : " << m_outputRecoNode << " already exists! " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string JetSmear::fetchCalibDir(const char *calibType)
{
  std::string calibName = calibType;
  return CDBInterface::instance()->getUrl(calibName);
}

// Grid matching (as in ppg08 dijetfinder::match_dijets_grid): build every
// truth-reco pair within dr_cut, sort by dR, and greedily accept pairs whose
// truth and reco jets are both still free. Each jet gets its nearest available
// partner. m_matched_jets holds (truth, reco) for every reco jet, in reco order;
// truth is nullptr for unmatched reco jets.
void JetSmear::match_reco_truth()
{
  struct candidate
  {
    float dR;
    size_t it;
    size_t ir;
  };
  std::vector<candidate> candidates;

  for (size_t it = 0; it < m_truth_jets.size(); it++)
  {
    for (size_t ir = 0; ir < m_calib_jets.size(); ir++)
    {
      float dR = fabs(getDR(m_calib_jets[ir], m_truth_jets[it]));
      if (dR < dr_cut)
      {
        candidates.push_back({dR, it, ir});
      }
    }
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const candidate &a, const candidate &b) { return a.dR < b.dR; });

  std::vector<bool> truth_used(m_truth_jets.size(), false);
  std::vector<Jet *> reco_match(m_calib_jets.size(), nullptr);

  for (const auto &c : candidates)
  {
    if (truth_used[c.it] || reco_match[c.ir])
    {
      continue;
    }
    reco_match[c.ir] = m_truth_jets[c.it];
    truth_used[c.it] = true;
  }

  m_matched_jets.reserve(m_calib_jets.size());
  for (size_t ir = 0; ir < m_calib_jets.size(); ir++)
  {
    m_matched_jets.emplace_back(reco_match[ir], m_calib_jets[ir]);
  }

  return;
}

double JetSmear::getDR(Jet *j1, Jet *j2)
{ 
  double dphi = fabs(j1->get_phi() - j2->get_phi());
  if (dphi > TMath::Pi())
    {
      dphi = 2*TMath::Pi() - dphi;
    }
  double dR = sqrt(TMath::Power(j1->get_eta() - j2->get_eta(), 2) + TMath::Power(dphi, 2));

  // Put in once calib jets have a vertex

  // double dR = dphi;
  // if (m_hasVertex)
  //   {
  //     dR = sqrt(TMath::Power(j1->get_eta() - j2->get_eta(), 2) + TMath::Power(dphi, 2));
  //   }
  return dR;
}

void JetSmear::smear_reco_jets()
{
  TAxis *xa = m_JerSmearBand->GetXaxis();
  const double x_lo = xa->GetBinCenter(1);
  const double x_hi = xa->GetBinCenter(xa->GetNbins());

  copy_container_info();

  // m_matched_jets is in the same order as m_calib_jets / m_raw_jets.
  for (size_t ijet = 0; ijet < m_matched_jets.size(); ijet++)
  {
    Jet *truth_jet = m_matched_jets[ijet].first;
    Jet *reco_jet = m_matched_jets[ijet].second;
    Jet *raw_jet = m_raw_jets[ijet];
    double scale = 1.;
    // Matched jets smear by the truth pT; unmatched jets (if enabled) by the reco pT.
    if (truth_jet || m_smear_unmatched_reco)
      {
	const double ref_pt = truth_jet ? truth_jet->get_pt() : reco_jet->get_pt();
	double lookup_pt = ref_pt;
	if (lookup_pt < x_lo) lookup_pt = x_lo;
	if (lookup_pt > x_hi) lookup_pt = x_hi;
	const double sigma_frac = m_JerSmearBand->Interpolate(lookup_pt);
	double smeared_pt = reco_jet->get_pt() + jetGaus(reco_jet->get_id()) * sigma_frac * ref_pt;
	// A negative smeared pT would flip the jet direction; floor it at 100 MeV.
	if (smeared_pt < 0)
	  {
	    smeared_pt = 0.1;
	  }
	scale = (reco_jet->get_pt() > 0) ? smeared_pt / reco_jet->get_pt() : 0.0;
      }
    Jet *out = m_output_jetcontainer->add_jet();
    copy_jet(reco_jet, raw_jet, out);
    out->set_px(reco_jet->get_px() * scale);
    out->set_py(reco_jet->get_py() * scale);
    out->set_pz(reco_jet->get_pz() * scale);
    out->set_e(out->get_e() * scale);
  }
}

// splitmix64 finalizer: maps a 64-bit value to a well-mixed 64-bit value.
uint64_t JetSmear::splitmix64(uint64_t x)
{
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

// Standard-normal draw for one jet, fully determined by (sample, run, event, jet id).
double JetSmear::jetGaus(int jet_id) const
{
  uint64_t h = splitmix64(static_cast<uint32_t>(m_sampleId));
  h = splitmix64(h ^ static_cast<uint32_t>(m_run));
  h = splitmix64(h ^ static_cast<uint32_t>(m_event));
  h = splitmix64(h ^ static_cast<uint32_t>(jet_id));

  // Two uniforms in (0,1) from the top 53 bits, then Box-Muller.
  const uint64_t h2 = splitmix64(h);
  const double u1 = (static_cast<double>(h >> 11) + 0.5) / 9007199254740992.0;
  const double u2 = (static_cast<double>(h2 >> 11) + 0.5) / 9007199254740992.0;
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

// Copy everything from a calib jet; the momentum and energy are scaled by the caller.
// JetCalib does not set E or the calo fractions on calib jets (they are NaN), so
// those are taken from the matching raw jet, with E scaled by calib pT / raw pT.
void JetSmear::copy_jet(Jet *in, Jet *raw, Jet *out)
{
  float e = in->get_e();
  float emcal_frac = in->get_emcal_frac();
  float ihcal_frac = in->get_ihcal_frac();
  float ohcal_frac = in->get_ohcal_frac();
  if (raw)
  {
    if (!std::isfinite(e) && raw->get_pt() > 0)
    {
      e = raw->get_e() * in->get_pt() / raw->get_pt();
    }
    if (!std::isfinite(emcal_frac)) emcal_frac = raw->get_emcal_frac();
    if (!std::isfinite(ihcal_frac)) ihcal_frac = raw->get_ihcal_frac();
    if (!std::isfinite(ohcal_frac)) ohcal_frac = raw->get_ohcal_frac();
  }

  out->set_id(in->get_id());
  out->set_px(in->get_px());
  out->set_py(in->get_py());
  out->set_pz(in->get_pz());
  out->set_e(e);
  out->set_isCalib(in->get_isCalib());
  out->set_emcal_frac(emcal_frac);
  out->set_ihcal_frac(ihcal_frac);
  out->set_ohcal_frac(ohcal_frac);
  out->insert_comp(in->get_comp_vec(), true);

  // Properties are stored by container index, so map input index -> output index.
  for (const auto &prop : m_calib_jetcontainer->property_indices())
  {
    out->set_property(m_output_jetcontainer->property_index(prop.first), in->get_property(prop.second));
  }
}

// Copy the container-level information (algorithm, R, sources, vertex, rho, properties).
void JetSmear::copy_container_info()
{
  m_output_jetcontainer->set_algo(m_calib_jetcontainer->get_algo());
  m_output_jetcontainer->set_jetpar_R(m_calib_jetcontainer->get_jetpar_R());
  for (auto src = m_calib_jetcontainer->begin_src(); src != m_calib_jetcontainer->end_src(); ++src)
  {
    m_output_jetcontainer->insert_src(*src);
  }
  m_output_jetcontainer->set_has_zvertex(m_calib_jetcontainer->has_zvertex());
  m_output_jetcontainer->set_vertex_type(m_calib_jetcontainer->get_vertex_type());
  m_output_jetcontainer->set_vertex_z(m_calib_jetcontainer->get_vertex_z());
  m_output_jetcontainer->set_rho_median(m_calib_jetcontainer->get_rho_median());
  for (const auto &prop : m_calib_jetcontainer->property_indices())
  {
    if (!m_output_jetcontainer->has_property(prop.first))
    {
      m_output_jetcontainer->add_property(prop.first);
    }
  }
}
