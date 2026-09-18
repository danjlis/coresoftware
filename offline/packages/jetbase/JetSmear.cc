#include "JetSmear.h"

#include "Jet.h"
#include "JetContainerv1.h"

#include <cdbobjects/CDBTF.h>  // for CDBTF1 (legacy method)

#include <ffamodules/CDBInterface.h>

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
#include <TH2.h>
#include <TRandom3.h>

#include <algorithm>  // for std::max, std::min
#include <cmath>
#include <iostream>  // for operator<<, basic_ostream
#include <string>

JetSmear::JetSmear(const std::string &name)
  : SubsysReco(name)
{

  m_rng = new TRandom3(0);

}

JetSmear::~JetSmear()
{
  if (m_JetSmearFile)
  {
    m_JetSmearFile->Close();
    delete m_JetSmearFile;
  }
}

int JetSmear::InitRun(PHCompositeNode *topNode)
{
  // Create calib jet node.
  int ret = CreateNodeTree(topNode);

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
    std::cout << "JetSmear::initSmearing() : No Jer Smear available! Will apply no smear." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_JetSmearFile = TFile::Open(calibFile.c_str(), "READ");
  if (!m_JetSmearFile || m_JetSmearFile->IsZombie())
  {
    std::cout << "JetSmear::initEMfracCalib() : Could not open calibration file " << calibFile << "!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
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
  const std::string histo_name = "h_jer_smear_" + r_suffix + pileup_suffix + jes_suffix + systematic_suffix;

  // First-pass map (required).
  m_JerSmearBand = dynamic_cast<TH1 *>(m_JetSmearFile->Get(histo_name.c_str()));
  if (!m_JerSmearBand)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const std::string reweight_name = "h_em_reweight" + r_suffix + pileup_suffix + jes_suffix;
  
  m_EMReweight = dynamic_cast<TH1 *>(m_JetSmearFile->Get(reweight_name.c_str()));
  if (!m_EMReweight)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}


int JetSmear::process_event(PHCompositeNode *topNode)
{
  m_event_reweight = 1.0;
  m_calib_jets.clear();
  m_truth_jets.clear();
  m_matched_jets.clear();
  
  // Get raw and calib jet nodes.
  m_reco_jetcontainer = findNode::getClass<JetContainer>(topNode, m_inputRecoNode);
  if (!m_reco_jetcontainer)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_calib_jetcontainer = findNode::getClass<JetContainer>(topNode, m_inputCalibNode);
  if (!m_calib_jetcontainer)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_truth_jetcontainer = findNode::getClass<JetContainer>(topNode, m_inputTruthNode);
  if (!m_truth_jetcontainer)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_output_jetcontainer = findNode::getClass<JetContainer>(topNode, m_outputRecoNode);

  if (!m_output_jetcontainer)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  float zvertex = 0;
  
  m_hasVertex = false;

  m_hasVertex = m_reco_jetcontainer->has_zvertex();
  if (m_hasVertex)
    {
      zvertex = m_reco_jetcontainer->get_vertex_z();
      if (!std::isfinite(zvertex))
      {
        m_hasVertex = false;
        zvertex = 0;
      }
    }

  // Apply calibration to jets and add to calib jet node.
  int ijet = 0;
  float max_reco_jet = 0;
  float max_reco_jet_emfrac = 0;
  m_calib_jets.reserve(m_calib_jetcontainer->size());
  for (auto jet : *m_calib_jetcontainer)
  {
    Jet *tempjet = m_reco_jetcontainer->get_jet(ijet);
    
    if (jet->get_pt() > max_reco_jet)
      {
	max_reco_jet_emfrac = tempjet->get_emcal_frac();
	max_reco_jet = jet->get_pt();
      }
	
    m_calib_jets.push_back(jet);
    ijet++;
  }

  int binLead = m_EMReweight->GetXaxis()->FindBin(max_reco_jet_emfrac);

  binLead = std::max(1, std::min(binLead, m_EMReweight->GetNbinsX()));

  m_event_reweight = m_EMReweight->GetBinContent(binLead);

  
  
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
  PHCompositeNode *antiktNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "ANTIKT"));
  if (!antiktNode)
  {
    std::cout << PHWHERE << "JetSmear::CreateNodeTree : ANTIKT node missing, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHCompositeNode *towerNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "TOWER"));
  if (!towerNode)
  {
    std::cout << PHWHERE << "TOWER node not found, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
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

void JetSmear::match_reco_truth()
{

  std::map<int,float> reco_dR;
  std::map<int,float> truth_dR;

  for (auto tjet : m_truth_jets)
    {
      auto truth_iter = std::find_if(m_matched_jets.begin(), m_matched_jets.end(), [=] (auto a) { return a.first == tjet;});
      float truth_closest_jet = dr_cut;
      float truth_closest_jet_pt = 0;
      bool already_matched = false;
      int truth_id = tjet->get_id();
      if (truth_iter != m_matched_jets.end())
	{
	  already_matched = true;
	  truth_closest_jet = truth_dR[truth_id];
	  truth_closest_jet_pt = truth_iter->second->get_pt();
	}

      
      for (auto jet : m_calib_jets)
	{
	  auto reco_iter = std::find_if(m_matched_jets.begin(), m_matched_jets.end(), [=] (auto a) { return a.second == jet;});
	  float reco_closest_jet = dr_cut;
	  float reco_closest_jet_pt = 0;
	  bool reco_already_matched = false;

	  float dR = fabs(getDR(jet, tjet));
	  bool overwrite = true;
	  int reco_id = jet->get_id();
	  if (reco_iter != m_matched_jets.end())
	    {
	      reco_already_matched = true;
	      reco_closest_jet = reco_dR[reco_id];
	      reco_closest_jet_pt = reco_iter->first->get_pt();
	      overwrite &= reco_closest_jet > dR && reco_closest_jet_pt < tjet->get_pt();
	    }
	  if (already_matched)
	    {
	      overwrite &= truth_closest_jet > dR && truth_closest_jet_pt < jet->get_pt();
	    }
	  if (dR < dr_cut && overwrite)
	    {
	      if (already_matched)
		{
		  auto truth_iter2 = std::find_if(m_matched_jets.begin(), m_matched_jets.end(), [=] (auto a) { return a.first == tjet;});
		  m_matched_jets.erase(truth_iter2);
		}
	      if (reco_already_matched)
		{
		  auto reco_iter2 = std::find_if(m_matched_jets.begin(), m_matched_jets.end(), [=] (auto a) { return a.second == jet;});
		  m_matched_jets.erase(reco_iter2);
		}

	      truth_dR[truth_id] = dR;
	      reco_dR[reco_id] = dR;
	      already_matched = true;
	      truth_closest_jet = dR;
	      truth_closest_jet_pt = jet->get_pt();

	      m_matched_jets.push_back(std::make_pair(tjet, jet));	
	      break;
	    }	  
	}
    }
  
  // collect the rest
  for (auto tjet : m_truth_jets)
    {
      auto truth_iter = std::find_if(m_matched_jets.begin(), m_matched_jets.end(), [=] (auto b) { return b.first == tjet;});

      if (truth_iter != m_matched_jets.end())
	{
	  continue;
	}

      int truth_id = tjet->get_id();
      for (auto jet : m_calib_jets)
	{
	  auto reco_iter = std::find_if(m_matched_jets.begin(), m_matched_jets.end(), [=] (auto b) { return b.second == jet;});

	  if (reco_iter != m_matched_jets.end())
	    {
	      continue;
	    }
	  int reco_id = jet->get_id();	  
	  float dR = fabs(getDR(jet, tjet));

	  if (dR < dr_cut)
	    {
	      truth_dR[truth_id] = dR;
	      reco_dR[reco_id] = dR;
	      m_matched_jets.push_back(std::make_pair(tjet, jet));	
	      break;
	    }	  
	}
    }
  for (auto jet : m_calib_jets)
    {
      auto reco_iter = std::find_if(m_matched_jets.begin(), m_matched_jets.end(), [=] (auto b) { return b.second == jet;});
      
      if (reco_iter != m_matched_jets.end())
	{
	  continue;
	}
      m_matched_jets.push_back(std::make_pair(nullptr, jet));	
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

  for (auto &matched : m_matched_jets)
  {
    Jet *reco_jet = matched.first;
    Jet *truth_jet = matched.second;
    double scale = 1.;
    if (truth_jet)
      {
	const double truth_pt = truth_jet->get_pt();
	double lookup_pt = truth_pt;
	if (lookup_pt < x_lo) lookup_pt = x_lo;
	if (lookup_pt > x_hi) lookup_pt = x_hi;
	const double sigma_frac = m_JerSmearBand->Interpolate(lookup_pt);
	const double smeared_pt = m_rng->Gaus(reco_jet->get_pt(), sigma_frac * truth_pt);
	scale = (reco_jet->get_pt() > 0) ? smeared_pt / reco_jet->get_pt() : 0.0;
      }
    Jet *out = m_output_jetcontainer->add_jet();
    out->set_px(reco_jet->get_px() * scale);
    out->set_py(reco_jet->get_py() * scale);
    out->set_pz(reco_jet->get_pz() * scale);
    out->set_e(reco_jet->get_e() * scale);
    out->set_reweight(m_event_reweight);
    out->set_id(reco_jet->get_id());

    
  }
}
