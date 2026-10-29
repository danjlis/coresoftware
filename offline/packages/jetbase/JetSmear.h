// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETBASE_JETSMEAR_H
#define JETBASE_JETSMEAR_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>
#include <vector>

class CDBHistos;
class PHCompositeNode;
class JetContainer;
class Jet;
class TF1;
class TH1;
class TH2;
class TFile;

class JetSmear : public SubsysReco
{
 public:
  explicit JetSmear(const std::string &name = "JetSmear");
  ~JetSmear() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int CreateNodeTree(PHCompositeNode *topNode);

  // Required: node names and jet radius have no defaults and must be set.
  void set_InputTruthNode(const std::string &inputNode) { m_inputTruthNode = inputNode; }
  void set_InputRecoNode(const std::string &inputNode) { m_inputRecoNode = inputNode; }
  void set_InputCalibNode(const std::string &inputNode) { m_inputCalibNode = inputNode; }
  void set_OutputRecoNode(const std::string &outputNode) { m_outputRecoNode = outputNode; }
  // Jet Radius
  void set_JetRadius(float radius) { jet_radius = radius; }

  // Method selection: true (default) = EMfrac method, false = legacy method.
  void set_UseEMfracCalib(bool use) { m_useEMfracCalib = use; }
  void set_UsePileupCorr(bool use) { m_usePileup = use; }
  void set_Systematic(int sys) { m_sys = sys; }
  // true = smear unmatched reco jets using their own reco pT; false (default) = leave them unsmeared.
  void set_SmearUnmatchedReco(bool smear) { m_smear_unmatched_reco = smear; }
  
  // Smearing seed: each jet's random draw is a hash of (sample id, run, event, jet id),
  // so it is reproducible and independent of job splitting. The systematic is not part of the
  // hash, so nominal/sysup/sysdown use the same draw per jet and only the width changes.
  // Set a distinct sample id for each simulation sample in a run.
  void set_SampleId(int id) { m_sampleId = id; }

  void set_CalibFile(const std::string &path) { m_calibFileOverride = path; }  // local file override (default: CDB)

  
 private:
  // Functions (common).
  static std::string fetchCalibDir(const char *calibType);
  bool m_hasVertex{false};
  // Matching
  double getDR(Jet *j1, Jet *j2);
  int initSmearing();
  void match_reco_truth();
  void smear_reco_jets();
  void copy_jet(Jet *in, Jet *raw, Jet *out);
  static uint64_t splitmix64(uint64_t x);
  double jetGaus(int jet_id) const;
  void copy_container_info();
  
  std::vector<Jet*> m_calib_jets;
  std::vector<Jet*> m_raw_jets;  // raw jet for each entry of m_calib_jets (nullptr if none)
  std::vector<Jet*> m_truth_jets;
  std::vector<std::pair<Jet*, Jet*>> m_matched_jets;
  
  // Input.
  std::string m_inputRecoNode;
  std::string m_inputCalibNode;
  std::string m_inputTruthNode;
  std::string m_outputRecoNode;

  JetContainer *m_reco_jetcontainer{nullptr};
  JetContainer *m_calib_jetcontainer{nullptr};
  JetContainer *m_truth_jetcontainer{nullptr};
  JetContainer *m_output_jetcontainer{nullptr};

  std::string m_calibFileOverride;              // EMfrac method: if set, read this file instead of the CDB payload
  float jet_radius{-1};                         // Jet radius (required, no default).
  float dr_cut{1.0};
  bool m_useEMfracCalib{true};                  // true = EMfrac method (default), false = legacy method
  bool m_usePileup{true};
  int m_sys{0};
  bool m_smear_unmatched_reco{false};
  int m_sampleId{0};
  int m_run{0};    // from EventHeader, per event
  int m_event{0};  // from EventHeader, per event

  TH1 *m_JerSmearBand{nullptr};  // owned clone; the CDB file is closed after loading

};

#endif  // JETBASE_JETSMEAR_H
