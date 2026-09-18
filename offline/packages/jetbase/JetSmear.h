// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETBASE_JETSMEAR_H
#define JETBASE_JETSMEAR_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class CDBHistos;
class PHCompositeNode;
class TRandom3;
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

  // Setters (both methods).
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
  
  void set_CalibFile(const std::string &path) { m_calibFileOverride = path; }  // local file override (default: CDB)

  
 private:
  // Functions (common).
  static std::string fetchCalibDir(const char *calibType);
  bool m_hasVertex{false};
  float m_event_reweight{0.0};
  // Matching
  double getDR(Jet *j1, Jet *j2);
  int initSmearing();
  void match_reco_truth();
  void smear_reco_jets();
  
  std::vector<Jet*> m_calib_jets;
  std::vector<Jet*> m_truth_jets;
  std::vector<std::pair<Jet*, Jet*>> m_matched_jets;
  
  // Input.
  std::string m_inputRecoNode{"AntiKt_Tower_r04"};
  std::string m_inputCalibNode{"AntiKt_Tower_r04_Calib"};
  std::string m_inputTruthNode{"AntiKt_Truth_r04"};
  std::string m_outputRecoNode{"AntiKt_Tower_r04_Calib_Smear"};

  JetContainer *m_reco_jetcontainer{nullptr};
  JetContainer *m_calib_jetcontainer{nullptr};
  JetContainer *m_truth_jetcontainer{nullptr};
  JetContainer *m_output_jetcontainer{nullptr};

  std::string m_calibFileOverride;              // EMfrac method: if set, read this file instead of the CDB payload
  float jet_radius{0.4};                        // Jet radius.
  float dr_cut{0.3};
  bool m_useEMfracCalib{true};                  // true = EMfrac method (default), false = legacy method
  bool m_useEMfracReweight{true};                  // true = EMfrac method (default), false = legacy method
  bool m_usePileup{true};
  int m_sys{0};
  TRandom3 *m_rng{nullptr};

  TFile *m_JetSmearFile{nullptr};
  TH1 *m_JerSmearBand{nullptr};
  TH1 *m_EMReweight{nullptr};   // h2d_jes_calib_r0R: truth-equivalent pT vs (reco pT, EMCal fraction)

};

#endif  // JETBASE_JETSMEAR_H
