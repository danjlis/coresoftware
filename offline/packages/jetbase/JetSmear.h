// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETBASE_JETSMEAR_H
#define JETBASE_JETSMEAR_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class CDBHistos;
class PHCompositeNode;
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
  void set_OutputRecoNode(const std::string &outputNode) { m_outputRecoNode = outputNode; }
  // Jet Radius
  void set_JetRadius(float radius) { jet_radius = radius; }

  // Method selection: true (default) = EMfrac method, false = legacy method.
  void set_UseEMfracCalib(bool use) { m_useEMfracCalib = use; }
  void set_UsePileupCorr(bool use) { m_usePileup = use; }
  void set_UseReweight(bool use) { m_useReweight = use; }

  void set_CalibFile(const std::string &path) { m_calibFileOverride = path; }  // local file override (default: CDB)

  
 private:
  // Functions (common).
  static std::string fetchCalibDir(const char *calibType);

  // Matching
  void match_reco_truth();
  std::vector<Jet*> m_reco_jets;
  std::vector<Jet*> m_truth_jets;
  
  // Input.
  std::string m_inputRecoNode{"AntiKt_Tower_r04_Calib"};
  std::string m_inputTruthNode{"AntiKt_Truth_r04"};
  std::string m_outputRecoNode{"AntiKt_Tower_r04_Calib_Smear"};
  std::string m_calibFileOverride;              // EMfrac method: if set, read this file instead of the CDB payload
  float jet_radius{0.4};                        // Jet radius.
  bool m_useEMfracCalib{true};                  // true = EMfrac method (default), false = legacy method
  bool m_useEMfracReweight{true};                  // true = EMfrac method (default), false = legacy method

  TFile *m_JetSmearFile{nullptr};
  TH1 *m_JerSmearBand{nullptr};
  TH3 *m_EMReweight{nullptr};   // h2d_jes_calib_r0R: truth-equivalent pT vs (reco pT, EMCal fraction)
  TH2 *m_ZetaCorrMap{nullptr};   // h2_zeta_corr_r0R: scale factor vs (signed z-vertex, f)
  TF1 *m_NozCorrFunc{nullptr};   // f_zeta_noz_corr_r0R: scale factor vs f (f at z = 0), no-vertex events
  bool m_warnedNoEmFrac{false};

  // Calibration objects (legacy method).
  CDBTF *m_LegacyCalibFile{nullptr};
  std::vector<std::vector<TF1 *>> m_JetSmearFunc{};
};

#endif  // JETBASE_JETSMEAR_H
