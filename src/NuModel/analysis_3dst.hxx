#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TRandom.h"
#include "TVectorD.h"
#include <iostream>
#include <vector>
#include "TMath.h"

using namespace std;
class MCMC {

public:

  TH2D* h2_trueE_trueY_nominal;
  TH2D* h2_recoE_recoY_nominal;

  TH2D* h2_nuwro_true;
  TH2D* h2_Ar_diff;

  TH2D* h2_trueE_trueY_var[20];
  TH2D* h2_recoE_recoY_var[20];

  TH2D* h2_trueE_trueY_nominal_CH;
  TH2D* h2_recoE_recoY_nominal_CH;
  TH2D* h2_nuwro_true_CH;
  TH2D* h2_Ar_diff_CH;

  TH2D* h2_trueE_trueY_var_CH[20];
  TH2D* h2_recoE_recoY_var_CH[20];

  TTree* genieAr;
  TTree* genieCH;
  TTree* nuwroAr;
  TTree* nuwroCH;

  Double_t wgt_MaCCQE[10], wgt_VecFFCCQEshape[10], wgt_MaCCRES[10], wgt_MvCCRES[10], wgt_AhtBY[10], wgt_BhtBY[10], wgt_CV1uBY[10],
         wgt_CV2uBY[10], wgt_FrElas_pi[10], wgt_FtInel_pi[10], wgt_FrAbs_pi[10], wgt_FormZone[10], wgt_FrPiProd_pi[10],
         wgt_MFP_N[10], wgt_FrCEx_N[10], wgt_CCQEEPauliSupViaKF[10], wgt_Mnv2p2hGaussEnhancement[10], wgt_MKSPP_ReWeight[10],
         wgt_E2p2h_A_nu[10], wgt_E2p2h_B_nu[10], wgt_E2p2h_A_nubar[10], wgt_E2p2h_B_nubar[10], wgt_BePRA_A[10], wgt_BePRA_B[10],
         wgt_BePRA_D[10], wgt_BePRA_E[10], wgt_C12ToAr40_2p2hScaling_nu[10], wgt_C12ToAr40_2p2hScaling_nubar[10],
         wgt_nuenuebar_xsec_ratio[10],  wgt_nuenumu_xsec_ratio[10], wgt_SPPLowQ2Suppression[10], wgt_FSILikeEAvailSmearing[10];

  Double_t Ev;
  Double_t YY;
  Double_t Ev_reco;
  Double_t Elep_reco, eP;
  Int_t isFD, isCC, nuPDG, LepPDG;
  Double_t LepE, eDepP, eDepN, eDepPip, eDepPim, eDepPi0, eDepOther;

  double finI,finII,finIII,finIIII,finJ,finJJ,finJJJ,finJJJJ,finK,finKK,finKKK,finKKKK;

  double Ev_CH, YY_CH;


  Double_t EvtVtx[3], StdHepP4[100][4];
  int StdHepPdg[100];
  double nuE, lepE;
  double pEnergy;
  double nuE_reco;

  void beauty();
  void process_Ar();
  void MCMC_master();
  double MCMC_processing(std::vector<double> currList);
  void applyToCH();
  void setIterationTime(int time){iterationTime = time;}
  int getIterationTime(){return iterationTime;}

  int iterationTime;

};
