{

  TFile f("/home/guang/work/DUNEPrismTools/FD_numode_OptimizedEngineeredNov2017Review_w_PPFX_fit_binning.root");
  TFile f2("/home/guang/work/DUNEPrismTools/FD_numode_OptimizedEngineeredNov2017Review_w_PPFX_fit_binning.root");
  TH1D* LBNF_nuebar_flux_Nom = (TH1D*)f.Get("LBNF_nuebar_flux_Nom");
  TH1D* LBNF_nue_flux_Nom = (TH1D*)f.Get("LBNF_nue_flux_Nom");
  TH1D* LBNF_numu_flux_Nom = (TH1D*)f.Get("LBNF_numu_flux_Nom");
  TH1D* var = (TH1D*)f2.Get("LBNF_numu_flux_Nom");
  TH1D* LBNF_numubar_flux_Nom = (TH1D*)f.Get("LBNF_numubar_flux_Nom");

  double stepDM = 1.8e-3/50.;
  double stepST = 0.5/50.;

  TFile* outf = TFile::Open("oscillated.root","recreate");

  for(int ii=0;ii<50;ii++){
    for(int jj=0;jj<50;jj++){

      double theta = TMath::ASin(sqrt(0.3 + stepST * ii));
      double dm2 = 1.8e-3 + stepDM * jj;
      double L = 1300;

      for(Int_t i=0;i< LBNF_numu_flux_Nom->GetNbinsX(); i++){
        double prob = 1 - TMath::Power(TMath::Sin(2*theta),2) * TMath::Power(TMath::Sin(1.267 * dm2 * L / LBNF_numu_flux_Nom->GetBinCenter(i+1)),2);
        var->SetBinContent(i+1, prob * LBNF_numu_flux_Nom->GetBinContent(i+1));
      }
      var->Write(Form("var_st%d_dm%d",ii,jj));
    }
  }

  new TCanvas();
  LBNF_numu_flux_Nom->Draw();

}
