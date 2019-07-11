#include "fluxRate.hxx"

using namespace std;

void fluxRate::doFluxRate(TString inputShift ){

  TFile f(inputShift);
  TH2D* norm = (TH2D*)f.Get("ND_nu_ppfx/LBNF_numu_flux_Nom");
  TH2D* targetDensity = (TH2D*)f.Get("ND_nu_TargetDensity_p1/LBNF_numu_flux");
  TH2D* BeamSigma = (TH2D*)f.Get("ND_nu_BeamSigma_p1/LBNF_numu_flux");
  TH2D* BeamOffsetX = (TH2D*)f.Get("ND_nu_BeamOffsetX_p1/LBNF_numu_flux");
  TH2D* BeamTheta = (TH2D*)f.Get("ND_nu_BeamTheta_p1/LBNF_numu_flux");
  TH2D* BeamThetaPhi =  (TH2D*)f.Get("ND_nu_BeamThetaPhi_p1/LBNF_numu_flux");
  TH2D* HC =  (TH2D*)f.Get("ND_nu_HC_p1/LBNF_numu_flux");
  TH2D* WL =  (TH2D*)f.Get("ND_nu_WL_p1/LBNF_numu_flux");
  TH2D* DPR =  (TH2D*)f.Get("ND_nu_DPR_p1/LBNF_numu_flux");
  TH2D* Horn1_XShift =  (TH2D*)f.Get("ND_nu_Horn1_XShift/LBNF_numu_flux");
  TH2D* Horn1_YShift =  (TH2D*)f.Get("ND_nu_Horn1_YShift/LBNF_numu_flux");
  TH2D* Horn1_X3mmShift =  (TH2D*)f.Get("ND_nu_Horn1_X3mmShift/LBNF_numu_flux");
  TH2D* Horn2_XShift =  (TH2D*)f.Get("ND_nu_Horn2_XShift/LBNF_numu_flux");
  TH2D* Horn2_YShift =  (TH2D*)f.Get("ND_nu_Horn2_YShift/LBNF_numu_flux");

  TH2D* fd_norm = (TH2D*)f.Get("FD_nu_ppfx/LBNF_numu_flux_Nom");
  TH2D* fd_targetDensity = (TH2D*)f.Get("FD_nu_TargetDensity_p1/LBNF_numu_flux");
  TH2D* fd_BeamOffsetX = (TH2D*)f.Get("FD_nu_BeamOffsetX_p1/LBNF_numu_flux");
  TH2D* fd_BeamTheta = (TH2D*)f.Get("FD_nu_BeamTheta_p1/LBNF_numu_flux");
  TH2D* fd_BeamThetaPhi =  (TH2D*)f.Get("FD_nu_BeamThetaPhi_p1/LBNF_numu_flux");
  TH2D* fd_HC =  (TH2D*)f.Get("FD_nu_HC_p1/LBNF_numu_flux");
  TH2D* fd_WL =  (TH2D*)f.Get("FD_nu_WL_p1/LBNF_numu_flux");
  TH2D* fd_DPR =  (TH2D*)f.Get("FD_nu_DPR_p1/LBNF_numu_flux");
  TH2D* fd_Horn1_XShift =  (TH2D*)f.Get("FD_nu_Horn1_XShift/LBNF_numu_flux");
  TH2D* fd_Horn1_YShift =  (TH2D*)f.Get("FD_nu_Horn1_YShift/LBNF_numu_flux");
  TH2D* fd_Horn1_X3mmShift =  (TH2D*)f.Get("FD_nu_Horn1_X3mmShift/LBNF_numu_flux");
  TH2D* fd_Horn2_XShift =  (TH2D*)f.Get("FD_nu_Horn2_XShift/LBNF_numu_flux");
  TH2D* fd_Horn2_YShift =  (TH2D*)f.Get("FD_nu_Horn2_YShift/LBNF_numu_flux");

  cout<<fd_norm->Integral()<<" "<<fd_targetDensity->Integral()<<" "<<fd_BeamOffsetX->Integral()<<" "<<fd_BeamTheta->Integral()<<" "<<fd_BeamThetaPhi->Integral()<<" "<<fd_HC->Integral()<<" "<<fd_WL->Integral()<<" "<<fd_DPR->Integral()<<" "<<fd_Horn1_XShift->Integral()<<" "<<fd_Horn1_YShift->Integral()<<" "<< fd_Horn1_X3mmShift->Integral()<<" "<< fd_Horn2_XShift->Integral()<<" "<<fd_Horn2_YShift->Integral()<<endl;


  TH1D* chg1[100];
  TH1D* chg2[100];
  TH1D* chg3[100];
  TH1D* chg4[100];
  TH1D* chg5[100];
  TH1D* chg6[100];
  TH1D* chg7[100];
  TH1D* chg8[100];
  TH1D* chg9[100];
  TH1D* chg10[100];
  TH1D* chg11[100];
  TH1D* chg12[100];
  TH1D* chg13[100];
  TH1D* usage1[100];
  TH1D* usage2[100];
  TH1D* usage3[100];
  TH1D* usage4[100];

  double binEdge[100];
  for(Int_t i=0;i<targetDensity->GetNbinsX();i++){
    binEdge[i] = targetDensity->GetXaxis()->GetBinLowEdge(i+1);
    cout<<binEdge[i]<<endl;
  }

  int nBins = targetDensity->GetNbinsX();
  cout<<" .. "<< nBins<<endl;

  for(int i =0; i< 100; i++){
    chg1[i] = new TH1D("","",nBins-1,binEdge);
    chg2[i] = new TH1D("","",nBins-1,binEdge);
    chg3[i] = new TH1D("","",nBins-1,binEdge);
    chg4[i] = new TH1D("","",nBins-1,binEdge);
    chg5[i] = new TH1D("","",nBins-1,binEdge);
    chg6[i] = new TH1D("","",nBins-1,binEdge);
    chg7[i] = new TH1D("","",nBins-1,binEdge);
    chg8[i] = new TH1D("","",nBins-1,binEdge);
    chg9[i] = new TH1D("","",nBins-1,binEdge);
    chg10[i] = new TH1D("","",nBins-1,binEdge);
    chg11[i] = new TH1D("","",nBins-1,binEdge);
    chg12[i] = new TH1D("","",nBins-1,binEdge);
    chg13[i] = new TH1D("","",nBins-1,binEdge);
    usage1[i] = new TH1D("","",nBins-1,binEdge);
    usage2[i] = new TH1D("","",nBins-1,binEdge);
    usage3[i] = new TH1D("","",nBins-1,binEdge);
    usage4[i] = new TH1D("","",nBins-1,binEdge);    
  }

  cout<<"3"<<endl;
  
  TH1D* cc[20];
  for(Int_t i=0;i<20;i++){
    cc[i] = new TH1D("","",70,0,35.);
  }  

  double weekNorm=0;
  for(int i=0;i<nBins-1;i++){
    weekNorm += norm->GetBinContent(i+1,1);
  }

  double weekScale = (7*37000*(1./(2.2*2.2*1.8))) / weekNorm;
  //double monthScale = (30*37000*(1./(2.2*2.2*1.8))) / weekNorm;

  cout<<5<<endl;

  for(int i=0; i<nBins-1; i++){
    for (int j=0;j<70; j++){
      cout<<11<<endl;
      chg1[j]-> SetBinContent(i+1, targetDensity->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      cout<<12<<endl;
      chg2[j]-> SetBinContent(i+1, BeamSigma->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg3[j]-> SetBinContent(i+1, BeamOffsetX->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      cout<<13<<endl;
      chg4[j]-> SetBinContent(i+1, BeamTheta->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg5[j] -> SetBinContent(i+1, BeamThetaPhi->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      cout<<14<<endl;
      chg6[j] -> SetBinContent(i+1, HC->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg7[j] -> SetBinContent(i+1, WL->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg8[j] -> SetBinContent(i+1, DPR->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg9[j] -> SetBinContent(i+1, Horn1_XShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg10[j] -> SetBinContent(i+1, Horn1_YShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg11[j] -> SetBinContent(i+1, Horn1_X3mmShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg12[j] -> SetBinContent(i+1, Horn2_XShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg13[j] -> SetBinContent(i+1, Horn2_YShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      cout<<9<<endl;
      usage1[j] -> SetBinContent(i+1, Horn1_XShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      usage2[j] -> SetBinContent(i+1, Horn1_YShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      usage3[j] -> SetBinContent(i+1, Horn2_XShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      usage4[j] -> SetBinContent(i+1, Horn2_YShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
    }
  }
  cout<<"6"<<endl;

  for(int i=0;i<70;i++){
    double tot1 = 0;
    for(int j=0;j<usage1[i]->GetNbinsX();j++){
      tot1 += usage1[i]->GetBinContent(j+1, i+1);
    }	  
    cc[0]-> SetBinContent(i+1, tot1*weekScale);
  }	  

  for(int i=0;i<70;i++){
    double tot1 = 0;
    for(int j=0;j<usage2[i]->GetNbinsX();j++){
      tot1 += usage2[i]->GetBinContent(j+1, i+1);
    }
    cc[1]-> SetBinContent(i+1, tot1*weekScale);
  }

  for(int i=0;i<70;i++){
    double tot1 = 0;
    for(int j=0;j<usage3[i]->GetNbinsX();j++){
      tot1 += usage3[i]->GetBinContent(j+1, i+1);
    }
    cc[2]-> SetBinContent(i+1, tot1*weekScale);
  }

  for(int i=0;i<70;i++){
    double tot1 = 0;
    for(int j=0;j<usage4[i]->GetNbinsX();j++){
      tot1 += usage4[i]->GetBinContent(j+1, i+1);
    }
    cc[3]-> SetBinContent(i+1, tot1*weekScale);
  }


  TCanvas* c1 = new TCanvas();
  c1->cd();
  cc[1]->SetTitle("integrated variations vs. off-axis position");
  cc[1]->GetXaxis()->SetTitle("off-axis location [m]");
  cc[1]->GetYaxis()->SetTitle("weekly scint. rate var. / ton");
  for(int i=0;i<4;i++){
    cc[i]->SetLineColor(i+1);
    cc[i]->SetLineWidth(3);
  }
  cc[1]->Draw("hist");
  cc[0]->Draw("same");
  cc[2]->Draw("same");
  cc[3]->Draw("same");
  TLegend* legend = new TLegend(0.58, 0.6, 0.9, 0.9);
  legend->AddEntry(cc[0],"Horn 1 X shift 0.5 mm","l");
  legend->AddEntry(cc[1],"Horn 1 Y shift 0.5 mm","l");
  legend->AddEntry(cc[2],"Horn 2 X shift 0.5 mm","l");
  legend->AddEntry(cc[3],"Horn 2 Y shift 0.5 mm","l");
  legend->Draw();  

}
