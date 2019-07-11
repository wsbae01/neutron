#include "TString.h"
#include "fluxChange.hxx"

#include <string>
#include <utility>
#include <vector>

#include "Riostream.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <limits>
#include <getopt.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TRandom.h"
#include <TMath.h>
#include "TSpline.h"

using namespace std;

void fluxChange::doFluxChange(TString inputShift, TString fileNDflux, TString fileFDflux, double cAmount, int st, int dm){

  TFile f(inputShift);
  TH2D* norm = (TH2D*)f.Get("ND_nu_ppfx/LBNF_numu_flux_Nom");
  TH2D* targetDensity = (TH2D*)f.Get("ND_nu_TargetDensity_p1/LBNF_numu_flux");
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
  cout<<"number of energy bins "<<norm->GetNbinsX()<<endl;
  cout<<"off-axis position bins "<<norm->GetNbinsY()<<endl;

  TH1D* td[100];
  TH1D* bo[100];
  TH1D* bt[100];
  TH1D* chg4[100]; 
  TH1D* chg5[100];
  TH1D* chg6[100];
  TH1D* chg7[100];
  TH1D* chg8[100];
  TH1D* chg9[100];
  TH1D* chg10[100];
  TH1D* chg11[100];
  TH1D* chg12[100];
  TH1D* diff[100];
  TH1D* usage1[100];
  TH1D* usage2[100];
  TH1D* usage3[100];
  TGraph* gra[100];

  for(int i =0; i< 100; i++){
    td[i] = new TH1D("","",12,0,6);
    bo[i] = new TH1D("","",12,0,6);
    bt[i] = new TH1D("","",12,0,6);
    chg4[i] = new TH1D("","",12,0,6);
    chg5[i] = new TH1D("","",12,0,6);
    chg6[i] = new TH1D("","",12,0,6);
    chg7[i] = new TH1D("","",12,0,6);
    chg8[i] = new TH1D("","",12,0,6);
    chg9[i] = new TH1D("","",12,0,6);
    chg10[i] = new TH1D("","",12,0,6);
    chg11[i] = new TH1D("","",12,0,6);
    chg12[i] = new TH1D("","",12,0,6);
    diff[i] = new TH1D("","",12,0,6);
    usage1[i] = new TH1D("","",12,0,6);
    usage2[i] = new TH1D("","",12,0,6);
    usage3[i] = new TH1D("","",12,0,6);
    gra[i] = new TGraph(12);
  }
  TH1D* onAxis = new TH1D("","",12,0,6);
  TH1D* onDiff = new TH1D("","",12,0,6);
  TGraph* onGra = new TGraph(12);

  //double count1=0;
  //double cAmount = 2;
  //int st = 20;
  //int dm = 20;

  cout<<"setting up the shifts "<<endl;
  for(int i=0; i<12; i++){
    for (int j=12;j<70; j++){
      td[j]-> SetBinContent(i+1, targetDensity->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      bo[j]-> SetBinContent(i+1, BeamOffsetX->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      bt[j]-> SetBinContent(i+1, BeamTheta->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg4[j] -> SetBinContent(i+1, BeamThetaPhi->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg5[j] -> SetBinContent(i+1, HC->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg6[j] -> SetBinContent(i+1, WL->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg7[j] -> SetBinContent(i+1, DPR->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg8[j] -> SetBinContent(i+1, Horn1_XShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg9[j] -> SetBinContent(i+1, Horn1_YShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg10[j] -> SetBinContent(i+1, Horn1_X3mmShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg11[j] -> SetBinContent(i+1, Horn2_XShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      chg12[j] -> SetBinContent(i+1, Horn2_YShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));

      usage1[j] -> SetBinContent(i+1, Horn1_X3mmShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      usage2[j] -> SetBinContent(i+1, Horn1_YShift->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
      usage3[j] -> SetBinContent(i+1, HC->GetBinContent(i+1,j+1)- norm->GetBinContent(i+1,j+1));
    }

    onAxis -> SetBinContent(i+1, fd_Horn1_X3mmShift->GetBinContent(i+1,1)- fd_norm->GetBinContent(i+1,1));
  }
    
  double minn[100] ;
  for(int i=0;i<100;i++){
    minn[i] = 1000.00;
  }

  double overallRate[100]={};
  double onAxisRate = 0;

  for(int ii=12;ii<70;ii++){
    double summ = 0;
    for(int m=0;m<12;m++){
      summ += cAmount* usage1[ii]->GetBinContent(m+1) ;
    }
    if(abs(summ) < minn[ii]){
      minn[ii] = abs(summ);
    }
    overallRate[ii] = summ/12.;
  }

  for(int m=0;m<12;m++){
    onAxisRate += cAmount* onAxis->GetBinContent(m+1) ;
  }
  onAxisRate /= 12.;

  for(int ii=12;ii<70;ii++){
    for(int m=0;m<12;m++){
      diff[ii]->SetBinContent(m+1, cAmount*usage1[ii]->GetBinContent(m+1) - overallRate[ii] );
      gra[ii]->SetPoint(m, usage1[ii]->GetBinCenter(m+1), cAmount*usage1[ii]->GetBinContent(m+1) - overallRate[ii] );
    }
  }

  for(int m=0;m<12;m++){
    onDiff ->SetBinContent(m+1, cAmount*onAxis->GetBinContent(m+1) - onAxisRate );
    onGra->SetPoint(m, onAxis->GetBinCenter(m+1), cAmount*onAxis->GetBinContent(m+1) - onAxisRate );
  }

  cout<<"getting files "<<endl;

  TFile f0(fileNDflux);
  TH2D* oriNumu = (TH2D*)f0.Get("LBNF_numu_flux_Nom");

  TFile ff0(fileNDflux);
  TH2D* oriNumu1 = (TH2D*)ff0.Get("LBNF_numu_flux_Nom");  

  TFile fff0(fileNDflux);
  TH2D* oriNumu2 = (TH2D*)fff0.Get("LBNF_numu_flux_Nom");
  oriNumu2 = oriNumu2;

  TFile ff(fileNDflux);
  TH2D* varNumu = (TH2D*) ff.Get("LBNF_numu_flux_Nom");

  TFile ffno(fileNDflux);
  TH2D* varNumu2 = (TH2D*) ffno.Get("LBNF_numu_flux_Nom");  

  TFile fdd(fileFDflux);
  TH1D* fdNom = (TH1D*)fdd.Get("LBNF_numu_flux_Nom");
 
  TFile fdd1(fileFDflux);
  TH1D* fRatio1 = (TH1D*)fdd1.Get("LBNF_numu_flux_Nom");  

  TFile fdd2(fileFDflux);
  TH1D* fRatio2 = (TH1D*)fdd2.Get("LBNF_numu_flux_Nom");

  TFile fdd3(fileFDflux);
  TH1D* fRatio3 = (TH1D*)fdd3.Get("LBNF_numu_flux_Nom");

  TFile fdd4(fileFDflux);
  TH1D* fRatio4 = (TH1D*)fdd4.Get("LBNF_numu_flux_Nom");

  TFile ffff(Form("output_st%d_dm%d.root",st,dm));
  TH1D* coef = (TH1D*) ffff.Get("Coeffs");

  TFile fd("oscillated.root");
  TH1D* varFD = (TH1D*) fd.Get(Form("var_st%d_dm%d",st,dm));

  //TFile* output = TFile::Open("shiftedHist.root","recreate");

  cout<<"applying variations "<<endl;

  for(int i=0;i<varNumu->GetNbinsX();i++){
    for(int j=0;j<varNumu->GetNbinsY();j++){
      if(j >=12 && j<70)
      {
        varNumu->SetBinContent(i+1,j+1, varNumu->GetBinContent(i+1,j+1) + gra[j]->Eval(varNumu->ProjectionX()->GetBinCenter(i+1) ) );
      }
      else{
        varNumu->SetBinContent(i+1,j+1,varNumu->GetBinContent(i+1,j+1));
      }
      oriNumu1 -> SetBinContent(i+1,j+1, varNumu->GetBinContent(i+1,j+1) - oriNumu->GetBinContent(i+1,j+1));
    }
  }

  TH1D* example1d[16];
  TH1D* exampleNom[16];
  for(int i =0;i<16;i++){
    example1d[i] = new TH1D("","",oriNumu1->GetNbinsX(), oriNumu1->ProjectionX()->GetBinLowEdge(1), oriNumu1->ProjectionX()->GetBinLowEdge(oriNumu1->GetNbinsX())+oriNumu1->ProjectionX()->GetBinWidth(oriNumu1->GetNbinsX()));
    exampleNom[i] = new TH1D("","",oriNumu1->GetNbinsX(), oriNumu1->ProjectionX()->GetBinLowEdge(1), oriNumu1->ProjectionX()->GetBinLowEdge(oriNumu1->GetNbinsX())+oriNumu1->ProjectionX()->GetBinWidth(oriNumu1->GetNbinsX()));
  }
  for(int ii=0;ii<16;ii++){
    for(int i =0 ;i<example1d[0]->GetNbinsX();i++){
      example1d[ii]->SetBinContent(i+1, varNumu->GetBinContent(i+1,ii+13));
      exampleNom[ii]->SetBinContent(i+1, oriNumu->GetBinContent(i+1,ii+13));
    }
  }

  TCanvas* c100 = new TCanvas();
  c100->Divide(4,4);
  for(int i=0;i<16;i++){
    c100->cd(i+1);
    example1d[i]->SetTitle(Form("off axis %f m",6.5+i*0.5));
    example1d[i]->SetLineColor(4);
    example1d[i]->GetXaxis()->SetTitle("Energy");
    example1d[i]->Draw();
    exampleNom[i]->SetLineColor(1);
    exampleNom[i]->Draw("same");
  }

  TH1D* onSpec = new TH1D("","",varFD->GetNbinsX(),0,10);
  TH1D* onSpec1 = new TH1D("","",varFD->GetNbinsX(),0,10);
  TH1D* onSpec2 = new TH1D("","",varFD->GetNbinsX(),0,10);
  TH1D* onSpec3 = new TH1D("","",varFD->GetNbinsX(),0,10);
  TH1D* onSpec4 = new TH1D("","",varFD->GetNbinsX(),0,10);

  for(int i=0;i<onSpec->GetNbinsX();i++){
    onSpec -> SetBinContent(i+1, varFD->GetBinContent(i+1) * 0.5 + (varFD->GetBinContent(i+1) + onGra->Eval(varNumu->ProjectionX()->GetBinCenter(i+1) )) * 0.5 );	  
    onSpec1 -> SetBinContent(i+1, varFD->GetBinContent(i+1) * 0.3 + (varFD->GetBinContent(i+1) + onGra->Eval(varNumu->ProjectionX()->GetBinCenter(i+1) )) * 0.7 );
    onSpec2 -> SetBinContent(i+1, varFD->GetBinContent(i+1) * 0.5 + (varFD->GetBinContent(i+1) + onGra->Eval(varNumu->ProjectionX()->GetBinCenter(i+1) )) * 0.5 );
    onSpec3 -> SetBinContent(i+1, varFD->GetBinContent(i+1) * 0.6 + (varFD->GetBinContent(i+1) + onGra->Eval(varNumu->ProjectionX()->GetBinCenter(i+1) )) * 0.4 );
    onSpec4 -> SetBinContent(i+1, varFD->GetBinContent(i+1) * 1 + (varFD->GetBinContent(i+1) + onGra->Eval(varNumu->ProjectionX()->GetBinCenter(i+1) )) * 0 );
    //cout<<onSpec->GetBinContent(i+1)<<endl;
  }
  // varNumu1:  24 - 70 : varNumu->SetBinContent(i+1,j+1, varNumu->GetBinContent(i+1,j+1)*(1 + gra[j]->Eval(varNumu->ProjectionX()->GetBinCenter(i+1) )) );
  // varNumu2:  12 - 70 : same
  // varNumu3:  12 - 70 : varNumu->SetBinContent(i+1,j+1,0);
  // varNumu4:  12 - 70 : varNumu->SetBinContent(i+1,j+1, varNumu->GetBinContent(i+1,j+1) *0.9 );

  varNumu->Write("varNumu",TObject::kOverwrite);
 
  TH1D* chSpec = new TH1D("","",varFD->GetNbinsX(),0,10);
  TH1D* chSpec2 = new TH1D("","",varFD->GetNbinsX(),0,10);

  for(int i=0;i<varNumu->GetNbinsX();i++){
    double tVal = 0;
    double tVal2 = 0;
    for(int j=0;j<varNumu->GetNbinsY();j++){  
      tVal += varNumu ->GetBinContent(i+1,j+1) * coef->GetBinContent(j+1);
      tVal2 += varNumu2 ->GetBinContent(i+1,j+1) * coef->GetBinContent(j+1);
    }
    chSpec ->SetBinContent ( i+1, tVal ); 
    chSpec2->SetBinContent ( i+1, tVal2 );
  }

  //cout<<"total number of bins "<<varNumu->GetNbinsX()<<" "<<varFD->GetNbinsX()<<" "<<fdNom->GetNbinsX()<<endl;
  for(int i=0;i<varNumu->GetNbinsX();i++){
    fRatio1->SetBinContent(i+1, (onSpec1->GetBinContent(i+1) - chSpec->GetBinContent(i+1) ) / fdNom->GetBinContent(i+1));
    fRatio2->SetBinContent(i+1, (onSpec2->GetBinContent(i+1) - chSpec->GetBinContent(i+1) ) / fdNom->GetBinContent(i+1));
    fRatio3->SetBinContent(i+1, (onSpec3->GetBinContent(i+1) - chSpec->GetBinContent(i+1) ) / fdNom->GetBinContent(i+1));
    fRatio4->SetBinContent(i+1, (onSpec4->GetBinContent(i+1) - chSpec->GetBinContent(i+1) ) / fdNom->GetBinContent(i+1));
  }  

  TGraph* grr1 = new TGraph(onSpec->GetNbinsX());
  TGraph* grr2 = new TGraph(chSpec->GetNbinsX());
  //TGraph* grr3 = new TGraph(chSpec2->GetNbinsX());

  for(int i=0; i<onSpec->GetNbinsX();i++){
    grr1->SetPoint(i, onSpec->GetBinCenter(i+1), onSpec->GetBinContent(i+1) / chSpec2->GetBinContent(i+1));
    grr2->SetPoint(i, onSpec->GetBinCenter(i+1), chSpec->GetBinContent(i+1) / chSpec2->GetBinContent(i+1));
  }

  TH1D* er1 = new TH1D("","",oriNumu1->GetNbinsX(), oriNumu1->ProjectionX()->GetBinLowEdge(1), oriNumu1->ProjectionX()->GetBinLowEdge(oriNumu1->GetNbinsX())+oriNumu1->ProjectionX()->GetBinWidth(oriNumu1->GetNbinsX()));
  TH1D* er2 = new TH1D("","",oriNumu1->GetNbinsX(), oriNumu1->ProjectionX()->GetBinLowEdge(1), oriNumu1->ProjectionX()->GetBinLowEdge(oriNumu1->GetNbinsX())+oriNumu1->ProjectionX()->GetBinWidth(oriNumu1->GetNbinsX()));
  TH1D* er3 = new TH1D("","",oriNumu1->GetNbinsX(), oriNumu1->ProjectionX()->GetBinLowEdge(1), oriNumu1->ProjectionX()->GetBinLowEdge(oriNumu1->GetNbinsX())+oriNumu1->ProjectionX()->GetBinWidth(oriNumu1->GetNbinsX()));

  TFile nonswap("FD_FHC_nonswap.root");
  TTree* t = (TTree*)nonswap.Get("caf");

  double Ev; //E_reco;
  t->SetBranchAddress("Ev",&Ev);

  double stepDM = 1.8e-3/50.;
  double stepST = 0.5/50.;
  double theta = TMath::ASin(sqrt(0.3 + stepST * st));
  double dm2 = 1.8e-3 + stepDM * dm;
  double L = 1300;

  cout<<"theta and dm2 "<<TMath::Sin(theta) * TMath::Sin(theta)<<" "<<dm2<<endl;

  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    double prob = 1 - TMath::Power(TMath::Sin(2*theta),2) * TMath::Power(TMath::Sin(1.267 * dm2 * L / Ev),2);
    er1->Fill(Ev, prob);
    er2->Fill(Ev, grr1->Eval(Ev)* prob);
    er3->Fill(Ev, grr2->Eval(Ev)* prob);    
  }

  new TCanvas();
  er1->SetLineColor(2);
  er1->GetXaxis()->SetTitle("Energy (GeV)");
  er1->GetYaxis()->SetTitle("Spectrum");
  er1->Draw("hist");
  er2->Draw("same");
  er2->SetLineColor(4);
  er3->SetLineColor(3);
  er3->Draw("same");
  TLegend* legend = new TLegend(0.58, 0.6, 0.9, 0.9);
  legend->AddEntry(er1,"nominal ND extrapolated FD spectrum","l");
  legend->AddEntry(er2,"shifted FD spectrum","l");
  legend->AddEntry(er3,"shifted ND extrapolated FD spectrum","l");
  legend->Draw();

  new TCanvas();
  fRatio1->GetYaxis()->SetTitle("(true - ND extrapolated) / unoscillated");
  fRatio1->GetXaxis()->SetTitle("Energy (GeV)");
  fRatio1->SetLineColor(3);
  fRatio1->SetLineWidth(3);
  fRatio2->SetLineColor(4);
  fRatio2->SetLineWidth(4);
  fRatio3->SetLineColor(2);
  fRatio3->SetLineWidth(3);
  fRatio4->SetLineColor(1);
  fRatio4->SetLineWidth(3);
  //fRatio1->Draw();
  fRatio2->Draw("");
  //fRatio3->Draw("same");
  //fRatio4->Draw("same");
  legend = new TLegend(0.58, 0.6, 0.9, 0.9);
  //legend->AddEntry(fRatio1,"30\% on-axis","l");
  legend->AddEntry(fRatio2,"50\% on-axis","l");
  //legend->AddEntry(fRatio3,"60\% on-axis","l");
  //legend->AddEntry(fRatio4,"70\% on-axis","l");
  legend->Draw();

  new TCanvas();
  onSpec->SetLineColor(4);
  chSpec->GetXaxis()->SetTitle("Energy (GeV)");
  chSpec->GetYaxis()->SetTitle("Flux");
  chSpec->Draw("hist");
  onSpec->Draw("same");
  chSpec2->SetLineColor(2);
  chSpec->SetLineColor(3);
  chSpec2->Draw("same");
  onSpec4->SetLineColor(1);
  //onSpec4->Draw("same");
  legend = new TLegend(0.58, 0.6, 0.9, 0.9);
  legend->AddEntry(chSpec,"shifted ND extrapolated FD flux","l");
  legend->AddEntry(chSpec2,"nominal ND extrapolated FD flux","l");
  legend->AddEntry(onSpec,"shifted FD flux","l");
  //legend->AddEntry(onSpec4,"nominal FD flux","l");
  legend->Draw();

  new TCanvas();
  varNumu->SetTitle("varied");
  varNumu->Draw("colz");

  new TCanvas();
  oriNumu->SetTitle("Original");
  oriNumu->Draw("colz");

  new TCanvas();
  oriNumu1->SetTitle("Difference");
  oriNumu1->Draw("colz");

  TCanvas* c0[10];
  for(int ii=0;ii<10;ii++){
    c0[ii] = new TCanvas(Form("offAxis_%03fm",6+ii*0.5),Form("offAxis_%03fm",6+ii*0.5));
    c0[ii]->Divide(3,4);
    c0[ii]->cd(1);
    td[12+ii]->SetTitle("Target density");
    td[12+ii]->Draw();
    c0[ii]->cd(2);
    bo[12+ii]->SetTitle("Beam offset X");
    bo[12+ii]->Draw();
    c0[ii]->cd(3);
    bt[12+ii]->SetTitle("Beam theta");
    bt[12+ii]->Draw();
    c0[ii]->cd(4);
    chg4[12+ii]->SetTitle("Beam theta phi");
    chg4[12+ii]->Draw();
    c0[ii]->cd(5);
    chg5[12+ii]->SetTitle("HC");
    chg5[12+ii]->Draw();
    c0[ii]->cd(6);
    chg6[12+ii]->SetTitle("WL");
    chg6[12+ii]->Draw();
    c0[ii]->cd(7);
    chg7[12+ii]->SetTitle("DPR");
    chg7[12+ii]->Draw();
    c0[ii]->cd(8);
    chg8[12+ii]->SetTitle("Horn1_XShift");
    chg8[12+ii]->Draw();
    c0[ii]->cd(9);
    chg9[12+ii]->SetTitle("Horn1_YShift");
    chg9[12+ii]->Draw();
    c0[ii]->cd(10);
    chg10[12+ii]->SetTitle("Horn1_X3mmShift");
    chg10[12+ii]->Draw();  
    c0[ii]->cd(11);
    chg11[12+ii]->SetTitle("Horn2_XShift");
    chg11[12+ii]->Draw();
    c0[ii]->cd(12);
    chg12[12+ii]->SetTitle("Horn2_YShift");
    chg12[12+ii]->Draw();
  }

  TCanvas* c1 = new TCanvas();
  c1->Divide(4,4);
  for(int i=0;i<16;i++){
    c1->cd(i+1);
    diff[12+i]->SetTitle(Form("off axis %f m",6+i*0.5));
    diff[12+i]->Draw();   
  }

  TCanvas* c2 = new TCanvas();
  c2->Divide(4,4);
  for(int i=0;i<16;i++){
    c2->cd(i+1);
    gra[12+i]->GetXaxis()->SetTitle("Energy (GeV)");
    gra[12+i]->GetYaxis()->SetTitle("flux shift");
    gra[12+i]->SetTitle(Form("off axis %f m",6+i*0.5));
    gra[12+i]->Draw();
  }


}
