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
#include "beamSpectrum.hxx"

using namespace std;

void beamSpectrum::doBeamSpectrum(double muonSmear, double pipmSmear , double neutronSmear, double protonSmear, double pi0Smear, int onlyMuon, int incNeutron, int fileLimit, TString filepath, int caseN ) {

  double neutronMass = 0.939565;
  double protonMass = 0.938272;
  double pipmMass = 0.13957;
  double pi0Mass = 0.13498;
  double muonMass = 0.105658;
  //double electronMass = 0.0005;

  double cAmount = 1;
  double cAmount2 = 1;
  double cAmount3 = 1;
  double cAmount4 = 1;

  TFile f("/dune/data/users/gyang/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root");
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

  std::cout<<targetDensity->Integral()<<" "<<BeamSigma->Integral()<<" "<<BeamOffsetX->Integral()<<" "<<BeamTheta->Integral()<<" "<<BeamThetaPhi->Integral()<<" "<<HC->Integral()<<" "<<WL->Integral()<<" "<<DPR->Integral()<<" "<<Horn1_XShift->Integral()<<" "<<Horn1_YShift->Integral()<<" "<<Horn1_X3mmShift->Integral()<<" "<<Horn2_XShift->Integral()<<" "<<Horn2_YShift->Integral()<<std::endl;

  TH1D* usage1[100];
  TH1D* usage2[100];
  TH1D* usage3[100];
  TH1D* usage4[100];
  TH1D* usage11[100];
  TH1D* usage12[100];
  TH1D* usage13[100];
  TH1D* usage14[100];
  //TGraph* gra[100];

  for(int i =0; i< 100; i++){
    usage1[i] = new TH1D("","",12,0,6);
    usage2[i] = new TH1D("","",12,0,6);
    usage3[i] = new TH1D("","",12,0,6);
    usage4[i] = new TH1D("","",12,0,6);
    usage11[i] = new TH1D("","",12,0,6);
    usage12[i] = new TH1D("","",12,0,6);
    usage13[i] = new TH1D("","",12,0,6);
    usage14[i] = new TH1D("","",12,0,6);

    //gra[i] = new TGraph(12);
  }

  for(int i=0; i<12; i++){
    for (int j=0;j<70; j++){
      usage1[j] -> SetBinContent(i+1, Horn1_XShift->GetBinContent(i+1,j+1) - norm->GetBinContent(i+1,j+1));
      usage2[j] -> SetBinContent(i+1, Horn1_YShift->GetBinContent(i+1,j+1) - norm->GetBinContent(i+1,j+1));
      usage3[j] -> SetBinContent(i+1, Horn2_XShift->GetBinContent(i+1,j+1) - norm->GetBinContent(i+1,j+1));
      usage4[j] -> SetBinContent(i+1, Horn2_YShift->GetBinContent(i+1,j+1) - norm->GetBinContent(i+1,j+1));
    }
  }

  double overallRate[100]={};
  double overallRate2[100]={};
  double overallRate3[100]={};
  double overallRate4[100]={};
  for(int ii=0;ii<70;ii++){
    double summ = 0;
    double summ2 = 0;
    double summ3 = 0;
    double summ4 = 0;
    for(int m=0;m<12;m++){
      summ += cAmount* usage1[ii]->GetBinContent(m+1) ;
      summ2 += cAmount2* usage2[ii]->GetBinContent(m+1) ;
      summ3 += cAmount3* usage3[ii]->GetBinContent(m+1) ;
      summ4 += cAmount4* usage4[ii]->GetBinContent(m+1) ;
    }
    //overallRate[ii] = summ/12.;
    //overallRate2[ii] = summ2/12.;
    //overallRate3[ii] = summ3/12.;
    //overallRate4[ii] = summ4/12.;
  }

  cout<<"overall rate "<<overallRate[0]<<" "<<overallRate2[0]<<endl;

  for(int ii=0;ii<70;ii++){
    for(int m=0;m<12;m++){
      usage11[ii]->SetBinContent(m+1, (norm->GetBinContent(m+1,ii+1) + cAmount*usage1[ii]->GetBinContent(m+1) - overallRate[ii])/norm->GetBinContent(m+1,ii+1) );
      usage12[ii]->SetBinContent(m+1, (norm->GetBinContent(m+1,ii+1) + cAmount2*usage2[ii]->GetBinContent(m+1) - overallRate2[ii])/norm->GetBinContent(m+1,ii+1) );
      usage13[ii]->SetBinContent(m+1, (norm->GetBinContent(m+1,ii+1) + cAmount3*usage3[ii]->GetBinContent(m+1) - overallRate3[ii])/norm->GetBinContent(m+1,ii+1) );
      usage14[ii]->SetBinContent(m+1, (norm->GetBinContent(m+1,ii+1) + cAmount4*usage4[ii]->GetBinContent(m+1) - overallRate4[ii])/norm->GetBinContent(m+1,ii+1) );
    }
  }


TH1F* h[10];

for (int i=0;i<10;i++)
{
  h[i] =  new TH1F(Form("shift%i",i),Form("shift%i",i),40,0,10);
}

h[0]->SetBinContent(1,0.994);
h[0]->SetBinContent(2,0.994);
h[0]->SetBinContent(3,0.99);
h[0]->SetBinContent(4,0.99);
h[0]->SetBinContent(5,0.991);
h[0]->SetBinContent(6,0.991);
h[0]->SetBinContent(7,0.993);
h[0]->SetBinContent(8,0.993);
h[0]->SetBinContent(9,0.993);
h[0]->SetBinContent(10,0.993);
h[0]->SetBinContent(11,0.994);
h[0]->SetBinContent(12,0.994);
h[0]->SetBinContent(13,0.994);
h[0]->SetBinContent(14,0.994);
h[0]->SetBinContent(15,0.996);
h[0]->SetBinContent(16,0.996);
h[0]->SetBinContent(17,0.997);
h[0]->SetBinContent(18,0.997);
h[0]->SetBinContent(19,1);
h[0]->SetBinContent(20,1);
h[0]->SetBinContent(21,1.012);
h[0]->SetBinContent(22,1.012);
h[0]->SetBinContent(23,1.025);
h[0]->SetBinContent(24,1.025);
h[0]->SetBinContent(25,1.02);
h[0]->SetBinContent(26,1.02);
h[0]->SetBinContent(27,1.019);
h[0]->SetBinContent(28,1.019);
h[0]->SetBinContent(29,1.021);
h[0]->SetBinContent(30,1.021);
h[0]->SetBinContent(31,1.017);
h[0]->SetBinContent(32,1.017);
h[0]->SetBinContent(33,1.023);
h[0]->SetBinContent(34,1.023);
h[0]->SetBinContent(35,1.017);
h[0]->SetBinContent(36,1.017);
h[0]->SetBinContent(37,1.017);
h[0]->SetBinContent(38,1.017);
h[0]->SetBinContent(39,1.016);
h[0]->SetBinContent(40,1.016);


h[1]->SetBinContent(1,1.002);
h[1]->SetBinContent(2,1.002);
h[1]->SetBinContent(3,1.001);
h[1]->SetBinContent(4,1.001);
h[1]->SetBinContent(5,0.998);
h[1]->SetBinContent(6,0.998);
h[1]->SetBinContent(7,0.998);
h[1]->SetBinContent(8,0.998);
h[1]->SetBinContent(9,1.002);
h[1]->SetBinContent(10,1.002);
h[1]->SetBinContent(11,0.994);
h[1]->SetBinContent(12,0.994);
h[1]->SetBinContent(13,0.983);
h[1]->SetBinContent(14,0.983);
h[1]->SetBinContent(15,0.976);
h[1]->SetBinContent(16,0.979);
h[1]->SetBinContent(17,0.973);
h[1]->SetBinContent(18,0.973);
h[1]->SetBinContent(19,0.968);
h[1]->SetBinContent(20,0.968);
h[1]->SetBinContent(21,0.952);
h[1]->SetBinContent(22,0.952);
h[1]->SetBinContent(23,0.974);
h[1]->SetBinContent(24,0.974);
h[1]->SetBinContent(25,0.989);
h[1]->SetBinContent(26,0.989);
h[1]->SetBinContent(27,0.986);
h[1]->SetBinContent(28,0.986);
h[1]->SetBinContent(29,0.987);
h[1]->SetBinContent(30,0.987);
h[1]->SetBinContent(31,0.993);
h[1]->SetBinContent(32,0.993);
h[1]->SetBinContent(33,1);
h[1]->SetBinContent(34,1);
h[1]->SetBinContent(35,1);
h[1]->SetBinContent(36,1);
h[1]->SetBinContent(37,0.998);
h[1]->SetBinContent(38,0.998);
h[1]->SetBinContent(39,1.005);
h[1]->SetBinContent(40,1.005);


h[2]->SetBinContent(1,0.995);
h[2]->SetBinContent(2,0.995);
h[2]->SetBinContent(3,0.99);
h[2]->SetBinContent(4,0.99);
h[2]->SetBinContent(5,0.989);
h[2]->SetBinContent(6,0.989);
h[2]->SetBinContent(7,0.996);
h[2]->SetBinContent(8,0.996);
h[2]->SetBinContent(9,1.001);
h[2]->SetBinContent(10,1.001);
h[2]->SetBinContent(11,0.998);
h[2]->SetBinContent(12,0.998);
h[2]->SetBinContent(13,0.987);
h[2]->SetBinContent(14,0.987);
h[2]->SetBinContent(15,0.985);
h[2]->SetBinContent(16,0.985);
h[2]->SetBinContent(17,0.985);
h[2]->SetBinContent(18,0.985);
h[2]->SetBinContent(19,0.993);
h[2]->SetBinContent(20,0.993);
h[2]->SetBinContent(21,1.005);
h[2]->SetBinContent(22,1.005);
h[2]->SetBinContent(23,1.027);
h[2]->SetBinContent(24,1.027);
h[2]->SetBinContent(25,1.033);
h[2]->SetBinContent(26,1.033);
h[2]->SetBinContent(27,1.045);
h[2]->SetBinContent(28,1.045);
h[2]->SetBinContent(29,1.046);
h[2]->SetBinContent(30,1.046);
h[2]->SetBinContent(31,1.047);
h[2]->SetBinContent(32,1.047);
h[2]->SetBinContent(33,1.053);
h[2]->SetBinContent(34,1.053);
h[2]->SetBinContent(36,1.045);
h[2]->SetBinContent(37,1.044);
h[2]->SetBinContent(38,1.044);
h[2]->SetBinContent(39,1.043);
h[2]->SetBinContent(40,1.043);


h[3]->SetBinContent(1,0.999);
h[3]->SetBinContent(2,0.999);
h[3]->SetBinContent(3,1);
h[3]->SetBinContent(4,1);
h[3]->SetBinContent(5,1);
h[3]->SetBinContent(6,1);
h[3]->SetBinContent(7,0.996);
h[3]->SetBinContent(8,0.996);
h[3]->SetBinContent(9,1.001);
h[3]->SetBinContent(10,1.001);
h[3]->SetBinContent(11,1.001);
h[3]->SetBinContent(12,1.001);
h[3]->SetBinContent(13,0.993);
h[3]->SetBinContent(14,0.993);
h[3]->SetBinContent(15,0.983);
h[3]->SetBinContent(16,0.983);
h[3]->SetBinContent(17,0.97);
h[3]->SetBinContent(18,0.97);
h[3]->SetBinContent(19,0.956);
h[3]->SetBinContent(20,0.956);
h[3]->SetBinContent(21,0.942);
h[3]->SetBinContent(22,0.942);
h[3]->SetBinContent(23,0.975);
h[3]->SetBinContent(24,0.975);
h[3]->SetBinContent(25,0.985);
h[3]->SetBinContent(26,0.985);
h[3]->SetBinContent(27,0.984);
h[3]->SetBinContent(28,0.984);
h[3]->SetBinContent(29,0.985);
h[3]->SetBinContent(30,0.985);
h[3]->SetBinContent(31,0.994);
h[3]->SetBinContent(32,0.994);
h[3]->SetBinContent(33,0.997);
h[3]->SetBinContent(34,0.997);
h[3]->SetBinContent(35,1);
h[3]->SetBinContent(36,1);
h[3]->SetBinContent(37,0.994);
h[3]->SetBinContent(38,0.994);
h[3]->SetBinContent(39,0.999);
h[3]->SetBinContent(40,0.999);


TGraph* graShift1 = new TGraph(12);
TGraph* graShift2 = new TGraph(12);
TGraph* graShift3 = new TGraph(12);
TGraph* graShift4 = new TGraph(12);

for(int m=0;m<12;m++){
  graShift1->SetPoint(m, usage11[0]->GetBinCenter(m+1), usage11[0]->GetBinContent(m+1) );
  graShift2->SetPoint(m, usage12[0]->GetBinCenter(m+1), usage12[0]->GetBinContent(m+1) );
  graShift3->SetPoint(m, usage13[0]->GetBinCenter(m+1), usage13[0]->GetBinContent(m+1) );
  graShift4->SetPoint(m, usage14[0]->GetBinCenter(m+1), usage14[0]->GetBinContent(m+1) );
  cout<<usage11[0]->GetBinContent(m+1)<<" "<<usage12[0]->GetBinContent(m+1)<<" "<<usage13[0]->GetBinContent(m+1)<<" "<<usage14[0]->GetBinContent(m+1)<<endl;
}
//TGraph graShift1(usage11[0]);
//TGraph graShift2(usage12[0]);


//graShift1 = (TH1F*)h[2];
//graShift2 = (TH1F*)h[3];

/*
  for(Int_t j=0;j< shift_numu3[0]->GetNbinsX();j++)
    shift_numu3[0]->SetBinContent(j+1, (shift_numu[0]->GetBinContent(j+1) + shift_numu2[0]->GetBinContent(j+1))/2.)  ;

  for(Int_t j=0;j< shift_numu3[1]->GetNbinsX();j++)
    shift_numu3[1]->SetBinContent(j+1, (shift_numu[1]->GetBinContent(j+1) + shift_numu2[1]->GetBinContent(j+1))/2.)  ;

  for(Int_t j=0;j< shift_numu3[2]->GetNbinsX();j++)
    shift_numu3[2]->SetBinContent(j+1, (shift_numu[2]->GetBinContent(j+1) + shift_numu2[2]->GetBinContent(j+1))/2.)  ;

  for(Int_t j=0;j< shift_numu3[3]->GetNbinsX();j++)
    shift_numu3[3]->SetBinContent(j+1, (shift_numu[3]->GetBinContent(j+1) + shift_numu2[3]->GetBinContent(j+1))/2.)  ;

  for(Int_t j=0;j< shift_numu3[4]->GetNbinsX();j++)
    shift_numu3[4]->SetBinContent(j+1, (shift_numu[4]->GetBinContent(j+1) + shift_numu2[4]->GetBinContent(j+1))/2.)  ;
  

  for(Int_t i=0;i<5;i++){
    shift_numu[i]->Divide(shift_numu3[i]);
    spline[i] = new TSpline3(shift_numu[i]);
  }
*/                 

   TFile *rootFile1 = new TFile(Form("Smearing-mu%0.2f-pipm%0.2f-n%0.2f-p%0.2f-pi0%0.2f-test_mu%d_case%d.root",muonSmear,pipmSmear,neutronSmear,protonSmear,pi0Smear, onlyMuon, caseN),"RECREATE");
 
   int num=0;
   int num1=0;
   string line;
   ifstream myfile(filepath+"test1.txt");
   while (getline(myfile, line)){
        ++num;
        ++num1;
        }
    num = fileLimit;
    std::cout << "Number of lines in text file: " << num1 <<endl;;
   if(num>num1) exit(1);
   const int n=num;
   //const int n1=num1;

   //TH1F* histo[n];
   //TH1F* histo1[n];
   //TH1F* histo22[n];
   TH1F* histo_true[n];
   TH1F* histo_reco_nominal[n];
   TH1F* histo_reco_shift1[n];
   TH1F* histo_reco_shift2[n];
   TH1F* histo_reco_shift3[n];
   TH1F* histo_reco_shift4[n];

   //TH1F* xpos[40];
   //TH1F* xposm[40];
   //TH1F* xposmm[40];
   TH1F* xpos_true[40];
   TH1F* xpos_reco_nominal[40];
   TH1F* xpos_reco_shift1[40];
   TH1F* xpos_reco_shift2[40];
   TH1F* xpos_reco_shift3[40];
   TH1F* xpos_reco_shift4[40];

   //TH1F* xposms[40]; 
   TH2F* histo2;

   //double EE[n];
   //double bins[n];
   //double binsm[n];
   //double binsmm[n];
   //double binsms[n];

   double bins_true[n];
   double bins_reco_nominal[n];
   double bins_reco_shift1[n];
   double bins_reco_shift2[n];
   double bins_reco_shift3[n];
   double bins_reco_shift4[n];

   int c=0;
   double max=0;
   double maxm=0;
   double min=0;
   double minm=0;
   //double maxms=0;
   //double minms=0;
   //double mean[40];
   //double rms[40];
   int totNFile = 0;
   string input1;
   ifstream infile;
   infile.open (filepath+"test1.txt");
   
  histo2 =  new TH2F("Energy vs Days","Energy vs Days",40,0,10,num,0,num);
  while(infile >> input1) // To get you all the lines.
  {          
    totNFile ++;
    if(totNFile> fileLimit) break; 
           
    cout << input1 <<endl;
    TFile input(Form("%s",input1.c_str()));
    TTree* h1 = (TTree*) input.Get("gRooTracker");
           
    Int_t StdHepN;
    Int_t StdHepPdg[1000];
    Int_t StdHepStatus[1000];
    Double_t vtx[4]={};
    Double_t ivtx[1000][4]={};
    Double_t imom[1000][4]={};
    //Bool_t flag=false;
    //Bool_t wflag=false;
    //Int_t rightSign=0,wrongSign=0;
    Int_t StdHepFm[50000];     // stdhep-like particle array: first mother
  
    h1->SetBranchAddress("StdHepPdg", &StdHepPdg);
    h1->SetBranchAddress("StdHepStatus", &StdHepStatus);
    h1->SetBranchAddress("EvtVtx", &vtx);
    h1->SetBranchAddress("StdHepX4", &ivtx);
    h1->SetBranchAddress("StdHepP4", &imom); 
    h1->SetBranchAddress("StdHepN", &StdHepN);
    h1->SetBranchAddress("StdHepFm",&StdHepFm);
  
    double inFV=0;
    double inFV1=0;
    double inFV2=0;
    int allEvt=0;
  
    //double vetoY=0.10;
    //double vetoX=0.10;
  
    rootFile1->cd();
    Int_t nentries = (Int_t) h1->GetEntries();
    //cout<< nentries <<endl;
    histo_true[c] =  new TH1F(Form("Energy Spectrum Day %i",c+1),Form("Energy Spectrum Day %i",c+1),40,0,10);
    histo_reco_nominal[c] =  new TH1F(Form("E_reco Spectrum Day %i",c+1),Form("E_reco Spectrum Day %i",c+1),40,0,10);

    histo_reco_shift1[c] =  new TH1F(Form("E_reco shift 1 Spectrum Day (M) %i",c+1),Form("E_reco shift 1 Spectrum Day (M) %i",c+1),40,0,10);
    histo_reco_shift2[c] =  new TH1F(Form("E_reco shift 2 Spectrum Day (M) %i",c+1),Form("E_reco shift 2 Spectrum Day (M) %i",c+1),40,0,10);
    histo_reco_shift3[c] =  new TH1F(Form("E_reco shift 3 Spectrum Day (M) %i",c+1),Form("E_reco shift 3 Spectrum Day (M) %i",c+1),40,0,10);
    histo_reco_shift4[c] =  new TH1F(Form("E_reco shift 4 Spectrum Day (M) %i",c+1),Form("E_reco shift 4 Spectrum Day (M) %i",c+1),40,0,10);

    //histo1s[c] =  new TH1F(Form("Muon Energy Spectrum Day %i - Smearing 4%%",c+1),Form("Muon Energy Spectrum Day %i - Smearing 4%%",c+1),40,0,10);

    for(Int_t i=0; i<nentries ; i++) 
    {
            
      h1->GetEntry(i);
      
      double muonMom = 0;
      double protonMom = 0;
      double pipmMom = 0;
      double pi0Mom = 0;
      //double neutronMom = 0;
      
      double nuE = 0;
      double Efsmuon = 0;
      double Efsproton = 0;
      double Efspipm = 0;
      double Efspi0 = 0;
      double Efsneutron = 0;

      
      if(vtx[0]>-1.1 && vtx[0]<1.1 && vtx[1]>-1.1 && vtx[1]<1.1 && vtx[2]<0.9 && vtx[2]>-0.9 ){

      inFV1++;
      for(int ip=0; ip<StdHepN; ip++) {
      inFV2++;
            
      
      if(StdHepPdg[ip]==14 && StdHepStatus[ip]== 0 && StdHepFm[ip]==-1){
        histo_true[c]->Fill(imom[ip][3]);
        histo2->Fill(imom[ip][3],c); 
        nuE=imom[ip][3];        
        inFV++;
      }
       

     if(StdHepPdg[ip]==13 && StdHepFm[ip]==0 && StdHepStatus[ip]==1){
       muonMom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
       if (muonMom > 0.05)
         Efsmuon=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,muonSmear)* muonMom,2) + TMath::Power(muonMass,2));
       //cout<<muonMom<<endl;
     }

    // take out neutron kinetic energy
    if(StdHepPdg[ip]==2112 && StdHepStatus[ip]==1){
      //neutronMom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
      //Efsneutron=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,neutronSmear)*neutronMom,2) + TMath::Power(neutronMass,2));
      // you lose all this energy
      Efsneutron = gRandom->Gaus(1,neutronSmear)*(imom[ip][3] - neutronMass);
    }
    
    if(StdHepPdg[ip]==2212 && StdHepStatus[ip]==1){
      protonMom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
      if (protonMom > 0.3 && gRandom->Uniform(0,1)<0.95)
        Efsproton=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,protonSmear)*protonMom,2) + TMath::Power(protonMass,2));
    }

    if(abs(StdHepPdg[ip])==211 && StdHepStatus[ip]==1){
      pipmMom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
      if (pipmMom > 0.1 && gRandom->Uniform(0,1)<0.9)
        Efspipm=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,pipmSmear)*pipmMom,2) + TMath::Power(pipmMass,2));
    }
    
    if(abs(StdHepPdg[ip])==111 && StdHepStatus[ip]==1){
      pi0Mom = TMath::Sqrt(imom[ip][0]*imom[ip][0]+imom[ip][1]*imom[ip][1]+imom[ip][2]*imom[ip][2]);
      Efspi0=TMath::Sqrt(TMath::Power(gRandom->Gaus(1,pi0Smear)*pi0Mom,2) + TMath::Power(pi0Mass,2));
    }
       
      //histo1s[c]->Fill(gRandom->Gaus(1,0.04)*imom[ip][3]);
      //cout<<"smearing: "<<gRandom->Gaus(1,0.04)*imom[ip][3]<<endl;
      //cout<<"energy: "<<imom[ip][3]<<endl;
      
    }
    if(nuE>0 && Efsmuon>0){
      if(onlyMuon==1){
        histo_reco_nominal[c]->Fill(Efsmuon );
        histo_reco_shift1[c]->Fill(Efsmuon , graShift1->Eval(nuE) );
        histo_reco_shift2[c]->Fill(Efsmuon , graShift2->Eval(nuE));
        histo_reco_shift3[c]->Fill(Efsmuon , graShift3->Eval(nuE) );
        histo_reco_shift4[c]->Fill(Efsmuon , graShift4->Eval(nuE));
      }
      else if(incNeutron==1){
        histo_reco_nominal[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 + Efsneutron - neutronMass);
        histo_reco_shift1[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 + Efsneutron - neutronMass, graShift1->Eval(nuE) );
        histo_reco_shift2[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 + Efsneutron - neutronMass, graShift2->Eval(nuE));
        histo_reco_shift3[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 + Efsneutron - neutronMass, graShift3->Eval(nuE) );
        histo_reco_shift4[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 + Efsneutron - neutronMass, graShift4->Eval(nuE));
      }
      else{
        histo_reco_nominal[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 - neutronMass);
        histo_reco_shift1[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 - neutronMass, graShift1->Eval(nuE) );
        histo_reco_shift2[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 - neutronMass, graShift2->Eval(nuE));
        histo_reco_shift3[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 - neutronMass, graShift3->Eval(nuE) );
        histo_reco_shift4[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 - neutronMass, graShift4->Eval(nuE));
      }
      //histo1[c]->Fill(Efsmuon + Efsproton + Efspipm + Efspi0 + Efsneutron -neutronMass );

      //histo22[c]->Fill( Efsmuon + Efsproton + Efspipm + Efspi0 + Efsneutron -neutronMass, 1.05 );
    }
      
    }
           
     allEvt ++;  
      
  }
           
  //histo[c]->Write();
  //histo1[c]->Write();
  //histo22[c]->Write();
  histo_true[c]->Write();
  histo_reco_nominal[c]->Write();
  histo_reco_shift1[c]->Write();
  histo_reco_shift2[c]->Write();
  histo_reco_shift3[c]->Write();
  histo_reco_shift4[c]->Write();

  c++;
  }
        
  histo2->Write();
  
  TH1F* xp_true =  new TH1F("Mean_Energy_Spectrum_per_Day","Mean_Energy_Spectrum_per_Day",40,0,10);
  TH1F* xp_reco_nominal =  new TH1F("Mean_E_reco_Energy_Spectrum_per_Day","Mean_E_reco_Energy_Spectrum_per_Day",40,0,10);
  TH1F* xp_reco_shift1 =  new TH1F("Mean_E_reco_shift1_Energy_Spectrum_per_Day(M)","Mean_E_reco_shift1_Energy_Spectrum_per_Day(M)",40,0,10);
  TH1F* xp_reco_shift2 =  new TH1F("Mean_E_reco_shift2_Energy_Spectrum_per_Day(M)","Mean_E_reco_shift2_Energy_Spectrum_per_Day(M)",40,0,10);
  TH1F* xp_reco_shift3 =  new TH1F("Mean_E_reco_shift3_Energy_Spectrum_per_Day(M)","Mean_E_reco_shift3_Energy_Spectrum_per_Day(M)",40,0,10);
  TH1F* xp_reco_shift4 =  new TH1F("Mean_E_reco_shift4_Energy_Spectrum_per_Day(M)","Mean_E_reco_shift4_Energy_Spectrum_per_Day(M)",40,0,10);

  //TH1F* xpms =  new TH1F("Mean Muon Energy Spectrum per Day - Smearing 4%","Mean Muon Energy Spectrum per Day - Smearing 4%",40,0,10);
  //TH1F* xp1 =  new TH1F("Mean RMS per Day","Mean RMS per Day",40,0,10);
  //TH1F* xp1m =  new TH1F("Mean E_reco RMS per Day","Mean E_reco RMS per Day",40,0,10);
  //TH1F* xp1ms =  new TH1F("Mean Muon RMS per Day - Smearing 4%","Mean Muon RMS per Day - Smearing 4%",40,0,10);
  
  for(int b=0;b<40;b++) 
  {
  
  //cout<<"N: "<<b<<endl;
             
  for(int p=0;p<num;p++)      
  {
  
  bins_true[p]=histo_true[p]->GetBinContent(b+1);
  bins_reco_nominal[p]=histo_reco_nominal[p]->GetBinContent(b+1);
  bins_reco_shift1[p]=histo_reco_shift1[p]->GetBinContent(b+1);
  bins_reco_shift2[p]=histo_reco_shift2[p]->GetBinContent(b+1);
  bins_reco_shift3[p]=histo_reco_shift3[p]->GetBinContent(b+1);
  bins_reco_shift4[p]=histo_reco_shift4[p]->GetBinContent(b+1);

  //binsms[p]=histo1s[p]->GetBinContent(b+1);
  //cout<<"Bin: "<<bins[p]<<endl;
  //cout<<"Bin 1: "<<binsm[p]<<endl;
  
  }
  
   std::vector<double> v(bins_true, bins_true + num);
   max = *max_element(v.begin(), v.end());
   min = *min_element(v.begin(), v.end());
   
   std::vector<double> v1(bins_reco_nominal, bins_reco_nominal + num);
   maxm = *max_element(v1.begin(), v1.end());
   minm = *min_element(v1.begin(), v1.end());
   
   //std::vector<double> v1s(binsms, binsms + num);
   //maxms = *max_element(v1s.begin(), v1s.end());
   //minms = *min_element(v1s.begin()+1, v1s.end());
   
//   cout<<"Max value: "<<max<<endl;
//   cout<<"Min value: "<<min<<endl;
  cout<<"________________________"<<endl;
   xpos_true[b] =  new TH1F(Form("EventsvsDays(%0.2fGeV)",0.25*(b+1)),Form("EventsvsDays(%0.2fGeV)",0.25*(b+1)),20,min,max+(max-min)/20);
   xpos_reco_nominal[b] =  new TH1F(Form("E_recoEventsvsDays (%0.2f GeV)",0.25*(b+1)),Form("E_recoEventsvsDays(%0.2fGeV)",0.25*(b+1)),40,minm*0.5,(maxm+(max-min)/40)*2);
   xpos_reco_shift1[b] =  new TH1F(Form("E_recoShift1EventsvsDays(%0.2fGeV)(M)",0.25*(b+1)),Form("E_recoShift1EventsvsDays(%0.2fGeV)(M)",0.25*(b+1)),40,minm*0.5,(maxm+(max-min)/40)*2);
   xpos_reco_shift2[b] =  new TH1F(Form("E_recoSsfift2EventsvsDays(%0.2fGeV)(M)",0.25*(b+1)),Form("E_recoShift2EventsvsDays(%0.2fGeV)(M)",0.25*(b+1)),40,minm*0.5,(maxm+(max-min)/40)*2);
   xpos_reco_shift3[b] =  new TH1F(Form("E_recoShift3EventsvsDays(%0.2fGeV)(M)",0.25*(b+1)),Form("E_recoShift3EventsvsDays(%0.2fGeV)(M)",0.25*(b+1)),40,minm*0.5,(maxm+(max-min)/40)*2);
   xpos_reco_shift4[b] =  new TH1F(Form("E_recoSsfift4EventsvsDays(%0.2fGeV)(M)",0.25*(b+1)),Form("E_recoShift4EventsvsDays(%0.2fGeV)(M)",0.25*(b+1)),40,minm*0.5,(maxm+(max-min)/40)*2);

   //xposms[b] =  new TH1F(Form("Muon Events vs Days (%0.2f GeV) - Smearing 4%%",0.25*(b+1)),Form("Muon Events vs Days (%0.2f GeV) - Smearing 4%%",0.25*(b+1)),20,minms,maxms);
   for(int r=0;r<num;r++) {
     xpos_true[b]->Fill(bins_true[r]); 
     xpos_reco_nominal[b]->Fill(bins_reco_nominal[r]),
     xpos_reco_shift1[b]->Fill(bins_reco_shift1[r]);
     xpos_reco_shift2[b]->Fill(bins_reco_shift2[r]);
     xpos_reco_shift3[b]->Fill(bins_reco_shift3[r]);
     xpos_reco_shift4[b]->Fill(bins_reco_shift4[r]);
   }
   xpos_true[b]->Write();
   xpos_reco_nominal[b]->Write();
   xpos_reco_shift1[b]->Write();
   xpos_reco_shift2[b]->Write();
   xpos_reco_shift3[b]->Write();
   xpos_reco_shift4[b]->Write();
   //xposms[b]->Write();
   //mean[b]=xpos[b]->GetMean();
   //rms[b]=xpos[b]->GetRMS();   
   
   xp_true->SetBinContent(b+1,xpos_true[b]->GetMean());
   xp_true->SetBinError(b+1,xpos_true[b]->GetRMS());

   xp_reco_nominal->SetBinContent(b+1,xpos_reco_nominal[b]->GetMean());
   xp_reco_nominal->SetBinError(b+1,xpos_reco_nominal[b]->GetRMS());

   xp_reco_shift1->SetBinContent(b+1,xpos_reco_shift1[b]->GetMean());
   xp_reco_shift1->SetBinError(b+1,xpos_reco_shift1[b]->GetRMS());

   xp_reco_shift2->SetBinContent(b+1,xpos_reco_shift2[b]->GetMean());
   xp_reco_shift2->SetBinError(b+1,xpos_reco_shift2[b]->GetRMS());

   xp_reco_shift3->SetBinContent(b+1,xpos_reco_shift3[b]->GetMean());
   xp_reco_shift3->SetBinError(b+1,xpos_reco_shift3[b]->GetRMS());

   xp_reco_shift4->SetBinContent(b+1,xpos_reco_shift4[b]->GetMean());
   xp_reco_shift4->SetBinError(b+1,xpos_reco_shift4[b]->GetRMS());

   //xpmm->SetBinContent(b+1,xposmm[b]->GetMean());
   //xpmm->SetBinError(b+1,xposmm[b]->GetMean());

   //xpms->SetBinContent(b+1,xposms[b]->GetMean());
   //xpms->SetBinError(b+1,xposms[b]->GetRMS());
      
   //xp1->SetBinContent(b+1,xpos[b]->GetRMS());
   //xp1m->SetBinContent(b+1,xposm[b]->GetRMS());
   //xp1ms->SetBinContent(b+1,xposms[b]->GetRMS());
  }

  //gPad->SaveAs("plot.pdf");
  xp_true->Write();
  xp_reco_nominal->Write();
  xp_reco_shift1->Write();
  xp_reco_shift2->Write();
  xp_reco_shift3->Write();
  xp_reco_shift4->Write();

  rootFile1->Close();
   
}
   


