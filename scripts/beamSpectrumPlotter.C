{

  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat("");

  gStyle->SetLabelFont(102,"");
  gStyle->SetLabelSize(0.06,"");
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.001,"x");
  gStyle->SetLabelOffset(0.01,"y");

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(0.9,"x");
  gStyle->SetTitleOffset(0.9,"y");

  gStyle->SetStripDecimals(kFALSE);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetStatW(0.35);
  gStyle->SetStatH(0.25);

  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(99);

  gStyle->SetHistLineWidth(2);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetFuncWidth(2);

  gStyle->SetStatFont(42);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0000);

  int nDay = 7;

  TFile f(Form("%s",gApplication->Argv(4)));
  TH1D* mean1 = (TH1D*)f.Get("Mean_E_reco_Energy_Spectrum_per_Day");

  cout<<"point 1"<<endl;
  TH1D* shift1;
  TH1D* shift2;
  TH1D* shift3;
  TH1D* shift4;
  shift1 = (TH1D*)f.Get("Mean_E_reco_shift1_Energy_Spectrum_per_Day(M)");
  shift2 = (TH1D*)f.Get("Mean_E_reco_shift2_Energy_Spectrum_per_Day(M)");
  shift3 = (TH1D*)f.Get("Mean_E_reco_shift3_Energy_Spectrum_per_Day(M)");
  shift4 = (TH1D*)f.Get("Mean_E_reco_shift4_Energy_Spectrum_per_Day(M)");
  /*
  for(Int_t i=0;i<mean1->GetNbinsX(); i++){
    mean1->SetBinError(i+1, mean1->GetBinError(i+1)* (1./TMath::Sqrt(nDay)) );
    shift1->SetBinError(i+1, shift1->GetBinError(i+1)* (1./TMath::Sqrt(nDay)) );
    shift2->SetBinError(i+1, shift2->GetBinError(i+1)* (1./TMath::Sqrt(nDay)) );  
  }
*/
  for(Int_t i=0;i<mean1->GetNbinsX(); i++){
    mean1->SetBinContent(i+1, mean1->GetBinContent(i+1)* nDay );
    shift1->SetBinContent(i+1, shift1->GetBinContent(i+1)* nDay );
    shift2->SetBinContent(i+1, shift2->GetBinContent(i+1)* nDay );
    shift3->SetBinContent(i+1, shift3->GetBinContent(i+1)* nDay );
    shift4->SetBinContent(i+1, shift4->GetBinContent(i+1)* nDay );
  }

  for(Int_t i=0;i<mean1->GetNbinsX(); i++){
    mean1->SetBinError(i+1, TMath::Sqrt(mean1->GetBinContent(i+1)) );
    shift1->SetBinError(i+1, TMath::Sqrt(shift1->GetBinContent(i+1)) );
    shift2->SetBinError(i+1, TMath::Sqrt(shift2->GetBinContent(i+1)) );
    shift3->SetBinError(i+1, TMath::Sqrt(shift3->GetBinContent(i+1)) );
    shift4->SetBinError(i+1, TMath::Sqrt(shift4->GetBinContent(i+1)) );
  }

  cout<<"point 2"<<endl;
  new TCanvas();
  mean1->SetTitle(Form("per %d day(s) spetrum comparison",nDay));
  mean1->SetLineColor(7);
  mean1->SetLineWidth(3);
  mean1->GetXaxis()->SetTitle("Energy (GeV)");
  mean1->GetYaxis()->SetTitle("Events per day per 0.25 GeV");
  mean1->Draw("hist e");
  cout<<"point 3"<<endl;

  shift1->SetLineColor(1);
  shift1->SetLineWidth(3);
  shift1->SetLineStyle(2);
  shift1->Draw("hist same");
  shift2->SetLineColor(2);
  shift2->SetLineWidth(3);
  shift2->SetLineStyle(2);
  shift2->Draw("hist same");
  shift3->SetLineColor(3);
  shift3->SetLineWidth(3);
  shift3->SetLineStyle(2);
  shift3->Draw("hist same");
  shift4->SetLineColor(4);
  shift4->SetLineWidth(3);
  shift4->SetLineStyle(2);
  shift4->Draw("hist same");

  cout<<"point 4"<<endl;
  TText *t1 = new TText(3.8, 700 , Form("smeared 4%% mu, 4%% pi+-, 10%% p, 10%% pi0, 20%% n") );
  t1->SetTextSize(0.035);
  t1->Draw();
  TText *t2 = new TText(3.8, 1000 , Form("Stat. Error with smearing") );
  t2->SetTextSize(0.035);
  t2->Draw();
/*  
  TLegend* legend = new TLegend(0.58, 0.62, 0.9, 0.9);
  legend->AddEntry(mean1,"nominal","l");
  legend->AddEntry(shift1,"target density","l");
  legend->AddEntry(shift2,"beam sigma","l");
  legend->AddEntry(shift3,"beam off set X","l");
  legend->AddEntry(shift4,"beam theta","l");
  legend->Draw();
*/
  TLegend* legend = new TLegend(0.58, 0.62, 0.9, 0.9);
  legend->AddEntry(shift1,"horn 1 X shift 0.5 mm","l");
  legend->AddEntry(shift2,"horn 1 Y shift 0.5 mm","l");
  legend->AddEntry(shift3,"horn 2 X shift 0.5 mm","l");
  legend->AddEntry(shift4,"horn 2 Y shift 0.5 mm","l");
  legend->Draw();

  TH1D* rshift1 = new TH1D("","",shift1->GetNbinsX(), 0, 10);
  TH1D* rshift2 = new TH1D("","",shift2->GetNbinsX(), 0, 10);
  TH1D* rshift3 = new TH1D("","",shift3->GetNbinsX(), 0, 10);
  TH1D* rshift4 = new TH1D("","",shift4->GetNbinsX(), 0, 10);

  for(Int_t i=0;i<mean1->GetNbinsX(); i++){
    rshift1->SetBinContent(i+1, (mean1->GetBinContent(i+1) - shift1->GetBinContent(i+1)) / mean1->GetBinError(i+1));
    rshift2->SetBinContent(i+1, (mean1->GetBinContent(i+1) - shift2->GetBinContent(i+1)) / mean1->GetBinError(i+1));
    rshift3->SetBinContent(i+1, (mean1->GetBinContent(i+1) - shift3->GetBinContent(i+1)) / mean1->GetBinError(i+1));
    rshift4->SetBinContent(i+1, (mean1->GetBinContent(i+1) - shift4->GetBinContent(i+1)) / mean1->GetBinError(i+1));
  }

  new TCanvas();
  rshift1->SetTitle("Stat. Error and detector effect (smearing + efficiency applied)");
  rshift1->SetLineColor(1);
  rshift1->SetLineWidth(3);
  rshift1->SetLineStyle(2);
  rshift1->Draw("hist");
  rshift2->SetLineColor(2);
  rshift2->SetLineWidth(3);
  rshift2->SetLineStyle(2);
  rshift2->Draw("hist same");
  rshift3->SetLineColor(3);
  rshift3->SetLineWidth(3);
  rshift3->SetLineStyle(2);
  rshift3->Draw("hist same");
  rshift4->SetLineColor(4);
  rshift4->SetLineWidth(3);
  rshift4->SetLineStyle(2);
  rshift4->Draw("hist same");
  rshift1->GetXaxis()->SetTitle("Energy (GeV)");
  rshift1->GetYaxis()->SetTitle("Pull (Significance)");
/*
  legend = new TLegend(0.58, 0.62, 0.9, 0.9);
  legend->AddEntry(rshift1,"target density","l");
  legend->AddEntry(rshift2,"beam sigma","l");
  legend->AddEntry(rshift3,"beam off set X","l");
  legend->AddEntry(rshift4,"beam theta","l");
  legend->Draw();
*/
/*
  legend = new TLegend(0.58, 0.62, 0.9, 0.9);
  legend->AddEntry(rshift1,"beam theta phi","l");
  legend->AddEntry(rshift2,"horn current","l");
  legend->AddEntry(rshift3,"water layer thickness","l");
  legend->AddEntry(rshift4,"decay pipe radius","l");
  legend->Draw();
*/
  legend = new TLegend(0.58, 0.62, 0.9, 0.9);
  legend->AddEntry(rshift1,"horn 1 X shift 0.5 mm","l");
  legend->AddEntry(rshift2,"horn 1 Y shift 0.5 mm","l");
  legend->AddEntry(rshift3,"horn 2 X shift 0.5 mm","l");
  legend->AddEntry(rshift4,"horn 2 Y shift 0.5 mm","l");
  legend->Draw();
  
  //int compare = 1;
  //if(compare == 1){
    TFile ff(Form("%s",gApplication->Argv(5)));
    TH1D* mean11 = (TH1D*)ff.Get("Mean_E_reco_Energy_Spectrum_per_Day");

    cout<<"point 3"<<endl;
    TH1D* shift11;
    TH1D* shift12;
    shift11 = (TH1D*)ff.Get("Mean_E_reco_shift1_Energy_Spectrum_per_Day(M)");
    shift12 = (TH1D*)ff.Get("Mean_E_reco_shift2_Energy_Spectrum_per_Day(M)");
/*
    for(Int_t i=0;i<mean11->GetNbinsX(); i++){
      mean11->SetBinError(i+1, mean11->GetBinError(i+1)* (1./TMath::Sqrt(nDay)) );
      shift11->SetBinError(i+1, shift11->GetBinError(i+1)* (1./TMath::Sqrt(nDay)) );
      shift12->SetBinError(i+1, shift12->GetBinError(i+1)* (1./TMath::Sqrt(nDay)) );
    }
*/  
    for(Int_t i=0;i<mean11->GetNbinsX(); i++){
      mean11->SetBinContent(i+1, mean11->GetBinContent(i+1)* nDay );
      shift11->SetBinContent(i+1, shift11->GetBinContent(i+1)* nDay );
      shift12->SetBinContent(i+1, shift12->GetBinContent(i+1)* nDay );
    }

    for(Int_t i=0;i<mean1->GetNbinsX(); i++){
      mean11->SetBinError(i+1, TMath::Sqrt(mean11->GetBinContent(i+1)) );
      shift11->SetBinError(i+1, TMath::Sqrt(shift11->GetBinContent(i+1)) );
      shift12->SetBinError(i+1, TMath::Sqrt(shift12->GetBinContent(i+1)) );
    }

    TH1D* rshift11 = new TH1D("","",shift11->GetNbinsX(), 0, 10);
    TH1D* rshift12 = new TH1D("","",shift12->GetNbinsX(), 0, 10);

    for(Int_t i=0;i<mean11->GetNbinsX(); i++){
      rshift11->SetBinContent(i+1, (mean11->GetBinContent(i+1) - shift11->GetBinContent(i+1)) / mean11->GetBinError(i+1));
      rshift12->SetBinContent(i+1, (mean11->GetBinContent(i+1) - shift12->GetBinContent(i+1)) / mean11->GetBinError(i+1));
    }

    new TCanvas();
    rshift11->SetTitle("Stat. Error and detector effect (smearing + efficiency applied)");
    rshift11->SetLineColor(2);
    rshift11->SetLineWidth(3);
    rshift11->SetLineStyle(1);
    rshift11->Draw("hist");
    rshift12->SetLineColor(3);
    rshift12->SetLineWidth(3);
    rshift12->SetLineStyle(1);
    rshift12->Draw("hist same");
    rshift1->Draw("hist same");
    rshift2->Draw("hist same");
    rshift11->GetXaxis()->SetTitle("Energy (GeV)");
    rshift11->GetYaxis()->SetTitle("Pull (Significance)");

    TF1* line1 = new TF1("","0",0,10);
    line1->SetLineWidth(3);
    line1->Draw("same");

    legend = new TLegend(0.58, 0.62, 0.9, 0.9);
    legend->AddEntry(rshift1,"Horn 1 X 3 mm neutrino reconstructed","l");
    legend->AddEntry(rshift2,"Horn 1 Y 0.5 mm neutrino reconstructed","l");
    legend->AddEntry(rshift11,"Horn 1 X 3 mm muon reconstructed","l");
    legend->AddEntry(rshift12,"Horn 1 Y 0.5 mm muon reconstructed","l");
    legend->Draw();

    //}
}
