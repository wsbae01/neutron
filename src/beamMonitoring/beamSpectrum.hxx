#ifndef BEAMSPECTRUM_SEEN
#define BEAMSPECTRUM_SEEN

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class beamSpectrum {

  public:

  double para1;

  void doBeamSpectrum(double muonSmear = 0.04, double pipmSmear = 0.04, double neutronSmear = 0.2, double protonSmear = 0.1, double pi0Smear = 0.1, int onlyMuon=0, int incNeutron=1, int fileLimit = 40, TString filepath = "/dune/app/users/gyang/genie-0.0-0.101/", int caseN = 3) ;

};

#endif
