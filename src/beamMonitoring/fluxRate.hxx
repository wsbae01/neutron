#ifndef FLUXRATE_SEEN
#define FLUXRATE_SEEN

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <memory>
#include <sstream>
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

class fluxRate {

  public:

  double para1;

  void doFluxRate(TString inputShift="DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root");

};

#endif
