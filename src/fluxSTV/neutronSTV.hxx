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
#include "TChain.h"
#include "TVector3.h"

using namespace std;

class neutronSTV {

  public:

  void doNeutronSTV(
    double totRate = 5e6,
    double STVcut1 = 0.10000,
    double STVcut2 = 0.05000,
    double STVcut3 = 0.02000,
    double muonThreshold = 0.1,
    double pionThreshold = 0.1,
    double neutronThreshold = 0.01,
    double protonThreshold = 0.3,
    double electronThreshold = 0,
    double muonSmear = 0.04,
    double pipSmear = 0.04,
    double pi0Smear = 0.1,
    double electronSmear = 0.1,
    double protonSmear = 0.1,
    double neutronSmear = 0.2 );

};
