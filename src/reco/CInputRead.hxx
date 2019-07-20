#ifndef CInputRead_hxx_seen
#define CInputRead_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include "CHit2D.hxx"
#include <iostream>
#include <vector>

#include "TH2F.h"
#include "TPad.h"
#include "TH1.h"
#include "TVector3.h"

class CInputRead{

public:

bool SortHits3DTrue(const TVector3& lhs,  const TVector3& rhs);

void CInputRead(TTree *tree,int& nevts,int& skip);
    

#endif
