#ifndef FLUXCHANGE_SEEN
#define FLUXCHANGE_SEEN

#include "TString.h"

class fluxChange {

  public:

  void doFluxChange (TString inputShift="DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root", TString fileNDflux="ND_numode_OptimizedEngineeredNov2017Review_fit_binning_wppfx_0_45m.root", TString fileFDflux="FD_numode_OptimizedEngineeredNov2017Review_w_PPFX_fit_binning.root", double cAmount=1., int st=20 , int dm=20 );

};

#endif
