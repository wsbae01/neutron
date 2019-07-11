#include "fluxChange.hxx"
#include "TString.h"

int main(int argc, char const *argv[]){
  fluxChange obj2;
  //obj2. doFluxChange();
  TString st1 = "DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root";
  TString st2 = "ND_numode_OptimizedEngineeredNov2017Review_fit_binning_wppfx_0_45m.root";
  TString st3 = "FD_numode_OptimizedEngineeredNov2017Review_w_PPFX_fit_binning.root";
  double amount = 1.;
  int st = 20;
  int dm = 20;
  obj2. doFluxChange( st1, st2, st3, amount, st, dm);
}
