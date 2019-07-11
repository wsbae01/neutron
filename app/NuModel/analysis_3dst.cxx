#include "analysis_3dst.hxx"

int main(){
  //beauty();
  MCMC mcmc;
  mcmc.setIterationTime(1000000);
  mcmc.process_Ar();
  mcmc.MCMC_master();
}
