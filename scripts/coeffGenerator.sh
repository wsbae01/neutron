for i in 0 10 20 30 40; do
  for j in 0 10 20 30 40; do
    /home/guang/work/DUNEPrismTools/build/app/flux_tools/dp_FluxLinearSolver_Standalone -N /home/guang/work/DUNEPrismTools/ND_numode_OptimizedEngineeredNov2017Review_fit_binning_wppfx_0_45m.root,LBNF_numu_flux_Nom -F /home/guang/work/DUNEPrismTools/oscillated.root,var_st${i}_dm${j} -o output_st${i}_dm${j}.root -M 4 
  done
done

