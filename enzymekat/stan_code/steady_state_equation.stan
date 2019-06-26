real[] steady_state_equation(real t,
                             real[] metabolites,
                             real[] params,
                             real[] known_reals,
                             int[] known_ints){
  for (m in 1:size(metabolites)){
    if (metabolites[m] < 0){
      reject("Metabolite ", m, " is ", metabolites[m], " but should be greater than zero.");
    }
  }
  return get_odes(get_fluxes(metabolites, params, known_reals));
}
