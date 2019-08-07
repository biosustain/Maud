parameters {
  real<lower=0> params[N_param];
}
transformed parameters {
  real metabolite_concentration[N_metabolite, N_experiment];
  real flux[N_reaction, N_experiment];
  for (e in 1:N_experiment){
    metabolite_concentration[,e] = integrate_ode_bdf(steady_state_equation,
                                                     initial_concentration,
                                                     initial_time,
                                                     {steady_time},
                                                     params,
                                                     known_reals[,e],
                                                     known_ints,
                                                     rel_tol, abs_tol, max_steps)[1];
    flux[,e] = get_fluxes(metabolite_concentration[,e], params, known_reals[,e]);
  }
}
