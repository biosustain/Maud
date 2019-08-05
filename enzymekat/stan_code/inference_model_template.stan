data {
  // dimensions
  int<lower=1> N_metabolite;
  int<lower=1> N_param;
  int<lower=1> N_reaction;
  int<lower=1> N_experiment;
  int<lower=1> N_known_real;
  int<lower=0> N_measurement_flux;
  int<lower=0> N_measurement_conc;
  int<lower=0> N_unbalanced_metabolites;
  // measurements
  int<lower=0,upper=N_metabolite> metabolite_ix[N_measurement_conc];
  int<lower=1,upper=N_experiment> experiment_ix_conc[N_measurement_conc];
  vector[N_measurement_conc] measurement_conc;
  vector<lower=0>[N_measurement_conc] measurement_scale_conc;
  int<lower=0,upper=N_reaction> reaction_ix[N_measurement_flux];
  int<lower=1,upper=N_experiment> experiment_ix_flux[N_measurement_flux];
  vector[N_measurement_flux] measurement_flux;
  vector<lower=0>[N_measurement_flux] measurement_scale_flux;
  int unbalanced_met_arr[N_unbalanced_metabolites, 2];
  // hardcoded
  real known_reals[N_known_real, N_experiment];
  vector[N_param] prior_location;
  vector[N_param] prior_scale;
  // ode stuff
  real initial_concentration[N_metabolite, N_experiment];
  real initial_time;
  real steady_time;
  real rel_tol;
  real abs_tol;
  int max_steps;
  // likelihood config
  int<lower=0,upper=1> LIKELIHOOD;
}
transformed data {
  int known_ints[0];
}
parameters {
  real<lower=0> params[N_param];
  real<lower=0> unbalanced_metabolite[N_unbalanced_metabolites];
}
transformed parameters {
  real metabolite_concentration[N_metabolite, N_experiment];
  real flux[N_reaction, N_experiment];
  real new_init_concentration [N_metabolite, N_experiment];
  new_init_concentration = initial_concentration;
  for (um in 1:N_unbalanced_metabolites){
    new_init_concentration[unbalanced_met_arr[um][1], unbalanced_met_arr[um][2]] = unbalanced_metabolite[um];
  }
  for (e in 1:N_experiment){
    metabolite_concentration[,e] = integrate_ode_bdf(steady_state_equation,
                                                     new_init_concentration[,e],
                                                     initial_time,
                                                     {steady_time},
                                                     params,
                                                     known_reals[,e],
                                                     known_ints,
                                                     rel_tol, abs_tol, max_steps)[1];
    flux[,e] = get_fluxes(metabolite_concentration[,e], params, known_reals[,e]);
  }
}
model {
  params ~ lognormal(prior_location, prior_scale);
  if (LIKELIHOOD == 1){
    vector[N_measurement_conc] conc_hat;
    vector[N_measurement_flux] flux_hat;
    for (mc in 1:N_measurement_conc){
      conc_hat[mc] = metabolite_concentration[metabolite_ix[mc], experiment_ix_conc[mc]];
    }
    for (mf in 1:N_measurement_flux){
      flux_hat[mf] = flux[reaction_ix[mf], experiment_ix_flux[mf]];
    }
    measurement_conc ~ lognormal(log(conc_hat), measurement_scale_conc);
    measurement_flux ~ normal(flux_hat, measurement_scale_flux);
  }
}
generated quantities {
  vector[N_measurement_flux] flux_pred;
  vector[N_measurement_conc] conc_pred;
  real metabolite_flux[N_metabolite, N_experiment];
  for (e in 1:N_experiment){
    metabolite_flux[, e] = get_odes(flux[, e]);
  }
  for (mc in 1:N_measurement_conc){
    conc_pred[mc] = lognormal_rng(log(metabolite_concentration[metabolite_ix[mc], experiment_ix_conc[mc]]),
                                  measurement_scale_conc[mc]);
  }
         for (mf in 1:N_measurement_flux){
           flux_pred[mf] = normal_rng(flux[reaction_ix[mf], experiment_ix_flux[mf]], measurement_scale_flux[mf]);
         }
}
