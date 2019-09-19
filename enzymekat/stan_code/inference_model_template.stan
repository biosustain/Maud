data {
  // dimensions
  int<lower=1> N_balanced;    // 'Balanced' metabolites must have constant concentration at steady state
  int<lower=1> N_unbalanced;  // 'Unbalanced' metabolites can have changing concentration at steady state
  int<lower=1> N_kinetic_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_enzyme;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_conc_measurement;
  // measurements
  int<lower=1,upper=N_experiment> experiment_yconc[N_conc_measurement];
  int<lower=1,upper=N_balanced+N_unbalanced> metabolite_yconc[N_conc_measurement];
  vector[N_conc_measurement] yconc;
  vector<lower=0>[N_conc_measurement] sigma_conc;
  int<lower=1,upper=N_experiment> experiment_yflux[N_flux_measurement];
  int<lower=1,upper=N_reaction> reaction_yflux[N_flux_measurement];
  vector[N_flux_measurement] yflux;
  vector<lower=0>[N_flux_measurement] sigma_flux;
  // hardcoded priors
  vector[N_kinetic_parameter] prior_loc_kinetic_parameter;
  vector<lower=0>[N_kinetic_parameter] prior_scale_kinetic_parameter;
  real prior_loc_unbalanced[N_unbalanced, N_experiment];
  real<lower=0> prior_scale_unbalanced[N_unbalanced, N_experiment];
  real prior_loc_enzyme[N_enzyme, N_experiment];
  real<lower=0> prior_scale_enzyme[N_enzyme, N_experiment];
  vector<lower=0>[N_balanced] balanced_guess;
  // algebra solver configuration
  real rel_tol;
  real f_tol;
  int max_steps;
  // likelihood configuration - set to 0 for priors-only mode
  int<lower=0,upper=1> LIKELIHOOD;
}
transformed data {
  real xr[0];
  int xi[0];
}
parameters {
  vector<lower=0>[N_kinetic_parameter] kinetic_parameter;
  vector<lower=0>[N_unbalanced] unbalanced[N_experiment];
  vector<lower=0>[N_enzyme] enzyme_concentration[N_experiment];
}
transformed parameters {
  vector<lower=0>[N_balanced+N_unbalanced] conc[N_experiment];
  vector[N_reaction] flux[N_experiment];
  for (e in 1:N_experiment){
    vector[N_unbalanced+N_enzyme+N_kinetic_parameter] theta = append_row(unbalanced[e], append_row(enzyme_concentration[e], kinetic_parameter));
    conc[e, { {{-balanced_codes|join(',')-}} }] = algebra_solver(steady_state_function, balanced_guess, theta, xr, xi, rel_tol, f_tol, max_steps);
    conc[e, { {{-unbalanced_codes|join(',')-}} }] = unbalanced[e];
    flux[e] = get_fluxes(to_array_1d(conc[e]), append_array(to_array_1d(enzyme_concentration[e]), to_array_1d(kinetic_parameter)));
  }
}
model {
  kinetic_parameter ~ lognormal(prior_loc_kinetic_parameter, prior_scale_kinetic_parameter);
  for (e in 1:N_experiment){
    unbalanced[e] ~ lognormal(prior_loc_unbalanced[,e], prior_scale_unbalanced[,e]);
    enzyme_concentration[e] ~ lognormal(prior_loc_enzyme[,e], prior_scale_enzyme[,e]);
  }
  if (LIKELIHOOD == 1){
    for (c in 1:N_conc_measurement){
      target += lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], metabolite_yconc[c]]), sigma_conc[c]);
    }
    for (f in 1:N_flux_measurement){
      target += normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
    }
  }
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_flux_measurement] yflux_sim;
  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], metabolite_yconc[c]]), sigma_conc[c]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}
