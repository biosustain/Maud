data {
  // dimensions
  int<lower=1> N_balanced;    // 'Balanced' metabolites must have constant concentration at steady state
  int<lower=1> N_unbalanced;  // 'Unbalanced' metabolites can have changing concentration at steady state
  int<lower=1> N_kinetic_parameters;
  int<lower=1> N_reaction;
  int<lower=1> N_enzyme;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=1> stoichiometric_rank;
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
  vector[N_enzyme] prior_loc_delta_g;
  vector<lower=0>[N_enzyme] prior_scale_delta_g;
  vector[N_kinetic_parameters] prior_loc_kinetic_parameters;
  vector<lower=0>[N_kinetic_parameters] prior_scale_kinetic_parameters;
  real prior_loc_unbalanced[N_unbalanced, N_experiment];
  real<lower=0> prior_scale_unbalanced[N_unbalanced, N_experiment];
  real prior_loc_enzyme[N_enzyme, N_experiment];
  real<lower=0> prior_scale_enzyme[N_enzyme, N_experiment];
  // network properties
  matrix[N_enzyme, stoichiometric_rank] delta_g_kernel;
  // configuration
  vector<lower=0>[N_balanced] as_guess;
  real rtol;
  real ftol;
  int steps;
  int<lower=0,upper=1> LIKELIHOOD;  // set to 0 for priors-only mode
}
transformed data {
  real xr[0];
  int xi[0];
  real minus_RT = - 0.008314 * 298.15;
}
parameters {
  vector[stoichiometric_rank] delta_g_basis_contribution;
  vector<lower=0>[N_kinetic_parameters] kinetic_parameters;
  vector<lower=0>[N_unbalanced] unbalanced[N_experiment];
  vector<lower=0>[N_enzyme] enzyme_concentration[N_experiment];
}
transformed parameters {
  vector<lower=0>[N_balanced+N_unbalanced] conc[N_experiment];
  vector[N_reaction] flux[N_experiment];
  vector[N_enzyme] delta_g = delta_g_kernel * delta_g_basis_contribution;  // linear transformation so no need for Jacobian adjustment
  for (e in 1:N_experiment){
    vector[N_enzyme] keq = exp(delta_g / minus_RT);
    vector[N_unbalanced+N_enzyme+N_enzyme+N_kinetic_parameters] theta;
    theta[1:N_unbalanced] = unbalanced[e];
    theta[N_unbalanced+1:N_unbalanced+N_enzyme] = enzyme_concentration[e];
    theta[N_unbalanced+N_enzyme+1:N_unbalanced+N_enzyme+N_enzyme] = keq;
    theta[N_unbalanced+N_enzyme+N_enzyme+1:] = kinetic_parameters;
    conc[e, { {{-balanced_codes|join(',')-}} }] = algebra_solver(steady_state_function, as_guess, theta, xr, xi, rtol, ftol, steps);
    conc[e, { {{-unbalanced_codes|join(',')-}} }] = unbalanced[e];
    flux[e] = get_fluxes(to_array_1d(conc[e]), to_array_1d(theta));
  }
}
model {
  kinetic_parameters ~ lognormal(log(prior_loc_kinetic_parameters), prior_scale_kinetic_parameters);
  delta_g ~ normal(prior_loc_delta_g, prior_scale_delta_g);
  for (e in 1:N_experiment){
    unbalanced[e] ~ lognormal(log(prior_loc_unbalanced[,e]), prior_scale_unbalanced[,e]);
    enzyme_concentration[e] ~ lognormal(log(prior_loc_enzyme[,e]), prior_scale_enzyme[,e]);
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
