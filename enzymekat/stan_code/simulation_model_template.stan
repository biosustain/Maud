data {
  // dimensions
  int<lower=1> N_balanced;    // 'Balanced' metabolites must have constant concentration at steady state
  int<lower=1> N_unbalanced;  // 'Unbalanced' metabolites can have changing concentration at steady state
  int<lower=1> N_kinetic_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_experiment;
  int<lower=1> N_known_real;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_conc_measurement;
  // measurements
  int<lower=1,upper=N_experiment> experiment_yconc[N_conc_measurement];
  int<lower=1,upper=N_balanced+N_unbalanced> metabolite_yconc[N_conc_measurement];
  vector<lower=0>[N_conc_measurement] sigma_conc;
  int<lower=1,upper=N_experiment> experiment_yflux[N_flux_measurement];
  int<lower=1,upper=N_reaction> reaction_yflux[N_flux_measurement];
  vector<lower=0>[N_flux_measurement] sigma_flux;
  // hardcoded numbers
  real xr[N_experiment, N_known_real];
  vector<lower=0>[N_balanced] balanced_guess;
  vector<lower=0>[N_kinetic_parameter] kinetic_parameter;
  vector<lower=0>[N_unbalanced] unbalanced[N_experiment];
  // ode configuration
  real rel_tol;
  real f_tol;
  int max_steps;
}
transformed data {
  int xi[0];
}
parameters {
  real z;
}
model {
  z ~ std_normal();
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_flux_measurement] yflux_sim;
  vector<lower=0>[N_balanced+N_unbalanced] conc[N_experiment];
  vector[N_reaction] flux[N_experiment];
  for (e in 1:N_experiment){
    vector[N_unbalanced+N_kinetic_parameter] theta = append_row(unbalanced[e], kinetic_parameter);
    conc[e, { {{-balanced_codes|join(',')-}} }] = algebra_solver(steady_state_system, balanced_guess, theta, xr[e], xi, rel_tol, f_tol, max_steps);
    conc[e, { {{-unbalanced_codes|join(',')-}} }] = unbalanced[e];
    flux[e] = get_fluxes(to_array_1d(conc[e]), to_array_1d(kinetic_parameter), xr[e]);
  }
  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], metabolite_yconc[c]]), sigma_conc[c]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}
