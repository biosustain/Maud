functions {
#include REPLACE_THIS_WORD
}
data {
  // dimensions
  int<lower=1> N_ode;          // number of ode metabolites
  int<lower=1> N_kinetic_parameter;        // total number of kinetic parameters
  int<lower=1> N_known_real;   // number of known reals
  int<lower=1> N_measurement;  // number of measurements of ode metabolites
  // measurements
  int<lower=1,upper=N_ode> measurement_ix[N_measurement];
  vector[N_measurement] measurement;
  // hardcoded priors
  vector[N_kinetic_parameter] prior_location;
  vector[N_kinetic_parameter] prior_scale;
  real<lower=0> sigma_measurement;
  real known_reals[N_known_real];
  // algebra solver config
  vector[N_ode] initial_guess;
  real rel_tol;
  real f_tol;
  int max_steps;
  // likelihood config
  int<lower=0,upper=1> LIKELIHOOD;
}
transformed data {
  int x_i[0];
}
parameters {
  vector<lower=0>[N_kinetic_parameter] kinetic_parameters;
}
transformed parameters {
  vector[N_ode] measurement_hat = algebra_solver(steady_state_equation,
                                                 initial_guess,
                                                 kinetic_parameters,
                                                 known_reals,
                                                 x_i,
                                                 rel_tol, f_tol, max_steps);
}
model {
  kinetic_parameters ~ lognormal(prior_location, prior_scale);
  if (LIKELIHOOD == 1){
    measurement ~ normal(measurement_hat[measurement_ix], sigma_measurement);
  }
}
generated quantities {
  vector[N_ode] measurement_pred;
  for (n in 1:N_ode)
    measurement_pred[n] = normal_rng(measurement_hat[n], sigma_measurement);
}
