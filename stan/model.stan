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
  real<lower=0> measurement_scale;
  real known_reals[N_known_real];
  // ode stuff
  real initial_state[N_ode];
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
  real<lower=0> kinetic_parameters[N_kinetic_parameter];
}
transformed parameters {
  real measurement_hat[N_ode] = integrate_ode_bdf(steady_state_equation,
                                                  initial_state,
                                                  initial_time,
                                                  {steady_time},
                                                  kinetic_parameters,
                                                  known_reals,
                                                  known_ints,
                                                  rel_tol, abs_tol, max_steps)[1] ;
}
model {
  kinetic_parameters ~ lognormal(prior_location, prior_scale);
  if (LIKELIHOOD == 1){
    measurement ~ lognormal(log(measurement_hat[measurement_ix]), measurement_scale);
  }
}
generated quantities {
  vector[N_ode] measurement_pred;
  for (n in 1:N_ode)
    measurement_pred[n] = lognormal_rng(log(measurement_hat[n]), measurement_scale);
}
