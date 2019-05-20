functions {
#include ../enzymekat/stan_code/big_k_rate_equations.stan
#include ../enzymekat/stan_code/haldane_relationships.stan
#include ../enzymekat/stan_code/ode_equations_yeast.stan
}
data {
  // dimensions
  int<lower=1> N_ode;          // number of ode metabolites
  int<lower=1> N_kinetic_parameter;        // total number of kinetic parameters
  int<lower=1> N_known_real;   // number of known reals
  int<lower=1> N_measurement;  // number of measurements of ode metabolites
  int<lower=1> N_reaction;
  // measurements
  int<lower=1,upper=N_ode> measurement_ix[N_measurement];
  vector[N_measurement] measurement;
  // hardcoded priors
  vector[N_kinetic_parameter] prior_location_kinetic;
  vector[N_kinetic_parameter] prior_scale_kinetic;
  vector[N_reaction] prior_location_thermodynamic;
  vector[N_reaction] prior_scale_thermodynamic;
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
  real thermodynamic_parameters[N_reaction];
  real<lower=0> kinetic_parameters[N_kinetic_parameter];
}
transformed parameters {
  real metabolite_concentration_hat[N_ode] = integrate_ode_bdf(
    steady_state_equation,
    initial_state,
    initial_time,
    {steady_time},
    append_array(thermodynamic_parameters, kinetic_parameters),
    known_reals,
    known_ints,
    rel_tol, abs_tol, max_steps
  )[1];
}
model {
  kinetic_parameters ~ lognormal(prior_location_kinetic, prior_scale_kinetic);
  thermodynamic_parameters ~ normal(prior_location_thermodynamic, prior_scale_thermodynamic);
  if (LIKELIHOOD == 1){
    measurement ~ lognormal(log(metabolite_concentration_hat[measurement_ix]), measurement_scale);
  }
}
generated quantities {
  vector[N_ode] measurement_pred;
  for (n in 1:N_ode)
    measurement_pred[n] = lognormal_rng(log(metabolite_concentration_hat[n]), measurement_scale);
}
