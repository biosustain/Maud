functions {
#include steady_state_equation.stan
}
data {
  // dimensions
  int<lower=1> N_ode;         // number of ode metabolites
  int<lower=1> N_derived;     // number of derived metabolites
  int<lower=1> N_flux;        // number of fluxes
  int<lower=1> M_ode;         // number of measurements of ode metabolites
  int<lower=1> M_derived;     // number of measurements of derived metabolites
  int<lower=1> M_flux;        // number of flux measurements
  int<lower=1> N_known_real;  // number of known reals
  int<lower=1> P;             // total number of parameters
  // measurements
  int<lower=1,upper=N_flux> measurement_ix_flux[M_flux];
  vector[M_flux] measurement_flux;
  int<lower=1,upper=N_ode> measurement_ix_ode[M_ode];
  vector[M_ode] measurement_ode;
  int<lower=1,upper=N_derived> measurement_ix_derived[M_derived];
  vector[M_derived] measurement_derived;
  // hardcoded priors
  vector[P] prior_location;
  vector[P] prior_scale;
  real<lower=0> sigma_metabolite;
  real<lower=0> sigma_flux;
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
  vector<lower=0>[P] kinetic_parameters;
}
transformed parameters {
  vector[N_ode] ode_hat = algebra_solver(steady_state_equation,
                                         initial_guess,
                                         kinetic_parameters,
                                         known_reals,
                                         x_i,
                                         rel_tol, f_tol, max_steps);
  real derived_hat[N_derived] = get_derived_metabolites(ode_hat, known_reals); 
  vector[N_flux] flux_hat = get_fluxes(ode_hat, kinetic_parameters, known_reals);
}
model {
  kinetic_parameters ~ lognormal(prior_location, prior_scale);
  if (LIKELIHOOD == 1){
    measurement_ode[measurement_ix_ode] ~ normal(ode_hat[measurement_ix_ode], sigma_metabolite);
    measurement_derived[measurement_ix_derived] ~ normal(derived_hat[measurement_ix_derived], sigma_metabolite);
    measurement_flux[measurement_ix_flux] ~ normal(flux_hat[measurement_ix_flux], sigma_flux);
  }
}
generated quantities {
  vector[N_ode] ode_pred;
  vector[N_derived] derived_pred;
  vector[N_flux] flux_pred;
  for (n in 1:N_ode)
    ode_pred[n] = normal_rng(ode_hat[n], sigma_metabolite);
  for (n in 1:N_derived)
    derived_pred[n] = normal_rng(derived_hat[n], sigma_metabolite);
  for (n in 1:N_flux)
    flux_pred[n] = normal_rng(flux_hat[n], sigma_flux);
}
