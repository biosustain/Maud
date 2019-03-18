functions {
#include steady_state_equation.stan
}
data {
  // dimensions
  int<lower=1> S;             // number of species
  int<lower=1> D;             // number of derived quantities
  int<lower=1> F;             // number of fluxes
  int<lower=1> SM;            // number of species measurements
  int<lower=1> DM;            // number of derived quantity measurements
  int<lower=1> FM;            // number of flux measurements
  int<lower=1> KR;            // number of known reals
  int<lower=1> P;             // total number of parameters
  int<lower=1> Q;             // total number of known quantities
  // measurements
  int[FM] flux_measurement_ix;
  vector[FM] flux_measurement;
  int[SM] species_measurement_ix;
  vector[SM] species_measurement;
  int[D] derived_quantity_measurement_ix;
  vector[DM] derived_quantity_measurement;
  // hardcoded priors
  vector[P] prior_location;
  vector[P] prior_scale;
  vector<lower=0>[SM] sigma_measurement;
  real known_reals[KR];
  // algebra solver config
  vector[S] initial_guess;
  real rel_tol;
  real f_tol;
  int max_steps;
  // likelihood config
  int<lower=0,upper=1> LIKELIHOOD;
}
parameters {
  real kinetic_parameters[P];
}
transformed parameters {
  int x_i[0];
  vector[S] species_hat = algebra_solver(steady_state_equation,
                                         initial_guess,
                                         kinetic_parameters,
                                         known_reals,
                                         x_i,
                                         rel_tol, f_tol, max_steps);
  vector[D] derived_quantity_hat = get_derived_quantities(species_hat, known_reals); 
}
model {
  kinetic_parameters ~ lognormal(prior_location, prior_scale);
  if (LIKELIHOOD == 1){
    species_measurement[species_measurement_ix]
      ~ normal(species_hat[species_measurement_ix], sigma_measurement);
    derived_quantity_measurement[derived_quantity_measurement_ix]
      ~ normal(derived_quantity_hat[derived_quantity_measurement_ix], sigma_measurement);
    flux_measurement[flux_measurement_ix] ~ normal(flux_hat[flux_measurement_ix], sigma_flux);
  }
}
generated quantities {
  vector[S] species_pred;
  vector[D] derived_quantity_pred;
  for (s in 1:S){
    species_pred[s] = normal_rng(species_hat[s], sigma_measurement);
  }
  for (dq in 1:D){
    derived_quantity_pred[dq] = normal_rng(derived_quantity_hat[dq], sigma_measurement);
  }
}
