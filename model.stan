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
  // measurements
  int flux_measurement_ix[FM];
  vector[FM] flux_measurement;
  int species_measurement_ix[SM];
  vector[SM] species_measurement;
  int derived_quantity_measurement_ix[DM];
  vector[DM] derived_quantity_measurement;
  // hardcoded priors
  vector[P] prior_location;
  vector[P] prior_scale;
  real<lower=0> sigma_measurement;
  real<lower=0> sigma_flux;
  real known_reals[KR];
  // algebra solver config
  vector[S] initial_guess;
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
  vector[S] species_hat = algebra_solver(steady_state_equation,
                                         initial_guess,
                                         kinetic_parameters,
                                         known_reals,
                                         x_i,
                                         rel_tol, f_tol, max_steps);
  real derived_quantity_hat[D] = get_derived_quantities(species_hat, known_reals); 
  vector[F] flux_hat = get_fluxes(species_hat, kinetic_parameters, known_reals);
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
  for (d in 1:D){
    derived_quantity_pred[d] = normal_rng(derived_quantity_hat[d], sigma_measurement);
  }
  for (f in 1:F){
    derived_quantity_pred[f] = normal_rng(derived_quantity_hat[f], sigma_flux);
  }
}
