functions {
#include steady_state_equations.stan
}
data {
  // dimensions
  int<lower=1> S;                    // number of species
  int<lower=1> R;                    // number of reactions
  int<lower=1> P;                    // total number of parameters
  int<lower=1> Q;                    // total number of known quantities
  int<lower=0> n_known_quantity[R];  // number of known quantities per reaction
  // measurements
  vector[S] species_measured;
  // hardcoded priors
  vector[P] prior_location;
  vector[P] prior_scale;
  real<lower=0> sigma_measurement;
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
  vector[S] species_hat = algebra_solver(steady_state_equations,
                                         initial_guess,
                                         kinetic_parameters,
                                         known_quantities,
                                         x_i,
                                         rel_tol, f_tol, max_steps);
}
model {
  kinetic_parameters ~ lognormal(prior_location, prior_scale);
  if (LIKELIHOOD == 1){
    species_measured ~ normal(species_hat, sigma_measurement);
  }
}
generated quantities {
  vector[S] species_pred;
  for (s in 1:S){
    species_pred[s] = normal_rng(species_hat[s], sigma_measurement);
  }
}
