functions {
#include steady_state_equations.stan
}
data {
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1> R;
  vector[S] species;
  vector[P] kinetic_parameters;
  real known_reals[R];
}
generated quantities {
  int known_ints[0];
  vector[S] dsdt = steady_state_equation(species,
                                         kinetic_parameters,
                                         known_reals,
                                         known_ints);
}
