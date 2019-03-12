functions {
#include steady_state_equations.stan
}
data {
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1> R;
  int<lower=1> I;
  vector[S] species;
  vector[P] kinetic_parameters;
  vector[R] known_reals;
  vector[I] known_ints;
}
generated quantities {
  vector[S] dsdt = steady_state_equation(species,
                                         kinetic_parameters,
                                         known_reals,
                                         known_ints);
}
