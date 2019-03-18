functions {
#include steady_state_equation.stan
  real[] ode(real t,        // time
             real[] s,      // state
             real[] theta,  // parameters
             real[] x_r,    // data (real)
             int[] x_i){   // data (integer)
    return to_array_1d(steady_state_equation(to_vector(s), to_vector(theta), x_r, x_i));
  }
}
data {
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1> R;
  int<lower=1> T;
  real species[S];
  real kinetic_parameters[P];
  real known_reals[R];
  real ts[T];
  real t0;
}
generated quantities {
  int known_ints[0];
  real species_sim[T+1,S]; 
  real derived_quantities_sim[T+1, 12];
  species_sim[1] = species;
  species_sim[2:T+1] = integrate_ode_rk45(ode, species, t0, ts, kinetic_parameters, known_reals, known_ints);
  for (t in 1:T+1){
    derived_quantities_sim[t] = get_derived_quantities(to_vector(species_sim[t]), known_reals);
  }
}
