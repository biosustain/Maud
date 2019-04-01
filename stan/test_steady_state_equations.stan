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
  int<lower=1> N_ode;
  int<lower=1> N_derived;
  int<lower=1> N_known_real;
  int<lower=1> P;
  int<lower=1> T;
  real initial_metabolite_ode[N_ode];
  real kinetic_parameters[P];
  real known_reals[N_known_real];
  real ts[T];
  real t0;
}
generated quantities {
  int known_ints[0];
  real metabolite_sim_ode[T+1,N_ode]; 
  real metabolite_sim_derived[T+1, N_derived];
  metabolite_sim_ode[1] = initial_metabolite_ode;
  metabolite_sim_ode[2:T+1] = integrate_ode_rk45(ode,
                                                 initial_metabolite_ode,
                                                 t0,
                                                 ts,
                                                 kinetic_parameters,
                                                 known_reals,
                                                 known_ints);
  for (t in 1:T+1){
    metabolite_sim_derived[t] = get_derived_metabolites(to_vector(metabolite_sim_ode[t]), known_reals);
  }
}
