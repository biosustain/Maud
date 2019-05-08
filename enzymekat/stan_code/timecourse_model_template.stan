functions {
#include REPLACE_THIS_WORD
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
  real rel_tol;
  real abs_tol;
  int max_num_steps;
}
generated quantities {
  int known_ints[0];
  real ode_metabolite_sim[T+1,N_ode]; 
  real derived_quantity_sim[T+1, N_derived];
  ode_metabolite_sim[1] = initial_metabolite_ode;
  ode_metabolite_sim[2:T+1] = integrate_ode_bdf(steady_state_equation,
                                                initial_metabolite_ode,
                                                t0,
                                                ts,
                                                kinetic_parameters,
                                                known_reals,
                                                known_ints,
                                                rel_tol, abs_tol, max_num_steps);
  for (t in 1:T+1){
    derived_quantity_sim[t] = get_derived_quantities(ode_metabolite_sim[t], known_reals);
  }
}
