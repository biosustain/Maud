functions {
//this path is relative to the cmdstan directory
#include ../stan/autogen/t_brucei.stan  
}
data {
  int N_ode;
  int N_kinetic_parameter;
  int N_known_real;
  vector[N_ode] initial_guess;
  vector[N_kinetic_parameter] log_kinetic_parameters;
  real known_reals[N_known_real];
  real rel_tol;
  real f_tol;
  int max_steps;
}
transformed data {
  int x_i[0];
}
generated quantities {
  vector<lower=0>[N_ode] measurement_hat = algebra_solver(steady_state_equation,
                                                          initial_guess,
                                                          exp(log_kinetic_parameters),
                                                          known_reals,
                                                          x_i,
                                                          rel_tol, f_tol, max_steps);
}
