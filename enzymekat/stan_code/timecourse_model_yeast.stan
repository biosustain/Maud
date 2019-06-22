functions {
#include big_k_rate_equations.stan
#include haldane_relationships.stan
#include ode_equations_yeast.stan
}
data {
  int<lower=1> N_ode;
  int<lower=1> N_known_real;
  int<lower=1> N_kinetic_parameter;
  int<lower=1> N_thermodynamic_parameter;
  int<lower=1> N_timepoint;
  int<lower=1> N_reaction;
  real initial_state[N_ode];
  real kinetic_parameters[N_kinetic_parameter];
  real thermodynamic_parameters[N_thermodynamic_parameter];
  real known_reals[N_known_real];
  real timepoints[N_timepoint-1];
  real initial_time;
  real rel_tol;
  real abs_tol;
  int max_steps;
}
parameters{
  real x;
}
model {
  x ~ normal(0, 1);
}
generated quantities {
  int known_ints[0];
  real flux_sim[N_timepoint, N_reaction+2];
  real ode_metabolite_sim[N_timepoint,N_ode]; 
  ode_metabolite_sim[1] = initial_state;
  ode_metabolite_sim[2:N_timepoint] = integrate_ode_bdf(steady_state_equation,
                                                        initial_state,
                                                        initial_time,
                                                        timepoints,
                                                        append_array(thermodynamic_parameters, kinetic_parameters),
                                                        known_reals,
                                                        known_ints,
                                                        rel_tol, abs_tol, max_steps);
  for (t in 1:N_timepoint){
    flux_sim[t] = get_fluxes(ode_metabolite_sim[t],
                             append_array(thermodynamic_parameters, kinetic_parameters),
                             known_reals);
  }
}
