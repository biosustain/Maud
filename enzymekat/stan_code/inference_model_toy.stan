functions {
#include big_k_rate_equations.stan
#include haldane_relationships.stan
#include ode_equations_toy.stan
}
data {
  // dimensions
  int<lower=1> N_metabolite;          // number of ode metabolites
  int<lower=1> N_kinetic_parameter;        // total number of kinetic parameters
  int<lower=1> N_thermodynamic_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_known_real;   // number of known reals
  int<lower=1> N_measurement; // number of measured fluxes
  // measurements
  int<lower=1,upper=N_metabolite> metabolite_ix[N_measurement];
  int<lower=1, upper=N_reaction> reaction_ix[N_measurement];
  int<lower=1,upper=N_experiment> experiment_ix[N_measurement];
  int<lower=0,upper=1> is_flux[N_measurement];
  vector[N_measurement] measurement;
  vector<lower=0>[N_measurement] measurement_scale;
  // hardcoded
  real known_reals[N_known_real, N_experiment];
  vector[N_kinetic_parameter] prior_location_kinetic;
  vector[N_kinetic_parameter] prior_scale_kinetic;
  vector[N_thermodynamic_parameter] prior_location_thermodynamic;
  vector[N_thermodynamic_parameter] prior_scale_thermodynamic;
  vector<lower=0>[N_metabolite_measurement] measurement_scale;
  vector<lower=0>[N_flux_measurement] flux_measurment_scale;
  // ode stuff
  real initial_concentrations[N_metabolite];
  real initial_time;
  real steady_time;
  real rel_tol;
  real abs_tol;
  int max_steps;
  // likelihood config
  int<lower=0,upper=1> LIKELIHOOD;
}
transformed data {
  int known_ints[0];
}
parameters {
  real thermodynamic_parameters[N_thermodynamic_parameter];
  real kinetic_parameters[N_kinetic_parameter];
}
transformed parameters {
  real metabolite_concentration[N_ode, N_experiment];
  real flux[N_reaction, N_experiment];
  for (e in 1:N_experiment){
    metabolite_concentration[,e] =
      integrate_ode_bdf(steady_state_equation,
                        initial_concentrations,
                        initial_time,
                        {steady_time},
                        append_array(thermodynamic_parameters, kinetic_parameters),
                        known_reals[,e],
                        known_ints,
                        rel_tol, abs_tol, max_steps)[1];
    flux[,e] =
      get_fluxes(metabolite_concentration[,e],
                 append_array(thermodynamic_parameters, kinetic_parameters),
                 known_reals);
  }
}
model {
  kinetic_parameters ~ lognormal(prior_location_kinetic, prior_scale_kinetic);
  thermodynamic_parameters ~ normal(prior_location_thermodynamic, prior_scale_thermodynamic);
  if (LIKELIHOOD == 1){
    for (m in 1:N_measurement){
      target += is_flux[m] ?
        normal_lpdf(measurement[m] | flux[reaction_ix[m], experiment_ix[m]], measurment_scale[m]):
        lognormal_lpdf(measurement[m] | log(metabolite_concentration[measurement_ix]), measurement_scale[m]);
    }
  }
}
generated quantities {
  vector[N_measurement] measurement_pred;
  for (m in 1:N_measurement){
    measurement_pred[m] = is_flux[m] ?
      normal_lpdf(measurement[m] | flux[reaction_ix[m], experiment_ix[m]], measurment_scale[m]):
      lognormal_lpdf(measurement[m] | log(metabolite_concentration[measurement_ix]), measurement_scale[m]);
  }
}
