functions {
#include ../enzymekat/stan_code/big_k_rate_equations.stan
#include ../enzymekat/stan_code/haldane_relationships.stan
#include ../enzymekat/stan_code/ode_equations_relative_metabolomics.stan
}
data {
  // dimensions
  int<lower=1> N_ode;          // number of ode metabolites
  int<lower=1> N_kinetic_parameter;        // total number of kinetic parameters
  int<lower=1> N_thermodynamic_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_known_real;   // number of known reals
  int<lower=1> N_measurement;  // number of measurements of ode metabolites
  int<lower=1> N_flux_measurement; // number of measured fluxes
  // measurements
  int<lower=1,upper=N_ode> measurement_ix[N_measurement];
  vector[N_measurement] measurement;
  int<lower=1, upper=N_reaction> flux_measurement_ix[N_flux_measurement];
  vector[N_flux_measurement] flux_measurment;
  // hardcoded priors
  vector[N_kinetic_parameter] prior_location_kinetic;
  vector[N_kinetic_parameter] prior_scale_kinetic;
  vector[N_reaction] prior_enzyme_concentrations;
  vector[N_reaction] prior_enzyme_scale;
  vector[N_thermodynamic_parameter] prior_location_thermodynamic;
  vector[N_thermodynamic_parameter] prior_scale_thermodynamic;
  vector<lower=0>[N_measurement] measurement_scale;
  vector<lower=0>[N_flux_measurement] flux_measurment_scale;
  real known_reals[N_known_real];
  // ode stuff
  real initial_state[N_ode];
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
  vector[N_kinetic_parameter] log_kinetic_parameters_z;
  real<lower=0> enzyme_concentrations[N_reaction];
  real<lower=0> metabolite_scaling[N_ode];
}
transformed parameters {
  real kinetic_parameters[N_kinetic_parameter] =
    to_array_1d(exp(prior_location_kinetic + log_kinetic_parameters_z .* prior_scale_kinetic));
  real ode_pars[N_kinetic_parameter+N_thermodynamic_parameter+N_reaction] = append_array(append_array(thermodynamic_parameters, kinetic_parameters), enzyme_concentrations);
  real metabolite_concentration_hat[N_ode] = integrate_ode_bdf(
    steady_state_equation,
    initial_state,
    initial_time,
    {steady_time},
    ode_pars,
    known_reals,
    known_ints,
    rel_tol, abs_tol, max_steps
  )[1];
  real flux_hat[N_reaction] = get_fluxes(metabolite_concentration_hat, ode_pars,
known_reals);
  real stability[N_ode] = fabs(get_odes(flux_hat));
  real abs_flux = sum(stability);
  real relative_metabolite_hat[N_ode];

  for (i in 1:N_ode){
    relative_metabolite_hat[i] = metabolite_concentration_hat[i] / metabolite_scaling[i];
  }

  if (abs_flux > 0.001){
    reject("steady state not achieved, sum of absolute fluxes was ", abs_flux);
  }
}
model {
  log_kinetic_parameters_z ~ normal(0, 1);
  enzyme_concentrations ~ normal(prior_enzyme_concentrations, prior_enzyme_scale);
  thermodynamic_parameters ~ normal(prior_location_thermodynamic, prior_scale_thermodynamic);
  if (LIKELIHOOD == 1){
    metabolite_scaling ~ lognormal(1, 1);
    flux_measurment ~ normal(flux_hat[flux_measurement_ix], flux_measurment_scale);
    measurement ~ lognormal(log(relative_metabolite_hat[measurement_ix]), measurement_scale);
  }
}
generated quantities {
  vector[N_measurement] measurement_pred;
  real flux[N_reaction] = get_fluxes(metabolite_concentration_hat,
                                     ode_pars,
                                     known_reals);
  for (n in 1:N_measurement)
    measurement_pred[n] = lognormal_rng(log(metabolite_concentration_hat[measurement_ix[n]]), measurement_scale[n]);
}
