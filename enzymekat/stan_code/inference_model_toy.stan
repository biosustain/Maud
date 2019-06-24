functions {
#include big_k_rate_equations.stan
#include haldane_relationships.stan
real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){
  real FBA_Keq = get_Keq(params[2], params[5], params[4]);
  real FBA_Kip = get_Kip_ordered_unibi(FBA_Keq, params[8], params[7], params[3], params[4]);
  real FBA_Kiq = get_Kiq_ordered_unibi(FBA_Keq, params[5], params[6], params[3], params[4]);
  real TDH_Keq = get_Keq(params[9], params[5], params[4]);
  real TPI_Keq = get_Keq(params[13], params[5], params[4]);
  return {
    irr_mass_action(metabolites[2], params[1]),
    ordered_unibi(metabolites[1], metabolites[2], metabolites[3], known_reals[1]*params[3], known_reals[1]*params[4], params[5], params[6], params[7], params[8], FBA_Kip, FBA_Kiq, FBA_Keq),
    uniuni(metabolites[3], metabolites[4], known_reals[2]*params[10], known_reals[2]*params[11], params[12], TDH_Keq),
    uniuni(metabolites[2], metabolites[3], known_reals[3]*params[14], known_reals[3]*params[15], params[16], TPI_Keq),
    irr_mass_action(metabolites[1], params[17]),
    irr_mass_action(metabolites[1], params[18])
  };
}
real[] get_odes(real[] fluxes){
  return {
    1*fluxes[1]-1*fluxes[2],
    1*fluxes[2]-1*fluxes[4]-1*fluxes[5],
    1*fluxes[2]-1*fluxes[3]+1*fluxes[4],
    1*fluxes[3]-1*fluxes[6]
  };
}
real[] steady_state_equation(real t,
                             real[] metabolites,
                             real[] params,
                             real[] known_reals,
                             int[] known_ints){
  for (m in 1:size(metabolites)){
    if (metabolites[m] < 0){
      reject("Metabolite ", m, " is ", metabolites[m], " but should be greater than zero");
    }
  }
  return get_odes(get_fluxes(metabolites, params, known_reals));
}
}

data {
  // dimensions
  int<lower=1> N_metabolite;          // number of ode metabolites
  int<lower=1> N_kinetic_parameter;        // total number of kinetic parameters
  int<lower=1> N_thermodynamic_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_experiment;
  int<lower=1> N_known_real;   // number of known reals
  int<lower=1> N_measurement; // number of measured fluxes
  // measurements
  int<lower=0,upper=N_metabolite> metabolite_ix[N_measurement];
  int<lower=0,upper=N_reaction> reaction_ix[N_measurement];
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
  // ode stuff
  real initial_concentration[N_metabolite];
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
  real<lower=0> kinetic_parameters[N_kinetic_parameter];
}
transformed parameters {
  real metabolite_concentration[N_metabolite, N_experiment];
  real flux[N_reaction, N_experiment];
  for (e in 1:N_experiment){
    metabolite_concentration[,e] =
      integrate_ode_bdf(steady_state_equation,
                        initial_concentration,
                        initial_time,
                        {steady_time},
                        append_array(thermodynamic_parameters, kinetic_parameters),
                        known_reals[,e],
                        known_ints,
                        rel_tol, abs_tol, max_steps)[1];
    flux[,e] =
      get_fluxes(metabolite_concentration[,e],
                 append_array(thermodynamic_parameters, kinetic_parameters),
                 known_reals[,e]);
  }
}
model {
  kinetic_parameters ~ lognormal(prior_location_kinetic, prior_scale_kinetic);
  thermodynamic_parameters ~ normal(prior_location_thermodynamic, prior_scale_thermodynamic);
  if (LIKELIHOOD == 1){
    for (m in 1:N_measurement){
      target += is_flux[m] ?
        normal_lpdf(measurement[m] | flux[reaction_ix[m], experiment_ix[m]], measurement_scale[m]):
        lognormal_lpdf(measurement[m] | log(metabolite_concentration[metabolite_ix, experiment_ix[m]]), measurement_scale[m]);
    }
  }
}
generated quantities {
  vector[N_measurement] measurement_pred;
  for (m in 1:N_measurement){
    measurement_pred[m] = is_flux[m] ?
      normal_lpdf(measurement[m] | flux[reaction_ix[m], experiment_ix[m]], measurement_scale[m]):
      lognormal_lpdf(measurement[m] | log(metabolite_concentration[metabolite_ix, experiment_ix[m]]), measurement_scale[m]);
  }
}
