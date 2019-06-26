functions {
#include big_k_rate_equations.stan
#include haldane_relationships.stan
real [] get_fluxes(real[] metabolites, real[] params, real[] known_reals){
  real r2_Keq = get_Keq(params[1], known_reals[5], known_reals[4]);
  real r2_Kip = get_Kip_ordered_unibi(r2_Keq, params[10], params[9], params[5], params[6]);
  real r2_Kiq = get_Kiq_ordered_unibi(r2_Keq, params[7], params[8], params[5], params[6]);
  real r3_Keq = get_Keq(params[2], known_reals[5], known_reals[4]);
  real r6_Keq = get_Keq(params[3], known_reals[5], known_reals[4]);
  return {
    fixed_flux(params[4]),
    ordered_unibi(metabolites[1], metabolites[2], metabolites[3], known_reals[1]*params[5], known_reals[1]*params[6], params[7], params[8], params[9], params[10], r2_Kip, r2_Kiq, r2_Keq),
    uniuni(metabolites[3], metabolites[4], known_reals[2]*params[11], known_reals[2]*params[12], params[13], r3_Keq),
    irr_mass_action(metabolites[4], params[14]),
    irr_mass_action(metabolites[2], params[15]),
    uniuni(metabolites[2], metabolites[3], known_reals[3]*params[16], known_reals[3]*params[17], params[18], r6_Keq)
  };
}
real[] get_odes(real[] fluxes){
  return {
    1*fluxes[1]-1*fluxes[2],
    1*fluxes[2]-1*fluxes[5]-1*fluxes[6],
    1*fluxes[2]-1*fluxes[3]+1*fluxes[6],
    1*fluxes[3]-1*fluxes[4]
  };
}
real[] steady_state_equation(real t,
                             real[] metabolites,
                             real[] params,
                             real[] known_reals,
                             int[] known_ints){
  for (m in 1:size(metabolites)){
    if (metabolites[m] < 0){
      reject("Metabolite ", m, " is ", metabolites[m], " but should be greater than zero.");
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
  int<lower=1> N_measurement_flux; // number of measured fluxes
  int<lower=1> N_measurement_conc; // number of measured concentrations
  // measurements
  int<lower=0,upper=N_metabolite> metabolite_ix[N_measurement_conc];
  int<lower=1,upper=N_experiment> experiment_ix_conc[N_measurement_conc];
  vector[N_measurement_conc] measurement_conc;
  vector<lower=0>[N_measurement_conc] measurement_scale_conc;
  int<lower=0,upper=N_reaction> reaction_ix[N_measurement_flux];
  int<lower=1,upper=N_experiment> experiment_ix_flux[N_measurement_flux];
  vector[N_measurement_flux] measurement_flux;
  vector<lower=0>[N_measurement_flux] measurement_scale_flux;
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
    vector[N_measurement_conc] conc_hat;
    vector[N_measurement_flux] flux_hat;
    for (mc in 1:N_measurement_conc){
      conc_hat[mc] = metabolite_concentration[metabolite_ix[mc], experiment_ix_conc[mc]];
    }
    for (mf in 1:N_measurement_flux){
      flux_hat[mf] = flux[reaction_ix[mf], experiment_ix_flux[mf]];
    }
    measurement_conc ~ lognormal(log(conc_hat), measurement_scale_conc);
    measurement_flux ~ normal(flux_hat, measurement_scale_flux);
  }
}
generated quantities {
  vector[N_measurement_flux] flux_pred;
  vector[N_measurement_conc] conc_pred;
  for (mc in 1:N_measurement_conc){
    conc_pred[mc] = lognormal_rng(log(metabolite_concentration[metabolite_ix[mc], experiment_ix_conc[mc]]),
                                  measurement_scale_conc[mc]);
  }
  for (mf in 1:N_measurement_flux){
    flux_pred[mf] = normal_rng(flux[reaction_ix[mf], experiment_ix_flux[mf]], measurement_scale_flux[mf]);
  }
}
