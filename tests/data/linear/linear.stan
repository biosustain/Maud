functions{
#include big_k_rate_equations.stan
#include partial_sums.stan
#include haldane_relationships.stan
#include allostery.stan
  vector get_fluxes(real[] m, real[] p){
  real empty_array[0];
  real Tr_r2 = p[3]*p[11]
                * (m[1]/p[9])^(-1*-1) 
                - p[3]*p[11]/p[6]
                * (p[9])^-1 
                * (p[10])^1 
                * (m[2]/p[10])^1 ;

real Dr_r2 = (1 + m[1]/p[9])^(-1*-1) 
                + (1 + m[2]/p[10])^1 
                - 1;


real Dr_reg_r2 = (m[2]/p[14]) ;


  real Tr_r1 = p[4]*p[17]
                * (p[1]/p[15])^(-1*-1) 
                - p[4]*p[17]/p[7]
                * (p[15])^-1 
                * (p[16])^1 
                * (m[1]/p[16])^1 ;

real Dr_r1 = (1 + p[1]/p[15])^(-1*-1) 
                + (1 + m[1]/p[16])^1 
                - 1;


real Dr_reg_r1 = 0;
  real Tr_r3 = p[5]*p[22]
                * (m[2]/p[20])^(-1*-1) 
                - p[5]*p[22]/p[8]
                * (p[20])^-1 
                * (p[21])^1 
                * (p[2]/p[21])^1 ;

real Dr_r3 = (1 + m[2]/p[20])^(-1*-1) 
                + (1 + p[2]/p[21])^1 
                - 1;


real Dr_reg_r3 = 0;
  real free_enzyme_ratio_r2 = get_free_enzyme_ratio_modular_rate_law(Tr_r2, Dr_r2, Dr_reg_r2);
  real free_enzyme_ratio_r1 = get_free_enzyme_ratio_modular_rate_law(Tr_r1, Dr_r1, Dr_reg_r1);
  return [
    modular_rate_law(Tr_r2, Dr_r2, Dr_reg_r2)*get_regulatory_effect(empty_array,{m[1]},free_enzyme_ratio_r2,empty_array,{p[12]},p[13]),
    modular_rate_law(Tr_r1, Dr_r1, Dr_reg_r1)*get_regulatory_effect({m[2]},empty_array,free_enzyme_ratio_r1,{p[18]},empty_array,p[19]),
    modular_rate_law(Tr_r3, Dr_r3, Dr_reg_r3)
  ]';
}
  real[] ode_func(real t, real[] m, real[] p, real[] xr, int[] xi){
  vector[3] fluxes = get_fluxes(m, p);
  return {
    -1*fluxes[1]+1*fluxes[2],
    1*fluxes[1]-1*fluxes[3]
  };
}
}
data {
  // dimensions
  int<lower=1> N_mic;         // Total number of metabolites in compartments
  int<lower=1> N_unbalanced;  // 'Unbalanced' metabolites can have changing concentration at steady state
  int<lower=1> N_kinetic_parameters;
  int<lower=1> N_reaction;
  int<lower=1> N_enzyme;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_enzyme_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=1> N_metabolite;  // NB metabolites in multiple compartments only count once here
  // measurements
  int<lower=1,upper=N_mic> unbalanced_mic_ix[N_unbalanced];
  int<lower=1,upper=N_mic> balanced_mic_ix[N_mic-N_unbalanced];
  int<lower=1,upper=N_experiment> experiment_yconc[N_conc_measurement];
  int<lower=1,upper=N_mic> mic_ix_yconc[N_conc_measurement];
  real yconc[N_conc_measurement];
  vector<lower=0>[N_conc_measurement] sigma_conc;
  int<lower=1,upper=N_experiment> experiment_yflux[N_flux_measurement];
  int<lower=1,upper=N_reaction> reaction_yflux[N_flux_measurement];
  real yflux[N_flux_measurement];
  vector<lower=0>[N_flux_measurement] sigma_flux;
  int<lower=1,upper=N_experiment> experiment_yenz[N_enzyme_measurement];
  int<lower=1,upper=N_enzyme> enzyme_yenz[N_enzyme_measurement];
  real yenz[N_enzyme_measurement];
  vector<lower=0>[N_enzyme_measurement] sigma_enz;
  // hardcoded priors
  vector[N_metabolite] prior_loc_formation_energy;
  vector<lower=0>[N_metabolite] prior_scale_formation_energy;
  vector[N_kinetic_parameters] prior_loc_kinetic_parameters;
  vector<lower=0>[N_kinetic_parameters] prior_scale_kinetic_parameters;
  real prior_loc_unbalanced[N_experiment, N_unbalanced];
  real<lower=0> prior_scale_unbalanced[N_experiment, N_unbalanced];
  real prior_loc_enzyme[N_experiment, N_enzyme];
  real<lower=0> prior_scale_enzyme[N_experiment, N_enzyme];
  // network properties
  matrix[N_mic, N_enzyme] stoichiometric_matrix;
  int<lower=1,upper=N_metabolite> metabolite_ix_stoichiometric_matrix[N_mic];
  matrix<lower=0,upper=1>[N_experiment, N_enzyme] knockout_enzymes;
  // configuration
  real<lower=0> conc_init[N_experiment, N_mic-N_unbalanced];
  real rtol;
  real ftol;
  int steps;
  int<lower=0,upper=1> LIKELIHOOD;  // set to 0 for priors-only mode
  real<lower=0> timepoint;
}
transformed data {
  real xr[0];
  int xi[0];
  real minus_RT = - 0.008314 * 298.15;

}
parameters {
  vector[N_metabolite] formation_energy;
  vector<lower=0>[N_kinetic_parameters] kinetic_parameters;
  matrix<lower=0>[N_experiment, N_enzyme] enzyme_concentration;
  matrix<lower=0>[N_experiment, N_unbalanced] conc_unbalanced;
}
transformed parameters {
  real initial_time = 0;
  vector<lower=0>[N_mic] conc[N_experiment];
  matrix[N_experiment, N_reaction] flux;
  vector[N_enzyme] delta_g = stoichiometric_matrix' * formation_energy[metabolite_ix_stoichiometric_matrix];
  for (e in 1:N_experiment){
    vector[N_enzyme] temp_enzyme_concentration = (enzyme_concentration[e, ] .* knockout_enzymes[e, ])';
    vector[N_enzyme] keq = exp(delta_g / minus_RT);
    vector[N_unbalanced+N_enzyme+N_enzyme+N_kinetic_parameters] theta = append_row(append_row(append_row(
      conc_unbalanced[e, ]', temp_enzyme_concentration), keq), kinetic_parameters);
    conc[e, balanced_mic_ix] = to_vector(integrate_ode_bdf(
                                    ode_func,
                                    conc_init[e,],
                                    initial_time,
                                    rep_array(timepoint, 1),
                                    to_array_1d(theta),
                                    xr,
                                    rep_array(0, 1),
                                    1e-8, 1e-12, 1e5
                                  )[1, ]); 
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e, ]';
    flux[e, ] = get_fluxes(to_array_1d(conc[e, balanced_mic_ix]), to_array_1d(theta))';
  }
}
model {
  kinetic_parameters ~ lognormal(log(prior_loc_kinetic_parameters), prior_scale_kinetic_parameters);
  formation_energy ~ normal(prior_loc_formation_energy, prior_scale_formation_energy);
  for (e in 1:N_experiment){
    conc_unbalanced[e] ~ lognormal(log(prior_loc_unbalanced[e]), prior_scale_unbalanced[e]);
    enzyme_concentration[e] ~ lognormal(log(prior_loc_enzyme[e]), prior_scale_enzyme[e]);
  }
  if (LIKELIHOOD == 1){
    target += reduce_sum(partial_sum_conc, yconc, 1,
                         conc, experiment_yconc, mic_ix_yconc, sigma_conc);
    target += reduce_sum(partial_sum_enz, yenz, 1,
                         enzyme_concentration, experiment_yenz, enzyme_yenz, sigma_enz);
    target += reduce_sum(partial_sum_flux, yflux, 1,
                         flux, experiment_yflux, reaction_yflux, sigma_flux);
  }
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_enzyme_measurement] yenz_sim;
  vector[N_flux_measurement] yflux_sim;
  vector[N_flux_measurement+N_conc_measurement] log_like;


  for (f in 1:N_flux_measurement){
    log_like[f] = normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }

  for (c in 1:N_conc_measurement){
    log_like[N_flux_measurement+c] = lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }


  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }
  for (ec in 1:N_enzyme_measurement){
    yenz_sim[ec] = lognormal_rng(log(enzyme_concentration[experiment_yenz[ec], enzyme_yenz[ec]]), sigma_enz[ec]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}
