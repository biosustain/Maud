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
  vector[N_conc_measurement] yconc;
  vector<lower=0>[N_conc_measurement] sigma_conc;
  int<lower=1,upper=N_experiment> experiment_yflux[N_flux_measurement];
  int<lower=1,upper=N_reaction> reaction_yflux[N_flux_measurement];
  vector[N_flux_measurement] yflux;
  vector<lower=0>[N_flux_measurement] sigma_flux;
  int<lower=1,upper=N_experiment> experiment_yenz[N_enzyme_measurement];
  int<lower=1,upper=N_enzyme> enzyme_yenz[N_enzyme_measurement];
  vector[N_enzyme_measurement] yenz;
  vector<lower=0>[N_enzyme_measurement] sigma_enz;
  // hardcoded priors
  vector[N_metabolite] prior_loc_formation_energy;
  vector<lower=0>[N_metabolite] prior_scale_formation_energy;
  vector[N_kinetic_parameters] prior_loc_kinetic_parameters;
  vector<lower=0>[N_kinetic_parameters] prior_scale_kinetic_parameters;
  real prior_loc_conc[N_experiment, N_mic];
  real<lower=0> prior_scale_conc[N_experiment, N_mic];
  real prior_loc_enzyme[N_experiment, N_enzyme];
  real<lower=0> prior_scale_enzyme[N_experiment, N_enzyme];
  // network properties
  matrix[N_mic, N_enzyme] stoichiometric_matrix;
  int<lower=1,upper=N_metabolite> metabolite_ix_stoichiometric_matrix[N_mic];
}
transformed data {
  real xr[0];
  int xi[0];
  real minus_RT = - 0.008314 * 298.15;

}
parameters {
  vector[N_metabolite] formation_energy;
  vector<lower=0>[N_kinetic_parameters] kinetic_parameters;
  vector<lower=0>[N_enzyme] enzyme_concentration[N_experiment];
  vector<lower=0>[N_mic] conc[N_experiment];
}
transformed parameters {
  real initial_time = 0;
  vector[N_reaction] flux[N_experiment];
  vector[N_enzyme] delta_g = stoichiometric_matrix' * formation_energy[metabolite_ix_stoichiometric_matrix];
  for (e in 1:N_experiment){ 
    vector[N_enzyme] keq = exp(delta_g / minus_RT);
    vector[N_unbalanced+N_enzyme+N_enzyme+N_kinetic_parameters] theta = append_row(append_row(append_row(
      conc[e, unbalanced_mic_ix], enzyme_concentration[e]), keq), kinetic_parameters);
    flux[e] = get_fluxes(to_array_1d(conc[e, balanced_mic_ix]), to_array_1d(theta));
  }
}
model {
  kinetic_parameters ~ lognormal(log(prior_loc_kinetic_parameters), prior_scale_kinetic_parameters);
  formation_energy ~ normal(prior_loc_formation_energy, prior_scale_formation_energy);
  for (e in 1:N_experiment){
    conc[e, ] ~ lognormal(log(prior_loc_conc[e]), prior_scale_conc[e]);
    enzyme_concentration[e] ~ lognormal(log(prior_loc_enzyme[e]), prior_scale_enzyme[e]);
  }
  for (c in 1:N_conc_measurement){
    target += lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }
  for (ec in 1:N_enzyme_measurement){
    target += lognormal_lpdf(yenz[ec] | log(enzyme_concentration[experiment_yenz[ec], enzyme_yenz[ec]]), sigma_enz[ec]);
  }
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_enzyme_measurement] yenz_sim;
  vector[N_flux_measurement] yflux_sim;


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
