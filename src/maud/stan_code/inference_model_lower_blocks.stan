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
  vector[N_metabolite] formation_energy;
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
  vector[N_enzyme] adjustment;
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
  vector<lower=0>[N_kinetic_parameters] kinetic_parameters;
  matrix<lower=0>[N_experiment, N_enzyme] enzyme_concentration;
  matrix<lower=0>[N_experiment, N_unbalanced] conc_unbalanced;
}
transformed parameters {
  real initial_time = 0;
  vector<lower=0>[N_mic] conc[N_experiment];
  matrix[N_experiment, N_reaction] flux;
  vector[N_enzyme] delta_g = stoichiometric_matrix' * formation_energy[metabolite_ix_stoichiometric_matrix] + adjustment;
  vector[N_enzyme] keq = exp(delta_g / minus_RT);

  for (e in 1:N_experiment){
    vector[N_enzyme] temp_enzyme_concentration = (enzyme_concentration[e, ] .* knockout_enzymes[e, ])';
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
                                    1e-10, 1e-12, 1e5
                                  )[1, ]); 
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e, ]';
    flux[e, ] = get_fluxes(to_array_1d(conc[e, balanced_mic_ix]), to_array_1d(theta))';
  }
}
model {
  kinetic_parameters ~ lognormal(log(prior_loc_kinetic_parameters), prior_scale_kinetic_parameters);
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
