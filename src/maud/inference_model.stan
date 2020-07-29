functions{
#include modular_rate_law.stan
#include dbalanced_dt.stan
#include partial_sums.stan
}
data {
  // dimensions
  int<lower=1> N_mic;
  int<lower=1> N_unbalanced;
  int<lower=1> N_metabolite;
  int<lower=1> N_km;
  int<lower=1> N_reaction;
  int<lower=1> N_enzyme;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_enzyme_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=0> N_competitive_inhibitor;
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
  vector[N_enzyme] prior_loc_kcat;
  vector<lower=0>[N_enzyme] prior_scale_kcat;
  vector[N_km] prior_loc_km;
  vector<lower=0>[N_km] prior_scale_km;
  vector[N_competitive_inhibitor] prior_loc_ki;
  vector<lower=0>[N_competitive_inhibitor] prior_scale_ki;
  real prior_loc_unbalanced[N_experiment, N_unbalanced];
  real<lower=0> prior_scale_unbalanced[N_experiment, N_unbalanced];
  real prior_loc_enzyme[N_experiment, N_enzyme];
  real<lower=0> prior_scale_enzyme[N_experiment, N_enzyme];
  // network properties
  matrix[N_mic, N_enzyme] S;
  int<lower=1,upper=N_metabolite> metabolite_ix_stoichiometric_matrix[N_mic];
  matrix<lower=0,upper=1>[N_experiment, N_enzyme] is_knockout;
  int<lower=0,upper=N_km> km_lookup[N_mic, N_enzyme];
  int<lower=0,upper=N_mic> n_ci[N_enzyme];
  int<lower=0,upper=N_mic> ci_ix[N_competitive_inhibitor];
  // configuration
  vector<lower=0>[N_mic] conc_init[N_experiment];
  real rel_tol;
  real abs_tol;
  int max_num_steps;
  int<lower=0,upper=1> LIKELIHOOD;  // set to 0 for priors-only mode
  real<lower=0> timepoint;
}
transformed data {
  real minus_RT = - 0.008314 * 298.15;
  real initial_time = 0;
  matrix[N_experiment, N_enzyme] knockout =
    rep_matrix(1, N_experiment, N_enzyme) - is_knockout;
}
parameters {
  vector[N_metabolite] formation_energy;
  vector<lower=0>[N_enzyme] kcat;
  vector<lower=0>[N_km] km;
  matrix<lower=0>[N_experiment, N_enzyme] enzyme;
  matrix<lower=0>[N_experiment, N_unbalanced] conc_unbalanced;
  vector<lower=0>[N_competitive_inhibitor] ki;
}
transformed parameters {
  vector<lower=0>[N_mic] conc[N_experiment];
  matrix[N_experiment, N_reaction] flux;
  vector[N_enzyme] delta_g = S' * formation_energy[metabolite_ix_stoichiometric_matrix];
  vector[N_enzyme] keq = exp(delta_g / minus_RT);
  for (e in 1:N_experiment){
    vector[N_enzyme] experiment_enzyme = enzyme[e]' .* knockout[e]';
    conc[e, balanced_mic_ix] = ode_bdf_tol(dbalanced_dt,
                                           conc_init[e, balanced_mic_ix],
                                           initial_time,
                                           rep_array(timepoint, 1),
                                           rel_tol, abs_tol, max_num_steps,
                                           conc_unbalanced[e]',
                                           balanced_mic_ix,
                                           unbalanced_mic_ix,
                                           experiment_enzyme,
                                           km,
                                           km_lookup,
                                           S,
                                           kcat,
                                           keq,
                                           ci_ix,
                                           n_ci,
                                           ki)[1]; 
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e]';
    flux[e] = get_flux(conc[e],
                       experiment_enzyme,
                       km,
                       km_lookup,
                       S,
                       kcat,
                       keq,
                       ci_ix,
                       n_ci,
                       ki)';
  }
  print(conc_unbalanced);
  print(conc);
  print(flux);
  print(km);
  print(ki);
  print(kcat);
  print(enzyme);
}
  model {
    target += lognormal_lpdf(kcat | log(prior_loc_kcat), prior_scale_kcat);
    target += lognormal_lpdf(km | log(prior_loc_km), prior_scale_km);
    target += lognormal_lpdf(ki | log(prior_loc_ki), prior_scale_ki);
    target += normal_lpdf(formation_energy |
                          prior_loc_formation_energy,
                          prior_scale_formation_energy);
    for (e in 1:N_experiment){
      target += lognormal_lpdf(conc_unbalanced[e] |
                               log(prior_loc_unbalanced[e]),
                               prior_scale_unbalanced[e]);
      target += lognormal_lpdf(enzyme[e] |
                               log(prior_loc_enzyme[e]),
                               prior_scale_enzyme[e]);
    }
    if (LIKELIHOOD == 1){
      target += reduce_sum(partial_sum_conc, yconc, 1,
                           conc, experiment_yconc, mic_ix_yconc, sigma_conc);
      target += reduce_sum(partial_sum_enz, yenz, 1,
                           enzyme, experiment_yenz, enzyme_yenz, sigma_enz);
      target += reduce_sum(partial_sum_flux, yflux, 1,
                           flux, experiment_yflux, reaction_yflux, sigma_flux);
    }
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_enzyme_measurement] yenz_sim;
  vector[N_flux_measurement] yflux_sim;
  vector[N_flux_measurement+N_conc_measurement] log_lik;
  for (f in 1:N_flux_measurement){
    log_lik[f] = normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
  for (c in 1:N_conc_measurement){
    log_lik[N_flux_measurement+c] = lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }
  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }
  for (ec in 1:N_enzyme_measurement){
    yenz_sim[ec] = lognormal_rng(log(enzyme[experiment_yenz[ec], enzyme_yenz[ec]]), sigma_enz[ec]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}
