functions{
#include allostery.stan
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
  int<lower=0> N_allosteric_inhibitor;
  int<lower=0> N_allosteric_activator;
  int<lower=0> N_allosteric_enzyme;
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
  vector[N_allosteric_inhibitor] prior_loc_diss_t;
  vector<lower=0>[N_allosteric_inhibitor] prior_scale_diss_t;
  vector[N_allosteric_activator] prior_loc_diss_r;
  vector<lower=0>[N_allosteric_activator] prior_scale_diss_r;
  vector<lower=0>[N_allosteric_enzyme] prior_loc_tc;
  vector<lower=0>[N_allosteric_enzyme] prior_scale_tc;
  real prior_loc_unbalanced[N_experiment, N_unbalanced];
  real<lower=0> prior_scale_unbalanced[N_experiment, N_unbalanced];
  real prior_loc_enzyme[N_experiment, N_enzyme];
  real<lower=0> prior_scale_enzyme[N_experiment, N_enzyme];
  // network properties
  matrix[N_mic, N_enzyme] S;
  matrix[N_reaction, N_enzyme] S_to_flux_map;
  int<lower=1,upper=N_metabolite> metabolite_ix_stoichiometric_matrix[N_mic];
  matrix<lower=0,upper=1>[N_experiment, N_enzyme] is_knockout;
  int<lower=0,upper=N_km> km_lookup[N_mic, N_enzyme];
  int<lower=0,upper=N_mic> n_ci[N_enzyme];
  int<lower=0,upper=N_mic> n_ai[N_enzyme];
  int<lower=0,upper=N_mic> n_aa[N_enzyme];
  int<lower=0,upper=N_mic> ci_ix[N_competitive_inhibitor];
  int<lower=0,upper=N_mic> ai_ix[N_allosteric_inhibitor];
  int<lower=0,upper=N_mic> aa_ix[N_allosteric_activator];
  int<lower=0> subunits[N_enzyme];
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
  vector<lower=0>[N_allosteric_inhibitor] dissociation_constant_t;
  vector<lower=0>[N_allosteric_activator] dissociation_constant_r;
  vector<lower=0>[N_allosteric_enzyme] transfer_constant;
}
transformed parameters {
  vector<lower=0>[N_mic] conc[N_experiment];
  matrix[N_experiment, N_reaction] flux;
  vector[N_enzyme] delta_g;
  vector[N_enzyme] keq;
  delta_g = S' * formation_energy[metabolite_ix_stoichiometric_matrix];
  keq = exp(delta_g / minus_RT);
  for (e in 1:N_experiment){
    vector[N_enzyme] experiment_enzyme = enzyme[e]' .* knockout[e]';
    vector[N_mic-N_unbalanced] conc_balanced[2] = ode_bdf_tol(dbalanced_dt,
                                                              conc_init[e, balanced_mic_ix],
                                                              initial_time,
                                                              {timepoint, timepoint + 10},
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
                                                              ai_ix,
                                                              aa_ix,
                                                              n_ci,
                                                              n_ai,
                                                              n_aa,
                                                              ki,
                                                              dissociation_constant_t,
                                                              dissociation_constant_r,
                                                              transfer_constant,
                                                              subunits);
    conc[e, balanced_mic_ix] = conc_balanced[1];
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e]';
    flux[e] = (S_to_flux_map * get_flux(conc[e],
                              experiment_enzyme,
                              km,
                              km_lookup,
                              S,
                              kcat,
                              keq,
                              ci_ix,
                              ai_ix,
                              aa_ix,
                              n_ci,
                              n_ai,
                              n_aa,
                              ki,
                              dissociation_constant_t,
                              dissociation_constant_r,
                              transfer_constant,
                              subunits))';
    if (squared_distance(conc_balanced[1], conc_balanced[2]) > 1){
      print("Balanced metabolite concentration at ", timepoint, " seconds is not steady.");
      print("Found ", conc_balanced[1], " at ", timepoint, " seconds and ",
            conc_balanced[2], " at ", timepoint + 10, " seconds.");
    }
    if (sum(conc[e]) > 100){
      print("Metabolite concentration is really high in experiment ", e, ".");
      print("metabolite concentration: ", conc[e]);
      print("kcat: ", kcat);
      print("km: ", kcat);
      print("keq: ", keq);
      print("flux: ", flux[e]);
      print("enzyme concentration: ", enzyme[e]);
      print("ki: ", ki);
      print("tense dissociation constants: ", dissociation_constant_t);
      print("relaxed dissociation constants: ", dissociation_constant_r);
      print("transfer constants: ", transfer_constant);
    }
  }
}
model {
  target += lognormal_lpdf(kcat | log(prior_loc_kcat), prior_scale_kcat);
  target += lognormal_lpdf(km | log(prior_loc_km), prior_scale_km);
  target += lognormal_lpdf(ki | log(prior_loc_ki), prior_scale_ki);
  target += lognormal_lpdf(ki | log(prior_loc_ki), prior_scale_ki);
  target += lognormal_lpdf(dissociation_constant_t |
                           log(prior_loc_diss_t), prior_scale_diss_t);
  target += lognormal_lpdf(dissociation_constant_r |
                           log(prior_loc_diss_r), prior_scale_diss_r);
  target += lognormal_lpdf(transfer_constant |
                           log(prior_loc_tc), prior_scale_tc);
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
