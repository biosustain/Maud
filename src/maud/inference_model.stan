functions{
#include allostery.stan
#include modular_rate_law.stan
#include drain_reaction.stan
#include dbalanced_dt.stan
#include partial_sums.stan
#include thermodynamics.stan
}
data {
  // dimensions
  int<lower=1> N_mic;
  int<lower=1> N_unbalanced;
  int<lower=1> N_metabolite;
  int<lower=1> N_km;
  int<lower=1> N_reaction;
  int<lower=1> N_enzyme;
  int<lower=0> N_drain;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=0> N_competitive_inhibitor;
  int<lower=0> N_allosteric_inhibitor;
  int<lower=0> N_allosteric_activator;
  int<lower=0> N_allosteric_enzyme;
  int<lower=0> N_enzyme_measurement;
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
  int<lower=0,upper=N_experiment> experiment_yenz[N_enzyme_measurement];
  int<lower=0,upper=N_enzyme> enzyme_yenz[N_enzyme_measurement];
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
  matrix[N_experiment, N_unbalanced] prior_loc_unbalanced;
  matrix<lower=0>[N_experiment, N_unbalanced] prior_scale_unbalanced;
  matrix[N_experiment, N_enzyme] prior_loc_enzyme;
  matrix<lower=0>[N_experiment, N_enzyme] prior_scale_enzyme;
  matrix[N_experiment, N_drain] prior_loc_drain;
  matrix<lower=0>[N_experiment, N_drain] prior_scale_drain;
  // network properties
  matrix[N_mic, N_enzyme] S_enz;
  matrix[N_mic, N_drain] S_drain;
  matrix[N_mic, N_drain+N_enzyme] S_full;
  int<lower=1,upper=N_metabolite> mic_to_met[N_mic];
  vector[N_enzyme] water_stoichiometry;
  matrix[N_reaction, N_enzyme] S_to_flux_map;
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
  real initial_time = 0;
  matrix[N_experiment, N_enzyme] knockout =
    rep_matrix(1, N_experiment, N_enzyme) - is_knockout;
}
parameters {
  vector[N_metabolite] formation_energy_z;
  matrix[N_experiment, N_drain] drain_z;
  vector[N_enzyme] log_kcat_z;
  vector[N_km] log_km_z;
  matrix[N_experiment, N_enzyme] log_enzyme_z;
  matrix[N_experiment, N_unbalanced] log_conc_unbalanced_z;
  vector[N_competitive_inhibitor] log_ki_z;
  vector[N_allosteric_inhibitor] log_dissociation_constant_t_z;
  vector[N_allosteric_activator] log_dissociation_constant_r_z;
  vector[N_allosteric_enzyme] log_transfer_constant_z;
}
transformed parameters {
  // rescale
  vector[N_metabolite] formation_energy =
    prior_loc_formation_energy + formation_energy_z .* prior_scale_formation_energy;
  vector[N_km] km = exp(log(prior_loc_km) + log_km_z .* prior_scale_km);
  vector[N_competitive_inhibitor] ki =
    exp(log(prior_loc_ki) + log_ki_z .* prior_scale_ki);
  vector[N_enzyme] kcat = exp(log(prior_loc_kcat) + log_kcat_z .* prior_scale_kcat);
  vector[N_allosteric_inhibitor] dissociation_constant_t =
    exp(log(prior_loc_diss_t) + log_dissociation_constant_t_z .* prior_scale_diss_t);
  vector[N_allosteric_activator] dissociation_constant_r =
    exp(log(prior_loc_diss_r) + log_dissociation_constant_r_z .* prior_scale_diss_r);
  vector[N_allosteric_enzyme] transfer_constant =
    exp(log(prior_loc_tc) + log_transfer_constant_z .* prior_scale_tc);
  matrix[N_experiment, N_drain] drain;
  matrix[N_experiment, N_enzyme] enzyme;
  matrix[N_experiment, N_unbalanced] conc_unbalanced;
  for (ex in 1:N_experiment){
    drain[ex] = prior_loc_drain[ex] + drain_z[ex] .* prior_scale_drain[ex];
    enzyme[ex] =
      exp(log(prior_loc_enzyme[ex]) + log_enzyme_z[ex] .* prior_scale_enzyme[ex]);
    conc_unbalanced[ex] =
      exp(log(prior_loc_unbalanced[ex])
          + log_conc_unbalanced_z[ex] .* prior_scale_unbalanced[ex]);
  }
  // transform
  vector<lower=0>[N_mic] conc[N_experiment];
  matrix[N_experiment, N_reaction] flux;
  vector[N_enzyme] keq = get_keq(S_enz,
                                 formation_energy,
                                 mic_to_met,
                                 water_stoichiometry);
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
                                                              S_enz,
                                                              S_drain,
                                                              S_full,
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
                                                              subunits,
                                                              drain[e,]');
    conc[e, balanced_mic_ix] = conc_balanced[1];
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e]';
    flux[e] = (S_to_flux_map * get_flux_enz(conc[e],
                              experiment_enzyme,
                              km,
                              km_lookup,
                              S_enz,
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
      print("km: ", km);
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
  target += std_normal_lpdf(log_kcat_z |);
  target += std_normal_lpdf(log_km_z |);
  target += std_normal_lpdf(log_ki_z |);
  target += std_normal_lpdf(log_dissociation_constant_t_z |);
  target += std_normal_lpdf(log_dissociation_constant_r_z |);
  target += std_normal_lpdf(log_transfer_constant_z |);
  target += std_normal_lpdf(formation_energy_z |);
  target += std_normal_lpdf(to_vector(log_conc_unbalanced_z) |);
  target += std_normal_lpdf(to_vector(log_enzyme_z) |);
  target += std_normal_lpdf(to_vector(drain_z) |);
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
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}
