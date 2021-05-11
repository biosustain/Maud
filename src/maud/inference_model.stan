functions{
#include allostery.stan
#include modular_rate_law.stan
#include drain_reaction.stan
#include dbalanced_dt.stan
#include partial_sums.stan
#include thermodynamics.stan
  vector unz_1d(vector[] priors, vector z){
    return priors[1] + priors[2] .* z;
  }
  vector unz_log_1d(vector[] priors, vector z){
    return exp(log(priors[1]) + priors[2] .* z);
  }
  vector[] unz_2d(vector[,] priors, vector[] z){
    array[size(z)] vector [rows(z[1])] out;
    for (ex in 1:size(z)){
      out[ex] = unz_1d(priors[:,ex], z[ex]);
    }
    return out;
  }
  vector[] unz_log_2d(vector[,] priors, vector[] z){
    array[size(z)] vector [rows(z[1])] out;
    for (ex in 1:size(z)){
      out[ex] = unz_log_1d(priors[:,ex], z[ex]);
    }
    return out;
  }
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
  int<lower=0> N_phosphorylation_enzymes;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=0> N_ki;
  int<lower=0> N_ai;
  int<lower=0> N_aa;
  int<lower=0> N_ae;
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
  array[2] vector[N_metabolite] priors_dgf;
  array[2] vector[N_enzyme] priors_kcat;
  array[2] vector[N_km] priors_km;
  array[2] vector[N_ki] priors_ki;
  array[2] vector[N_ai] priors_diss_t;
  array[2] vector[N_aa] priors_diss_r;
  array[2] vector[N_ae] priors_transfer_constant;
  array[2] vector[N_phosphorylation_enzymes] priors_kcat_phos;
  array[2, N_experiment] vector[N_phosphorylation_enzymes] priors_conc_phos;
  array[2, N_experiment] vector[N_unbalanced] priors_conc_unbalanced;
  array[2, N_experiment] vector[N_enzyme] priors_conc_enzyme;
  array[2, N_experiment] vector[N_drain] priors_drain;
  // network properties
  matrix[N_mic, N_enzyme] S_enz;
  matrix[N_mic, N_drain] S_drain;
  matrix[N_mic, N_drain+N_enzyme] S_full;
  int<lower=1,upper=N_metabolite> mic_to_met[N_mic];
  vector[N_enzyme] water_stoichiometry;
  matrix[N_reaction, N_enzyme] S_to_flux_map;
  matrix<lower=0,upper=1>[N_experiment, N_enzyme] is_knockout;
  matrix<lower=0,upper=1>[N_experiment, N_phosphorylation_enzymes] is_phos_knockout;
  matrix[N_phosphorylation_enzymes, N_enzyme] S_phos_act;
  matrix[N_phosphorylation_enzymes, N_enzyme] S_phos_inh;
  int<lower=0,upper=N_km> km_lookup[N_mic, N_enzyme];
  int<lower=0,upper=N_mic> n_ci[N_enzyme];
  int<lower=0,upper=N_mic> n_ai[N_enzyme];
  int<lower=0,upper=N_mic> n_aa[N_enzyme];
  int<lower=0,upper=N_mic> ci_ix[N_ki];
  int<lower=0,upper=N_mic> ai_ix[N_ai];
  int<lower=0,upper=N_mic> aa_ix[N_aa];
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
  matrix[N_experiment, N_phosphorylation_enzymes] phos_knockout =
    rep_matrix(1, N_experiment, N_phosphorylation_enzymes) - is_phos_knockout;
}
parameters {
  vector[N_metabolite] dgf_z;
  array[N_experiment] vector[N_drain] drain_z;
  vector[N_enzyme] log_kcat_z;
  vector[N_km] log_km_z;
  vector[N_phosphorylation_enzymes] log_kcat_phos_z;
  array[N_experiment] vector[N_enzyme] log_conc_enzyme_z;
  array[N_experiment] vector[N_phosphorylation_enzymes] log_conc_phos_z;
  array[N_experiment] vector[N_unbalanced] log_conc_unbalanced_z;
  vector[N_ki] log_ki_z;
  vector[N_ai] log_diss_t_z;
  vector[N_aa] log_diss_r_z;
  vector[N_ae] log_transfer_constant_z;
}
transformed parameters {
  // rescale
  vector[N_metabolite] dgf = unz_1d(priors_dgf, dgf_z);
  vector[N_km] km = unz_log_1d(priors_km, log_km_z);
  vector[N_ki] ki = unz_log_1d(priors_ki, log_ki_z);
  vector[N_enzyme] kcat = unz_log_1d(priors_kcat, log_kcat_z);
  vector[N_ai] diss_t = unz_log_1d(priors_diss_t, log_diss_t_z);
  vector[N_aa] diss_r = unz_log_1d(priors_diss_r, log_diss_r_z);
  vector[N_ae] transfer_constant = unz_log_1d(priors_transfer_constant, log_transfer_constant_z);
  vector[N_phosphorylation_enzymes] kcat_phos = unz_log_1d(priors_kcat_phos, log_kcat_phos_z);
  array[N_experiment] vector[N_drain] drain = unz_2d(priors_drain, drain_z);
  array[N_experiment] vector[N_enzyme] conc_enzyme = unz_log_2d(priors_conc_enzyme, log_conc_enzyme_z);
  array[N_experiment] vector[N_unbalanced] conc_unbalanced = unz_log_2d(priors_conc_unbalanced, log_conc_unbalanced_z);
  array[N_experiment] vector[N_phosphorylation_enzymes] conc_phos = unz_log_2d(priors_conc_phos, log_conc_phos_z);
  // transform
  array[N_experiment] vector<lower=0>[N_mic] conc;
  array[N_experiment] vector[N_reaction] flux;
  vector[N_enzyme] keq = get_keq(S_enz, dgf, mic_to_met, water_stoichiometry);
  for (e in 1:N_experiment){
    vector[N_enzyme] conc_enzyme_experiment = conc_enzyme[e] .* knockout[e]';
    vector[N_phosphorylation_enzymes] conc_phos_experiment = conc_phos[e] .* phos_knockout[e]';
    vector[N_mic-N_unbalanced] conc_balanced[2] = ode_bdf_tol(dbalanced_dt,
                                                              conc_init[e, balanced_mic_ix],
                                                              initial_time,
                                                              {timepoint, timepoint + 10},
                                                              rel_tol, abs_tol, max_num_steps,
                                                              conc_unbalanced[e],
                                                              balanced_mic_ix,
                                                              unbalanced_mic_ix,
                                                              conc_enzyme_experiment,
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
                                                              diss_t,
                                                              diss_r,
                                                              transfer_constant,
                                                              subunits,
                                                              conc_phos_experiment,
                                                              kcat_phos,
                                                              S_phos_act,
                                                              S_phos_inh,
                                                              drain[e]);
    conc[e, balanced_mic_ix] = conc_balanced[1];
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e];
    flux[e] = (S_to_flux_map * get_flux_enz(conc[e],
                              conc_enzyme_experiment,
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
                              diss_t,
                              diss_r,
                              transfer_constant,
                              subunits,
                              conc_phos_experiment,
                              kcat_phos,
                              S_phos_act,
                              S_phos_inh));
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
      print("enzyme concentration: ", conc_enzyme_experiment);
      print("ki: ", ki);
      print("tense dissociation constants: ", diss_t);
      print("relaxed dissociation constants: ", diss_r);
      print("transfer constants: ", transfer_constant);
    }
  }
}
model {
  target += std_normal_lpdf(log_kcat_z |);
  target += std_normal_lpdf(log_km_z |);
  target += std_normal_lpdf(log_ki_z |);
  target += std_normal_lpdf(log_diss_t_z |);
  target += std_normal_lpdf(log_diss_r_z |);
  target += std_normal_lpdf(log_transfer_constant_z |);
  target += std_normal_lpdf(dgf_z |);
  target += std_normal_lpdf(log_kcat_phos_z |);
  for (ex in 1:N_experiment){
    target += std_normal_lpdf(log_conc_unbalanced_z[ex] |);
    target += std_normal_lpdf(log_conc_enzyme_z[ex] |);
    target += std_normal_lpdf(log_conc_phos_z[ex] |);
    target += std_normal_lpdf(drain_z[ex] |);
  }
  if (LIKELIHOOD == 1){
    target += reduce_sum(partial_sum_conc, yconc, 1,
                         conc, experiment_yconc, mic_ix_yconc, sigma_conc);
    target += reduce_sum(partial_sum_enz, yenz, 1,
                         conc_enzyme, experiment_yenz, enzyme_yenz, sigma_enz);
    target += reduce_sum(partial_sum_flux, yflux, 1,
                         flux, experiment_yflux, reaction_yflux, sigma_flux);
  }
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_flux_measurement] yflux_sim;
  vector[N_conc_measurement] log_lik_conc;
  vector[N_flux_measurement] log_lik_flux;
  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
    log_lik_conc[c] = lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
    log_lik_flux[f] = normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}
