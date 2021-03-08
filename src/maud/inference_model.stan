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
  matrix[N_metabolite, 2] fe_priors;
  matrix[N_enzyme, 2] kcat_priors;
  matrix[N_km, 2] km_priors;
  matrix[N_ki, 2] ki_priors;
  matrix[N_ai, 2] diss_t_priors;
  matrix[N_aa, 2] diss_r_priors;
  matrix[N_ae, 2] tc_priors;
  matrix[N_phosphorylation_enzymes, 2] phos_kcat_priors;
  matrix[N_phosphorylation_enzymes, 2] phos_conc_priors[N_experiment];
  matrix[N_unbalanced, 2] unbalanced_priors[N_experiment];
  matrix[N_enzyme, 2] enzyme_priors[N_experiment];
  matrix[N_drain, 2] drain_priors[N_experiment];
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
  vector[N_metabolite] fe_z;
  matrix[N_experiment, N_drain] drain_z;
  vector[N_enzyme] log_kcat_z;
  vector[N_km] log_km_z;
  vector[N_phosphorylation_enzymes] log_phos_kcat_z;
  matrix[N_experiment, N_enzyme] log_enzyme_z;
  matrix[N_experiment, N_phosphorylation_enzymes] log_phos_conc_z;
  matrix[N_experiment, N_unbalanced] log_conc_unbalanced_z;
  vector[N_ki] log_ki_z;
  vector[N_ai] log_dt_z;
  vector[N_aa] log_dr_z;
  vector[N_ae] log_tc_z;
}
transformed parameters {
  // rescale
  vector[N_metabolite] formation_energy = fe_priors[,1] + fe_z .* fe_priors[,2];
  vector[N_km] km = exp(log(km_priors[,1]) + log_km_z .* km_priors[,2]);
  vector[N_ki] ki = exp(log(ki_priors[,1]) + log_ki_z .* ki_priors[,2]);
  vector[N_enzyme] kcat = exp(log(kcat_priors[,1]) + log_kcat_z .* kcat_priors[,2]);
  vector[N_ai] diss_t = exp(log(diss_t_priors[,1]) + log_dt_z .* diss_t_priors[,2]);
  vector[N_aa] diss_r = exp(log(diss_r_priors[,1]) + log_dt_z .* diss_r_priors[,2]);
  vector[N_ae] transfer_constant =
    exp(log(tc_priors[,1]) + log_tc_z .* tc_priors[,2]);
  vector[N_phosphorylation_enzymes] phos_enzyme_kcat =
    exp(log(phos_kcat_priors[,1]) + log_phos_kcat_z .* phos_kcat_priors[,2]);
  matrix[N_experiment, N_drain] drain;
  matrix[N_experiment, N_enzyme] enzyme;
  matrix[N_experiment, N_unbalanced] conc_unbalanced;
  matrix[N_experiment, N_phosphorylation_enzymes] phos_enzyme_conc;
  for (ex in 1:N_experiment){
    drain[ex] = drain_priors[ex][,1]' + drain_z[ex] .* drain_priors[ex][,2]';
    enzyme[ex] =
      exp(log(enzyme_priors[ex][,1])' + log_enzyme_z[ex] .* enzyme_priors[ex][,2]');
    conc_unbalanced[ex] =
      exp(log(unbalanced_priors[ex][,1])'
          + log_conc_unbalanced_z[ex] .* unbalanced_priors[ex][,2]');
    phos_enzyme_conc[ex] = 
      exp(log(phos_conc_priors[ex][,1])' + log_phos_conc_z[ex] .* phos_conc_priors[ex][,2]');
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
    vector[N_phosphorylation_enzymes] experiment_phos_conc = 
      phos_enzyme_conc[e]' .* phos_knockout[e]';
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
                                                              diss_t,
                                                              diss_r,
                                                              transfer_constant,
                                                              subunits,
                                                              experiment_phos_conc,
                                                              phos_enzyme_kcat,
                                                              S_phos_act,
                                                              S_phos_inh,
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
                              diss_t,
                              diss_r,
                              transfer_constant,
                              subunits,
                              phos_enzyme_conc[e]',
                              phos_enzyme_kcat,
                              S_phos_act,
                              S_phos_inh))';
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
  target += std_normal_lpdf(log_dt_z |);
  target += std_normal_lpdf(log_dr_z |);
  target += std_normal_lpdf(log_tc_z |);
  target += std_normal_lpdf(fe_z |);
  target += std_normal_lpdf(log_phos_kcat_z |);
  target += std_normal_lpdf(to_vector(log_conc_unbalanced_z) |);
  target += std_normal_lpdf(to_vector(log_enzyme_z) |);
  target += std_normal_lpdf(to_vector(log_phos_conc_z) |);
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
