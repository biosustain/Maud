functions{
#include functions.stan
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
  int<lower=1> N_edge;
  int<lower=0> N_phosphorylation_enzymes;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=0> N_enzyme_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=0> N_ki;
  int<lower=0> N_ai;
  int<lower=0> N_aa;
  int<lower=0> N_ae;
  int<lower=0> N_pa;
  int<lower=0> N_pi;
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
  vector[N_metabolite] prior_loc_dgf;
  cov_matrix[N_metabolite] prior_cov_dgf;
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
  matrix[N_mic, N_edge] S;
  int<lower=1,upper=3> edge_type[N_edge];  // 1 = reversible modular rate law, 2 = drain
  int<lower=0,upper=N_enzyme> edge_to_enzyme[N_edge];  // 0 if drain
  int<lower=0,upper=N_drain> edge_to_drain[N_edge];  // 0 if enzyme
  int<lower=0,upper=N_reaction> edge_to_reaction[N_edge];
  int<lower=1,upper=N_metabolite> mic_to_met[N_mic];
  vector[N_edge] water_stoichiometry;
  matrix<lower=0,upper=1>[N_experiment, N_enzyme] is_knockout;
  matrix<lower=0,upper=1>[N_experiment, N_phosphorylation_enzymes] is_phos_knockout;
  int<lower=0,upper=N_km> km_lookup[N_mic, N_edge];
  int<lower=0,upper=N_mic> n_ci[N_edge];
  int<lower=0,upper=N_mic> n_ai[N_edge];
  int<lower=0,upper=N_mic> n_aa[N_edge];
  int<lower=0,upper=N_phosphorylation_enzymes> n_pa[N_edge];
  int<lower=0,upper=N_phosphorylation_enzymes> n_pi[N_edge];
  int<lower=0,upper=N_mic> ix_ci[N_ki];
  int<lower=0,upper=N_mic> ix_ai[N_ai];
  int<lower=0,upper=N_mic> ix_aa[N_aa];
  int<lower=0,upper=N_phosphorylation_enzymes> ix_pa[N_pa];
  int<lower=0,upper=N_phosphorylation_enzymes> ix_pi[N_pi];
  int<lower=1> subunits[N_enzyme];
  // configuration
  vector<lower=0>[N_mic] conc_init[N_experiment];
  real rel_tol_forward; 
  vector[N_mic - N_unbalanced] abs_tol_forward;
  real rel_tol_backward; 
  vector[N_mic - N_unbalanced] abs_tol_backward; 
  real rel_tol_quadrature;
  real abs_tol_quadrature;
  int max_num_steps;
  int num_steps_between_checkpoints;
  int interpolation_polynomial;
  int solver_forward;
  int solver_backward;
  int<lower=0,upper=1> LIKELIHOOD;  // set to 0 for priors-only mode
  real<lower=0> timepoint;
  int<lower=0,upper=1> reject_non_steady;
}
transformed data {
  real initial_time = 0;
  matrix[N_experiment, N_enzyme] knockout = rep_matrix(1, N_experiment, N_enzyme) - is_knockout;
  matrix[N_experiment, N_phosphorylation_enzymes] phos_knockout =
    rep_matrix(1, N_experiment, N_phosphorylation_enzymes) - is_phos_knockout;
}
parameters {
  vector[N_metabolite] dgf;
  vector[N_enzyme] log_kcat_z;
  vector[N_km] log_km_z;
  vector[N_phosphorylation_enzymes] log_kcat_phos_z;
  vector[N_ki] log_ki_z;
  vector[N_ai] log_diss_t_z;
  vector[N_aa] log_diss_r_z;
  vector[N_ae] log_transfer_constant_z;
  array[N_experiment] vector[N_drain] drain_z;
  array[N_experiment] vector[N_enzyme] log_conc_enzyme_z;
  array[N_experiment] vector[N_phosphorylation_enzymes] log_conc_phos_z;
  array[N_experiment] vector[N_unbalanced] log_conc_unbalanced_z;
}
transformed parameters {
  // rescale
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
  vector[N_edge] keq = get_keq(S, dgf, mic_to_met, water_stoichiometry);
  for (e in 1:N_experiment){
    flux[e] = rep_vector(0, N_reaction);
    real timepoints[2] = {timepoint, timepoint + 10};
    vector[N_enzyme] conc_enzyme_experiment = conc_enzyme[e] .* knockout[e]';
    vector[N_phosphorylation_enzymes] conc_phos_experiment = conc_phos[e] .* phos_knockout[e]';
    vector[N_mic-N_unbalanced] conc_balanced[2] =
      ode_adjoint_tol_ctl(dbalanced_dt,
                  conc_init[e, balanced_mic_ix],
                  initial_time,
                  timepoints,
                  rel_tol_forward, 
                  abs_tol_forward,
                  rel_tol_backward, 
                  abs_tol_backward, 
                  rel_tol_quadrature,
                  abs_tol_quadrature,
                  max_num_steps,
                  num_steps_between_checkpoints,
                  interpolation_polynomial,
                  solver_forward,
                  solver_backward,
                  conc_unbalanced[e],
                  balanced_mic_ix,
                  unbalanced_mic_ix,
                  conc_enzyme_experiment,
                  km,
                  drain[e],
                  km_lookup,
                  S,
                  edge_type,
                  edge_to_drain,
                  edge_to_enzyme,
                  kcat,
                  keq,
                  ix_ci,
                  ix_ai,
                  ix_aa,
                  ix_pa,
                  ix_pi,
                  n_ci,
                  n_ai,
                  n_aa,
                  n_pa,
                  n_pi,
                  ki,
                  diss_t,
                  diss_r,
                  transfer_constant,
                  subunits,
                  kcat_phos,
                  conc_phos_experiment);
    conc[e, balanced_mic_ix] = conc_balanced[1];
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e];
    {
    vector[N_edge] flux_edge = get_flux(conc[e],
                                        conc_enzyme_experiment,
                                        km,
                                        drain[e],
                                        km_lookup,
                                        S,
                                        edge_type,
                                        edge_to_drain,
                                        edge_to_enzyme,
                                        kcat,
                                        keq,
                                        ix_ci,
                                        ix_ai,
                                        ix_aa,
                                        ix_pa,
                                        ix_pi,
                                        n_ci,
                                        n_ai,
                                        n_aa,
                                        n_pa,
                                        n_pi,
                                        ki,
                                        diss_t,
                                        diss_r,
                                        transfer_constant,
                                        subunits,
                                        kcat_phos,
                                        conc_phos_experiment);
    for (j in 1:N_edge)
      flux[e, edge_to_reaction[j]] += flux_edge[j];
    }
    if (reject_non_steady == 1 && check_steady_state(conc_balanced,
                                                     e,
                                                     flux[e],
                                                     conc_init[e],
                                                     timepoints,
                                                     conc_unbalanced[e],
                                                     conc_enzyme_experiment,
                                                     km,
                                                     drain[e],
                                                     kcat,
                                                     keq,
                                                     ki,
                                                     diss_t,
                                                     diss_r,
                                                     transfer_constant,
                                                     kcat_phos,
                                                     conc_phos_experiment) == 0) {
      reject("Non-steady state in experiment ", e);
    }
  }
}
model {
  log_kcat_z ~ std_normal();
  log_km_z ~ std_normal();
  log_ki_z ~ std_normal();
  log_diss_t_z ~ std_normal();
  log_diss_r_z ~ std_normal();
  log_transfer_constant_z ~ std_normal();
  dgf ~ multi_normal(prior_loc_dgf, prior_cov_dgf);
  log_kcat_phos_z ~ std_normal();
  for (ex in 1:N_experiment){
    log_conc_unbalanced_z[ex] ~ std_normal();
    log_conc_enzyme_z[ex] ~ std_normal();
    log_conc_phos_z[ex] ~ std_normal();
    drain_z[ex] ~ std_normal();
  }
  if (LIKELIHOOD == 1){
    for (c in 1:N_conc_measurement)
      yconc[c] ~ lognormal(log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
    for (e in 1:N_enzyme_measurement)
      yenz[e] ~ lognormal(log(conc_enzyme[experiment_yenz[e], enzyme_yenz[e]]), sigma_enz[e]);
    for (f in 1:N_flux_measurement)
      yflux[f] ~ normal(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
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
