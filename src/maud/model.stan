#include functions.stan
data {
  // dimensions
  int<lower=1> N_mic;
  int<lower=1> N_edge_sub;
  int<lower=1> N_edge_prod;
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
  int<lower=0> N_ci;
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
  array[2] vector[N_ci] priors_ki;
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
  int<lower=0,upper=N_ae> edge_to_tc[N_edge];  // 0 if non-allosteric
  int<lower=0,upper=N_drain> edge_to_drain[N_edge];  // 0 if enzyme
  int<lower=0,upper=N_reaction> edge_to_reaction[N_edge];
  int<lower=0,upper=N_km> km_lookup[N_mic, N_edge];
  int<lower=0,upper=N_ci> ki_lookup[N_mic, N_edge];
  int<lower=0,upper=N_ai> dt_lookup[N_mic, N_edge];
  int<lower=0,upper=N_aa> dr_lookup[N_mic, N_edge];
  array[N_edge_sub] int sub_by_edge_long;
  array[N_edge,2] int sub_by_edge_bounds;
  array[N_edge_prod] int prod_by_edge_long;
  array[N_edge, 2] int prod_by_edge_bounds;
  array[N_ci] int ci_by_edge_long;
  array[N_edge, 2] int ci_by_edge_bounds;
  array[N_ai] int ai_ix_long;
  array[N_edge, 2] int ai_ix_bounds;
  array[N_aa] int aa_ix_long;
  array[N_edge, 2] int aa_ix_bounds;
  array[N_pa] int pa_ix_long;
  array[N_edge, 2] int pa_ix_bounds;
  array[N_pi] int pi_ix_long;
  array[N_edge, 2] int pi_ix_bounds;
  int<lower=1,upper=N_metabolite> mic_to_met[N_mic];
  vector[N_edge] water_stoichiometry;
  matrix<lower=0,upper=1>[N_experiment, N_enzyme] is_knockout;
  matrix<lower=0,upper=1>[N_experiment, N_phosphorylation_enzymes] is_phos_knockout;
  vector<lower=1>[N_enzyme] subunits;
  // configuration
  vector<lower=0>[N_mic] conc_init[N_experiment];
  real rel_tol; 
  real abs_tol;
  int max_num_steps;
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
  vector[N_ci] log_ki_z;
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
  vector[N_ci] ki = unz_log_1d(priors_ki, log_ki_z);
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
  vector[N_edge] dgrs = get_dgrs(S, dgf, mic_to_met, water_stoichiometry);
  for (e in 1:N_experiment){
    flux[e] = rep_vector(0, N_reaction);
    vector[N_enzyme] conc_enzyme_experiment = conc_enzyme[e] .* knockout[e]';
    vector[N_phosphorylation_enzymes] conc_phos_experiment = conc_phos[e] .* phos_knockout[e]';
    vector[N_mic-N_unbalanced] conc_balanced[1] =
      ode_bdf_tol(dbalanced_dt,
                  conc_init[e, balanced_mic_ix],
                  initial_time,
                  {timepoint},
                  rel_tol, 
                  abs_tol,
                  max_num_steps,
                  conc_unbalanced[e],
                  balanced_mic_ix,
                  unbalanced_mic_ix,
                  conc_enzyme_experiment,
                  dgrs,
                  kcat,
                  km,
                  ki,
                  transfer_constant,
                  diss_t,
                  diss_r,
                  kcat_phos,
                  conc_phos_experiment,
                  drain[e],
                  S,
                  subunits,
                  edge_type,
                  edge_to_enzyme,
                  edge_to_tc,
                  edge_to_drain,
                  km_lookup,
                  ki_lookup,
                  dt_lookup,
                  dr_lookup,
                  sub_by_edge_long,
                  sub_by_edge_bounds,
                  prod_by_edge_long,
                  prod_by_edge_bounds,
                  ci_by_edge_long,
                  ci_by_edge_bounds,
                  ai_ix_long,
                  ai_ix_bounds,
                  aa_ix_long,
                  aa_ix_bounds,
                  pa_ix_long,
                  pa_ix_bounds,
                  pi_ix_long,
                  pi_ix_bounds);
    conc[e, balanced_mic_ix] = conc_balanced[1];
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e];
    {
    vector[N_edge] edge_flux = get_edge_flux(conc[e],
                                             conc_enzyme_experiment,
                                             dgrs,
                                             kcat,
                                             km,
                                             ki,
                                             transfer_constant,
                                             diss_t,
                                             diss_r,
                                             kcat_phos,
                                             conc_phos_experiment,
                                             drain[e],
                                             S,
                                             subunits,
                                             edge_type,
                                             edge_to_enzyme,
                                             edge_to_tc,
                                             edge_to_drain,
                                             km_lookup,
                                             ki_lookup,
                                             dt_lookup,
                                             dr_lookup,
                                             sub_by_edge_long,
                                             sub_by_edge_bounds,
                                             prod_by_edge_long,
                                             prod_by_edge_bounds,
                                             ci_by_edge_long,
                                             ci_by_edge_bounds,
                                             ai_ix_long,
                                             ai_ix_bounds,
                                             aa_ix_long,
                                             aa_ix_bounds,
                                             pa_ix_long,
                                             pa_ix_bounds,
                                             pi_ix_long,
                                             pi_ix_bounds);
    for (j in 1:N_edge)
      flux[e, edge_to_reaction[j]] += edge_flux[j];
    if (reject_non_steady == 1 &&
        check_steady_state((S * edge_flux)[balanced_mic_ix], conc_balanced[1]) == 0){
        print("Non-steady state in experiment ", e);
        print("Balanced metabolite concentration", conc_balanced[1]);
        print("flux: ", flux);
        print("conc_init: ", conc_init);
        print("conc_unbalanced: ", conc_unbalanced[e]);
        print("conc_enzyme_experiment: ", conc_enzyme_experiment);
        print("km: ", km);
        print("drain: ", drain[e]);
        print("kcat: ", kcat);
        print("dgrs: ", dgrs);
        print("ki: ", ki);
        print("diss_t: ", diss_t);
        print("diss_r: ", diss_r);
        print("transfer_constant: ", transfer_constant);
        print("kcat_phos: ", kcat_phos);
        print("conc_phos_experiment: ", conc_phos_experiment);
        reject("Rejecting");
      }
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
  array[N_experiment] vector[N_edge] free_enzyme_ratio;
  array[N_experiment] vector[N_edge] saturation;
  array[N_experiment] vector[N_edge] allostery;
  array[N_experiment] vector[N_edge] phosphorylation;
  array[N_experiment] vector[N_edge] reversibility;
  vector[N_edge] keq = get_keq(S, dgf, mic_to_met, water_stoichiometry);
  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
    log_lik_conc[c] = lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], mic_ix_yconc[c]]), sigma_conc[c]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
    log_lik_flux[f] = normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
  for (e in 1:N_experiment){
    free_enzyme_ratio[e] = get_free_enzyme_ratio(conc[e],
                                                 S,
                                                 km,
                                                 ki,
                                                 edge_type,
                                                 km_lookup,
                                                 ki_lookup,
                                                 sub_by_edge_long,
                                                 sub_by_edge_bounds,
                                                 prod_by_edge_long,
                                                 prod_by_edge_bounds,
                                                 ci_by_edge_long,
                                                 ci_by_edge_bounds);
    saturation[e] = get_saturation(conc[e],
                                   km,
                                   free_enzyme_ratio[e],
                                   km_lookup,
                                   sub_by_edge_long,
                                   sub_by_edge_bounds,
                                   edge_type);
    allostery[e] = get_allostery(conc[e],
                                 free_enzyme_ratio[e],
                                 transfer_constant,
                                 diss_t,
                                 diss_r,
                                 subunits,
                                 dt_lookup,
                                 dr_lookup,
                                 edge_to_tc,
                                 ai_ix_long,
                                 ai_ix_bounds,
                                 aa_ix_long,
                                 aa_ix_bounds);
    phosphorylation[e] = get_phosphorylation(kcat_phos,
                                             conc_phos[e] .* phos_knockout[e]',
                                             pa_ix_long,
                                             pa_ix_bounds,
                                             pi_ix_long,
                                             pi_ix_bounds,
                                             subunits);

    reversibility[e] = get_reversibility(dgrs, S, conc[e], edge_type);
  }
}
