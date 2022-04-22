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
  int<lower=0> N_allostery;
  int<lower=0> N_allosteric_enzyme;
  int<lower=0> N_phosphorylation;
  int<lower=0> N_competitive_inhibition;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=0> N_enzyme_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=0> N_enzyme_knockout;
  int<lower=0> N_phosphorylation_knockout;
  // measurements
  array[N_conc_measurement] int<lower=1,upper=N_experiment> experiment_yconc;
  array[N_conc_measurement] int<lower=1,upper=N_mic> mic_ix_yconc;
  array[N_conc_measurement] real yconc;
  vector<lower=0>[N_conc_measurement] sigma_conc;
  array[N_flux_measurement] int<lower=1,upper=N_experiment> experiment_yflux;
  array[N_flux_measurement] int<lower=1,upper=N_reaction> reaction_yflux;
  array[N_flux_measurement] real yflux;
  vector<lower=0>[N_flux_measurement] sigma_flux;
  array[N_enzyme_measurement] int<lower=0,upper=N_experiment> experiment_yenz;
  array[N_enzyme_measurement] int<lower=0,upper=N_enzyme> enzyme_yenz;
  array[N_enzyme_measurement] real yenz;
  vector<lower=0>[N_enzyme_measurement] sigma_enz;
  // hardcoded priors
  vector[N_metabolite] prior_loc_dgf;
  cov_matrix[N_metabolite] prior_cov_dgf;
  array[2] vector[N_enzyme] priors_kcat;
  array[2] vector[N_km] priors_km;
  array[2] vector[N_competitive_inhibition] priors_ki;
  array[2] vector[N_allostery] priors_dissociation_constant;
  array[2] vector[N_allosteric_enzyme] priors_transfer_constant;
  array[2] vector[N_phosphorylation] priors_kcat_phos;
  array[2, N_experiment] vector[N_phosphorylation] priors_conc_phos;
  array[2, N_experiment] vector[N_unbalanced] priors_conc_unbalanced;
  array[2, N_experiment] vector[N_enzyme] priors_conc_enzyme;
  array[2, N_experiment] vector[N_drain] priors_drain;
  // network properties
  matrix[N_mic, N_edge] S;
  array[N_mic-N_unbalanced] int<lower=1,upper=N_mic> balanced_mic_ix;
  array[N_unbalanced] int<lower=1,upper=N_mic> unbalanced_mic_ix;
  array[N_competitive_inhibition] int<lower=1,upper=N_mic> ci_mic_ix;
  array[N_edge] int<lower=1,upper=3> edge_type;  // 1 = reversible modular rate law, 2 = drain
  array[N_edge] int<lower=0,upper=N_enzyme> edge_to_enzyme;  // 0 if drain
  array[N_edge] int<lower=0,upper=N_allostery> edge_to_tc;  // 0 if non-allosteric
  array[N_edge] int<lower=0,upper=N_drain> edge_to_drain;  // 0 if enzyme
  array[N_edge] int<lower=0,upper=N_reaction> edge_to_reaction;
  array[N_allostery] int<lower=1,upper=2> allostery_type;
  array[N_allostery] int<lower=1,upper=N_mic> allostery_mic;
  array[N_phosphorylation] int<lower=1,upper=2> phosphorylation_type;
  array[N_edge_sub] int sub_by_edge_long;
  array[N_edge,2] int sub_by_edge_bounds;
  array[N_edge_prod] int prod_by_edge_long;
  array[N_edge, 2] int prod_by_edge_bounds;
  array[N_edge_sub] int sub_km_ix_by_edge_long;
  array[N_edge, 2] int sub_km_ix_by_edge_bounds;
  array[N_edge_prod] int prod_km_ix_by_edge_long;
  array[N_edge,2] int prod_km_ix_by_edge_bounds;
  array[N_competitive_inhibition] int ci_ix_long;
  array[N_edge, 2] int ci_ix_bounds;
  array[N_allostery] int allostery_ix_long;
  array[N_edge, 2] int allostery_ix_bounds;
  array[N_phosphorylation] int phosphorylation_ix_long;
  array[N_edge, 2] int phosphorylation_ix_bounds;
  array[N_mic] int<lower=1,upper=N_metabolite> mic_to_met;
  vector[N_edge] water_stoichiometry;
  array[N_enzyme_knockout] int<lower=0,upper=N_enzyme> enzyme_knockout_long;
  array[N_experiment, 2] int enzyme_knockout_bounds;
  array[N_phosphorylation_knockout] int<lower=0,upper=N_phosphorylation> phosphorylation_knockout_long;
  array[N_experiment, 2] int phosphorylation_knockout_bounds;
  vector<lower=1>[N_enzyme] subunits;
  // configuration
  array[N_experiment] vector<lower=0>[N_mic] conc_init;
  real rel_tol; 
  real abs_tol;
  real steady_state_threshold_abs;
  real steady_state_threshold_rel;
  int max_num_steps;
  int<lower=0,upper=1> likelihood;  // set to 0 for priors-only mode
  real<lower=0> timepoint;
  int<lower=0,upper=1> reject_non_steady;
}
transformed data {
  real initial_time = 0;
  matrix[N_metabolite, N_metabolite] prior_cov_dgf_chol = cholesky_decompose(prior_cov_dgf);
}
parameters {
  vector[N_metabolite] dgf;
  vector[N_enzyme] log_kcat_z;
  vector[N_km] log_km_z;
  vector[N_phosphorylation] log_kcat_phos_z;
  vector[N_competitive_inhibition] log_ki_z;
  vector[N_allostery] log_dissociation_constant_z;
  vector[N_allosteric_enzyme] log_transfer_constant_z;
  array[N_experiment] vector[N_drain] drain_z;
  array[N_experiment] vector[N_enzyme] log_conc_enzyme_z;
  array[N_experiment] vector[N_phosphorylation] log_conc_phos_z;
  array[N_experiment] vector[N_unbalanced] log_conc_unbalanced_z;
}
transformed parameters {
  // rescale
  vector[N_km] km = unz_log_1d(priors_km, log_km_z);
  vector[N_competitive_inhibition] ki = unz_log_1d(priors_ki, log_ki_z);
  vector[N_enzyme] kcat = unz_log_1d(priors_kcat, log_kcat_z);
  vector[N_allostery] dissociation_constant = unz_log_1d(priors_dissociation_constant, log_dissociation_constant_z);
  vector[N_allosteric_enzyme] transfer_constant = unz_log_1d(priors_transfer_constant, log_transfer_constant_z);
  vector[N_phosphorylation] kcat_phos = unz_log_1d(priors_kcat_phos, log_kcat_phos_z);
  array[N_experiment] vector[N_drain] drain = unz_2d(priors_drain, drain_z);
  array[N_experiment] vector[N_enzyme] conc_enzyme = unz_log_2d(priors_conc_enzyme, log_conc_enzyme_z);
  array[N_experiment] vector[N_unbalanced] conc_unbalanced = unz_log_2d(priors_conc_unbalanced, log_conc_unbalanced_z);
  array[N_experiment] vector[N_phosphorylation] conc_phos = unz_log_2d(priors_conc_phos, log_conc_phos_z);
  // transform
  array[N_experiment] vector<lower=0>[N_mic] conc;
  array[N_experiment] vector[N_reaction] flux;
  vector[N_edge] dgrs = get_dgrs(S, dgf, mic_to_met, water_stoichiometry);
  for (e in 1:N_experiment){
    flux[e] = rep_vector(0, N_reaction);
    vector[N_enzyme] conc_enzyme_experiment = conc_enzyme[e];
    vector[N_phosphorylation] conc_phos_experiment = conc_phos[e];
    array[1] vector[N_mic-N_unbalanced] conc_balanced;
    int N_eko_experiment = measure_ragged(enzyme_knockout_bounds, e);
    int N_pko_experiment = measure_ragged(phosphorylation_knockout_bounds, e);
    if (N_eko_experiment > 0){
      array[N_eko_experiment] int eko_experiment =
        extract_ragged(enzyme_knockout_long, enzyme_knockout_bounds, e);
      conc_enzyme_experiment[eko_experiment] = rep_vector(0, N_eko_experiment);
    }
    if (N_pko_experiment > 0){
      array[N_pko_experiment] int pko_experiment =
        extract_ragged(phosphorylation_knockout_long, phosphorylation_knockout_bounds, e);
      conc_phos_experiment[pko_experiment] = rep_vector(0, N_pko_experiment);
    }
    conc_balanced =
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
                  dissociation_constant,
                  kcat_phos,
                  conc_phos_experiment,
                  drain[e],
                  S,
                  subunits,
                  edge_type,
                  edge_to_enzyme,
                  edge_to_drain,
                  ci_mic_ix,
                  sub_km_ix_by_edge_long,
                  sub_km_ix_by_edge_bounds,
                  prod_km_ix_by_edge_long,
                  prod_km_ix_by_edge_bounds,
                  sub_by_edge_long,
                  sub_by_edge_bounds,
                  prod_by_edge_long,
                  prod_by_edge_bounds,
                  ci_ix_long,
                  ci_ix_bounds,
                  allostery_ix_long,
                  allostery_ix_bounds,
                  allostery_type,
                  allostery_mic,
                  edge_to_tc,
                  phosphorylation_ix_long,
                  phosphorylation_ix_bounds,
                  phosphorylation_type);
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
                                             dissociation_constant,
                                             kcat_phos,
                                             conc_phos_experiment,
                                             drain[e],
                                             S,
                                             subunits,
                                             edge_type,
                                             edge_to_enzyme,
                                             edge_to_drain,
                                             ci_mic_ix,
                                             sub_km_ix_by_edge_long,
                                             sub_km_ix_by_edge_bounds,
                                             prod_km_ix_by_edge_long,
                                             prod_km_ix_by_edge_bounds,
                                             sub_by_edge_long,
                                             sub_by_edge_bounds,
                                             prod_by_edge_long,
                                             prod_by_edge_bounds,
                                             ci_ix_long,
                                             ci_ix_bounds,
                                             allostery_ix_long,
                                             allostery_ix_bounds,
                                             allostery_type,
                                             allostery_mic,
                                             edge_to_tc,
                                             phosphorylation_ix_long,
                                             phosphorylation_ix_bounds,
                                             phosphorylation_type);
    for (j in 1:N_edge)
      flux[e, edge_to_reaction[j]] += edge_flux[j];
    if (reject_non_steady == 1 &&
        check_steady_state((S * edge_flux)[balanced_mic_ix], conc_balanced[1], steady_state_threshold_abs, steady_state_threshold_rel) == 0){
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
        print("dissociation_constant: ", dissociation_constant);
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
  log_dissociation_constant_z ~ std_normal();
  log_transfer_constant_z ~ std_normal();
  dgf ~ multi_normal_cholesky(prior_loc_dgf, prior_cov_dgf_chol);
  log_kcat_phos_z ~ std_normal();
  for (ex in 1:N_experiment){
    log_conc_unbalanced_z[ex] ~ std_normal();
    log_conc_enzyme_z[ex] ~ std_normal();
    log_conc_phos_z[ex] ~ std_normal();
    drain_z[ex] ~ std_normal();
  }
  if (likelihood == 1){
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
                                                 ci_mic_ix,
                                                 sub_km_ix_by_edge_long,
                                                 sub_km_ix_by_edge_bounds,
                                                 prod_km_ix_by_edge_long,
                                                 prod_km_ix_by_edge_bounds,
                                                 sub_by_edge_long,
                                                 sub_by_edge_bounds,
                                                 prod_by_edge_long,
                                                 prod_by_edge_bounds,
                                                 ci_ix_long,
                                                 ci_ix_bounds);
    saturation[e] = get_saturation(conc[e],
                                   km,
                                   free_enzyme_ratio[e],
                                   sub_km_ix_by_edge_long,
                                   sub_km_ix_by_edge_bounds,
                                   sub_by_edge_long,
                                   sub_by_edge_bounds,
                                   edge_type);
    allostery[e] = get_allostery(conc[e],
                                 free_enzyme_ratio[e],
                                 transfer_constant,
                                 dissociation_constant,
                                 subunits,
                                 allostery_ix_long,
                                 allostery_ix_bounds,
                                 allostery_type,
                                 allostery_mic,
                                 edge_to_tc);
    phosphorylation[e] = get_phosphorylation(kcat_phos,
                                             conc_phos[e],
                                             phosphorylation_ix_long,
                                             phosphorylation_ix_bounds,
                                             phosphorylation_type,
                                             subunits);
    reversibility[e] = get_reversibility(dgrs, S, conc[e], edge_type);
  }
}
