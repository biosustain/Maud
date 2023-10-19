#include functions.stan
data {
  // network properties
  int<lower=1> N_mic;
  int<lower=1> N_edge_sub;
  int<lower=1> N_edge_prod;
  int<lower=1> N_unbalanced;
  int<lower=1> N_metabolite;
  int<lower=1> N_km;
  int<lower=1> N_sub_km;
  int<lower=1> N_prod_km;
  int<lower=1> N_reaction;
  int<lower=1> N_enzyme;
  int<lower=0> N_drain;
  int<lower=1> N_edge;
  int<lower=0> N_allostery;
  int<lower=0> N_allosteric_enzyme;
  int<lower=0> N_phosphorylation;
  int<lower=0> N_pme; // phosphorylation modifying enzyme
  int<lower=0> N_competitive_inhibition;
  int<lower=0> N_dgf_fixed;
  matrix[N_mic, N_edge] S;
  array[N_mic - N_unbalanced] int<lower=1, upper=N_mic> balanced_mic_ix;
  array[N_unbalanced] int<lower=1, upper=N_mic> unbalanced_mic_ix;
  array[N_competitive_inhibition] int<lower=1, upper=N_mic> ci_mic_ix;
  array[N_edge] int<lower=1, upper=3> edge_type; // 1 = reversible modular rate law, 2 = drain
  array[N_edge] int<lower=0, upper=N_enzyme> edge_to_enzyme; // 0 if drain
  array[N_edge] int<lower=0, upper=N_allostery> edge_to_tc; // 0 if non-allosteric
  array[N_edge] int<lower=0, upper=N_drain> edge_to_drain; // 0 if enzyme
  array[N_edge] int<lower=0, upper=N_reaction> edge_to_reaction;
  array[N_allostery] int<lower=1, upper=2> allostery_type;
  array[N_allostery] int<lower=1, upper=N_mic> allostery_mic;
  array[N_phosphorylation] int<lower=1, upper=2> phosphorylation_type;
  array[N_phosphorylation] int<lower=1, upper=N_pme> phosphorylation_pme;
  array[N_edge_sub] int sub_by_edge_long;
  array[N_edge, 2] int sub_by_edge_bounds;
  array[N_edge_prod] int prod_by_edge_long;
  array[N_edge, 2] int prod_by_edge_bounds;
  array[N_sub_km] int sub_km_ix_by_edge_long;
  array[N_edge, 2] int sub_km_ix_by_edge_bounds;
  array[N_prod_km] int prod_km_ix_by_edge_long;
  array[N_edge, 2] int prod_km_ix_by_edge_bounds;
  array[N_competitive_inhibition] int ci_ix_long;
  array[N_edge, 2] int ci_ix_bounds;
  array[N_allostery] int allostery_ix_long;
  array[N_edge, 2] int allostery_ix_bounds;
  array[N_phosphorylation] int phosphorylation_ix_long;
  array[N_edge, 2] int phosphorylation_ix_bounds;
  array[N_mic] int<lower=1, upper=N_metabolite> mic_to_met;
  vector[N_edge] water_stoichiometry;
  vector[N_edge] transported_charge;
  vector<lower=1>[N_enzyme] subunits;
  vector[N_dgf_fixed] dgf_fixed;
  array[N_dgf_fixed] int ix_dgf_fixed;
  array[N_metabolite-N_dgf_fixed] int ix_dgf_free;
  // experiment properties
  int<lower=1> N_experiment_train;
  int<lower=0> N_flux_measurement_train;
  int<lower=0> N_enzyme_measurement_train;
  int<lower=0> N_conc_measurement_train;
  int<lower=0> N_enzyme_knockout_train;
  int<lower=0> N_pme_knockout_train;
  array[N_enzyme_knockout_train] int<lower=0, upper=N_enzyme> enzyme_knockout_train_long;
  array[N_experiment_train, 2] int enzyme_knockout_train_bounds;
  array[N_pme_knockout_train] int<lower=0, upper=N_pme> pme_knockout_train_long;
  array[N_experiment_train, 2] int pme_knockout_train_bounds;
  vector[N_experiment_train] temperature_train;
  array[N_conc_measurement_train] int<lower=1, upper=N_experiment_train> experiment_yconc_train;
  array[N_conc_measurement_train] int<lower=1, upper=N_mic> mic_ix_yconc_train;
  array[N_conc_measurement_train] real yconc_train;
  vector<lower=0>[N_conc_measurement_train] sigma_yconc_train;
  array[N_flux_measurement_train] int<lower=1, upper=N_experiment_train> experiment_yflux_train;
  array[N_flux_measurement_train] int<lower=1, upper=N_reaction> reaction_yflux_train;
  array[N_flux_measurement_train] real yflux_train;
  vector<lower=0>[N_flux_measurement_train] sigma_yflux_train;
  array[N_enzyme_measurement_train] int<lower=0, upper=N_experiment_train> experiment_yenz_train;
  array[N_enzyme_measurement_train] int<lower=0, upper=N_enzyme> enzyme_yenz_train;
  array[N_enzyme_measurement_train] real yenz_train;
  vector<lower=0>[N_enzyme_measurement_train] sigma_yenz_train;
  // hardcoded priors (system)
  vector[N_metabolite] prior_loc_dgf;
  cov_matrix[N_metabolite] prior_cov_dgf;
  array[2] vector[N_enzyme] priors_kcat;
  array[2] vector[N_km] priors_km;
  array[2] vector[N_competitive_inhibition] priors_ki;
  array[2] vector[N_allostery] priors_dissociation_constant;
  array[2] vector[N_allosteric_enzyme] priors_transfer_constant;
  array[2] vector[N_pme] priors_kcat_pme;
  // hardcoded priors (boundary conditions)
  array[2] vector[N_experiment_train] priors_psi_train;
  array[2, N_experiment_train] vector[N_pme] priors_conc_pme_train;
  array[2, N_experiment_train] vector[N_unbalanced] priors_conc_unbalanced_train;
  array[2, N_experiment_train] vector[N_enzyme] priors_conc_enzyme_train;
  array[2, N_experiment_train] vector[N_drain] priors_drain_train;
  // configuration
  array[N_experiment_train] vector<lower=0>[N_mic - N_unbalanced] conc_init;
  real rel_tol_ode;
  real abs_tol_ode;
  int max_num_steps_ode;
  real rel_tol_alg;
  real abs_tol_alg;
  int max_num_steps_alg;
  real steady_state_threshold_abs;
  real steady_state_threshold_rel;
  real steady_state_penalty_rel;
  int<lower=0, upper=1> likelihood; // set to 0 for priors-only mode
  real drain_small_conc_corrector;
  int<lower=0, upper=1> penalize_non_steady;
}
transformed data {
  real initial_time = 0;
  complex complex_step = 1e-14i;
  matrix[N_metabolite-N_dgf_fixed, N_metabolite-N_dgf_fixed] prior_cov_dgf_free;
  for (mi in 1:N_metabolite-N_dgf_fixed){
    for (mj in 1:N_metabolite-N_dgf_fixed)
      prior_cov_dgf_free[mi,mj] = prior_cov_dgf[ix_dgf_free[mi],ix_dgf_free[mj]];
  }
  matrix[N_metabolite-N_dgf_fixed, N_metabolite-N_dgf_fixed]
    prior_cov_dgf_free_chol = cholesky_decompose(prior_cov_dgf_free);
  int N_balanced = N_mic - N_unbalanced;
}
parameters {
  vector[N_metabolite-N_dgf_fixed] dgf_free;
  vector[N_enzyme] log_kcat_z;
  vector[N_km] log_km_z;
  vector[N_pme] log_kcat_pme_z;
  vector[N_competitive_inhibition] log_ki_z;
  vector[N_allostery] log_dissociation_constant_z;
  vector[N_allosteric_enzyme] log_transfer_constant_z;
  vector[N_experiment_train] psi_train_z;
  array[N_experiment_train] vector[N_drain] drain_train_z;
  array[N_experiment_train] vector[N_enzyme] log_conc_enzyme_train_z;
  array[N_experiment_train] vector[N_pme] log_conc_pme_train_z;
  array[N_experiment_train] vector[N_unbalanced] log_conc_unbalanced_train_z;
}
transformed parameters {
  // unfix
  vector[N_metabolite] dgf = unfix(dgf_free, dgf_fixed, ix_dgf_free, ix_dgf_fixed);
  // rescale
  vector[N_km] km = unz_log_1d(priors_km, log_km_z);
  vector[N_competitive_inhibition] ki = unz_log_1d(priors_ki, log_ki_z);
  vector[N_enzyme] kcat = unz_log_1d(priors_kcat, log_kcat_z);
  vector[N_allostery] dissociation_constant = unz_log_1d(priors_dissociation_constant,
                                                         log_dissociation_constant_z);
  vector[N_allosteric_enzyme] transfer_constant = unz_log_1d(priors_transfer_constant,
                                                             log_transfer_constant_z);
  vector[N_pme] kcat_pme = unz_log_1d(priors_kcat_pme, log_kcat_pme_z);
  vector[N_experiment_train] psi_train = unz_1d(priors_psi_train,
                                                psi_train_z);
  array[N_experiment_train] vector[N_drain] drain_train = unz_2d(priors_drain_train,
                                                                 drain_train_z);
  array[N_experiment_train] vector[N_enzyme] conc_enzyme_train = unz_log_2d(priors_conc_enzyme_train,
                                                                    log_conc_enzyme_train_z);
  array[N_experiment_train] vector[N_unbalanced] conc_unbalanced_train = unz_log_2d(priors_conc_unbalanced_train,
                                                                    log_conc_unbalanced_train_z);
  array[N_experiment_train] vector[N_pme] conc_pme_train = unz_log_2d(priors_conc_pme_train,
                                                                    log_conc_pme_train_z);
  // transform
  array[N_experiment_train] vector<lower=0>[N_mic] conc_train;
  array[N_experiment_train] vector[N_reaction] flux_train;
  array[N_experiment_train] vector[N_edge] dgr_train;
  matrix[N_experiment_train, N_mic - N_unbalanced] steady_dev;
  for (e in 1 : N_experiment_train) {
    dgr_train[e] = get_dgr(S, dgf, temperature_train[e], mic_to_met,
                           water_stoichiometry, transported_charge,
                           psi_train[e]);
    flux_train[e] = rep_vector(0, N_reaction);
    vector[N_enzyme] conc_enzyme_experiment = conc_enzyme_train[e];
    vector[N_pme] conc_pme_experiment = conc_pme_train[e];
    vector[N_mic - N_unbalanced] conc_balanced_experiment;
    int N_eko_experiment = measure_ragged(enzyme_knockout_train_bounds, e);
    int N_pko_experiment = measure_ragged(pme_knockout_train_bounds, e);
    if (N_eko_experiment > 0) {
      array[N_eko_experiment] int eko_experiment = extract_ragged(enzyme_knockout_train_long,
                                                                  enzyme_knockout_train_bounds,
                                                                  e);
      conc_enzyme_experiment[eko_experiment] = rep_vector(0,
                                                          N_eko_experiment);
    }
    if (N_pko_experiment > 0) {
      array[N_pko_experiment] int pko_experiment = extract_ragged(pme_knockout_train_long,
                                                                  pme_knockout_train_bounds,
                                                                  e);
      conc_pme_experiment[pko_experiment] = rep_vector(0, N_pko_experiment);
    }
    conc_balanced_experiment = solve_newton_tol(maud_ae_system,
                                               conc_init[e],
                                               rel_tol_alg, abs_tol_alg,
                                               max_num_steps_alg,
                                               rel_tol_ode, abs_tol_ode,
                                               max_num_steps_ode,
                                               conc_unbalanced_train[e],
                                               balanced_mic_ix,
                                               unbalanced_mic_ix,
                                               conc_enzyme_experiment,
                                               dgr_train[e], kcat, km, ki,
                                               transfer_constant,
                                               dissociation_constant, kcat_pme,
                                               conc_pme_experiment,
                                               drain_train[e],
                                               temperature_train[e],
                                               drain_small_conc_corrector, S,
                                               subunits, edge_type,
                                               edge_to_enzyme, edge_to_drain,
                                               ci_mic_ix, sub_km_ix_by_edge_long,
                                               sub_km_ix_by_edge_bounds,
                                               prod_km_ix_by_edge_long,
                                               prod_km_ix_by_edge_bounds,
                                               sub_by_edge_long,
                                               sub_by_edge_bounds,
                                               prod_by_edge_long,
                                               prod_by_edge_bounds, ci_ix_long,
                                               ci_ix_bounds, allostery_ix_long,
                                               allostery_ix_bounds,
                                               allostery_type, allostery_mic,
                                               edge_to_tc,
                                               phosphorylation_ix_long,
                                               phosphorylation_ix_bounds,
                                               phosphorylation_type,
                                               phosphorylation_pme);
    conc_train[e, balanced_mic_ix] = conc_balanced_experiment;
    conc_train[e, unbalanced_mic_ix] = conc_unbalanced_train[e];
    {
      vector[N_edge] edge_flux = get_edge_flux(conc_train[e],
                                               conc_enzyme_experiment,
                                               dgr_train[e], kcat, km, ki,
                                               transfer_constant,
                                               dissociation_constant,
                                               kcat_pme, conc_pme_experiment,
                                               drain_train[e],
                                               temperature_train[e],
                                               drain_small_conc_corrector, S,
                                               subunits, edge_type,
                                               edge_to_enzyme, edge_to_drain,
                                               ci_mic_ix,
                                               sub_km_ix_by_edge_long,
                                               sub_km_ix_by_edge_bounds,
                                               prod_km_ix_by_edge_long,
                                               prod_km_ix_by_edge_bounds,
                                               sub_by_edge_long,
                                               sub_by_edge_bounds,
                                               prod_by_edge_long,
                                               prod_by_edge_bounds,
                                               ci_ix_long, ci_ix_bounds,
                                               allostery_ix_long,
                                               allostery_ix_bounds,
                                               allostery_type, allostery_mic,
                                               edge_to_tc,
                                               phosphorylation_ix_long,
                                               phosphorylation_ix_bounds,
                                               phosphorylation_type,
                                               phosphorylation_pme);
      steady_dev[e] = (S * edge_flux)[balanced_mic_ix]';
      for (j in 1 : N_edge) {
        flux_train[e, edge_to_reaction[j]] += edge_flux[j];
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
  dgf_free ~ multi_normal_cholesky(prior_loc_dgf[ix_dgf_free],
                                   prior_cov_dgf_free_chol);
  log_kcat_pme_z ~ std_normal();
  for (ex in 1 : N_experiment_train) {
    log_conc_unbalanced_train_z[ex] ~ std_normal();
    log_conc_enzyme_train_z[ex] ~ std_normal();
    log_conc_pme_train_z[ex] ~ std_normal();
    drain_train_z[ex] ~ std_normal();
    psi_train_z[ex] ~ std_normal();
  }
  if (likelihood == 1) {
    for (c in 1 : N_conc_measurement_train) {
      yconc_train[c] ~ lognormal(log(conc_train[experiment_yconc_train[c], mic_ix_yconc_train[c]]),
                                 sigma_yconc_train[c]);
    }
    for (e in 1 : N_enzyme_measurement_train) {
      yenz_train[e] ~ lognormal(log(conc_enzyme_train[experiment_yenz_train[e], enzyme_yenz_train[e]]),
                                sigma_yenz_train[e]);
    }
    for (f in 1 : N_flux_measurement_train) {
      yflux_train[f] ~ normal(flux_train[experiment_yflux_train[f], reaction_yflux_train[f]],
                              sigma_yflux_train[f]);
    }
    if (penalize_non_steady == 1) {
      for (xpt in 1 : N_experiment_train) {
        steady_dev[xpt] ~ normal(0.0,
                                 conc_train[xpt, balanced_mic_ix]
                                 * steady_state_penalty_rel);
      }
    }
  }
}
generated quantities {
  vector[N_conc_measurement_train] yrep_conc_train;
  vector[N_flux_measurement_train] yrep_flux_train;
  vector[N_conc_measurement_train] llik_conc_train;
  vector[N_flux_measurement_train] llik_flux_train;
  array[N_experiment_train] matrix[N_mic - N_unbalanced, N_edge] concentration_control_matrix;
  array[N_experiment_train] matrix[N_edge, N_edge] flux_control_matrix;
  array[N_experiment_train] matrix[N_edge, N_balanced] elasticity;
  array[N_experiment_train] matrix[N_edge, N_enzyme] sensitivity;
  array[N_experiment_train] matrix[N_edge, N_enzyme] flux_response_coefficient;
  array[N_experiment_train] matrix[N_balanced, N_enzyme] concentration_response_coefficient;
  array[N_experiment_train] vector[N_edge] free_enzyme_ratio_train;
  array[N_experiment_train] vector[N_edge] saturation_train;
  array[N_experiment_train] vector[N_edge] allostery_train;
  array[N_experiment_train] vector[N_edge] phosphorylation_train;
  array[N_experiment_train] vector[N_edge] reversibility_train;
  array[N_experiment_train] vector[N_edge] keq_train;
  // Simulating measurements from the posterior predictive distribution
  for (c in 1 : N_conc_measurement_train) {
    yrep_conc_train[c] = lognormal_rng(log(conc_train[experiment_yconc_train[c], mic_ix_yconc_train[c]]),
                                       sigma_yconc_train[c]);
    llik_conc_train[c] = lognormal_lpdf(yconc_train[c] | log(conc_train[experiment_yconc_train[c], mic_ix_yconc_train[c]]), sigma_yconc_train[c]);
  }
  for (f in 1 : N_flux_measurement_train) {
    yrep_flux_train[f] = normal_rng(flux_train[experiment_yflux_train[f], reaction_yflux_train[f]],
                                    sigma_yflux_train[f]);
    llik_flux_train[f] = normal_lpdf(yflux_train[f] | flux_train[experiment_yflux_train[f], reaction_yflux_train[f]], sigma_yflux_train[f]);
  }
  // Calculating regulatory decomposition
  for (e in 1 : N_experiment_train) {
    keq_train[e] = get_keq(S, dgf, temperature_train[e], mic_to_met,
                           water_stoichiometry, transported_charge,
                           psi_train[e]);
    free_enzyme_ratio_train[e] = get_free_enzyme_ratio(conc_train[e], S, km,
                                                       ki, edge_type,
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
    saturation_train[e] = get_saturation(conc_train[e], km,
                                         free_enzyme_ratio_train[e],
                                         sub_km_ix_by_edge_long,
                                         sub_km_ix_by_edge_bounds,
                                         sub_by_edge_long,
                                         sub_by_edge_bounds, edge_type);
    allostery_train[e] = get_allostery(conc_train[e],
                                       free_enzyme_ratio_train[e],
                                       transfer_constant,
                                       dissociation_constant, subunits,
                                       allostery_ix_long,
                                       allostery_ix_bounds, allostery_type,
                                       allostery_mic, edge_to_tc);
    phosphorylation_train[e] = get_phosphorylation(kcat_pme,
                                                   conc_pme_train[e],
                                                   phosphorylation_ix_long,
                                                   phosphorylation_ix_bounds,
                                                   phosphorylation_type,
                                                   phosphorylation_pme,
                                                   subunits);
    reversibility_train[e] = get_reversibility(dgr_train[e],
                                               temperature_train[e], S,
                                               conc_train[e], edge_type);
  }
  // Calculating control coefficients
  for (e in 1 : N_experiment_train) {
    vector[N_pme] conc_pme_experiment = conc_pme_train[e];
    for (bal_met in 1 : N_balanced) {
      int met_idx = balanced_mic_ix[bal_met];
      complex_vector[N_mic] complex_step_conc = conc_train[e];
      complex_step_conc[met_idx] = complex_step_conc[met_idx] + complex_step;
      complex_vector[N_edge] edge_flux_met = get_complex_edge_flux_metabolite(complex_step_conc,
                                                                    conc_enzyme_train[e],
                                                                    dgr_train[e],
                                                                    kcat, km,
                                                                    ki,
                                                                    transfer_constant,
                                                                    dissociation_constant,
                                                                    kcat_pme,
                                                                    conc_pme_experiment,
                                                                    drain_train[e],
                                                                    temperature_train[e],
                                                                    drain_small_conc_corrector,
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
                                                                    phosphorylation_type,
                                                                    phosphorylation_pme);
      elasticity[e][ : , bal_met] = get_flux_jacobian(edge_flux_met,
                                                      complex_step);
    }
    for (enz in 1 : N_enzyme) {
      complex_vector[N_enzyme] numerical_step_enzyme = conc_enzyme_train[e];
      numerical_step_enzyme[enz] = conc_enzyme_train[e][enz] + complex_step;
      complex_vector[N_edge] edge_flux_enz = get_complex_edge_flux_enzyme(conc_train[e],
                                                                    numerical_step_enzyme,
                                                                    dgr_train[e],
                                                                    kcat, km,
                                                                    ki,
                                                                    transfer_constant,
                                                                    dissociation_constant,
                                                                    kcat_pme,
                                                                    conc_pme_experiment,
                                                                    drain_train[e],
                                                                    temperature_train[e],
                                                                    drain_small_conc_corrector,
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
                                                                    phosphorylation_type,
                                                                    phosphorylation_pme);
      sensitivity[e][ : , enz] = get_flux_jacobian(edge_flux_enz,
                                                   complex_step);
    }
    concentration_control_matrix[e] = get_concentration_control_matrix(S,
                                                                    elasticity[e],
                                                                    balanced_mic_ix,
                                                                    N_edge,
                                                                    N_balanced);
    flux_control_matrix[e] = get_flux_control_matrix(S, elasticity[e],
                                                     balanced_mic_ix, N_edge,
                                                     N_balanced);
    flux_response_coefficient[e] = flux_control_matrix[e] * sensitivity[e];
    concentration_response_coefficient[e] = concentration_control_matrix[e]
                                            * sensitivity[e];
  }
}
