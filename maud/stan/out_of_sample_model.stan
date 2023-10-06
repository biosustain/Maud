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
  // experiment properties
  int<lower=1> N_experiment_test;
  int<lower=0> N_enzyme_knockout_test;
  int<lower=0> N_pme_knockout_test;
  array[N_enzyme_knockout_test] int<lower=0, upper=N_enzyme> enzyme_knockout_test_long;
  array[N_experiment_test, 2] int enzyme_knockout_test_bounds;
  array[N_pme_knockout_test] int<lower=0, upper=N_pme> pme_knockout_test_long;
  array[N_experiment_test, 2] int pme_knockout_test_bounds;
  vector[N_experiment_test] temperature_test;
  // hardcoded priors
  array[2, N_experiment_test] vector[N_pme] priors_conc_phos_test;
  array[2, N_experiment_test] vector[N_unbalanced] priors_conc_unbalanced_test;
  array[2, N_experiment_test] vector[N_enzyme] priors_conc_enzyme_test;
  array[2, N_experiment_test] vector[N_drain] priors_drain_test;
  // configuration
  array[N_experiment_test] vector<lower=0>[N_mic - N_unbalanced] conc_init;
  real rel_tol;
  real abs_tol;
  int max_num_steps;
  int<lower=0, upper=1> likelihood; // set to 0 for priors-only mode
  real drain_small_conc_corrector;
  real<lower=0> timepoint;
}
transformed data {
  real initial_time = 0;
}
parameters {
  vector[N_km] km;
  vector[N_competitive_inhibition] ki;
  vector[N_enzyme] kcat;
  vector[N_allostery] dissociation_constant;
  vector[N_allosteric_enzyme] transfer_constant;
  vector[N_pme] kcat_pme;
  array[N_experiment_test] vector[N_edge] dgr_test;
}
generated quantities {
  array[N_experiment_test] vector<lower=0>[N_mic] conc_test;
  array[N_experiment_test] vector[N_reaction] flux_test;
  array[N_experiment_test] vector[N_pme] conc_pme_test;
  array[N_experiment_test] vector[N_unbalanced] conc_unbalanced_test;
  array[N_experiment_test] vector[N_enzyme] conc_enzyme_test;
  array[N_experiment_test] vector[N_drain] drain_test;
  array[N_experiment_test] vector[N_edge] free_enzyme_ratio_test;
  array[N_experiment_test] vector[N_edge] saturation_test;
  array[N_experiment_test] vector[N_edge] allostery_test;
  array[N_experiment_test] vector[N_edge] phosphorylation_test;
  array[N_experiment_test] vector[N_edge] reversibility_test;
  // Sampling experiment boundary conditions from priors
  for (e in 1 : N_experiment_test) {
    drain_test[e] = to_vector(normal_rng(priors_drain_test[1, e],
                                         priors_drain_test[2, e]));
    conc_pme_test[e] = to_vector(lognormal_rng(log(priors_conc_phos_test[1, e]),
                                               priors_conc_phos_test[2, e]));
    conc_unbalanced_test[e] = to_vector(lognormal_rng(log(priors_conc_unbalanced_test[1, e]),
                                                      priors_conc_unbalanced_test[2, e]));
    conc_enzyme_test[e] = to_vector(lognormal_rng(log(priors_conc_enzyme_test[1, e]),
                                                  priors_conc_enzyme_test[2, e]));
  }
  // Simulation of experiments
  for (e in 1 : N_experiment_test) {
    flux_test[e] = rep_vector(0, N_reaction);
    vector[N_enzyme] conc_enzyme_experiment = conc_enzyme_test[e];
    vector[N_pme] conc_pme_experiment = conc_pme_test[e];
    array[1] vector[N_mic - N_unbalanced] conc_balanced_experiment;
    int N_eko_experiment = measure_ragged(enzyme_knockout_test_bounds, e);
    int N_pko_experiment = measure_ragged(pme_knockout_test_bounds, e);
    if (N_eko_experiment > 0) {
      array[N_eko_experiment] int eko_experiment = extract_ragged(enzyme_knockout_test_long,
                                                                  enzyme_knockout_test_bounds,
                                                                  e);
      conc_enzyme_experiment[eko_experiment] = rep_vector(0,
                                                          N_eko_experiment);
    }
    if (N_pko_experiment > 0) {
      array[N_pko_experiment] int pko_experiment = extract_ragged(pme_knockout_test_long,
                                                                  pme_knockout_test_bounds,
                                                                  e);
      conc_pme_experiment[pko_experiment] = rep_vector(0, N_pko_experiment);
    }
    conc_balanced_experiment = ode_bdf_tol(dbalanced_dt,
                                           conc_init[e, balanced_mic_ix],
                                           initial_time, {timepoint},
                                           rel_tol, abs_tol, max_num_steps,
                                           conc_unbalanced_test[e],
                                           balanced_mic_ix,
                                           unbalanced_mic_ix,
                                           conc_enzyme_experiment,
                                           dgr_test[e], kcat, km, ki,
                                           transfer_constant,
                                           dissociation_constant, kcat_pme,
                                           conc_pme_experiment,
                                           drain_test[e],
                                           temperature_test[e],
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
    conc_test[e, balanced_mic_ix] = conc_balanced_experiment[1];
    conc_test[e, unbalanced_mic_ix] = conc_unbalanced_test[e];
    vector[N_edge] edge_flux = get_edge_flux(conc_test[e],
                                             conc_enzyme_experiment,
                                             dgr_test[e], kcat, km, ki,
                                             transfer_constant,
                                             dissociation_constant, kcat_pme,
                                             conc_pme_experiment,
                                             drain_test[e],
                                             temperature_test[e],
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
                                             prod_by_edge_bounds, ci_ix_long,
                                             ci_ix_bounds, allostery_ix_long,
                                             allostery_ix_bounds,
                                             allostery_type, allostery_mic,
                                             edge_to_tc,
                                             phosphorylation_ix_long,
                                             phosphorylation_ix_bounds,
                                             phosphorylation_type,
                                             phosphorylation_pme);
    for (j in 1 : N_edge) {
      flux_test[e, edge_to_reaction[j]] += edge_flux[j];
    }
  }
  for (e in 1 : N_experiment_test) {
    free_enzyme_ratio_test[e] = get_free_enzyme_ratio(conc_test[e], S, km,
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
    saturation_test[e] = get_saturation(conc_test[e], km,
                                        free_enzyme_ratio_test[e],
                                        sub_km_ix_by_edge_long,
                                        sub_km_ix_by_edge_bounds,
                                        sub_by_edge_long, sub_by_edge_bounds,
                                        edge_type);
    allostery_test[e] = get_allostery(conc_test[e],
                                      free_enzyme_ratio_test[e],
                                      transfer_constant,
                                      dissociation_constant, subunits,
                                      allostery_ix_long, allostery_ix_bounds,
                                      allostery_type, allostery_mic,
                                      edge_to_tc);
    phosphorylation_test[e] = get_phosphorylation(kcat_pme, conc_pme_test[e],
                                                  phosphorylation_ix_long,
                                                  phosphorylation_ix_bounds,
                                                  phosphorylation_type,
                                                  phosphorylation_pme,
                                                  subunits);
    reversibility_test[e] = get_reversibility(dgr_test[e],
                                              temperature_test[e], S,
                                              conc_test[e], edge_type);
  }
}
