#include functions.stan
data {
  // dimensions
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
  int<lower=0> N_pme;  // phosphorylation modifying enzyme
  int<lower=0> N_competitive_inhibition;
  int<lower=1> N_experiment;
  int<lower=0> N_enzyme_knockout;
  int<lower=0> N_pme_knockout;
  // hardcoded priors
  array[2, N_experiment] vector[N_pme] priors_conc_phos;
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
  array[N_phosphorylation] int<lower=1,upper=N_pme> phosphorylation_pme;
  array[N_edge_sub] int sub_by_edge_long;
  array[N_edge,2] int sub_by_edge_bounds;
  array[N_edge_prod] int prod_by_edge_long;
  array[N_edge, 2] int prod_by_edge_bounds;
  array[N_sub_km] int sub_km_ix_by_edge_long;
  array[N_edge, 2] int sub_km_ix_by_edge_bounds;
  array[N_prod_km] int prod_km_ix_by_edge_long;
  array[N_edge,2] int prod_km_ix_by_edge_bounds;
  array[N_competitive_inhibition] int ci_ix_long;
  array[N_edge, 2] int ci_ix_bounds;
  array[N_allostery] int allostery_ix_long;
  array[N_edge, 2] int allostery_ix_bounds;
  array[N_phosphorylation] int phosphorylation_ix_long;
  array[N_edge, 2] int phosphorylation_ix_bounds;
  array[N_mic] int<lower=1,upper=N_metabolite> mic_to_met;
  vector[N_edge] water_stoichiometry;
  vector[N_edge] transported_charge;
  array[N_enzyme_knockout] int<lower=0,upper=N_enzyme> enzyme_knockout_long;
  array[N_experiment, 2] int enzyme_knockout_bounds;
  array[N_pme_knockout] int<lower=0,upper=N_pme> pme_knockout_long;
  array[N_experiment, 2] int pme_knockout_bounds;
  vector<lower=1>[N_enzyme] subunits;
  vector[N_experiment] temperature;
  // configuration
  array[N_experiment] vector<lower=0>[N_mic] conc_init;
  real rel_tol; 
  real abs_tol;
  int max_num_steps;
  int<lower=0,upper=1> likelihood;  // set to 0 for priors-only mode
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
  array[N_experiment] vector[N_edge] dgrs;
}

generated quantities {
  array[N_experiment] vector<lower=0>[N_mic] conc;
  array[N_experiment] vector[N_reaction] flux;
  array[N_experiment] vector[N_pme] conc_pme;
  array[N_experiment] vector[N_unbalanced] conc_unbalanced;
  array[N_experiment] vector[N_enzyme] conc_enzyme;
  array[N_experiment] vector[N_drain] drain;
  array[N_experiment] vector[N_edge] free_enzyme_ratio;
  array[N_experiment] vector[N_edge] saturation;
  array[N_experiment] vector[N_edge] allostery;
  array[N_experiment] vector[N_edge] phosphorylation;
  array[N_experiment] vector[N_edge] reversibility;
  // Sampling experiment boundary conditions from priors
  for (e in 1:N_experiment){
    drain[e] = to_vector(normal_rng(priors_drain[1,e], priors_drain[2,e]));
    conc_pme[e] = to_vector(lognormal_rng(log(priors_conc_pme[1,e]), priors_conc_pme[2,e]));
    conc_unbalanced[e] = to_vector(lognormal_rng(log(priors_conc_unbalanced[1,e]), priors_conc_unbalanced[2,e]));
    conc_enzyme[e] = to_vector(lognormal_rng(log(priors_conc_enzyme[1,e]), priors_conc_enzyme[2,e]));
  }
  // Simulation of experiments
  for (e in 1:N_experiment){
    flux[e] = rep_vector(0, N_reaction);
    vector[N_enzyme] conc_enzyme_experiment = conc_enzyme[e];
    vector[N_pme] conc_pme_experiment = conc_pme[e];
    array[1] vector[N_mic-N_unbalanced] conc_balanced;
    int N_eko_experiment = measure_ragged(enzyme_knockout_bounds, e);
    int N_pko_experiment = measure_ragged(pme_knockout_bounds, e);
    if (N_eko_experiment > 0){
      array[N_eko_experiment] int eko_experiment =
        extract_ragged(enzyme_knockout_long, enzyme_knockout_bounds, e);
      conc_enzyme_experiment[eko_experiment] = rep_vector(0, N_eko_experiment);
    }
    if (N_pko_experiment > 0){
      array[N_pko_experiment] int pko_experiment =
        extract_ragged(pme_knockout_long, pme_knockout_bounds, e);
      conc_pme_experiment[pko_experiment] = rep_vector(0, N_pko_experiment);
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
                  dgrs[e],
                  kcat,
                  km,
                  ki,
                  transfer_constant,
                  dissociation_constant,
                  kcat_pme,
                  conc_pme_experiment,
                  drain[e],
                  temperature[e],
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
    conc[e, balanced_mic_ix] = conc_balanced[1];
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e];
    vector[N_edge] edge_flux = get_edge_flux(conc[e],
                                             conc_enzyme_experiment,
                                             dgrs[e],
                                             kcat,
                                             km,
                                             ki,
                                             transfer_constant,
                                             dissociation_constant,
                                             kcat_pme,
                                             conc_pme_experiment,
                                             drain[e],
                                             temperature[e],
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
    for (j in 1:N_edge){
      flux[e, edge_to_reaction[j]] += edge_flux[j];
    }
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
    phosphorylation[e] = get_phosphorylation(kcat_pme,
                                             conc_pme[e],
                                             phosphorylation_ix_long,
                                             phosphorylation_ix_bounds,
                                             phosphorylation_type,
                                             phosphorylation_pme,
                                             subunits);
    reversibility[e] = get_reversibility(dgrs[e], temperature[e], S, conc[e], edge_type);
  }
}
