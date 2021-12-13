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
  // hardcoded priors
  array[2, N_experiment] vector[N_phosphorylation_enzymes] priors_conc_phos;
  array[2, N_experiment] vector[N_unbalanced] priors_conc_unbalanced;
  array[2, N_experiment] vector[N_enzyme] priors_conc_enzyme;
  array[2, N_experiment] vector[N_drain] priors_drain;
  // network properties
  int<lower=1,upper=N_mic> unbalanced_mic_ix[N_unbalanced];
  int<lower=1,upper=N_mic> balanced_mic_ix[N_mic-N_unbalanced];
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
}
transformed data {
  real initial_time = 0;
  matrix[N_experiment, N_enzyme] knockout = rep_matrix(1, N_experiment, N_enzyme) - is_knockout;
  matrix[N_experiment, N_phosphorylation_enzymes] phos_knockout =
    rep_matrix(1, N_experiment, N_phosphorylation_enzymes) - is_phos_knockout;
}
parameters {
  vector[N_km] km;
  vector[N_ci] ki;
  vector[N_enzyme] kcat;
  vector[N_ai] diss_t;
  vector[N_aa] diss_r;
  vector[N_ae] transfer_constant;
  vector[N_phosphorylation_enzymes] kcat_phos;
  vector[N_edge] dgrs;
}

generated quantities {
  array[N_experiment] vector<lower=0>[N_mic] conc;
  array[N_experiment] vector[N_reaction] flux;
  array[N_experiment] vector[N_phosphorylation_enzymes] conc_phos;
  array[N_experiment] vector[N_unbalanced] conc_unbalanced;
  array[N_experiment] vector[N_enzyme] conc_enzyme;
  array[N_experiment] vector[N_drain] drain;
  // Sampling experiment boundary conditions from priors
  for (e in 1:N_experiment){
    drain[e] = to_vector(normal_rng(priors_drain[1,e], priors_drain[2,e]));
    conc_phos[e] = to_vector(lognormal_rng(log(priors_conc_phos[1,e]), priors_conc_phos[2,e]));
    conc_unbalanced[e] = to_vector(lognormal_rng(log(priors_conc_unbalanced[1,e]), priors_conc_unbalanced[2,e]));
    conc_enzyme[e] = to_vector(lognormal_rng(log(priors_conc_enzyme[1,e]), priors_conc_enzyme[2,e]));
  }
  // Simulation of experiments
  for (e in 1:N_experiment){
    flux[e] = rep_vector(0, N_reaction);
    real timepoints[2] = {timepoint, timepoint + 10};
    vector[N_enzyme] conc_enzyme_experiment = conc_enzyme[e] .* knockout[e]';
    vector[N_phosphorylation_enzymes] conc_phos_experiment = conc_phos[e] .* phos_knockout[e]';
    vector[N_mic-N_unbalanced] conc_balanced[2] = ode_bdf_tol(dbalanced_dt,
                  conc_init[e, balanced_mic_ix],
                  initial_time,
                  timepoints,
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
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e,:];
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
    for (j in 1:N_edge){
      flux[e, edge_to_reaction[j]] += edge_flux[j];
    }
  }
}
