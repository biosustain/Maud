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
  int<lower=0> N_ki;
  int<lower=0> N_ai;
  int<lower=0> N_aa;
  int<lower=0> N_ae;
  int<lower=0> N_pa;
  int<lower=0> N_pi;
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
  int<lower=1,upper=N_mic> unbalanced_mic_ix[N_unbalanced];
  int<lower=1,upper=N_mic> balanced_mic_ix[N_mic-N_unbalanced];
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
}
transformed data {
  real initial_time = 0;
  matrix[N_experiment, N_enzyme] knockout = rep_matrix(1, N_experiment, N_enzyme) - is_knockout;
  matrix[N_experiment, N_phosphorylation_enzymes] phos_knockout =
    rep_matrix(1, N_experiment, N_phosphorylation_enzymes) - is_phos_knockout;
}
parameters {
  vector[N_km] km;
  vector[N_ki] ki;
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
                  conc_unbalanced[e,:],
                  balanced_mic_ix,
                  unbalanced_mic_ix,
                  conc_enzyme_experiment,
                  km,
                  drain[e,:],
                  km_lookup,
                  S,
                  edge_type,
                  edge_to_drain,
                  edge_to_enzyme,
                  kcat,
                  dgrs,
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
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e,:];
    {
    vector[N_edge] flux_edge = get_flux(conc[e],
                                        conc_enzyme_experiment,
                                        km,
                                        drain[e,:],
                                        km_lookup,
                                        S,
                                        edge_type,
                                        edge_to_drain,
                                        edge_to_enzyme,
                                        kcat,
                                        dgrs,
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
  }
}
