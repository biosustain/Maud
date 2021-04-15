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
  int<lower=1> N_samples;
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
  // hardcoded priors
  vector[N_metabolite] formation_energy;
  vector[N_enzyme] kcat;
  vector[N_km] km;
  vector[N_ki] ki;
  vector[N_ai] diss_t;
  vector[N_aa] diss_r;
  vector[N_ae] tc;
  vector[N_phosphorylation_enzymes] phos_kcat;
  array[N_experiment] vector[N_phosphorylation_enzymes] phos_conc;
  array[N_experiment] vector[N_unbalanced] unbalanced;
  array[N_experiment] vector[N_enzyme] enzyme;
  array[N_experiment] vector[N_drain] drain;
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
  real<lower=0> timepoint;
}
transformed data {
  real initial_time = 0;
  matrix[N_experiment, N_enzyme] knockout =
    rep_matrix(1, N_experiment, N_enzyme) - is_knockout;
  matrix[N_experiment, N_phosphorylation_enzymes] phos_knockout =
    rep_matrix(1, N_experiment, N_phosphorylation_enzymes) - is_phos_knockout;
}
transformed parameters {
}
model {
}
generated quantities {
  array[N_experiment] vector<lower=0>[N_mic] conc;
  array[N_experiment] vector[N_reaction] flux;
  vector[N_enzyme] keq = get_keq(S_enz,
                                 formation_energy,
                                 mic_to_met,
                                 water_stoichiometry);
  for (e in 1:N_experiment){
    vector[N_enzyme] experiment_enzyme = enzyme[e] .* knockout[e]';
    vector[N_phosphorylation_enzymes] experiment_phos_conc = 
      phos_enzyme_conc[e] .* phos_knockout[e]';
    vector[N_mic-N_unbalanced] conc_balanced = ode_bdf_tol(dbalanced_dt,
                                                              conc_init[e, balanced_mic_ix],
                                                              initial_time,
                                                              {timepoint},
                                                              rel_tol, abs_tol, max_num_steps,
                                                              conc_unbalanced[e],
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
                                                              drain[e]);
    conc[e, balanced_mic_ix] = conc_balanced;
    conc[e, unbalanced_mic_ix] = conc_unbalanced[e];
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
                              phos_enzyme_conc[e],
                              phos_enzyme_kcat,
                              S_phos_act,
                              S_phos_inh));
  }
