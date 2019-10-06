functions{
#include big_k_rate_equations.stan
#include haldane_relationships.stan
#include allostery.stan
  vector get_fluxes(real[] m, real[] p){
  real empty_array[0];
  real free_enzyme_ratio_r2 = get_free_enzyme_ratio_uniuni(m[3],m[4],p[2]*p[9],p[2]*p[10],p[11],p[8]);
  return [
    uniuni(m[1],m[3],p[1]*p[5],p[1]*p[6],p[7],p[4]),
    uniuni(m[3],m[4],p[2]*p[9],p[2]*p[10],p[11],p[8])*get_regulatory_effect(empty_array,{m[3]},free_enzyme_ratio_r2,empty_array,{p[13]},p[12]),
    uniuni(m[4],m[2],p[3]*p[15],p[3]*p[16],p[17],p[14])
  ]';
}
  real[] ode_func(real t, real[] m, real[] p, real[] xr, int[] xi){
  vector[3] fluxes = get_fluxes(m, p);
  return {
    0,
    0,
    1*fluxes[1]-1*fluxes[2],
    1*fluxes[2]-1*fluxes[3]
  };
}
  vector steady_state_function(vector balanced, vector theta, real[] xr, int[] xi){
  int N_unbalanced = 2;
  int N_balanced = 2;
  real initial_time = 0;
  real time_step = 0.05;
  real conc[4];
  real balanced_new[2];
  conc[{3,4}] = to_array_1d(balanced);
  conc[{1,2}] = to_array_1d(theta[1:2]);
  balanced_new = integrate_ode_bdf(
    ode_func,
    conc,
    initial_time,
    rep_array(time_step, 1),
    to_array_1d(theta[N_unbalanced+1:]),
    xr,
    rep_array(0, 1),
    1e-8, 1e-8, 1e5
  )[1, {3,4}]; 
  return to_vector(balanced_new) - balanced;
}
}
data {
  // dimensions
  int<lower=1> N_balanced;    // 'Balanced' metabolites must have constant concentration at steady state
  int<lower=1> N_unbalanced;  // 'Unbalanced' metabolites can have changing concentration at steady state
  int<lower=1> N_kinetic_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_enzyme;
  int<lower=1> N_experiment;
  int<lower=1> N_flux_measurement;
  int<lower=1> N_conc_measurement;
  int<lower=1> stoichiometric_rank;
  // measurements
  int<lower=1,upper=N_experiment> experiment_yconc[N_conc_measurement];
  int<lower=1,upper=N_balanced+N_unbalanced> metabolite_yconc[N_conc_measurement];
  vector[N_conc_measurement] yconc;
  vector<lower=0>[N_conc_measurement] sigma_conc;
  int<lower=1,upper=N_experiment> experiment_yflux[N_flux_measurement];
  int<lower=1,upper=N_reaction> reaction_yflux[N_flux_measurement];
  vector[N_flux_measurement] yflux;
  vector<lower=0>[N_flux_measurement] sigma_flux;
  // hardcoded priors
  vector[N_kinetic_parameter] prior_loc_kinetic_parameter;
  vector<lower=0>[N_kinetic_parameter] prior_scale_kinetic_parameter;
  real prior_loc_unbalanced[N_unbalanced, N_experiment];
  real<lower=0> prior_scale_unbalanced[N_unbalanced, N_experiment];
  real prior_loc_enzyme[N_enzyme, N_experiment];
  real<lower=0> prior_scale_enzyme[N_enzyme, N_experiment];
  vector<lower=0>[N_balanced] balanced_guess;
  // network properties
  matrix[N_enzyme, stoichiometric_rank] ln_equilibrium_basis;
  // algebra solver configuration
  real rel_tol;
  real f_tol;
  int max_steps;
  // likelihood configuration - set to 0 for priors-only mode
  int<lower=0,upper=1> LIKELIHOOD;
}
transformed data {
  real xr[0];
  int xi[0];
  int<lower=0> Keq_pos[N_enzyme] = {1,5,11};
}
parameters {
  vector<lower=0>[N_kinetic_parameter] kinetic_parameter;
  vector<lower=0>[N_unbalanced] unbalanced[N_experiment];
  vector<lower=0>[N_enzyme] enzyme_concentration[N_experiment];
  matrix[stoichiometric_rank, 1] basis_contribution;
}
transformed parameters {
  vector<lower=0>[N_balanced+N_unbalanced] conc[N_experiment];
  vector[N_reaction] flux[N_experiment];
  vector[N_kinetic_parameter] updated_kinetic_parameters = kinetic_parameter;
  matrix[1, N_enzyme] Keq = exp(ln_equilibrium_basis*basis_contribution)';
  for (k in 1:N_enzyme){
    updated_kinetic_parameters[Keq_pos[k]] = Keq[1, k];
  }
  for (e in 1:N_experiment){
    vector[N_unbalanced+N_enzyme+N_kinetic_parameter] theta = append_row(unbalanced[e], append_row(enzyme_concentration[e], updated_kinetic_parameters));
    conc[e, {3,4}] = algebra_solver(steady_state_function, balanced_guess, theta, xr, xi, rel_tol, f_tol, max_steps);
    conc[e, {1,2}] = unbalanced[e];
    flux[e] = get_fluxes(to_array_1d(conc[e]), append_array(to_array_1d(enzyme_concentration[e]), to_array_1d(kinetic_parameter)));
  }
}
model {
  kinetic_parameter ~ lognormal(log(prior_loc_kinetic_parameter), prior_scale_kinetic_parameter);
  for (e in 1:N_experiment){
    unbalanced[e] ~ lognormal(log(prior_loc_unbalanced[,e]), prior_scale_unbalanced[,e]);
    enzyme_concentration[e] ~ lognormal(log(prior_loc_enzyme[,e]), prior_scale_enzyme[,e]);
  }
  if (LIKELIHOOD == 1){
    for (c in 1:N_conc_measurement){
      target += lognormal_lpdf(yconc[c] | log(conc[experiment_yconc[c], metabolite_yconc[c]]), sigma_conc[c]);
    }
    for (f in 1:N_flux_measurement){
      target += normal_lpdf(yflux[f] | flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
    }
  }
}
generated quantities {
  vector[N_conc_measurement] yconc_sim;
  vector[N_flux_measurement] yflux_sim;
  for (c in 1:N_conc_measurement){
    yconc_sim[c] = lognormal_rng(log(conc[experiment_yconc[c], metabolite_yconc[c]]), sigma_conc[c]);
  }
  for (f in 1:N_flux_measurement){
    yflux_sim[f] = normal_rng(flux[experiment_yflux[f], reaction_yflux[f]], sigma_flux[f]);
  }
}
