data {
  // dimensions
  int<lower=1> N_balanced;    // 'Balanced' metabolites must have constant concentration at steady state
  int<lower=1> N_unbalanced;  // 'Unbalanced' metabolites can have changing concentration at steady state
  int<lower=1> N_training_experiments;  // Number of training conditions
  int<lower=1> N_holdout_experiments;   // Number of holdout conditions
  int<lower=1> N_kinetic_parameter;
  int<lower=1> N_reaction;
  int<lower=1> N_experiment;
  int<lower=1> N_known_real;
  int<lower=1> N_flux_measurement_training;
  int<lower=1> N_flux_measurement_holdout;
  int<lower=1> N_concentration_measurement_training;
  int<lower=1> N_concentration_measurement_holdout;
  // what experiments are used for training and which are held out
  int<lower=1> training_experiments[N_training_experiments];
  int<lower=1> holdout_experiments[N_holdout_experiments];
  // position of balanced and unbalanced metabolites in overall metabolite array
  int<lower=1,upper=N_balanced+N_unbalanced> pos_balanced[N_balanced];
  int<lower=1,upper=N_balanced+N_unbalanced> pos_unbalanced[N_unbalanced];
  // which training and holdout measurements are associated with each experiments
  int<lower=1,upper=N_experiment> ix_training_concentration_experiment[N_concentration_measurement_training];
  int<lower=1,upper=N_experiment> ix_holdout_concentration_experiment[N_concentration_measurement_holdout];
  int<lower=1,upper=N_experiment> ix_training_flux_experiment[N_flux_measurement_training];
  int<lower=1,upper=N_experiment> ix_holdout_flux_experiment[N_flux_measurement_holdout];
  // inter-group indexing for metabolites and reactions for holdout and training datasets
  int<lower=1,upper=N_training_experiments> ig_training_concentration_experiment[N_concentration_measurement_training];
  int<lower=1,upper=N_holdout_experiments> ig_holdout_concentration_experiment[N_concentration_measurement_holdout];
  int<lower=1,upper=N_training_experiments> ig_training_flux_experiment[N_flux_measurement_training];
  int<lower=1,upper=N_holdout_experiments> ig_holdout_flux_experiment[N_flux_measurement_holdout];
  // which measurments are associated with holdout and training metabolites and reactions
  int<lower=1,upper=N_balanced+N_unbalanced> ix_training_concentration_measurement[N_concentration_measurement_training];
  int<lower=1,upper=N_balanced+N_unbalanced> ix_holdout_concentration_measurement[N_concentration_measurement_holdout];
  int<lower=1,upper=N_reaction> ix_training_flux_measurement[N_flux_measurement_training];
  int<lower=1,upper=N_reaction> ix_holdout_flux_measurement[N_flux_measurement_holdout];
  // measurements
  vector[N_flux_measurement_training] flux_measurement_training;
  vector[N_concentration_measurement_training] concentration_measurement_training;
  vector[N_flux_measurement_holdout] flux_measurement_holdout;
  vector[N_concentration_measurement_holdout] concentration_measurement_holdout;
  // measurement scales, i.e. how much the measured values deviate from the true values
  vector<lower=0>[N_flux_measurement_training] flux_measurement_scale_training;
  vector<lower=0>[N_flux_measurement_holdout] flux_measurement_scale_holdout;
  vector<lower=0>[N_concentration_measurement_training] concentration_measurement_scale_training;
  vector<lower=0>[N_concentration_measurement_holdout] concentration_measurement_scale_holdout;
  // hardcoded priors
  vector[N_kinetic_parameter] prior_location_kinetic_parameter;
  vector<lower=0>[N_kinetic_parameter] prior_scale_kinetic_parameter;
  real prior_location_unbalanced[N_unbalanced, N_experiment];
  real<lower=0> prior_scale_unbalanced[N_unbalanced, N_experiment];
  // other hardcoded numbers
  real known_reals[N_known_real, N_experiment];
  // ode configuration
  real initial_time;
  real steady_time;
  real rel_tol;
  real abs_tol;
  int max_steps;
  // likelihood configuration - set to 0 for priors-only mode
  int<lower=0,upper=1> LIKELIHOOD;
}
transformed data {
  int known_ints[0];
}
parameters {
  real<lower=0> kinetic_parameter[N_kinetic_parameter];
  real<lower=0> concentration_unbalanced[N_unbalanced, N_experiment];
}
transformed parameters {
  real concentration_training[N_balanced+N_unbalanced, N_training_experiments];
  real flux_training[N_reaction, N_training_experiments];
  for (e in 1:N_training_experiments){
    int te = training_experiments[e];
    real initial_concentration[N_balanced+N_unbalanced];
    initial_concentration[pos_balanced] = rep_array(1.0, N_balanced);
    initial_concentration[pos_unbalanced] = concentration_unbalanced[,te];
    concentration_training[,e] = integrate_ode_bdf(steady_state_equation,
                                          initial_concentration,
                                          initial_time,
                                          {steady_time},
                                          kinetic_parameter,
                                          known_reals[,te],
                                          known_ints,
                                          rel_tol, abs_tol, max_steps)[1];
    flux_training[,e] = get_fluxes(concentration_training[,e], kinetic_parameter, known_reals[,te]);
  }
}
model {
  kinetic_parameter ~ lognormal(log(prior_location_kinetic_parameter), prior_scale_kinetic_parameter);
  for (e in 1:N_experiment){
    prior_location_unbalanced[,e] ~ lognormal(log(concentration_unbalanced[,e]), prior_scale_unbalanced[,e]);
  }
  if (LIKELIHOOD == 1){
    real concentration_hat[N_concentration_measurement_training];
    real flux_hat[N_flux_measurement_training];
    for (mc in 1:N_concentration_measurement_training){
      concentration_hat[mc] = concentration_training[ix_training_concentration_measurement[mc],
      ig_training_concentration_experiment[mc]];
    }
    for (mf in 1:N_flux_measurement_training){
      flux_hat[mf] = flux_training[ix_training_flux_measurement[mf], ig_training_flux_experiment[mf]];
    }
    concentration_measurement_training ~ lognormal(log(concentration_hat), concentration_measurement_scale_training);
    flux_measurement_training ~ normal(flux_hat, flux_measurement_scale_training);
  }
}
generated quantities {
  vector[N_concentration_measurement_training] simulated_concentration_measurement;
  vector[N_concentration_measurement_holdout+N_flux_measurement_holdout] log_like;
  real initial_concentration[N_balanced+N_unbalanced];
  real concentration_holdout[N_balanced+N_unbalanced, N_holdout_experiments];
  real concentration_hat_holdout[N_concentration_measurement_holdout];
  real flux_hat_holdout[N_flux_measurement_holdout];
  real flux_holdout[N_reaction, N_holdout_experiments];
  vector[N_flux_measurement_training] simulated_flux_measurement;
  real balanced_metabolite_rate_of_change[N_balanced, N_training_experiments];
  for (e in 1:N_training_experiments){
    balanced_metabolite_rate_of_change[, e] = get_odes(flux_training[, e])[pos_balanced];
  }
  for (mc in 1:N_concentration_measurement_training){
    simulated_concentration_measurement[mc] =
      lognormal_rng(log(concentration_training[ix_training_concentration_measurement[mc],                               ig_training_concentration_experiment[mc]]),
                                  concentration_measurement_scale_training[mc]);
  }
  for (mf in 1:N_flux_measurement_training){
    simulated_flux_measurement[mf] =
      normal_rng(flux_training[ix_training_flux_measurement[mf],
                      ig_training_flux_experiment[mf]],
                      flux_measurement_scale_training[mf]);
  }
  // Cross validation experiments
  for (e in 1:N_holdout_experiments){
    int he = holdout_experiments[e];
    initial_concentration[pos_balanced] = rep_array(1.0, N_balanced);
    initial_concentration[pos_unbalanced] = concentration_unbalanced[,he];
    concentration_holdout[,e] = integrate_ode_bdf(steady_state_equation,
                                          initial_concentration,
                                          initial_time,
                                          {steady_time},
                                          kinetic_parameter,
                                          known_reals[,he],
                                          known_ints,
                                          rel_tol, abs_tol, max_steps)[1];
    flux_holdout[,e] = get_fluxes(concentration_holdout[,e], kinetic_parameter, known_reals[,he]);
  }

  for (mc in 1:N_concentration_measurement_holdout){
    concentration_hat_holdout[mc] = concentration_holdout[ix_holdout_concentration_measurement[mc],
                          ig_holdout_concentration_experiment[mc]];

    log_like[mc] = lognormal_lpdf(concentration_measurement_holdout[mc]|concentration_hat_holdout[mc],
                                  concentration_measurement_scale_holdout[mc]);
  }

  for (mf in 1:N_flux_measurement_holdout){
    flux_hat_holdout[mf] = flux_holdout[ix_holdout_flux_measurement[mf],
                                ig_holdout_flux_experiment[mf]];

    log_like[mf+N_concentration_measurement_holdout] = normal_lpdf(flux_measurement_holdout[mf]|flux_hat_holdout[mf],                                               flux_measurement_scale_holdout[mf]);
  }
}
