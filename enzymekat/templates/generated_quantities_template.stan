generated quantities {
  vector[N_measurement_flux] flux_pred;
  vector[N_measurement_conc] conc_pred;
  real metabolite_flux[N_metabolite, N_experiment];
  for (e in 1:N_experiment){
    metabolite_flux[, e] = get_odes(flux[, e]);
  }
  for (mc in 1:N_measurement_conc){
    conc_pred[mc] = lognormal_rng(log(metabolite_concentration[metabolite_ix[mc], experiment_ix_conc[mc]]),
                                  measurement_scale_conc[mc]);
  }
         for (mf in 1:N_measurement_flux){
           flux_pred[mf] = normal_rng(flux[reaction_ix[mf], experiment_ix_flux[mf]], measurement_scale_flux[mf]);
         }
}
