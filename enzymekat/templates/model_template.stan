model {
  params ~ lognormal(prior_location, prior_scale);
  if (LIKELIHOOD == 1){
    vector[N_measurement_conc] conc_hat;
    vector[N_measurement_flux] flux_hat;
    for (mc in 1:N_measurement_conc){
      conc_hat[mc] = metabolite_concentration[metabolite_ix[mc], experiment_ix_conc[mc]];
    }
    for (mf in 1:N_measurement_flux){
      flux_hat[mf] = flux[reaction_ix[mf], experiment_ix_flux[mf]];
    }
    measurement_conc ~ lognormal(log(conc_hat), measurement_scale_conc);
    measurement_flux ~ normal(flux_hat, measurement_scale_flux);
  }
}
