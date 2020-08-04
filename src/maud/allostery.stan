real get_allostery(vector activator_concentration,    // metabolite
                   vector inhibitor_concentration,    // metabolite
                   real free_enzyme_ratio,            // derived from rate equation
                   vector dissociation_constant_r,    // parameter
                   vector dissociation_constant_t,    // parameter
                   real transfer_constant){           // parameter
  if ((rows(activator_concentration) == 0) && (rows(inhibitor_concentration) == 0)){
    return 1;
  }
  else {
    real Q_num = 1 + sum(inhibitor_concentration ./ dissociation_constant_t);
    real Q_denom = 1 + sum(activator_concentration ./ dissociation_constant_r);
    return inv(1 + transfer_constant * free_enzyme_ratio * Q_num / Q_denom);
  }
}
