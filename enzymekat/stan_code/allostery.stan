real get_regulatory_effect(real[] activator_concentration,    // metabolite
                           real[] inhibitor_concentration,    // metabolite
                           real free_enzyme_ratio,            // derived from rate equation
                           real[] dissociation_constant_r,    // parameter
                           real[] dissociation_constant_t,    // parameter
                           real transfer_constant){           // parameter
  real Q_num = size(inhibitor_concentration) == 0 ? 1 :
    1 + sum(to_vector(inhibitor_concentration) ./ to_vector(dissociation_constant_t));
  real Q_denom = size(activator_concentration) == 0 ? 1
    : 1 + sum(to_vector(activator_concentration) ./ to_vector(dissociation_constant_r));
  real Q = transfer_constant * free_enzyme_ratio * Q_num / Q_denom;
  return inv(1 + Q); 
}
real get_free_enzyme_ratio_uniuni(real A, real P, real V1, real V2, real Ka, real Keq){
  real Kp = get_Kp_uniuni(V1, Ka, V2, Keq);
  return 1 / (1 + A/Ka + P/Kp);
}
real get_free_enzyme_ratio_ordered_unibi(real A, real P, real Q,
                                         real V1, real V2,
                                         real Ka, real Kp, real Kq,
                                         real Kia, real Kip, real Kiq,
                                         real Keq){
  real num = Ka *V2 + Kq*V1*P/Keq;
  real denom =
    Ka*V2
    + V2*A
    + Kq*V1*P/Keq
    + Kp*V1*Q/Keq
    + V1*P*Q/Keq
    + V2*A*P/Kip;
  return num / denom;
}
