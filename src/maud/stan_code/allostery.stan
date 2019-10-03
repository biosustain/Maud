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
real get_free_enzyme_ratio_ordered_bibi(real A, real B, real P, real Q,
                   real V1, real V2,
                   real Ka, real Kb, real Kp, real Kq,
                   real Kia, real Kib, real Kip, real Kiq,
                   real Keq){
  real num = Kia*Kb*V2 + Kq*V1*P/Keq + Ka*V2*B;
  real denom =
    Kia*Kb*V2
    + Kb*V2*A
    + Ka*V2*B
    + V2*A*B
    + Kq*V1*P/Keq
    + Kp*V1*Q/Keq
    + V1*P*Q/Keq
    + Kq*V1*A*P/(Kia*Keq)
    + Ka*V2*B*Q/Kiq
    + V2*A*B*P/Kip
    + V1*B*P*Q/(Kib*Keq);
  return num / denom;
}
real get_free_enzyme_ratio_ordered_terbi(real A, real B, real C, real P, real Q,
                   real V1, real V2,
                   real Ka, real Kb, real Kc, real Kp, real Kq,
                   real Kia, real Kib, real Kic, real Kip, real Kiq,
                   real Keq){
  real num = Kq*V1*P/Keq + Ka*V2*B*C + Kia*Kb*V2*C + Kia*Kib*Kc*V2;
  real denom = 
      Kia*Kib*Kc*V2
    + Kib*Kc*V2*A
    + Kia*Kb*V2*C
    + Kc*V2*A*B
    + Kb*V2*A*C
    + Ka*V2*B*C
    + V2*A*B*C
    + Kp*V1*Q/Keq
    + Kq*V1*P/Keq
    + V1*P*Q/Keq
    + Kq*V1*A*P/(Kia*Keq)
    + Kia*Kb*V2*C*Q/Kiq
    + Kq*V1*A*B*P/(Kia*Kib*Keq)
    + Ka*V2*B*C*Q/Kiq
    + Ka*Kic*V2*B*P*Q/(Kip*Kiq)
    + Kia*Kb*V2*C*P*Q/(Kip*Kiq)
    + Kq*V1*A*B*C*P/(Kia*Kib*Kic*Keq)
    + Ka*V2*B*C*P*Q/(Kip*Kiq);
  return num / denom;
}

