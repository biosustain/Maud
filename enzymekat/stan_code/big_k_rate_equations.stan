real uniuni(real V1, real V2, real A, real P, real Ka, real Kp){
  return (V1*A/Ka-V2*P/Kp)/(1+A/Ka+P/Kp);
}
real   ordered_bibi(real Kib, real Kip, real A, real B, real Kiq, real Kp, real P, real Kb, real Kq, real Kia, real Q, real Ka, real V1, real V2){
  return  (Kib*Kip*(A*B*Kiq*Kp*V1 - Kb*Kia*P*Q*V2))/(A*Kib*(B*Kiq*Kp*(Kip + P) + Kb*Kip*(Kiq*Kp + Kq*P)) + Kip*(B*Kb*Kia*P*Q + B*Ka*Kib*Kp*(Kiq + Q) + Kb*Kia*Kib*(Kiq*Kp + Kq*P + (Kp + P)*Q)));
}
real   ordered_unibi(real Kip, real Ka, real Kiq, real Kp, real A, real P, real Kia, real Kq, real Q, real V1, real V2){
  return  (Kip*(A*Kiq*Kp*V1 - Kia*P*Q*V2))/(Ka*Kip*Kiq*Kp + A*Kiq*Kp*(Kip + P) + Kia*Kip*(Kq*P + (Kp + P)*Q));
}

