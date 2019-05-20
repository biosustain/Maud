real uniuni(real A, real P, real V1, real V2, real Ka, real Kp){
  return (V1*A/Ka-V2*P/Kp)/(1+A/Ka+P/Kp);
}
real ordered_bibi(real A, real B, real P, real Q,
                  real V1, real V2,
                  real Ka, real Kb, real Kp, real Kq,
                  real Kia, real Kib, real Kip, real Kiq){
  return (Kib*Kip*(A*B*Kiq*Kp*V1 - Kb*Kia*P*Q*V2))/(A*Kib*(B*Kiq*Kp*(Kip + P) + Kb*Kip*(Kiq*Kp + Kq*P)) + Kip*(B*Kb*Kia*P*Q + B*Ka*Kib*Kp*(Kiq + Q) + Kb*Kia*Kib*(Kiq*Kp + Kq*P + (Kp + P)*Q)));
}
real ordered_unibi(real A, real P, real Q,
                   real V1, real V2,
                   real Ka, real Kp, real Kq,
                   real Kia, real Kip, real Kiq){
  return (Kip*(A*Kiq*Kp*V1 - Kia*P*Q*V2))/(Ka*Kip*Kiq*Kp + A*Kiq*Kp*(Kip + P) + Kia*Kip*(Kq*P + (Kp + P)*Q));
}

