real ordered_bibi(real A, real Kib, real Kiq, real B, real Kia, real Kip, real P, real Kb, real Kq, real Q, real Ka, real Kp, real V1, real V2){
  return (A*V1 - P*Q*V2)/(A*(Kip + P) + Kip*(Ka + Kq*P + (Kp + P)*Q));
}
real ordered_unibi(real A, real Kip, real P, real Ka, real Kq, real Kp, real Q, real V1, real V2){
  return (A*V1 - P*Q*V2)/(A*(Kip + P) + Kip*(Ka + Kq*P + (Kp + P)*Q));
}
real uniuni(real V1, real V2, real A, real P, real Ka, real Kp){
  return (V1*A/Ka-V2*P/Kp)/(1+A/Ka+P/Kp);
}
