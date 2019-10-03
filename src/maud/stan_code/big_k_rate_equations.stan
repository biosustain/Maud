/* 

   These are mostly copied from this paper:

  - Cleland, W. W. (1963). The kinetics of enzyme-catalyzed reactions with two or
  more substrates or products: I. Nomenclature and rate equations. Biochimica
  et Biophysica Acta (BBA) - Specialized Section on Enzymological Subjects,
  67(), 104–137. http://dx.doi.org/10.1016/0926-6569(63)90211-6

 */

real uniuni(real A, real P, real V1, real V2, real Ka, real Keq){
  /*
    NB this way of formulating the mechanism already takes into account the
    relevant Haldane relationship!
  */
  real num = V1*V2*(A-P/Keq);
  real denom = Ka*V2 + V2*A + V1*P/Keq;
  return num / denom;
}
real ordered_bibi(real A, real B, real P, real Q,
                  real V1, real V2,
                  real Ka, real Kb, real Kp, real Kq,
                  real Kia, real Kib, real Kip, real Kiq,
                  real Keq){
  real num = V1*V2*(A*B - P*Q/Keq);
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
real ordered_unibi(real A, real P, real Q,
                   real V1, real V2,
                   real Ka, real Kp, real Kq,
                   real Kia, real Kip, real Kiq,
                   real Keq){
  real num = V1*V2*(A-P*Q/Keq); 
  real denom =
    Ka*V2
    + V2*A
    + Kq*V1*P/Keq
    + Kp*V1*Q/Keq
    + V1*P*Q/Keq
    + V2*A*P/Kip;
  return num / denom;
}
real pingpong(real A, real B, real P, real Q,
              real V1, real V2,
              real Ka, real Kb, real Kp, real Kq,
              real Kia, real Kib, real Kip, real Kiq,
              real Keq){
  real num = V1*V2*(A*B-P*Q/Keq);
  real denom =
      Kb*V2*A
    + Ka*V2*B
    + V2*A*B
    + Kq*V1*P/Keq
    + Kp*V1*Q/Keq
    + V1*P*Q/Keq
    + Kq*V1*A*P/(Kia*Keq)
    + Ka*V2*B*Q/Kiq;
  return num / denom;
}
real ordered_terbi(real A, real B, real C, real P, real Q,
                   real V1, real V2,
                   real Ka, real Kb, real Kc, real Kp, real Kq,
                   real Kia, real Kib, real Kic, real Kip, real Kiq,
                   real Keq){
  real num = V1*V2*(A*B*C-P*Q/Keq);
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
  return num/denom;
}
real modular_rate_law(real Tr, real Dr){
  return(Tr/(Dr));
}
real irr_mass_action(real A, real V1){
  return(A*V1);
}
real fixed_flux(real f){
  return f;
}
