/* 
Functions for getting kinetic parameters from other kinetic parameters plus a
Keq parameter. Copied from this paper:

- Cleland, W. W. (1963). The kinetics of enzyme-catalyzed reactions with two or
  more substrates or products: I. Nomenclature and rate equations. Biochimica
  et Biophysica Acta (BBA) - Specialized Section on Enzymological Subjects,
  67(), 104–137. http://dx.doi.org/10.1016/0926-6569(63)90211-6

*/

real get_Kp_uniuni(real V1, real Ka, real V2, real Keq){
  return Keq * V2 * Ka / V1;
}
real get_Kip_ordered_unibi(real Keq, real Kia, real Kq, real V1, real V2){
  return Keq*Kia*V2 / (Kq*V1);
}
real get_Kiq_ordered_unibi(real Keq, real Ka, real Kp, real V1, real V2){
  return Keq*V2*Ka / (V1*Kp);
}
real get_Kip_ordered_bibi(real Keq, real Ka, real Kq, real Kib, real V1, real V2){
  return  Ka*Keq*Kib*V2^2 / (Kq*V1^2);
}
real get_Kia_ordered_bibi(real Keq, real Kb, real Kp, real Kiq, real V1, real V2){
  return V1*Kp*Kiq / (V2*Keq*Kb);
}
real get_Kp_ping_pong(real Keq, real Ka, real Kb, real Kq, real V1, real V2){
  return (V2/V1)^2*Keq*Ka*Kb/Kq;
}
real get_Kp_ordered_terbi(real Keq, real Kc, real Kia, real Kib, real Kiq, real V1, real V2){
  return Keq*V2*Kia*Kib*Kc/(V1*Kiq);
}
real get_Kip_ordered_terbi(real Keq, real Kia, real Kib, real Kic, real Kiq){
  return Keq*Kia*Kib*Kic/Kiq;
}

/* Function for getting Keq from delta g, temperature and gas constant*/
real get_Keq(real delta_g, real T, real R){
  return exp(delta_g / (-R * T));
}
