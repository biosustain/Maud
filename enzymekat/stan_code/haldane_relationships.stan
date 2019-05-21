/* Functions for getting kinetic parameters from other kinetic parameters plus
   a Keq parameter. */

real get_Kp_uniuni(real V1, real V2, real Keq, real Ka){
  return Keq*V2*Ka/V1;
}
real get_Kip_ordered_unibi(real Keq, real Kia, real Kq, real V1, real V2){
  return (Keq*Kia*V2)/(Kq*V1);
}
real get_Kip_ordered_bibi(real Keq, real Ka, real Kq, real Kib, real Kiq, real V1, real V2){
  return  (Ka*Keq*Kib*V2^2)/(Kq*V1^2);
}
real get_Keq(real delta_g, real T, real R){
/* Function for getting Keq from delta g, temperature and gas constant*/
  return exp(delta_g / (-R * T));
}
