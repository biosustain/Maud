real get_Tr(vector metabolite,
            vector km,
            vector stoichiometry,
            real kcat,
            real keq){
  real plus_product = 1;
  real minus_product = 1;
  real k_minus = (kcat / keq);
  for (m in 1:size(metabolite)){
    real multiplicand = (metabolite[m] / km[m]) ^ abs(stoichiometry[m]);
    k_minus *= km[m] ^ stoichiometry[m];
    if (stoichiometry[m] < 0)
      plus_product *= multiplicand;
    else
      minus_product *= multiplicand;
  }
  return kcat * plus_product - k_minus * minus_product;
}

real get_Dr_common_rate_law(vector metabolite, vector km, vector stoichiometry){
  real psi_plus = 1;
  real psi_minus = 1;
  for (m in 1:size(metabolite)){
    real multiplicand = (1 + metabolite[m] / km[m]) ^ abs(stoichiometry[m]);
    if (stoichiometry[m] < 0)
      psi_plus *= multiplicand;
    else
      psi_minus *= multiplicand;
  }
  return psi_plus + psi_minus - 1;
}

real get_Dr_reg(){
  return 0;
}

real modular_rate_law(vector metabolite,
                      vector km,
                      vector stoichiometry,
                      real kcat,
                      real keq,
                      real enz){
  real Tr = get_Tr(metabolite, km, stoichiometry, kcat, keq);
  real Dr = get_Dr_common_rate_law(metabolite, km, stoichiometry);
  real Dr_reg = get_Dr_reg();
  return enz * Tr / (Dr + Dr_reg);
}
