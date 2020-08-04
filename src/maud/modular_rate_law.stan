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

real get_Dr_reg(vector conc_ci, vector ki){
  if (rows(conc_ci) > 0){
    return 0;
  }
  else {
    return prod(conc_ci ./ ki);
  }
}

real get_free_enzyme_ratio(vector metabolite,
                           vector km,
                           vector stoichiometry,
                           vector conc_ci,
                           vector ki){
  real Dr = get_Dr_common_rate_law(metabolite, km, stoichiometry);
  real Dr_reg = get_Dr_reg(conc_ci, ki);
  return 1 / (Dr + Dr_reg);
}

real modular_rate_law(vector metabolite,
		      vector km,
		      vector stoichiometry,
		      real kcat,
		      real keq,
		      real enz,
		      vector conc_ci,
		      vector ki){
  real free_enzyme_ratio = get_free_enzyme_ratio(metabolite,
                                                 km,
                                                 stoichiometry,
                                                 conc_ci,
                                                 ki);
  real Tr = get_Tr(metabolite, km, stoichiometry, kcat, keq);
  return enz * free_enzyme_ratio * Tr;
}

