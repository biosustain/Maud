real[] get_fluxes(real[] ode_metabolites,
                  real[] params,
                  real[] known_reals){
  // unpacking
  real NA = known_reals[1];
  real sum_AXP = known_reals[2];
  real sum_NAD = known_reals[3];
  real sum_UXP = known_reals[4];
  real volume = known_reals[5];
  real ENO1 = known_reals[6];
  real ENO2 = known_reals[7];
  real FBA = known_reals[8];
  real GPM = known_reals[9];
  real PGK = known_reals[10];
  real TDH1 = known_reals[11];
  real TDH2 = known_reals[12];
  real TDH3 = known_reals[13];
  real TPI = known_reals[14];
  real influx_fbp = known_reals[15];
  real outflux_pep = known_reals[16];
  real NAD = known_reals[17];
  real NADH = known_reals[18];
  real ATP = known_reals[19];
  real ADP = known_reals[20];
  real temperature = known_reals[21];
  real gas_constant = known_reals[22];
  real BPG = ode_metabolites[1];
  real DHAP = ode_metabolites[2];
  real F16bP = ode_metabolites[3];
  real GAP = ode_metabolites[4];
  real P2G = ode_metabolites[5];
  real P3G = ode_metabolites[6];
  real PEP = ode_metabolites[7];

  // thermodynamic parameters
  real ENO_delta_g = params[1];
  real FBA_delta_g = params[2];
  real GPM_delta_g = params[3];
  real PGK_delta_g = params[4];
  real TDH_delta_g = params[5];
  real TPI_delta_g = params[6];

  // kinetic parameters
  real ENO1_Kcat1 = params[7];
  real ENO1_Kcat2 = params[8];
  real ENO1_Ka = params[9];
  real ENO2_Kcat1 = params[10];
  real ENO2_Kcat2 = params[11];
  real ENO2_Ka = params[12];
  real FBA_Kcat1 = params[13];
  real FBA_Kcat2 = params[14];
  real FBA_Kia = params[15];
  real FBA_Kiq = params[16];
  real FBA_Ka = params[17];
  real FBA_Kp = params[18];
  real FBA_Kq = params[19];
  real GPM_Kcat1 = params[20];
  real GPM_Kcat2 = params[21];
  real GPM_Ka = params[22];
  real PGK_Kcat1 = params[23];
  real PGK_Kcat2 = params[24];
  real PGK_Kia = params[25];
  real PGK_Kib = params[26];
  real PGK_Kiq = params[27];
  real PGK_Ka = params[28];
  real PGK_Kb = params[29];
  real PGK_Kp = params[30];
  real PGK_Kq = params[31];
  real TDH1_Kcat1 = params[32];
  real TDH1_Kcat2 = params[33];
  real TDH1_Kia = params[34];
  real TDH1_Kib = params[35];
  real TDH1_Kiq = params[36];
  real TDH1_Ka = params[37];
  real TDH1_Kb = params[38];
  real TDH1_Kp = params[39];
  real TDH1_Kq = params[40];
  real TDH2_Kcat1 = params[41];
  real TDH2_Kcat2 = params[42];
  real TDH2_Kib = params[43];
  real TDH2_Kiq = params[44];
  real TDH2_Kia = params[45];
  real TDH2_Ka = params[46];
  real TDH2_Kb = params[47];
  real TDH2_Kp = params[48];
  real TDH2_Kq = params[49];
  real TDH3_Kcat1 = params[50];
  real TDH3_Kcat2 = params[51];
  real TDH3_Kia = params[52];
  real TDH3_Kib = params[53];
  real TDH3_Kiq = params[54];
  real TDH3_Ka = params[55];
  real TDH3_Kb = params[56];
  real TDH3_Kp = params[57];
  real TDH3_Kq = params[58];
  real TPI_Kcat1 = params[59];
  real TPI_Kcat2 = params[60];
  real TPI_Ka = params[61];
 
  // derived numbers
  real sum_PXG = P2G + P3G;
  real energy_charge = (ATP + ADP / 2) / sum_AXP;

  // haldane_relationships
  real FBA_Keq = get_Keq(FBA_delta_g, temperature, gas_constant);
  real TPI_Keq = get_Keq(TPI_delta_g, temperature, gas_constant);
  real TDH_Keq = get_Keq(TDH_delta_g, temperature, gas_constant);
  real PGK_Keq = get_Keq(PGK_delta_g, temperature, gas_constant);
  real GPM_Keq = get_Keq(GPM_delta_g, temperature, gas_constant);
  real ENO_Keq = get_Keq(ENO_delta_g, temperature, gas_constant);
  real FBA_Kip = get_Kip_ordered_unibi(FBA_Keq, FBA_Kia, FBA_Kq, FBA_Kcat1, FBA_Kcat2);
  real TPI_Kp = get_Kp_uniuni(TPI_Kcat1, TPI_Kcat2, TPI_Keq, TPI_Ka);
  real TDH1_Kip = get_Kip_ordered_bibi(TDH_Keq, TDH1_Ka, TDH1_Kq, TDH1_Kib, TDH1_Kiq, TDH1_Kcat1, TDH1_Kcat2);
  real TDH2_Kip = get_Kip_ordered_bibi(TDH_Keq, TDH2_Ka, TDH2_Kq, TDH2_Kib, TDH2_Kiq, TDH2_Kcat1, TDH2_Kcat2);
  real TDH3_Kip = get_Kip_ordered_bibi(TDH_Keq, TDH3_Ka, TDH3_Kq, TDH3_Kib, TDH3_Kiq, TDH3_Kcat1, TDH3_Kcat2);
  real PGK_Kip = get_Kip_ordered_bibi(PGK_Keq, PGK_Ka, PGK_Kq, PGK_Kib, PGK_Kiq, PGK_Kcat1, PGK_Kcat2);
  real GPM_Kp = get_Kp_uniuni(GPM_Kcat1, GPM_Kcat2, GPM_Keq, ENO1_Ka);
  real ENO1_Kp = get_Kp_uniuni(ENO1_Kcat1, ENO1_Kcat2, ENO_Keq, ENO1_Ka);
  real ENO2_Kp = get_Kp_uniuni(ENO2_Kcat1, ENO2_Kcat2, ENO_Keq, ENO2_Ka);

  // fluxes
  real flux_ENO1 = uniuni(P2G, PEP, ENO1 * ENO1_Kcat1, ENO1 * ENO1_Kcat2, ENO1_Ka, ENO1_Kp);
  real flux_ENO2 = uniuni(P2G, PEP, ENO2 * ENO2_Kcat1, ENO2 * ENO2_Kcat2, ENO2_Ka, ENO2_Kp);
  real flux_FBA = ordered_unibi(F16bP, GAP, DHAP,
                                FBA * FBA_Kcat1, FBA * FBA_Kcat2,
                                FBA_Ka, FBA_Kp, FBA_Kq,
                                FBA_Kia, FBA_Kip, FBA_Kiq);
  real flux_GPM = uniuni(P3G, P2G, GPM * GPM_Kcat1, GPM * GPM_Kcat2, GPM_Ka, GPM_Kp);
  real flux_PGK = ordered_bibi(BPG, ADP, P3G, ATP,
                               PGK * PGK_Kcat1, PGK * PGK_Kcat2,
                               PGK_Ka, PGK_Kb, PGK_Kp, PGK_Kq,
                               PGK_Kia, PGK_Kib, PGK_Kip, PGK_Kiq);
  real flux_TDH1 = ordered_bibi(GAP, NAD, BPG, NADH,
                                TDH1 * TDH1_Kcat1, TDH1 * TDH1_Kcat2,
                                TDH1_Ka, TDH1_Kb, TDH1_Kp, TDH1_Kq,
                                TDH1_Kia, TDH1_Kib, TDH1_Kip, TDH1_Kiq);
  real flux_TDH2 = ordered_bibi(GAP, NAD, BPG, NADH,
                                TDH2 * TDH2_Kcat1, TDH2 * TDH2_Kcat2,
                                TDH2_Ka, TDH2_Kb, TDH2_Kp, TDH2_Kq,
                                TDH2_Kia, TDH2_Kib, TDH2_Kip, TDH2_Kiq);
  real flux_TDH3 = ordered_bibi(GAP, NAD, BPG, NADH,
                                TDH3 * TDH3_Kcat1, TDH3 * TDH3_Kcat2,
                                TDH3_Ka, TDH3_Kb, TDH3_Kp, TDH3_Kq,
                                TDH3_Kia, TDH3_Kib, TDH3_Kip, TDH3_Kiq);
  real flux_TPI = uniuni(DHAP, GAP, TPI * TPI_Kcat1, TPI * TPI_Kcat2, TPI_Ka, TPI_Kp);
  
  return {influx_fbp, flux_ENO1, flux_ENO2, flux_FBA, flux_GPM, flux_PGK,
          flux_TDH1, flux_TDH2, flux_TDH3, flux_TPI, outflux_pep};
}

real[] get_odes(real[] fluxes){
  real influx_fbp = fluxes[1];
  real ENO1 = fluxes[2];
  real ENO2 = fluxes[3];
  real FBA = fluxes[4];
  real GPM = fluxes[5];
  real PGK = fluxes[6];
  real TDH1 = fluxes[7];
  real TDH2 = fluxes[8];
  real TDH3 = fluxes[9];
  real TPI = fluxes[10];
  real outflux_pep = fluxes[11];
  return {
    TDH1+TDH3+TDH2-PGK,   // BPG
    FBA-TPI,   // DHAP
    influx_fbp - FBA,   // F16bP
    FBA+TPI-TDH1-TDH3-TDH2,   // GAP
    GPM-ENO1-ENO2,   // P2G
    PGK-GPM,   // P3G
    ENO1+ENO2 - outflux_pep   // PEP
  };
}

real[] steady_state_equation(real t,
                             real[] ode_metabolites,
                             real[] params,
                             real[] known_reals,
                             int[] known_ints){
  for (m in 1:size(ode_metabolites)){
    if (ode_metabolites[m] < 0){
      print(ode_metabolites);
      print(get_fluxes(ode_metabolites, params, known_reals));
      reject("Metabolite ", m, " is ", ode_metabolites[m], " but should be greater than zero");
    }
  }

  return get_odes(get_fluxes(ode_metabolites, params, known_reals));
}
