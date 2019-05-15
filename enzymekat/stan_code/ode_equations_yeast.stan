real[] get_fluxes(real[] ode_metabolites,
                  real[] parameters,
                  real[] known_reals){
  real cell = known_reals[1];
  real Keq_ADH = known_reals[2];
  real Keq_ENO = known_reals[3];
  real Keq_HXK = known_reals[4];
  real Keq_PYK = known_reals[5];
  real Keq_TDH = known_reals[6];
  real NA = known_reals[7];
  real sum_AXP = known_reals[8];
  real sum_NAD = known_reals[9];
  real sum_UXP = known_reals[10];
  real volume = known_reals[11];
  real ENO1 = known_reals[12];
  real ENO2 = known_reals[13];
  real FBA1 = known_reals[14];
  real GPM1 = known_reals[15];
  real PGK1 = known_reals[16];
  real TDH1 = known_reals[17];
  real TDH2 = known_reals[18];
  real TDH3 = known_reals[19];
  real TPI1 = known_reals[20];
  real influx_fbp = known_reals[21];
  real outflux_pep = known_reals[22];
  real NAD = known_reals[23];
  real NADH = known_reals[24];
  real ATP = known_reals[25];
  real ADP = known_reals[26];
  real BPG = ode_metabolites[1];
  real DHAP = ode_metabolites[2];
  real F16bP = ode_metabolites[3];
  real GAP = ode_metabolites[4];
  real P2G = ode_metabolites[5];
  real P3G = ode_metabolites[6];
  real PEP = ode_metabolites[7];
  real FBA_delta_g = parameters[1];
  real TPI_delta_g = parameters[6];
  real TDH1_delta_g = parameters[10];
  real TDH2_delta_g = parameters[16];
  real TDH3_delta_g = parameters[22];
  real PGK_delta_g = parameters[28];
  real GPM_delta_g = parameters[34];
  real ENO1_delta_g = parameters[38];
  real ENO2_delta_g = parameters[42];
  real FBA_Ka = parameters[2];
  real FBA_Kp = parameters[3];
  real FBA_V1 = parameters[4];
  real FBA_V2 = parameters[5];
  real TPI_V1 = parameters[7];
  real TPI_V2 = parameters[8];
  real TPI_Ka = parameters[9];
  real TDH1_Kb = parameters[11];
  real TDH1_Kia = parameters[12];
  real TDH1_Kp = parameters[13];
  real TDH1_V1 = parameters[14];
  real TDH1_V2 = parameters[15];
  real TDH2_Kb = parameters[17];
  real TDH2_Kia = parameters[18];
  real TDH2_Kp = parameters[19];
  real TDH2_V1 = parameters[20];
  real TDH2_V2 = parameters[21];
  real TDH3_Kb = parameters[23];
  real TDH3_Kia = parameters[24];
  real TDH3_Kp = parameters[25];
  real TDH3_V1 = parameters[26];
  real TDH3_V2 = parameters[27];
  real PGK_Kb = parameters[29];
  real PGK_Kia = parameters[30];
  real PGK_Kp = parameters[31];
  real PGK_V1 = parameters[32];
  real PGK_V2 = parameters[33];
  real GPM_V1 = parameters[35];
  real GPM_V2 = parameters[36];
  real GPM_Ka = parameters[37];
  real ENO1_V1 = parameters[39];
  real ENO1_V2 = parameters[40];
  real ENO1_Ka = parameters[41];
  real ENO2_V1 = parameters[43];
  real ENO2_V2 = parameters[44];
  real ENO2_Ka = parameters[45];

  // derived numbers
  real sum_PXG = P2G + P3G;
  real energy_charge = (ATP + ADP / 2) / sum_AXP;

  // haldane_relationships
  real FBA_Keq = get_Keq(FBA_delta_g, T, R);
  real TPI_Keq = get_Keq(TPI_delta_g, T, R);
  real TDH1_Keq = get_Keq(TDH1_delta_g, T, R);
  real TDH2_Keq = get_Keq(TDH2_delta_g, T, R);
  real TDH3_Keq = get_Keq(TDH3_delta_g, T, R);
  real PGK_Keq = get_Keq(PGK_delta_g, T, R);
  real GPM_Keq = get_Keq(GPM_delta_g, T, R);
  real ENO1_Keq = get_Keq(ENO1_delta_g, T, R);
  real ENO2_Keq = get_Keq(ENO2_delta_g, T, R);
  real FBA_Kip = get_Kip_ordered_unibi(FBA_Keq, FBA_Kia, FBA_Kq, FBA_V1, FBA_V2);
  real TPI_Kp = get_Kp_uniuni(TPI_V1, TPI_V2, TPI_Keq, TPI_Ka);
  real TDH1_Kip = get_Kip_ordered_bibi(TDH1_Kb, TDH1_Keq, TDH1_Kia, TDH1_Kp, TDH1_V1, TDH1_V2);
  real TDH2_Kip = get_Kip_ordered_bibi(TDH2_Kb, TDH2_Keq, TDH2_Kia, TDH2_Kp, TDH2_V1, TDH2_V2);
  real TDH3_Kip = get_Kip_ordered_bibi(TDH3_Kb, TDH3_Keq, TDH3_Kia, TDH3_Kp, TDH3_V1, TDH3_V2);
  real PGK_Kip = get_Kip_ordered_bibi(PGK_Kb, PGK_Keq, PGK_Kia, PGK_Kp, PGK_V1, PGK_V2);
  real GPM_Kp = get_Kp_uniuni(GPM_V1, GPM_V2, GPM_Keq, ENO1_Ka);
  real ENO1_Kp = get_Kp_uniuni(ENO1_V1, ENO1_V2, ENO1_Keq, ENO1_Ka);
  real ENO2_Kp = get_Kp_uniuni(ENO2_V1, ENO2_V2, ENO2_Keq, ENO2_Ka);

  return {
    influx_fbp,
    uniuni(ENO1_V1, ENO1_V2, ENO1_A, ENO1_P, ENO1_Ka, ENO1_Kp), // ENO1
    uniuni(ENO2_V1, ENO2_V2, ENO2_A, ENO2_P, ENO2_Ka, ENO2_Kp), // ENO2
    ordered_unibi(FBA_A, FBA_Kip, FBA_Ka, FBA_Kq, FBA_Kp, FBA_Q, FBA_V1, FBA_V2), // FBA
    uniuni(GPM_V1, GPM_V2, GPM_A, GPM_P, GPM_Kq, GPM_Kp),  // GPM
    ordered_bibi(PGK_Kib, PGK_Kiq, PGK_B, PGK_Kia, PGK_Kip, PGK_P, PGK_Kb, PGK_Kq, PGK_Q, PGK_Ka, PGK_Kp, PGK_V1, PGK_V2),  // PGK
    ordered_bibi(TDH1_Kib, TDH1_Kiq, TDH1_B, TDH1_Kia, TDH1_Kip, TDH1_P, TDH1_Kb, TDH1_Kq, TDH1_Q, TDH1_Ka, TDH1_Kp, TDH1_V1, TDH1_V2),  // TDH1
    ordered_bibi(TDH2_Kib, TDH2_Kiq, TDH2_B, TDH2_Kia, TDH2_Kip, TDH2_P, TDH2_Kb, TDH2_Kq, TDH2_Q, TDH2_Ka, TDH2_Kp, TDH2_V1, TDH2_V2),  // TDH2
    ordered_bibi(TDH3_Kib, TDH3_Kiq, TDH3_B, TDH3_Kia, TDH3_Kip, TDH3_P, TDH3_Kb, TDH3_Kq, TDH3_Q, TDH3_Ka, TDH3_Kp, TDH3_V1, TDH3_V2),  // TDH3
    uniuni(TPI_V1, TPI_V2, TPI_A, TPI_P, TPI_Kq, TPI_Kp),  // TPI
    outflux_pep
  };
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
                             real[] parameters,
                             real[] known_reals,
                             int[] known_ints){
  for (m in 1:size(ode_metabolites)){
    if (ode_metabolites[m] < 0){
      reject("Metabolite ", m, " is ", ode_metabolites[m], " but should be greater than zero");
    }
  }

  return get_odes(get_fluxes(ode_metabolites, parameters, known_reals));
}
