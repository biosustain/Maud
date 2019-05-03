real[] get_derived_quantities(real[] ode_metabolites, real[] known_reals){
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
  real ADP = ode_metabolites[1];
  real ATP = ode_metabolites[2];
  real BPG = ode_metabolites[3];
  real DHAP = ode_metabolites[4];
  real F16bP = ode_metabolites[5];
  real GAP = ode_metabolites[6];
  real NAD = ode_metabolites[7];
  real P2G = ode_metabolites[8];
  real P3G = ode_metabolites[9];
  real PEP = ode_metabolites[10];
  real sum_PXG = P2G + P3G;
  real energy_charge = (ATP + ADP / 2) / sum_AXP;
  real NADH = max({sum_NAD - NAD, 0});
  return {sum_PXG, energy_charge, NADH};
}

real[] get_fluxes(real[] ode_metabolites,
                    real[] kinetic_parameters,
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
  real ADP = ode_metabolites[1];
  real ATP = ode_metabolites[2];
  real BPG = ode_metabolites[3];
  real DHAP = ode_metabolites[4];
  real F16bP = ode_metabolites[5];
  real GAP = ode_metabolites[6];
  real NAD = ode_metabolites[7];
  real P2G = ode_metabolites[8];
  real P3G = ode_metabolites[9];
  real PEP = ode_metabolites[10];
  real ENO1_Kp2g = kinetic_parameters[1];
  real ENO1_Kpep = kinetic_parameters[2];
  real ENO1_kcat = kinetic_parameters[3];
  real FBA_Kdhap = kinetic_parameters[4];
  real FBA_Keq = kinetic_parameters[5];
  real FBA_Kf16bp = kinetic_parameters[6];
  real FBA_Kgap = kinetic_parameters[7];
  real FBA_Kigap = kinetic_parameters[8];
  real FBA_kcat = kinetic_parameters[9];
  real GPM_Keq = kinetic_parameters[10];
  real GPM_Kp2g = kinetic_parameters[11];
  real GPM_Kp3g = kinetic_parameters[12];
  real GPM_kcat = kinetic_parameters[13];
  real PGK_Kadp = kinetic_parameters[14];
  real PGK_Katp = kinetic_parameters[15];
  real PGK_Kbpg = kinetic_parameters[16];
  real PGK_Keq = kinetic_parameters[17];
  real PGK_Kp3g = kinetic_parameters[18];
  real PGK_kcat = kinetic_parameters[19];
  real PGK_nHadp = kinetic_parameters[20];
  real TDH1_Kbpg = kinetic_parameters[21];
  real TDH1_Kgap = kinetic_parameters[22];
  real TDH1_Knad = kinetic_parameters[23];
  real TDH1_Knadh = kinetic_parameters[24];
  real TDH1_kcat = kinetic_parameters[25];
  real TPI_Kdhap = kinetic_parameters[26];
  real TPI_Keq = kinetic_parameters[27];
  real TPI_Kgap = kinetic_parameters[28];
  real TPI_Kigap = kinetic_parameters[29];
  real TPI_kcat = kinetic_parameters[30];
  real TDH3_Kbpg = kinetic_parameters[31];
  real TDH3_Kgap = kinetic_parameters[32];
  real TDH3_Knad = kinetic_parameters[33];
  real TDH3_Knadh = kinetic_parameters[34];
  real TDH3_kcat = kinetic_parameters[35];
  real TDH2_Kbpg = kinetic_parameters[36];
  real TDH2_Kgap = kinetic_parameters[37];
  real TDH2_Knad = kinetic_parameters[38];
  real TDH2_Knadh = kinetic_parameters[39];
  real TDH2_kcat = kinetic_parameters[40];
  real ENO2_Kp2g = kinetic_parameters[41];
  real ENO2_Kpep = kinetic_parameters[42];
  real ENO2_kcat = kinetic_parameters[43];
  real derived_quantities[3] = get_derived_quantities(ode_metabolites, known_reals);
  real sum_PXG = derived_quantities[1];
  real energy_charge = derived_quantities[2];
  real NADH = derived_quantities[3];
  return {cell*reversible_michaelis_menten(P2G,PEP,ENO1*ENO1_kcat,ENO1_Kp2g,ENO1_Kpep,Keq_ENO),
      cell*ordered_uni_bi(F16bP,DHAP,GAP,FBA1*FBA_kcat,FBA_Kf16bp,FBA_Kdhap,FBA_Kgap,FBA_Kigap,FBA_Keq),
      cell*reversible_michaelis_menten(P3G,P2G,GPM1*GPM_kcat,GPM_Kp3g,GPM_Kp2g,GPM_Keq),
      cell*phosphoglycerate_kinase_kinetics(BPG,ADP,P3G,ATP,PGK1*PGK_kcat,PGK_Kbpg,PGK_Kadp,PGK_Kp3g,PGK_Katp,PGK_Keq,PGK_nHadp),
      cell*two_noncompeting_couples(GAP,NAD,BPG,NADH,TDH1*TDH1_kcat,TDH1_Kgap,TDH1_Knad,TDH1_Kbpg,TDH1_Knadh,Keq_TDH),
      cell*triphosphate_isomerase_kinetics(DHAP,GAP,TPI1*TPI_kcat,TPI_Kdhap,TPI_Kgap,TPI_Kigap, TPI_Keq),
      cell*two_noncompeting_couples(GAP,NAD,BPG,NADH,TDH3*TDH3_kcat,TDH3_Kgap,TDH3_Knad,TDH3_Kbpg,TDH3_Knadh,Keq_TDH),
      cell*two_noncompeting_couples(GAP,NAD,BPG,NADH,TDH1*TDH2_kcat,TDH2_Kgap,TDH2_Knad,TDH2_Kbpg,TDH2_Knadh,Keq_TDH),
      cell*reversible_michaelis_menten(P2G,PEP,ENO2*ENO2_kcat,ENO2_Kp2g,ENO2_Kpep,Keq_ENO)};
}

real[] get_odes(real[] fluxes){
  real ENO1 = fluxes[1];
  real FBA = fluxes[2];
  real GPM = fluxes[3];
  real PGK = fluxes[4];
  real TDH1 = fluxes[5];
  real TPI = fluxes[6];
  real TDH3 = fluxes[7];
  real TDH2 = fluxes[8];
  real ENO2 = fluxes[9];
  return {
    (-1.0*PGK)/1.0,  // ADP
      (1.0*PGK)/1.0,  // ATP
      (1.0*TDH1+1.0*TDH3+1.0*TDH2-1.0*PGK)/1.0,  // BPG
      (1.0*FBA-1.0*TPI)/1.0,  // DHAP
      (-1.0*FBA)/1.0,  // F16bP
      (1.0*FBA+1.0*TPI-1.0*TDH1-1.0*TDH3-1.0*TDH2)/1.0,  // GAP
      (-1.0*TDH1-1.0*TDH3-1.0*TDH2)/1.0,  // NAD
      (1.0*GPM-1.0*ENO1-1.0*ENO2)/1.0,  // P2G
      (1.0*PGK-1.0*GPM)/1.0,  // P3G
      (1.0*ENO1+1.0*ENO2)/1.0  // PEP
      };
}

real[] steady_state_equation(real t,
                             real[] ode_metabolites,
                             real[] kinetic_parameters,
                             real[] known_reals,
                             int[] known_ints){
  for (m in 1:size(ode_metabolites)){
    if (ode_metabolites[m] < 0){
      reject("Metabolite ", m, " is ", ode_metabolites[m], " but should be greater than zero");
    }
  }

  return get_odes(get_fluxes(ode_metabolites, kinetic_parameters, known_reals));
}
