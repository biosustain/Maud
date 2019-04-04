real Function_for_Glucose_transport(real GlcE,real GlcI,real K1Glc,real Vm1,real Vt,real afac,real tot_cell){
  return tot_cell / Vt * Vm1 * (GlcE - GlcI) / (K1Glc + GlcE + GlcI + afac * GlcE * GlcI / K1Glc);
}

real Function_for_Hexokinase(real ADPg,real ATPg,real Glc6P,real GlcI,real K2ADPg,real K2ATPg,real K2Glc6P,real K2GlcI,real Vm2,real Vt,real tot_cell){
  return tot_cell / Vt * Vm2 * GlcI * ATPg / (K2ATPg * K2GlcI * (1 + Glc6P / K2Glc6P + GlcI / K2GlcI) * (1 + ATPg / K2ATPg + ADPg / K2ADPg));
}

real Function_for_Glucose_phosphate_isomerase(real Fru6P,real Glc6P,real K3Fru6P,real K3Glc6P,real Vm3,real Vt,real glycosome,real tot_cell){
  return tot_cell / Vt * Vm3 * (Glc6P / K3Glc6P - Fru6P / K3Fru6P) / (1 + Glc6P / K3Glc6P + Fru6P / K3Fru6P) / glycosome;
}

real Function_for_Phosphofructokinase(real ATPg,real Fru16BP,real Fru6P,real K4ATPg,real K4Fru6P,real K4i1Fru16BP,real K4i2Fru16BP,real Vm4,real Vt,real glycosome,real tot_cell){
  return tot_cell / Vt * K4i1Fru16BP * Vm4 * Fru6P * ATPg / (K4ATPg * K4Fru6P * (K4i1Fru16BP + Fru16BP) * (1 + Fru16BP / K4i2Fru16BP + Fru6P / K4Fru6P) * (1 + ATPg / K4ATPg)) / glycosome;
}

real Function_for_Aldolase(real ADPg,real ATPg,real DHAPg,real Fru16BP,real GAP,real K5DHAP,real K5GAP,real K5GAPi,real Vm5f,real Vm5r,real Vt,real sumAg,real tot_cell){
  return tot_cell / Vt * (Vm5f * Fru16BP / (0.009 * (1 + ATPg / 0.68 + ADPg / 1.51 + (sumAg - (ATPg + ADPg)) / 3.65)) - Vm5r * GAP * DHAPg / (K5DHAP * K5GAP)) / (1 + GAP / K5GAP + DHAPg / K5DHAP + GAP * DHAPg / (K5DHAP * K5GAP) + Fru16BP / (0.009 * (1 + ATPg / 0.68 + ADPg / 1.51 + (sumAg - (ATPg + ADPg)) / 3.65)) + Fru16BP * GAP / (K5GAPi * 0.009 * (1 + ATPg / 0.68 + ADPg / 1.51 + (sumAg - (ATPg + ADPg)) / 3.65)));
}

real Function_for_Triosephosphate_isomerase(real DHAPg,real GAP,real K6DHAPg,real K6GAP,real TPIact,real Vm6,real Vt,real tot_cell){
  return tot_cell / Vt * TPIact * Vm6 * (DHAPg / K6DHAPg - 5.7 * GAP / K6GAP) / (1 + GAP / K6GAP + DHAPg / K6DHAPg);
}

real Function_for_Glyceraldehyde_3_phosphate_dehydrogenase(real BPGA13,real GAP,real K7BPGA13,real K7GAP,real K7NAD,real K7NADH,real NAD,real NADH,real Vm7,real Vm7f,real Vm7r,real Vt,real glycosome,real tot_cell){
  return tot_cell / Vt * Vm7 * (Vm7f * (GAP * (NAD / K7GAP / K7NAD) - Vm7r / Vm7f * (BPGA13 * NADH / K7BPGA13 / K7NADH)) / ((1 + GAP / K7GAP + BPGA13 / K7BPGA13) * (1 + NAD / K7NAD + NADH / K7NADH))) / glycosome;
}

real Function_for_Glycerol_3_phosphate_dehydrogenase(real DHAPg,real Gly3Pg,real K8DHAPg,real K8Gly3Pg,real K8NAD,real K8NADH,real NAD,real NADH,real Vm8,real Vm8f,real Vm8r,real Vt,real tot_cell){
  return tot_cell / Vt * Vm8 * Vm8f * (NADH * DHAPg / (K8DHAPg * K8NADH) - Vm8r * NAD * Gly3Pg / (K8Gly3Pg * K8NAD * Vm8f)) / ((1 + NAD / K8NAD + NADH / K8NADH) * (1 + DHAPg / K8DHAPg + Gly3Pg / K8Gly3Pg));
}

real Function_for_Phosphoglycerate_kinase(real ADPg,real ATPg,real BPGA13,real K11ADPg,real K11ATPg,real K11BPGA13,real K11PGA3,real PGAg,real Vm11,real Vm11f,real Vm11r,real Vt,real tot_cell){
  return tot_cell / Vt * Vm11 * Vm11f * (-Vm11r * PGAg * ATPg / (K11ATPg * K11PGA3 * Vm11f) + BPGA13 * ADPg / (K11ADPg * K11BPGA13)) / ((1 + BPGA13 / K11BPGA13 + PGAg / K11PGA3) * (1 + ATPg / K11ATPg + ADPg / K11ADPg));
}

real Function_for_Pyruvate_transport(real K10Pyr,real Pyr,real Vm10,real Vt,real tot_cell){
  return tot_cell / Vt * Vm10 * Pyr / K10Pyr / (1 + Pyr / K10Pyr);
}

real Function_for_Pyruvate_kinase(real ADPc,real ATPc,real K12ADP,real PEPc,real Vm12,real Vt,real n12,real tot_cell){
  return tot_cell / Vt * Vm12 * pow(PEPc / (0.34 * (1 + ADPc / 0.57 + ATPc / 0.64)), n12) * ADPc / K12ADP / ((1 + pow(PEPc / (0.34 * (1 + ADPc / 0.57 + ATPc / 0.64)), n12)) * (1 + ADPc / K12ADP));
}

real Function_for_ATPase(real ADPc,real ATPc,real K13,real Vt,real cytosol,real tot_cell){
  return tot_cell / Vt * K13 * ATPc / ADPc / cytosol;
}

real Function_for_Glycerol_kinase(real ADPg,real ATPg,real Gly,real Gly3Pg,real K14ADPg,real K14ATPg,real K14Gly,real K14Gly3Pg,real Vm14,real Vm14f,real Vm14r,real Vt,real tot_cell){
  return tot_cell / Vt * Vm14 * (Vm14f * ADPg * Gly3Pg / (K14ADPg * K14Gly3Pg) - Gly * Vm14r * ATPg / (K14ATPg * K14Gly)) / ((1 + Gly / K14Gly + Gly3Pg / K14Gly3Pg) * (1 + ATPg / K14ATPg + ADPg / K14ADPg));
}

real Function_for_Glycerol_3_phosphate_oxidase(real Gly3Pc,real K9Gly3Pc,real Vm9,real Vt,real tot_cell){
  return tot_cell / Vt * Vm9 * Gly3Pc / (K9Gly3Pc * 1 + Gly3Pc);
}

real[] get_derived_quantities(vector ode_metabolites, real[] known_reals){
  real tot_cell = known_reals[1];
  real glycosome = known_reals[2];
  real cytosol = known_reals[3];
  real extracellular = known_reals[4];
  real Vt = known_reals[5];
  real TPIact = known_reals[6];
  real sumAc = known_reals[7];
  real sumAg = known_reals[8];
  real KeqAK = known_reals[9];
  real Keq_anti = known_reals[10];
  real sumc4 = known_reals[11];
  real sumc5 = known_reals[12];
  real Keq_PGM = known_reals[13];
  real Keq_ENO = known_reals[14];
  real PyrE = known_reals[15];
  real Gly = known_reals[16];
  real GlcE = known_reals[17];
  real GlcI = ode_metabolites[1];
  real Pg = ode_metabolites[2];
  real Glc6P = ode_metabolites[3];
  real Fru6P = ode_metabolites[4];
  real Fru16BP = ode_metabolites[5];
  real DHAP = ode_metabolites[6];
  real GAP = ode_metabolites[7];
  real NAD = ode_metabolites[8];
  real BPGA13 = ode_metabolites[9];
  real NADH = ode_metabolites[10];
  real Pyr = ode_metabolites[11];
  real Nb = ode_metabolites[12];
  real Pc = ode_metabolites[13];
  real Vc = cytosol * Vt / tot_cell;
  real Vg = glycosome * Vt / tot_cell;
  real ATPc = (Pc * (1 - 4 * KeqAK) - sumAc + pow(pow(sumAc - (1 - 4 * KeqAK) * Pc, 2) + 4 * (1 - 4 * KeqAK) * KeqAK * pow(Pc, 2), 0.5)) / (2 * (1 - 4 * KeqAK));
  real ADPc = Pc - 2 * ATPc;
  real ATPg = (Pg * (1 - 4 * KeqAK) - sumAg + pow(pow(sumAg - (1 - 4 * KeqAK) * Pg, 2) + 4 * (1 - 4 * KeqAK) * KeqAK * pow(Pg, 2), 0.5)) / (2 * (1 - 4 * KeqAK));
  real ADPg = Pg - 2 * ATPg;
  real DHAPc = sumc5 * (1 + Vc / Vg) * DHAP / (sumc4 + sumc5 * Vc / Vg - (BPGA13 + 2 * Fru16BP + Fru6P + GAP + Glc6P + Pg));
  real DHAPg = (DHAP * Vt - DHAPc * Vc) / Vg;
  real Gly3Pc = sumc5 - DHAPc;
  real Gly3Pg = Gly3Pc * DHAPg / (Keq_anti * DHAPc);
  real Gly3P = (Gly3Pc * cytosol + Gly3Pg * glycosome) / tot_cell;
  real PGAg = Nb * (1 + Vc / Vg) / (1 + (1 + Keq_PGM + Keq_PGM * Keq_ENO) * Vc / Vg);
  real PEPc = Keq_ENO * Keq_PGM * PGAg;
  return {Vc, Vg, ATPc, ADPc, ATPg, ADPg, DHAPc, DHAPg, Gly3Pc, Gly3Pg, Gly3P, PGAg, PEPc};
}

vector get_kinetics(vector ode_metabolites,
                    vector kinetic_parameters,
                    real[] known_reals){
  real tot_cell = known_reals[1];
  real glycosome = known_reals[2];
  real cytosol = known_reals[3];
  real extracellular = known_reals[4];
  real Vt = known_reals[5];
  real TPIact = known_reals[6];
  real sumAc = known_reals[7];
  real sumAg = known_reals[8];
  real KeqAK = known_reals[9];
  real Keq_anti = known_reals[10];
  real sumc4 = known_reals[11];
  real sumc5 = known_reals[12];
  real Keq_PGM = known_reals[13];
  real Keq_ENO = known_reals[14];
  real PyrE = known_reals[15];
  real Gly = known_reals[16];
  real GlcE = known_reals[17];
  real GlcI = ode_metabolites[1];
  real Pg = ode_metabolites[2];
  real Glc6P = ode_metabolites[3];
  real Fru6P = ode_metabolites[4];
  real Fru16BP = ode_metabolites[5];
  real DHAP = ode_metabolites[6];
  real GAP = ode_metabolites[7];
  real NAD = ode_metabolites[8];
  real BPGA13 = ode_metabolites[9];
  real NADH = ode_metabolites[10];
  real Pyr = ode_metabolites[11];
  real Nb = ode_metabolites[12];
  real Pc = ode_metabolites[13];
  real K1Glc = kinetic_parameters[1];
  real Vm1 = kinetic_parameters[2];
  real afac = kinetic_parameters[3];
  real K2ADPg = kinetic_parameters[4];
  real K2ATPg = kinetic_parameters[5];
  real K2Glc6P = kinetic_parameters[6];
  real K2GlcI = kinetic_parameters[7];
  real Vm2 = kinetic_parameters[8];
  real K3Fru6P = kinetic_parameters[9];
  real K3Glc6P = kinetic_parameters[10];
  real Vm3 = kinetic_parameters[11];
  real K4ATPg = kinetic_parameters[12];
  real K4Fru6P = kinetic_parameters[13];
  real K4i1Fru16BP = kinetic_parameters[14];
  real K4i2Fru16BP = kinetic_parameters[15];
  real Vm4 = kinetic_parameters[16];
  real K5DHAP = kinetic_parameters[17];
  real K5GAP = kinetic_parameters[18];
  real K5GAPi = kinetic_parameters[19];
  real Vm5f = kinetic_parameters[20];
  real Vm5r = kinetic_parameters[21];
  real K6DHAPg = kinetic_parameters[22];
  real K6GAP = kinetic_parameters[23];
  real Vm6 = kinetic_parameters[24];
  real K7BPGA13 = kinetic_parameters[25];
  real K7GAP = kinetic_parameters[26];
  real K7NAD = kinetic_parameters[27];
  real K7NADH = kinetic_parameters[28];
  real Vm7 = kinetic_parameters[29];
  real Vm7f = kinetic_parameters[30];
  real Vm7r = kinetic_parameters[31];
  real K8DHAPg = kinetic_parameters[32];
  real K8Gly3Pg = kinetic_parameters[33];
  real K8NAD = kinetic_parameters[34];
  real K8NADH = kinetic_parameters[35];
  real Vm8 = kinetic_parameters[36];
  real Vm8f = kinetic_parameters[37];
  real Vm8r = kinetic_parameters[38];
  real K9Gly3Pc = kinetic_parameters[39];
  real Vm9 = kinetic_parameters[40];
  real K10Pyr = kinetic_parameters[41];
  real Vm10 = kinetic_parameters[42];
  real K11ADPg = kinetic_parameters[43];
  real K11ATPg = kinetic_parameters[44];
  real K11BPGA13 = kinetic_parameters[45];
  real K11PGA3 = kinetic_parameters[46];
  real Vm11 = kinetic_parameters[47];
  real Vm11f = kinetic_parameters[48];
  real Vm11r = kinetic_parameters[49];
  real K12ADP = kinetic_parameters[50];
  real Vm12 = kinetic_parameters[51];
  real n12 = kinetic_parameters[52];
  real K13 = kinetic_parameters[53];
  real K14ADPg = kinetic_parameters[54];
  real K14ATPg = kinetic_parameters[55];
  real K14Gly = kinetic_parameters[56];
  real K14Gly3Pg = kinetic_parameters[57];
  real Vm14 = kinetic_parameters[58];
  real Vm14f = kinetic_parameters[59];
  real Vm14r = kinetic_parameters[60];
  real derived_quantities[13] = get_derived_quantities(ode_metabolites, known_reals);
  real Vc = derived_quantities[1];
  real Vg = derived_quantities[2];
  real ATPc = derived_quantities[3];
  real ADPc = derived_quantities[4];
  real ATPg = derived_quantities[5];
  real ADPg = derived_quantities[6];
  real DHAPc = derived_quantities[7];
  real DHAPg = derived_quantities[8];
  real Gly3Pc = derived_quantities[9];
  real Gly3Pg = derived_quantities[10];
  real Gly3P = derived_quantities[11];
  real PGAg = derived_quantities[12];
  real PEPc = derived_quantities[13];
  real vGlcTr = Function_for_Glucose_transport(GlcE, GlcI, K1Glc, Vm1, Vt, afac, tot_cell);
  real vHK = Function_for_Hexokinase(ADPg, ATPg, Glc6P, GlcI, K2ADPg, K2ATPg, K2Glc6P, K2GlcI, Vm2, Vt, tot_cell);
  real vPGI = glycosome * Function_for_Glucose_phosphate_isomerase(Fru6P, Glc6P, K3Fru6P, K3Glc6P, Vm3, Vt, glycosome, tot_cell);
  real vPFK = glycosome * Function_for_Phosphofructokinase(ATPg, Fru16BP, Fru6P, K4ATPg, K4Fru6P, K4i1Fru16BP, K4i2Fru16BP, Vm4, Vt, glycosome, tot_cell);
  real vALD = Function_for_Aldolase(ADPg, ATPg, DHAPg, Fru16BP, GAP, K5DHAP, K5GAP, K5GAPi, Vm5f, Vm5r, Vt, sumAg, tot_cell);
  real vTPI = Function_for_Triosephosphate_isomerase(DHAPg, GAP, K6DHAPg, K6GAP, TPIact, Vm6, Vt, tot_cell);
  real vGAPdh = glycosome * Function_for_Glyceraldehyde_3_phosphate_dehydrogenase(BPGA13, GAP, K7BPGA13, K7GAP, K7NAD, K7NADH, NAD, NADH, Vm7, Vm7f, Vm7r, Vt, glycosome, tot_cell);
  real vGDH = Function_for_Glycerol_3_phosphate_dehydrogenase(DHAPg, Gly3Pg, K8DHAPg, K8Gly3Pg, K8NAD, K8NADH, NAD, NADH, Vm8, Vm8f, Vm8r, Vt, tot_cell);
  real vGPO = Function_for_Glycerol_3_phosphate_oxidase(Gly3Pc, K9Gly3Pc, Vm9, Vt, tot_cell);
  real vPyrTr = Function_for_Pyruvate_transport(K10Pyr, Pyr, Vm10, Vt, tot_cell);
  real vPGK = Function_for_Phosphoglycerate_kinase(ADPg, ATPg, BPGA13, K11ADPg, K11ATPg, K11BPGA13, K11PGA3, PGAg, Vm11, Vm11f, Vm11r, Vt, tot_cell);
  real vPK = Function_for_Pyruvate_kinase(ADPc, ATPc, K12ADP, PEPc, Vm12, Vt, n12, tot_cell);
  real vATPase = cytosol * Function_for_ATPase(ADPc, ATPc, K13, Vt, cytosol, tot_cell);
  real vGlyK = Function_for_Glycerol_kinase(ADPg, ATPg, Gly, Gly3Pg, K14ADPg, K14ATPg, K14Gly, K14Gly3Pg, Vm14, Vm14f, Vm14r, Vt, tot_cell);
  return [vGlcTr, vHK, vPGI, vPFK, vALD, vTPI, vGAPdh, vGDH, vGPO, vPyrTr, vPGK, vPK, vATPase, vGlyK]';
}

vector get_odes(vector fluxes){
  real vGlcTr = fluxes[1];
  real vHK = fluxes[2];
  real vPGI = fluxes[3];
  real vPFK = fluxes[4];
  real vALD = fluxes[5];
  real vTPI = fluxes[6];
  real vGAPdh = fluxes[7];
  real vGDH = fluxes[8];
  real vGPO = fluxes[9];
  real vPyrTr = fluxes[10];
  real vPGK = fluxes[11];
  real vPK = fluxes[12];
  real vATPase = fluxes[13];
  real vGlyK = fluxes[14];
  return [
    (1.0*vGlcTr-1.0*vHK)/5.7,  // GlcI
    (1.0*vPGK+1.0*vGlyK-1.0*vHK-1.0*vPFK)/0.2446,  // Pg
    (1.0*vHK-1.0*vPGI)/0.2446,  // Glc6P
    (1.0*vPGI-1.0*vPFK)/0.2446,  // Fru6P
    (1.0*vPFK-1.0*vALD)/0.2446,  // Fru16BP
    (1.0*vALD+1.0*vGPO-1.0*vTPI-1.0*vGDH)/5.7,  // DHAP
    (1.0*vALD+1.0*vTPI-1.0*vGAPdh)/0.2446,  // GAP
    (1.0*vGDH-1.0*vGAPdh)/0.2446,  // NAD
    (1.0*vGAPdh-1.0*vPGK)/0.2446,  // BPGA13
    (1.0*vGAPdh-1.0*vGDH)/0.2446,  // NADH
    (1.0*vPK-1.0*vPyrTr)/5.4554,  // Pyr
    (1.0*vPGK-1.0*vPK)/5.7,  // Nb
    (1.0*vPK-1.0*vATPase)/5.4554  // Pc
  ]';
}

vector steady_state_equation(vector ode_metabolites,
                             vector kinetic_parameters,
                             real[] known_reals,
                             int[] known_ints){
    return get_odes(get_kinetics(ode_metabolites, kinetic_parameters, known_reals));
}