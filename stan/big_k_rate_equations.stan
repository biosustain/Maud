// from the yeast file...

real two_noncompeting_couples(real S1, real S2, real P1, real P2, real V, real KS1, real KS2, real KP1, real KP2, real Keq){
  return V*(S1*S2/(KS1*KS2)-P1*P2/(KS1*KS2*Keq))/((1+S1/KS1+P1/KP1)*(1+S2/KS2+P2/KP2));
}

real reversible_michaelis_menten(real S, real P, real V, real KS, real KP, real Keq){
  return V*(S/KS-P/(KS*Keq))/(1+S/KS+P/KP);
}

real ordered_uni_bi(real S, real P1, real P2, real V, real KS, real KP1, real KP2, real KiP2, real Keq){
  return V*(S/KS-P1*P2/(KS*Keq))/(1+S/KS+P1/KP1+P2/KP2+S*P2/(KS*KiP2)+P1*P2/(KP1*KP2));
}

real phosphoglycerate_kinase_kinetics(real S1, real S2, real P1, real P2, real V, real KS1, real KS2, real KP1, real KP2, real Keq, real nS2){
  return V*(S2/KS2)^(nS2-1)*(S1*S2/(KS1*KS2)-P1*P2/(KS1*KS2*Keq))/((1+S1/KS1+P1/KP1)*(1+(S2/KS2)^nS2+P2/KP2));
}

real triphosphate_isomerase_kinetics(real S, real P, real V, real KS, real KP, real KiP, real Keq){
  return V/KS*(S-P/Keq)/(1+S/KS+P/KP*(1+(P/KiP)^4));
}

// From the t brucei file...
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
