#include kinetics.stan
vector steady_state_equation(vector species, vector kinetic_parameters, real[] known_reals, int[] known_ints){
  // unpack species...
  real ATP = species[1];
  real GAP = species[2];
  real F6P = species[3];
  real PEP = species[4];
  real P = species[5];
  real DAP = species[6];
  real eiia = species[7];
  real GLCp = species[8];
  real PGA2 = species[9];
  real ei = species[10];
  real PGA3 = species[11];
  real G6P = species[12];
  real PGN = species[13];
  real hpr = species[14];
  real BPG = species[15];
  real eiicb = species[16];
  real ADP = species[17];
  real tal = species[18];
  real NAD = species[19];
  real AKG = species[20];
  real OAA = species[21];
  real R5P = species[22];
  real RU5P = species[23];
  real S7P = species[24];
  real SUCCOA = species[25];
  real X5P = species[26];

  // unpack parameters...
  // PGI
  real Keq_PGI = kinetic_parameters[1];
  real KmF6P_PGI = kinetic_parameters[2];
  real KmG6P_PGI = kinetic_parameters[3];
  real KmPEP_PGI = kinetic_parameters[4];
  real Vmax_PGI = kinetic_parameters[5];
  real KmPGN_PGI = kinetic_parameters[6];

  // PFK
  real KefrADP_PFK = kinetic_parameters[7];
  real KefrPEP_PFK = kinetic_parameters[8];
  real KeftADP_PFK = kinetic_parameters[9];
  real KeftPEP_PFK = kinetic_parameters[10];
  real Keq_PFK = kinetic_parameters[11];
  real KirADP_PFK = kinetic_parameters[12];
  real KirATP_PFK = kinetic_parameters[13];
  real KirF6P_PFK = kinetic_parameters[14];
  real KirFDP_PFK = kinetic_parameters[15];
  real KitADP_PFK = kinetic_parameters[16];
  real KitATP_PFK = kinetic_parameters[17];
  real KitF6P_PFK = kinetic_parameters[18];
  real KitFDP_PFK = kinetic_parameters[19];
  real KmrADP_PFK = kinetic_parameters[20];
  real KmrATPMg_PFK = kinetic_parameters[21];
  real KmrF6P_PFK = kinetic_parameters[22];
  real KmrFDP_PFK = kinetic_parameters[23];
  real KmtADP_PFK = kinetic_parameters[24];
  real KmtATPMg_PFK = kinetic_parameters[25];
  real KmtF6P_PFK = kinetic_parameters[26];
  real KmtFDP_PFK = kinetic_parameters[27];
  real L0_PFK = kinetic_parameters[28];
  real Vmax_PFK = kinetic_parameters[29];
  real Wr_PFK = kinetic_parameters[30];
  real Wt_PFK = kinetic_parameters[31];
  real n_PFK = kinetic_parameters[32];

  // FBA
  real Keq_FBA = kinetic_parameters[33];
  real KmDAP_FBA = kinetic_parameters[34];
  real KmFDP_FBA = kinetic_parameters[35];
  real KmGAP_FBA = kinetic_parameters[36];
  real KmPEP_FBA = kinetic_parameters[37];
  real Vmax_FBA = kinetic_parameters[38];

  // TPI
  real Keq_TPI = kinetic_parameters[39];
  real KmDAP_TPI = kinetic_parameters[40];
  real KmGAP_TPI = kinetic_parameters[41];
  real Vmax_TPI = kinetic_parameters[42];

  // GDH
  real Keq_GDH = kinetic_parameters[43];
  real KmBPG_GDH = kinetic_parameters[44];
  real KmGAP_GDH = kinetic_parameters[45];
  real KmNAD_GDH = kinetic_parameters[46];
  real KmNADH_GDH = kinetic_parameters[47];
  real KmP_GDH = kinetic_parameters[48];
  real Vmax_GDH = kinetic_parameters[49];

  // PGK
  real Keq_PGK = kinetic_parameters[50];
  real KmADPMg_PGK = kinetic_parameters[51];
  real KmATPMg_PGK = kinetic_parameters[52];
  real KmBPG_PGK = kinetic_parameters[53];
  real KmPGA3_PGK = kinetic_parameters[54];
  real Vmax_PGK = kinetic_parameters[55];

  // GPM
  real Keq_GPM = kinetic_parameters[56];
  real KmPGA2_GPM = kinetic_parameters[57];
  real KmPGA1_GPM = kinetic_parameters[58];
  real Vmax_GPM = kinetic_parameters[59];

  // ENO
  real Keq_ENO = kinetic_parameters[60];
  real KmPEP_ENO = kinetic_parameters[61];
  real KmPGA2_ENO = kinetic_parameters[62];
  real Vmax_ENO = kinetic_parameters[63];

  // PYK
  real KefrFDP_PYK = kinetic_parameters[64];
  real KefrG6P_PYK = kinetic_parameters[65];
  real KefrGL6P_PYK = kinetic_parameters[66];
  real KefrR5P_PYK = kinetic_parameters[67];
  real KefrRU5P_PYK = kinetic_parameters[68];
  real KefrS7P_PYK = kinetic_parameters[69];
  real KefrX5P_PYK = kinetic_parameters[70];
  real KeftATP_PYK = kinetic_parameters[71];
  real KeftSUCCOA_PYK = kinetic_parameters[72];
  real KirADP_PYK = kinetic_parameters[73];
  real KirATP_PYK = kinetic_parameters[74];
  real KirPEP_PYK = kinetic_parameters[75];
  real KirPYR_PYK = kinetic_parameters[76];
  real KirPyrATP_PYK = kinetic_parameters[77];
  real KitADP_PYK = kinetic_parameters[78];
  real KitATP_PYK = kinetic_parameters[79];
  real KitPEP_PYK = kinetic_parameters[80];
  real KitPYR_PYK = kinetic_parameters[81];
  real KitPyrATP_PYK = kinetic_parameters[82];
  real KmrADPMg_PYK = kinetic_parameters[83];
  real KmrPEP_PYK = kinetic_parameters[84];
  real KmtADPMg_PYK = kinetic_parameters[85];
  real KmtPEP_PYK = kinetic_parameters[86];
  real L0_PYK = kinetic_parameters[87];
  real Vmax_PYK = kinetic_parameters[88];
  real n_PYK = kinetic_parameters[89];

  // F6P_GAP_TAL
  real Keq_F6P_GAP_TAL = kinetic_parameters[90];
  real kcat_F6P_GAP_TAL = kinetic_parameters[91];
  
  // FBP
  real KirAMP_FBP = kinetic_parameters[92];
  real KirAMPFDP_FBP = kinetic_parameters[93];
  real KirF6P_FBP = kinetic_parameters[94];
  real KirF6PMg_FBP = kinetic_parameters[95];
  real KirFDP_FBP = kinetic_parameters[96];
  real KirFDPMg_FBP = kinetic_parameters[97];
  real KirFDPMgMg_FBP = kinetic_parameters[98];
  real KirP_FBP = kinetic_parameters[99];
  real KirPF6P_FBP = kinetic_parameters[100];
  real KirPF6PMg_FBP = kinetic_parameters[101];
  real KirPMg_FBP = kinetic_parameters[102];
  real KitAMP_FBP = kinetic_parameters[103];
  real KitAMPFDP_FBP = kinetic_parameters[104];
  real KitF6P_FBP = kinetic_parameters[105];
  real KitF6PMg_FBP = kinetic_parameters[106];
  real KitFDP_FBP = kinetic_parameters[107];
  real KitFDPMg_FBP = kinetic_parameters[108];
  real KitFDPMgMg_FBP = kinetic_parameters[109];
  real KitP_FBP = kinetic_parameters[110];
  real KitPF6P_FBP = kinetic_parameters[111];
  real KitPF6PMg_FBP = kinetic_parameters[112];
  real KitPMg_FBP = kinetic_parameters[113];
  real KmrFDP_FBP = kinetic_parameters[114];
  real KmrMg_FBP = kinetic_parameters[115];
  real KmtFDP_FBP = kinetic_parameters[116];
  real KmtMg_FBP = kinetic_parameters[117];
  real L0_FBP = kinetic_parameters[118];
  real Vmax_FBP = kinetic_parameters[119];
  real n_FBP = kinetic_parameters[120];

  // PPS
  real KdAMP_PPS = kinetic_parameters[121];
  real KdATPMgPPS_PPS = kinetic_parameters[122];
  real KdMg_PPS = kinetic_parameters[123];
  real KdP_PPS = kinetic_parameters[124];
  real KdPEP_PPS = kinetic_parameters[125];
  real KdPYR_PPS = kinetic_parameters[126];
  real KefADP_PPS = kinetic_parameters[127];
  real KefAKG_PPS = kinetic_parameters[128];
  real KefATP_PPS = kinetic_parameters[129];
  real KefOAA_PPS = kinetic_parameters[130];
  real Keq_PPS = kinetic_parameters[131];
  real KmAMP_PPS = kinetic_parameters[132];
  real KmATPMg_PPS = kinetic_parameters[133];
  real KmP_PPS = kinetic_parameters[134];
  real KmPEP_PPS = kinetic_parameters[135];
  real KmPYR_PPS = kinetic_parameters[136];
  real Vmax_PPS = kinetic_parameters[137];
  real W_PPS = kinetic_parameters[138];
  real alpha_PPS = kinetic_parameters[139];

  // PTS_0
  real KmPEP_PTS_0 = kinetic_parameters[140];
  real KmPYR_PTS_0 = kinetic_parameters[141];
  real kF_PTS_0 = kinetic_parameters[142];
  real kR_PTS_0 = kinetic_parameters[143];

  // PTS_1
  real k1_PTS_1 = kinetic_parameters[144];
  real k2_PTS_1 = kinetic_parameters[145];

  // PTS_2
  real k1_PTS_2 = kinetic_parameters[146];
  real k2_PTS_2 = kinetic_parameters[147];

  // PTS_3
  real k1_PTS_3 = kinetic_parameters[148];
  real k2_PTS_3 = kinetic_parameters[149];

  // PTS_4
  real KmG6P_PTS_4 = kinetic_parameters[150];
  real KmGLC_PTS_4 = kinetic_parameters[151];
  real kF_PTS_4 = kinetic_parameters[152];
  real kR_PTS_4 = kinetic_parameters[153];

  // ATP_MAINTENANCE
  real Vmax_ATP_MAINTENANCE = kinetic_parameters[154];
  real Keq_ATP_MAINTENANCE = kinetic_parameters[155];

  // XCH_GLC
  real Vmax_XCH_GLC = kinetic_parameters[156];
  real Km_XCH_GLC = kinetic_parameters[157];

  // GL6P_HYDROLYSIS
  real KGl6Phydrol_GL6P_HYDROLYSIS = kinetic_parameters[158];
  real KeqGl6Phydrol_GL6P_HYDROLYSIS = kinetic_parameters[159];
  
  // flux equations
  real PGI = Vmax*(G6P_c-F6P_c/Keq)/KmG6P/(1+F6P_c/KmF6P+G6P_c/KmG6P+PEP_c/KmPEP+PGN_c/KmPGN);
  real PFK = Vmax_1*n*(MgATP_c*F6P_c-MgADP_c*FDP_c/Keq_1)/(KirF6P*KmrATPMg)/(1+KmrFDP/KirFDP*(MgADP_c/KmrADP)+KmrF6P/KirF6P*(MgATP_c/KmrATPMg)+KmrFDP/KirFDP*(MgADP_c/KmrADP)*(F6P_c/KirF6P)+MgATP_c/KmrATPMg*(F6P_c/KirF6P)+MgADP_c/KirADP*(MgATP_c/KmrATPMg)*(F6P_c/KirF6P)+(1+(ATP_c-MgATP_c)/KirATP)*(F6P_c/KirF6P)+FDP_c/KirFDP+MgADP_c/KmrADP*(FDP_c/KirFDP)+KmrF6P/KirF6P*(MgATP_c/KmrATPMg)*(FDP_c/KirFDP)+Wr*(KmrF6P/KirF6P)*(MgADP_c/KirADP)*(MgATP_c/KmrATPMg)*(FDP_c/KmrFDP))/(1+L0*((1+KmtFDP/KitFDP*(MgADP_c/KmtADP)+KmtF6P/KitF6P*(MgATP_c/KmtATPMg)+KmtFDP/KitFDP*(MgADP_c/KmtADP)*(F6P_c/KitF6P)+MgATP_c/KmtATPMg*(F6P_c/KitF6P)+MgADP_c/KitADP*(MgATP_c/KmtATPMg)*(F6P_c/KitF6P)+(1+(ATP_c-MgATP_c)/KitATP)*(F6P_c/KitF6P)+FDP_c/KitFDP+MgADP_c/KmtADP*(FDP_c/KitFDP)+KmtF6P/KitF6P*(MgATP_c/KmtATPMg)*(FDP_c/KitFDP)+Wt*(KmtF6P/KitF6P)*(MgADP_c/KitADP)*(MgATP_c/KmtATPMg)*(FDP_c/KmtFDP))*(1+MgADP_c/KeftADP+PEP_c/KeftPEP+MgADP_c/KeftADP*(PEP_c/KeftPEP))/((1+KmrFDP/KirFDP*(MgADP_c/KmrADP)+KmrF6P*MgATP_c/(KirF6P*KmrATPMg)+KmrFDP/KirFDP*(MgADP_c/KmrADP)*(F6P_c/KirF6P)+MgATP_c/KmrATPMg*(F6P_c/KirF6P)+MgADP_c/KirADP*(MgATP_c/KmrATPMg)*(F6P_c/KirF6P)+(1+(ATP_c-MgATP_c)/KirATP)*(F6P_c/KirF6P)+FDP_c/KirFDP+MgADP_c/KmrADP*(FDP_c/KirFDP)+KmrF6P/KirF6P*(MgATP_c/KmrATPMg)*(FDP_c/KirFDP)+Wr*(KmrF6P/KirF6P)*(MgADP_c/KirADP)*(MgATP_c/KmrATPMg)*(FDP_c/KmrFDP))*(1+MgADP_c/KefrADP+PEP_c/KefrPEP+MgADP_c/KefrADP*(PEP_c/KefrPEP))))^n);
  real FBA = Vmax_2*(FDP_c-DAP_c*GAP_c/Keq_2)/KmFDP/(1+FDP_c/KmFDP+DAP_c/KmDAP+DAP_c/KmDAP*(GAP_c/KmGAP)+PEP_c/KmPEP_1);
  real TPI = Vmax_3*(DAP_c-GAP_c/Keq_3)/KmDAP_1/(1+DAP_c/KmDAP_1+GAP_c/KmGAP_1);
  real GDH = Vmax_4*(P_c*GAP_c*NAD_c-BPG_c*NADH_c/Keq_4)/(KmP*KmGAP_2*KmNAD)/((1+P_c/KmP)*(1+GAP_c/KmGAP_2)*(1+NAD_c/KmNAD)+(1+BPG_c/KmBPG)*(1+NADH_c/KmNADH)-1);
  real PGK =
  real GPM =
  real ENO =
  real PYK =
  real F6P_GAP_TAL = 
  real FBP = 
  real PPS =
  real PTS_0 =
  real PTS_1 =
  real PTS_2 =
  real PTS_3 =
  real PTS_4 =
  real ATP_MAINTENANCE = 
  real XCH_GLC = 
  real GL6P_HYDROLYSIS = 
  
  // work out rate of change of each species concentration
  vector[N_reactant] dsdt = [PGK + PYK - PPS - ATP_MAINTENANCE,  // ATP
                             FBA + TPI - GDH - F6P,              // GAP
                             PGI - PFK + F6P_GAP_TAL,            // F6P
                             ENO - PYK + PPS - PTS_0,            // PEP
                             -GDH + FBP + PPS + ATP_MAINTENANCE, // P
                             FBA - TPI, // DAP
                             -PTS_2 + PTS_3, // eiia
                             -PTS_4 + XCH_GLC, // GLCp
                             GPM - ENO, // PGA2
                             -PTS_0 + PTS_1, // ei
                             PGK - ENO, // PGA3
                             -PGI + PTS_4, // G6P
                             GL6P_HYDROLYSIS, // PGN
                             -PTS_1 + PTS_2, // hppr
                             GDH - PGK, // BPG
                             -PTS_3 + PTS_4, // eiicb
                             PFK - PGK - PYK + ADP_MAINTENANCE, // ADP
                             F6P_GL6P, // tal
                             , // NAD
                             , // PYR
                             , // eiicbP
                             , // hprP
                             , // eiiaP
                             , // talC3
                             , // GK6P
                             , // AMP
                             , // FDP
                             , // NADH
                             , // GLCx
                             ,] // eiP
d/dt(ATP) = -PFK*cell_cytoplasm+PGK*cell_cytoplasm+PYK*cell_cytoplasm-PPS*cell_cytoplasm-ATP_MAINTENANCE_1*cell_cytoplasm		;  

d/dt(GAP) = FunctionForFBA*cell_cytoplasm+FunctionForTPI*cell_cytoplasm-FunctionForGDH*cell_cytoplasm-FunctionForF6P_GAP_TAL*cell_cytoplasm		;  
d/dt(F6P) = FunctionForPGI_1*cell_cytoplasm-FunctionForPFK*cell_cytoplasm+FunctionForF6P_GAP_TAL*cell_cytoplasm+FunctionForFBP*cell_cytoplasm		;  
d/dt(PEP) = FunctionForENO*cell_cytoplasm-FunctionForPYK*cell_cytoplasm+FunctionForPPS*cell_cytoplasm-FunctionForPTS_0*cell_cytoplasm		;  
d/dt(P) = -FunctionForGDH*cell_cytoplasm+FunctionForFBP*cell_cytoplasm+FunctionForPPS*cell_cytoplasm+FunctionForATP_MAINTENANCE_1*cell_cytoplasm		;  
d/dt(DAP) = FunctionForFBA*cell_cytoplasm-FunctionForTPI*cell_cytoplasm		;  
d/dt(eiia) = -MassAction_reversible__1*cell_cytoplasm+MassAction_reversible__2*cell_cytoplasm		;  
d/dt(GLCp) = -FunctionForPTS_4_1+FunctionForXCH_RMM_1		;  
d/dt(PGA2) = FunctionForGPM*cell_cytoplasm-FunctionForENO*cell_cytoplasm		;  
d/dt(ei) = -FunctionForPTS_0*cell_cytoplasm+MassAction_reversible_*cell_cytoplasm		;  
d/dt(PGA3) = FunctionForPGK*cell_cytoplasm-FunctionForGPM*cell_cytoplasm		;  
d/dt(G6P) = -FunctionForPGI_1*cell_cytoplasm+FunctionForPTS_4_1		;  
d/dt(PGN) = FunctionForGL6P_HYDROLYSIS_1*cell_cytoplasm		;  
d/dt(hpr) = -MassAction_reversible_*cell_cytoplasm+MassAction_reversible__1*cell_cytoplasm		;  
d/dt(BPG) = FunctionForGDH*cell_cytoplasm-FunctionForPGK*cell_cytoplasm		;  
d/dt(eiicb) = -MassAction_reversible__2*cell_cytoplasm+FunctionForPTS_4_1		;  
d/dt(ADP) = FunctionForPFK*cell_cytoplasm-FunctionForPGK*cell_cytoplasm-FunctionForPYK*cell_cytoplasm+FunctionForATP_MAINTENANCE_1*cell_cytoplasm		;  
d/dt(tal) = FunctionForF6P_GAP_TAL*cell_cytoplasm		;  
d/dt(NAD) = -FunctionForGDH*cell_cytoplasm		;  
d/dt(AKG) = 0		;  
d/dt(OAA) = 0		;  
d/dt(R5P) = 0		;  
d/dt(RU5P) = 0		;  
d/dt(S7P) = 0		;  
d/dt(SUCCOA) = 0		;  
d/dt(X5P) = 0		;  


  return dsdt;

}
