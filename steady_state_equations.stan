#include kinetics.stan
vector steady_state_equation(vector species, vector kinetic_parameters, real[] known_reals, int[] known_ints){
  // fixed numbers
  real MG = 1;
  real KdADPMg = 1.27771;
  real KdATPMg = 0.0847634;
  real KdFDPMg = 5.81;
  real KmICIT_ACN = 9.31352;
  real KmCIT_ACN = 0.0628882;
  real KmACO_ACN = 0.02001;
  real KeqNDH = 27.6193;

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
  real Keq = kinetic_parameters[1];
  real KmF6P = kinetic_parameters[2];
  real KmG6P = kinetic_parameters[3];
  real KmPEP = kinetic_parameters[4];
  real Vmax = kinetic_parameters[5];
  real KmPGN = kinetic_parameters[6];
  real KefrADP = kinetic_parameters[7];
  real KefrPEP = kinetic_parameters[8];
  real KeftADP = kinetic_parameters[9];
  real KeftPEP = kinetic_parameters[10];
  real Keq_1 = kinetic_parameters[11];
  real KirADP = kinetic_parameters[12];
  real KirATP = kinetic_parameters[13];
  real KirF6P = kinetic_parameters[14];
  real KirFDP = kinetic_parameters[15];
  real KitADP = kinetic_parameters[16];
  real KitATP = kinetic_parameters[17];
  real KitF6P = kinetic_parameters[18];
  real KitFDP = kinetic_parameters[19];
  real KmrADP = kinetic_parameters[20];
  real KmrATPMg = kinetic_parameters[21];
  real KmrF6P = kinetic_parameters[22];
  real KmrFDP = kinetic_parameters[23];
  real KmtADP = kinetic_parameters[24];
  real KmtATPMg = kinetic_parameters[25];
  real KmtF6P = kinetic_parameters[26];
  real KmtFDP = kinetic_parameters[27];
  real L0 = kinetic_parameters[28];
  real Vmax_1 = kinetic_parameters[29];
  real Wr = kinetic_parameters[30];
  real Wt = kinetic_parameters[31];
  real n = kinetic_parameters[32];
  real Keq_2 = kinetic_parameters[33];
  real KmDAP = kinetic_parameters[34];
  real KmFDP = kinetic_parameters[35];
  real KmGAP = kinetic_parameters[36];
  real KmPEP_1 = kinetic_parameters[37];
  real Vmax_2 = kinetic_parameters[38];
  real Keq_3 = kinetic_parameters[39];
  real KmDAP_1 = kinetic_parameters[40];
  real KmGAP_1 = kinetic_parameters[41];
  real Vmax_3 = kinetic_parameters[42];
  real Keq_4 = kinetic_parameters[43];
  real KmBPG = kinetic_parameters[44];
  real KmGAP_2 = kinetic_parameters[45];
  real KmNAD = kinetic_parameters[46];
  real KmNADH = kinetic_parameters[47];
  real KmP = kinetic_parameters[48];
  real Vmax_4 = kinetic_parameters[49];
  real Keq_5 = kinetic_parameters[50];
  real KmADPMg = kinetic_parameters[51];
  real KmATPMg = kinetic_parameters[52];
  real KmBPG_1 = kinetic_parameters[53];
  real KmPGA3 = kinetic_parameters[54];
  real Vmax_5 = kinetic_parameters[55];
  real Keq_6 = kinetic_parameters[56];
  real KmPGA2 = kinetic_parameters[57];
  real KmPGA3_1 = kinetic_parameters[58];
  real Vmax_6 = kinetic_parameters[59];
  real Keq_7 = kinetic_parameters[60];
  real KmPEP_2 = kinetic_parameters[61];
  real KmPGA2_1 = kinetic_parameters[62];
  real Vmax_7 = kinetic_parameters[63];
  real KefrFDP = kinetic_parameters[64];
  real KefrG6P = kinetic_parameters[65];
  real KefrGL6P = kinetic_parameters[66];
  real KefrR5P = kinetic_parameters[67];
  real KefrRU5P = kinetic_parameters[68];
  real KefrS7P = kinetic_parameters[69];
  real KefrX5P = kinetic_parameters[70];
  real KeftATP = kinetic_parameters[71];
  real KeftSUCCOA = kinetic_parameters[72];
  real KirADP_1 = kinetic_parameters[73];
  real KirATP_1 = kinetic_parameters[74];
  real KirPEP = kinetic_parameters[75];
  real KirPYR = kinetic_parameters[76];
  real KirPyrATP = kinetic_parameters[77];
  real KitADP_1 = kinetic_parameters[78];
  real KitATP_1 = kinetic_parameters[79];
  real KitPEP = kinetic_parameters[80];
  real KitPYR = kinetic_parameters[81];
  real KitPyrATP = kinetic_parameters[82];
  real KmrADPMg = kinetic_parameters[83];
  real KmrPEP = kinetic_parameters[84];
  real KmtADPMg = kinetic_parameters[85];
  real KmtPEP = kinetic_parameters[86];
  real L0_1 = kinetic_parameters[87];
  real Vmax_8 = kinetic_parameters[88];
  real n_1 = kinetic_parameters[89];
  real KirAMP = kinetic_parameters[90];
  real KirAMPFDP = kinetic_parameters[91];
  real KirF6P_1 = kinetic_parameters[92];
  real KirF6PMg = kinetic_parameters[93];
  real KirFDP_1 = kinetic_parameters[94];
  real KirFDPMg = kinetic_parameters[95];
  real KirFDPMgMg = kinetic_parameters[96];
  real KirP = kinetic_parameters[97];
  real KirPF6P = kinetic_parameters[98];
  real KirPF6PMg = kinetic_parameters[99];
  real KirPMg = kinetic_parameters[100];
  real KitAMP = kinetic_parameters[101];
  real KitAMPFDP = kinetic_parameters[102];
  real KitF6P_1 = kinetic_parameters[103];
  real KitF6PMg = kinetic_parameters[104];
  real KitFDP_1 = kinetic_parameters[105];
  real KitFDPMg = kinetic_parameters[106];
  real KitFDPMgMg = kinetic_parameters[107];
  real KitP = kinetic_parameters[108];
  real KitPF6P = kinetic_parameters[109];
  real KitPF6PMg = kinetic_parameters[110];
  real KitPMg = kinetic_parameters[111];
  real KmrFDP_1 = kinetic_parameters[112];
  real KmrMg = kinetic_parameters[113];
  real KmtFDP_1 = kinetic_parameters[114];
  real KmtMg = kinetic_parameters[115];
  real L0_2 = kinetic_parameters[116];
  real Vmax_9 = kinetic_parameters[117];
  real n_2 = kinetic_parameters[118];
  real KdAMP = kinetic_parameters[119];
  real KdATPMgPPS = kinetic_parameters[120];
  real KdMg = kinetic_parameters[121];
  real KdP = kinetic_parameters[122];
  real KdPEP = kinetic_parameters[123];
  real KdPYR = kinetic_parameters[124];
  real KefADP = kinetic_parameters[125];
  real KefAKG = kinetic_parameters[126];
  real KefATP = kinetic_parameters[127];
  real KefOAA = kinetic_parameters[128];
  real Keq_8 = kinetic_parameters[129];
  real KmAMP = kinetic_parameters[130];
  real KmATPMg_1 = kinetic_parameters[131];
  real KmP_1 = kinetic_parameters[132];
  real KmPEP_3 = kinetic_parameters[133];
  real KmPYR = kinetic_parameters[134];
  real Vmax_10 = kinetic_parameters[135];
  real W = kinetic_parameters[136];
  real alpha = kinetic_parameters[137];
  real KmPEP_4 = kinetic_parameters[139];
  real KmPYR_1 = kinetic_parameters[140];
  real kF = kinetic_parameters[141];
  real kR = kinetic_parameters[142];
  real k1 = kinetic_parameters[143];
  real k2 = kinetic_parameters[144];
  real k1_1 = kinetic_parameters[145];
  real k2_1 = kinetic_parameters[146];
  real k1_2 = kinetic_parameters[147];
  real k2_2 = kinetic_parameters[148];
  real KmG6P_1 = kinetic_parameters[149];
  real KmGLC = kinetic_parameters[150];
  real kF_1 = kinetic_parameters[151];
  real kR_1 = kinetic_parameters[152];
  real Vmax_11 = kinetic_parameters[153];
  real Keq_9 = kinetic_parameters[154];
  real Vmax_12 = kinetic_parameters[155];
  real Km = kinetic_parameters[156];
 
  // flux equations
  real flux_PGI_1 = Vmax*(G6P-F6P/Keq)/KmG6P/(1+F6P/KmF6P+G6P/KmG6P+PEP/KmPEP+PGN/KmPGN);
  real flux_PFK = Vmax_1*n*(MgATP*F6P-MgADP*FDP/Keq_1)/(KirF6P*KmrATPMg)/(1+KmrFDP/KirFDP*(MgADP/KmrADP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)+KmrFDP/KirFDP*(MgADP/KmrADP)*(F6P/KirF6P)+MgATP/KmrATPMg*(F6P/KirF6P)+MgADP/KirADP*(MgATP/KmrATPMg)*(F6P/KirF6P)+(1+(ATP-MgATP)/KirATP)*(F6P/KirF6P)+FDP/KirFDP+MgADP/KmrADP*(FDP/KirFDP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)*(FDP/KirFDP)+Wr*(KmrF6P/KirF6P)*(MgADP/KirADP)*(MgATP/KmrATPMg)*(FDP/KmrFDP))/(1+L0*((1+KmtFDP/KitFDP*(MgADP/KmtADP)+KmtF6P/KitF6P*(MgATP/KmtATPMg)+KmtFDP/KitFDP*(MgADP/KmtADP)*(F6P/KitF6P)+MgATP/KmtATPMg*(F6P/KitF6P)+MgADP/KitADP*(MgATP/KmtATPMg)*(F6P/KitF6P)+(1+(ATP-MgATP)/KitATP)*(F6P/KitF6P)+FDP/KitFDP+MgADP/KmtADP*(FDP/KitFDP)+KmtF6P/KitF6P*(MgATP/KmtATPMg)*(FDP/KitFDP)+Wt*(KmtF6P/KitF6P)*(MgADP/KitADP)*(MgATP/KmtATPMg)*(FDP/KmtFDP))*(1+MgADP/KeftADP+PEP/KeftPEP+MgADP/KeftADP*(PEP/KeftPEP))/((1+KmrFDP/KirFDP*(MgADP/KmrADP)+KmrF6P*MgATP/(KirF6P*KmrATPMg)+KmrFDP/KirFDP*(MgADP/KmrADP)*(F6P/KirF6P)+MgATP/KmrATPMg*(F6P/KirF6P)+MgADP/KirADP*(MgATP/KmrATPMg)*(F6P/KirF6P)+(1+(ATP-MgATP)/KirATP)*(F6P/KirF6P)+FDP/KirFDP+MgADP/KmrADP*(FDP/KirFDP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)*(FDP/KirFDP)+Wr*(KmrF6P/KirF6P)*(MgADP/KirADP)*(MgATP/KmrATPMg)*(FDP/KmrFDP))*(1+MgADP/KefrADP+PEP/KefrPEP+MgADP/KefrADP*(PEP/KefrPEP))))^n);
  real flux_FBA = Vmax_2*(FDP-DAP*GAP/Keq_2)/KmFDP/(1+FDP/KmFDP+DAP/KmDAP+DAP/KmDAP*(GAP/KmGAP)+PEP/KmPEP_1);
  real flux_TPI = Vmax_3*(DAP-GAP/Keq_3)/KmDAP_1/(1+DAP/KmDAP_1+GAP/KmGAP_1);
  real flux_GDH = Vmax_4*(P*GAP*NAD-BPG*NADH/Keq_4)/(KmP*KmGAP_2*KmNAD)/((1+P/KmP)*(1+GAP/KmGAP_2)*(1+NAD/KmNAD)+(1+BPG/KmBPG)*(1+NADH/KmNADH)-1);
  real flux_PGK = Vmax_5*(MgADP*BPG-MgATP*PGA3/Keq_5)/(KmADPMg*KmBPG_1)/(1+MgADP/KmADPMg+BPG/KmBPG_1+MgADP/KmADPMg*BPG/KmBPG_1+MgATP/KmATPMg+PGA3/KmPGA3+MgATP/KmATPMg*PGA3/KmPGA3);
  real flux_GPM = Vmax_6*(PGA3-PGA2/Keq_6)/KmPGA3_1/(1+PGA3/KmPGA3_1+PGA2/KmPGA2);
  real flux_ENO = Vmax_7*(PGA2-PEP/Keq_7)/KmPGA2_1/(1+PGA2/KmPGA2_1+PEP/KmPEP_2);
  real flux_PYK = Vmax_8*n_1*PEP*MgADP/(KirPEP*KmrADPMg)/(1+KmrPEP/KirPEP*(MgADP/KmrADPMg)+MgATP/KirATP_1+MgADP/KmrADPMg*(PEP/KirPEP)+KmrADPMg/KmrADPMg*(1+(ADP-MgADP)/KirADP_1)*(PEP/KirPEP)+PYR/KirPYR+MgATP/KirPyrATP*(PYR/KirPYR))/(1+L0_1*((1+KmtPEP/KitPEP*(MgADP/KmtADPMg)+MgATP/KitATP_1+MgADP*PEP/(KitPEP*KmtADPMg)+(1+(ADP-MgADP)/KitADP_1)*(PEP/KitPEP)+PYR/KitPYR+MgATP/KitPyrATP*(PYR/KitPYR))*(1+SUCCOA/KeftSUCCOA+MgATP*SUCCOA/(KeftATP*KeftSUCCOA))/((1+KmrPEP/KirPEP*(MgADP/KmrADPMg)+MgATP/KirATP_1+MgADP/KmrADPMg*(PEP/KirPEP)+(1+(ADP-MgADP)/KirADP_1)*(PEP/KirPEP)+PYR/KirPYR+MgATP/KirPyrATP*(PYR/KirPYR))*(1+FDP/KefrFDP+G6P/KefrG6P+GL6P/KefrGL6P+R5P/KefrR5P+RU5P/KefrRU5P+S7P/KefrS7P+X5P/KefrX5P)))^n_1);
  real flux_FBP = Vmax_9*n_2*MgFDP/KirFDPMg/(1+KmrFDP_1/KirFDP_1*(MG/KmrMg)+P/KirP+P/KirP*(MG/KirPMg)+F6P/KirF6P_1+F6P/KirF6P_1*(MG/KirF6PMg)+P/KirP*(F6P/KirPF6P)+P/KirP*(F6P/KirPF6P)*(MG/KirPF6PMg)+(FDP-MgFDP)/KirFDP_1+KdFDPMg/KmrMg*(MgFDP/KirFDP_1)+AMP/KirAMP+MgFDP/KirFDPMg+MgFDP/KirFDPMg*(MG/KirFDPMgMg)+AMP/KirAMP*((FDP-MgFDP)/KirAMPFDP))/(1+L0_2*((1+KmtFDP_1/KitFDP_1*(MG/KmtMg)+P/KitP+P/KitP*(MG/KitPMg)+F6P/KitF6P_1+F6P/KitF6P_1*(MG/KitF6PMg)+P/KitP*(F6P/KitPF6P)+P/KitP*(F6P/KitPF6P)*(MG/KitPF6PMg)+(FDP-MgFDP)/KitFDP_1+KdFDPMg/KmtMg*(MgFDP/KitFDP_1)+AMP/KitAMP+MgFDP/KitFDPMg+MgFDP/KitFDPMg*(MG/KitFDPMgMg)+AMP/KitAMP*((FDP-MgFDP)/KitAMPFDP))/(1+KmrFDP_1/KirFDP_1*(MG/KmrMg)+P/KirP+P/KirP*(MG/KirPMg)+F6P/KirF6P_1+F6P/KirF6P_1*(MG/KirF6PMg)+P/KirP*(F6P/KirPF6P)+P/KirP*(F6P/KirPF6P)*(MG/KirPF6PMg)+(FDP-MgFDP)/KirFDP_1+KdFDPMg/KmrMg*(MgFDP/KirFDP_1)+AMP/KirAMP+MgFDP/KirFDPMg+MgFDP/KirFDPMg*(MG/KirFDPMgMg)+AMP/KirAMP*((FDP-MgFDP)/KirAMPFDP)))^n_2)
    real flux_PPS = Vmax_10*(MgATP*PYR-AMP*PEP*P*MG/Keq_8)/(KmATPMg_1*KmPYR)/(MgATP/KmATPMg_1+alpha*(P/KdP)*(MgATP/KmATPMg_1)+alpha*(AMP/KdAMP)*(MgATP/KmATPMg_1)+alpha*(P/KdP)*(AMP/KdAMP)*(MgATP/KmATPMg_1)+alpha*(MG/KdMg)*(P/KmP_1)*(AMP/KdAMP)*(MgATP/KdATPMgPPS)/(W*(1+MG/KdMg))+MgATP/KmATPMg_1*(AKG/KefAKG)+(1+MG/KdMg)*(AKG/KefAKG)*(PEP/KmPEP_3)/W+MgATP/KmATPMg_1*(OAA/KefOAA)+(1+MG/KdMg)*(OAA/KefOAA)*(PEP/KmPEP_3)/W+MG/KdMg*(P/KmP_1)*(AMP/KdAMP)/W+alpha*(P/KdP)*(AMP/KdAMP)*(PEP/KmPEP_3)/W+alpha*(MG/KdMg)*(P/KmP_1)*(AMP/KdAMP)*(PEP/KmPEP_3)/W+alpha*(1+MG/KdMg)*(KmAMP/KdAMP*(P/KmP_1)*(PEP/KmPEP_3)+AMP/KdAMP*(PEP/KmPEP_3))/W+(1+MG/KdMg)*(PYR/KmPYR)+MgATP/KmATPMg_1*(PYR/KmPYR)+KdADPMg/KdMg*(P/KmP_1)*(MgADP/KefADP)*(AMP/KdAMP)/(W*(1+MG/KdMg))+(ADP-MgADP)/KefADP*(PYR/KmPYR)+KdATPMg/KdMg*(P/KmP_1)*(AMP/KdAMP)*(MgATP/KefATP)/(W*(1+MG/KdMg))+(ATP-MgATP)/KefATP*(PYR/KmPYR)+(1+MG/KdMg)*(PEP/KmPEP_3)/W+alpha*(1+MG/KdMg)*(PEP/KdPEP)*(PYR/KmPYR)+(1+MG/KdMg)*(PYR/KdPYR)*(PEP/KmPEP_3)/W);
  real flux_PTS_0 = kF*ei*PEP^2/(KmPEP_4^2+PEP^2)-kR*eiP*PYR^2/(KmPYR_1^2+PYR^2);
  real flux_PTS_1 = kF*ei*PEP^2/(KmPEP_4^2+PEP^2)-kR*eiP*PYR^2/(KmPYR_1^2+PYR^2);
  real flux_PTS_2 = k1_1*eiia*hprP-k2_1*eiiaP*hpr;
  real flux_PTS_3 = k1_2*eiicb*eiiaP-k2_2*eiicbP*eiia;
  real flux_PTS_4 = *(kF_1*eiicbP*GLCp/(KmGLC+GLCp)-kR_1*eiicb*G6P/(KmG6P_1+G6P));
  real flux_ATP_MAINTENANCE = Vmax_11*(ATP-ADP*P/Keq_9);
  real flux_XCH_RMM_1 = Vmax_12*(GLCx/Km-GLCp/Km)/(1+GLCx/Km+GLCp/Km);
  
  // work out rate of change of each species concentration
  vector[N_reactant] dsdt = [-flux_PFK*+flux_PGK*+flux_PYK*-flux_PPS*-flux_ATP_MAINTENANCE_1*, // ATP
                             flux_ENO*-flux_PYK*+flux_PPS*-flux_PTS_0*, // PEP
                             -flux_GDH*+flux_FBP*+flux_PPS*+flux_ATP_MAINTENANCE_1*, // P
                             flux_FBA*+flux_TPI*-flux_GDH*, // GAP
                             flux_PGI_1*-flux_PFK*+flux_FBP*, // F6P
                             flux_FBA*-flux_TPI*, // DAP
                             -PTS_2*+PTS_3*, // eiia
                             -flux_PTS_4_1+flux_XCH_RMM_1, // GLCp
                             flux_GPM*-flux_ENO*, // PGA2
                             -flux_PTS_0*+MassAction_reversible_*, // ei
                             flux_PGK*-flux_GPM*, // PGA3
                             -PTS_3*+flux_PTS_4_1, // eiicb
                             flux_PFK*-flux_FBA*-flux_FBP*, // FDP
                             -MassAction_reversible_*+PTS_2*, // hpr
                             flux_PFK*-flux_PGK*-flux_PYK*+flux_ATP_MAINTENANCE_1*, // ADP
                             -flux_PGI_1*+flux_PTS_4_1, // G6P
                             flux_GDH*] // NADH
  return dsdt;

}
