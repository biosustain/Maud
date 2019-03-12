real mass_action_irrev(real S, real k){
  return k * S;
}

real mass_action_rev(real S, real P, real k_fwd, real k_rev){
  return k_fwd * S - k_rev * P;
}

real uni_uni(real S, real P, real Vf, real Keq, real Kmp, real Kms){
  return Vf * (S - P / Keq) / (S + Kms * (1 + P/Kmp));
}

real uncompetitive_inhibition_irrev(real S, real I, real Km, real V, real Ki){
  return V * S / (Km + S * (1 + I/Ki));
}

real uncompetitive_inhibition_rev(real S, real P, real I, real Kms, real Kmp, real Vf, real Vr, real Ki){
  return (Vf * S/Kms - Vr * P/Kmp) / (1 + (S/Kms + P/Kmp) * (1 + I/Ki));
}

real substrate_inhibition_irrev(real S, real Km, real V, real Ki){
  return V * S / (Km + S + Km * pow(S/Ki, 2));
}

real substrate_inhibition_rev(real S, real P, real Kms, real Kmp, real Vf, real Vr, real Ki){
  return (Vf * S/Kms - Vr * P/Kmp) / (1 + S/Kms + P/Kmp + pow(S/Ki, 2));
}

real substrate_activation_irrev(real S, real V, real Ksc, real Ksa){
  return V * pow(S/Ksa, 2) / (1 + S/Ksc + S/Ksa + pow(S/Ksa, 2));
}

real specific_activation_irrev(real S, real A, real Kms, real V, real Ka){
  return V * S * A / (Kms * Ka + A * (Kms + S));
}

real specific_activation_rev(real S, real P, real A, real Kms, real Kmp, real Vf, real Vr, real Ka){
  return A * (Vf * S/Kms - Vr * P/Kmp) / (Ka + A * (1 + S/Kms + P/Kmp));
}

real michaelis_menten_rev(real S, real P, real Kms, real Kmp, real Vf, real Vr){
  return (Vf * S/Kms - Vr * P/Kmp) / (1 + S/Kms + P/Kmp);
}

real michaelis_menten_irrev(real S, real Km, real V){
  return V * S / (Km + S);
}

real ordered_bi_bi(real S1,
                   real S2,
                   real P1,
                   real P2,
                   real Keq,
                   real Vf,
                   real Vr,
                   real KmS1,
                   real KmS2,
                   real KmP1,
                   real KmP2,
                   real KiS1,
                   real KiS2,
                   real KiP
                   ){
  real numerator = Vf * (S1 * S2 - (P1 * P2 / Keq));
  real denominator = (
                      S1 * S2 * (1 + P1/KiP)
                      + KmS1 * S2
                      + KmS2 * (S1 + KiS1)
                      + Vf/(Vr*Keq) * (KmP2 * P1 * (1 + S1/KiS1))
                      + P2 * (Kmp * (1 + (KmS1 * S2)/(KiS1 * KmS2)) + P1 * (1 + S2/KiS2))
                      );
  return numerator / denominator;
}

real flux_atp_maintenance(real ATP, real ADP, real P, real Keq, real Vmax){
  return Vmax * (ATP - ADP * P/Keq);
}

real flux_eno(real S, real P, real Vmax, real Keq, real KmS, real KmP){
  return (Vmax * (S - P/Keq) / KmS) / (1 + S/KmS + P/KmP);
}

real flux_F6P_GAP_TAL(real GAP, real talC3, real FCP, real tal, real Keq, real Kcat){
  return Kcat * (GAP * talC3 - F6P * tal/Keq);
}

real flux_FBA(real DAP real FDP, real GAP, real Keq, real KmDAP, real KmFDP, real KmGAP, real KmPEP, real PEP, real Vmax){
  /* PEP is a 'modifier', not sure what this means...  */
  return (Vmax * (FDP - (DAP * GAP/Keq))/ KmFDP) / (1 + FDP/KmFDP + 2 * DAP/KmDAP + GAP/KmGAP + PEP/KmPEP);
}

real flux_FBP(real AMP, real F6P, real FDP, real KdFDPMg, real KirAMP, real KirAMPFDP,
              real KirF6P, real KirF6PMg, real KirFDP, real KirFDPMg, real KirFDPMgMg,
              real KirP, real KirPF6P, real KirPF6PMg, real KirPMg, real KitAMP,
              real KitAMPFDP, real KitF6P, real KitF6PMg, real KitFDP, real KitFDPMg,
              real KitFDPMgMg, real KitP, real KitPF6P, real KitPF6PMg, real KitPMg,
              real KmrMg, real KmtFDP, real KmtMg, real L0, real MG, real MgFDP, real P,
              real Vmax, real n){
  real num_num = Vmax * n * MgFDP/KirFDPMg;
  real num_denom = (1
                    + KmrFDP/KirFDP
                    + MG/KmrMg
                    + P/KirP
                    + (P/KirP) * (MG/KirPMG)
                    + F6P/KirF6P
                    + (F6P/KirF6P) * (MG/KirF6PMg)
                    + (P/KirP) * (F6P/KirPF6P) * (MG/KirPF6PMg)
                    + (FDP - MgFDP) / KirFDP
                    + (KdFDPMg/KmrMg) * (MgFDP/KirFDP)
                    + AMP/KirAMP
                    + MgFDP/KirFDPMg
                    + (MgFDP/KirFDPMg) * (MG/KirFDPMgMg)
                    + (AMP/KirAMP) * (FDP - MgFDP) / KirAMPFDP);
  real denom_num = (1
                    + (KmtFDP/KitFDP) * (MG/KmtMg)
                    + P/KitP
                    + (P/KitP) * (MG/KitPMg)
                    + F6P/KitF6P
                    + (F6P/KitF6P) * (MG/KitF6PMg)
                    + (P/KitP) * (F6P/KitPF6P)
                    + (P/KitP) * (F6P/KitPF6P) * (MG/KitPF6PMg)
                    + (FDP - MgFDP) / KitFDP
                    + (KdFDPMg/KmtMg) * (MgFDP/KitFDP)
                    + AMP/KitAMP
                    + MgFDP/KitFDPMg
                    + (MgFDP/KitFDPMg) * (MG/KitFDPMgMg)
                    + (AMP/KitAMP) * (FDP - MgFDP) / KitAMPFDP);
  real denom_denom = (1
                      + (KmrFDP/KirFDP) * (MG/KmrMg)
                      + P/KirP
                      + (P/KirP) * (MG/KirPMg)
                      + F6P/KirF6P
                      + (F6P/KirF6P) * (MG/KirF6PMg)
                      + (P/KirP) * (F6P/KirPF6P)
                      + (P/KirP) * (F6P/KirPF6P) * (MG/KirPF6PMg)
                      + (FDP - MgFDP) / KirFDP
                      + (KdFDPMg/KmrMg) * (MgFDP/KirFDP)
                      + AMP/KirAMP
                      + MgFDP/KirFDPMg
                      + (MgFDP/KirFDPMg) * (MG/KirFDPMgMg)
                      + (AMP/KirAMP) * (FDP - MgFDP) / KirAMPFDP);
  real numerator = num_num / num_denom;
  real denominator = 1 + L0 * pow(denom_num / denom_denom, n);
  return numerator / denominator;
}

real flux_GDH(real BPG, real GAP, real Keq, real KmBPG, real KmGAP, real KmNAD,
              real KmNADH, real KmP, real NAD, real NADH, real P, real Vmax){
  real numerator = Vmax * (P * GAP * NAD - (BPG * NADH/Keq));
  real denominator = (-1
                      + (1 + P/KmP) * (1 + GAP/KmGAP) * (1 + NAD/KmNAD)
                      + (1 + BPG/KmBPG) * (1 + NADH/KmNADH));
  return numerator / denominator;
}

real flux_GL6P_HYDROLYSIS_1(real GL6P, real KGI6Phydrol, real KeqGI6Phydrol, real PGN){
  return KGI6Phydrol * (GL6P - PGN/KeqGI6Phydrol);
}

real flux_GPM(real PGA3, real PGA2, real Keq, real KmPGA2, real KmPGA3, real Vmax){
  real numerator = Vmax * (PGA3 - PGA2/Keq) / KmPGA3;
  real denominator = 1 + PGA3/KmPGA3 + PGA2/KmPGA2;
  return numerator / denominator;
}

real flux_PFK(real ATP, real F6P, real FDP, real KefrADP, real KefrPEP, real KeftADP,
              real KeftPEP, real Keq, real KirADP, real KirATP, real KirF6P, real KirFDP,
              real KitADP, real KitATP, real KitF6P, real KitFDP, real KmrADP, real KmrATPMg,
              real KmrF6P, real KmrFDP, real KmtADP, real KmtATPMg, real KmtF6P, real KmtFDP,
              real L0, real MgADP, real MgATP, real PEP, real Vmax, real Wr, real Wt, real n){
  real numnum =  Vmax * n * (MgATP * F6P - MgADP * FDP/Keq)/(KirF6P*KmrATPMg);
  real numdenom = (1
                   + KmrFDP/KirFDP*(MgADP/KmrADP)
                   + KmrF6P/KirF6P*(MgATP/KmrATPMg)
                   + KmrFDP/KirFDP*(MgADP/KmrADP)*(F6P/KirF6P)
                   + MgATP/KmrATPMg*(F6P/KirF6P)
                   + MgADP/KirADP*(MgATP/KmrATPMg)*(F6P/KirF6P)
                   + (1+(ATP-MgATP)/KirATP)*(F6P/KirF6P)
                   + FDP/KirFDP
                   + MgADP/KmrADP*(FDP/KirFDP)
                   + KmrF6P/KirF6P*(MgATP/KmrATPMg)*(FDP/KirFDP)
                   + Wr*(KmrF6P/KirF6P)*(MgADP/KirADP)*(MgATP/KmrATPMg)*(FDP/KmrFDP));
  real numerator = numnum / numdenom;
  real denomnum = (1
                   +KmtFDP/KitFDP*(MgADP/KmtADP)
                   + KmtF6P/KitF6P*(MgATP/KmtATPMg)
                   + KmtFDP/KitFDP*(MgADP/KmtADP)*(F6P/KitF6P)
                   + MgATP/KmtATPMg*(F6P/KitF6P)
                   + MgADP/KitADP*(MgATP/KmtATPMg)*(F6P/KitF6P)
                   + (1+(ATP-MgATP)/KitATP)*(F6P/KitF6P)
                   + FDP/KitFDP
                   + MgADP/KmtADP*(FDP/KitFDP)
                   + KmtF6P/KitF6P*(MgATP/KmtATPMg)*(FDP/KitFDP)
                   + Wt*(KmtF6P/KitF6P)*(MgADP/KitADP)*(MgATP/KmtATPMg)*(FDP/KmtFDP))
    * (1+MgADP/KeftADP+PEP/KeftPEP+MgADP/KeftADP*(PEP/KeftPEP));
  real denomdenom = (1
                     + KmrFDP/KirFDP*(MgADP/KmrADP)
                     + KmrF6P*MgATP/(KirF6P*KmrATPMg)
                     + KmrFDP/KirFDP*(MgADP/KmrADP)*(F6P/KirF6P)
                     + MgATP/KmrATPMg*(F6P/KirF6P)
                     + MgADP/KirADP*(MgATP/KmrATPMg)*(F6P/KirF6P)
                     + (1+(ATP-MgATP)/KirATP)*(F6P/KirF6P)
                     + FDP/KirFDP
                     + MgADP/KmrADP*(FDP/KirFDP)
                     + KmrF6P/KirF6P*(MgATP/KmrATPMg)*(FDP/KirFDP)
                     + Wr*(KmrF6P/KirF6P)*(MgADP/KirADP)*(MgATP/KmrATPMg)*(FDP/KmrFDP))
    * (1+MgADP/KefrADP+PEP/KefrPEP+MgADP/KefrADP*(PEP/KefrPEP));
  real denominator = 1 + L0 * pow(denomnum / denomdenom, n);

  return numerator / denominator;
}

real flux_PGI_1(real G6P, real F6P, real Keq, real KmF6P, real KmG6P,
                real KmPEP, real KmPGN, real PEP, real PGN, real Vmax){
  real numerator = Vmax*(G6P-F6P/Keq)/KmG6P;
  real denominator = (1+F6P/KmF6P+G6P/KmG6P+PEP/KmPEP+PGN/KmPGN);
  return numerator / denominator;
}

real flux_PGK(real BPG, real Keq, real KmADPMg, real KmATPMg, real KmBPG,
              real KmPGA3, real MgADP, real MgATP, real PGA3, real Vmax){
  real numerator = Vmax*(MgADP*BPG-MgATP*PGA3/Keq)/(KmADPMg*KmBPG);
  real denominator = (1
                      + MgADP/KmADPMg
                      + BPG/KmBPG
                      + MgADP/KmADPMg*BPG/KmBPG
                      + MgATP/KmATPMg
                      + PGA3/KmPGA3
                      + MgATP/KmATPMg*PGA3/KmPGA3);
  return numerator / denominator;
}

real flux_PPS(real ADP, real AKG, real AMP, real ATP, real KdADPMg, real KdAMP,
              real KdATPMg, real KdATPMgPPS, real KdMg, real KdP, real KdPEP,
              real KdPYR, real KefADP, real KefAKG, real KefATP, real KefOAA,
              real Keq, real KmAMP, real KmATPMg, real KmP, real KmPEP,
              real KmPYR, real MG, real MgADP, real MgATP, real OAA, real P,
              real PEP, real PYR, real Vmax, real W, real alpha){
  real numerator = Vmax*(MgATP*PYR-AMP*PEP*P*MG/Keq)/(KmATPMg*KmPYR);
  real denominator = (MgATP/KmATPMg
                      + alpha*(P/KdP)*(MgATP/KmATPMg)
                      + alpha*(AMP/KdAMP)*(MgATP/KmATPMg)
                      + alpha*(P/KdP)*(AMP/KdAMP)*(MgATP/KmATPMg)
                      + alpha*(MG/KdMg)*(P/KmP)*(AMP/KdAMP)*(MgATP/KdATPMgPPS)/(W*(1+MG/KdMg))
                      + MgATP/KmATPMg*(AKG/KefAKG)
                      + (1+MG/KdMg)*(AKG/KefAKG)*(PEP/KmPEP)/W
                      + MgATP/KmATPMg*(OAA/KefOAA)
                      + (1+MG/KdMg)*(OAA/KefOAA)*(PEP/KmPEP)/W
                      + MG/KdMg*(P/KmP)*(AMP/KdAMP)/W
                      + alpha*(P/KdP)*(AMP/KdAMP)*(PEP/KmPEP)/W
                      + alpha*(MG/KdMg)*(P/KmP)*(AMP/KdAMP)*(PEP/KmPEP)/W
                      + alpha*(1+MG/KdMg)*(KmAMP/KdAMP*(P/KmP)*(PEP/KmPEP)+AMP/KdAMP*(PEP/KmPEP))/W
                      + (1+MG/KdMg)*(PYR/KmPYR)
                      + MgATP/KmATPMg*(PYR/KmPYR)
                      + KdADPMg/KdMg*(P/KmP)*(MgADP/KefADP)*(AMP/KdAMP)/(W*(1+MG/KdMg))
                      + (ADP-MgADP)/KefADP*(PYR/KmPYR)
                      + KdATPMg/KdMg*(P/KmP)*(AMP/KdAMP)*(MgATP/KefATP)/(W*(1+MG/KdMg))
                      + (ATP-MgATP)/KefATP*(PYR/KmPYR)
                      + (1+MG/KdMg)*(PEP/KmPEP)/W
                      + alpha*(1+MG/KdMg)*(PEP/KdPEP)*(PYR/KmPYR)
                      + (1+MG/KdMg)*(PYR/KdPYR)*(PEP/KmPEP)/W);
  return numerator / denominator;
}

real flux_PTS_0(real KmPEP, real KmPYR, real PEP, real PYR, real ei, real eiP,
                real kF, real kR){
  return (kF*ei*pow(PEP, 2)/(pow(KmPEP, 2)+pow(PEP,2))
          - kR*eiP*pow(PYR,2)/(pow(KmPYR, 2)+pow(PYR,2)));
}

real flux_PTS_4_1(real G6P, real GLCp, real KmG6P, real KmGLC, real cell,
                  real eiicb, real eiicbP, real kF, real kR){
  return cell * (kF*eiicbP*GLCp/(KmGLC+GLCp)
                 -kR*eiicb*G6P/(KmG6P+G6P));
}

real flux_PYK(real ADP, real FDP, real G6P, real GL6P, real KefrFDP,
              real KefrG6P, real KefrGL6P, real KefrR5P, real KefrRU5P,
              real KefrS7P, real KefrX5P, real KeftATP, real KeftSUCCOA,
              real KirADP, real KirATP, real KirPEP, real KirPYR,
              real KirPyrATP, real KitADP, real KitATP, real KitPEP, real KitPYR,
              real KitPyrATP, real KmrADPMg, real KmrPEP, real KmtADPMg,
              real KmtPEP, real L0, real MgADP, real MgATP, real PEP, real PYR,
              real r5P, real RU5P, real S7P, real SUCCOA, real Vmax, real X5P,
              real n){
  real numnum = Vmax*n*PEP*MgADP/(KirPEP*KmrADPMg);
  real numdenom = (1
                   + KmrPEP/KirPEP*(MgADP/KmrADPMg)
                   + MgATP/KirATP
                   + MgADP/KmrADPMg*(PEP/KirPEP)
                   + KmrADPMg/KmrADPMg*(1+(ADP-MgADP)/KirADP)*(PEP/KirPEP)
                   + PYR/KirPYR
                   + MgATP/KirPyrATP*(PYR/KirPYR));
  real numerator = numnum / numdenom;
  real denomnum = ((1
                    + KmtPEP/KitPEP*(MgADP/KmtADPMg)
                    + MgATP/KitATP
                    + MgADP*PEP/(KitPEP*KmtADPMg)
                    + (1+(ADP-MgADP)/KitADP)*(PEP/KitPEP)
                    + PYR/KitPYR
                    + MgATP/KitPyrATP*(PYR/KitPYR))
                   * (1
                      + SUCCOA/KeftSUCCOA
                      + MgATP*SUCCOA/(KeftATP*KeftSUCCOA)));
  real denomdenom = ((1
                      + KmrPEP/KirPEP*(MgADP/KmrADPMg)
                      + MgATP/KirATP
                      + MgADP/KmrADPMg*(PEP/KirPEP)
                      + (1+(ADP-MgADP)/KirADP)*(PEP/KirPEP)
                      + PYR/KirPYR
                      + MgATP/KirPyrATP*(PYR/KirPYR))
                     * (1
                        + FDP/KefrFDP
                        + G6P/KefrG6P
                        + GL6P/KefrGL6P
                        + R5P/KefrR5P
                        + RU5P/KefrRU5P
                        + S7P/KefrS7P
                        + X5P/KefrX5P));
  real denominator = 1 + L0 * pow(denomnum / denomnum, n);
  return numerator / denominator;
}

real flux_TPI(real DAP, real GAP, real Keq, real KmDAP, real KmGAP, real Vmax){
  return Vmax*(DAP-GAP/Keq)/KmDAP/(1+DAP/KmDAP+GAP/KmGAP);
}

real flux_XCH_RMM_1(real GLCp, real GLCx, real Km, real Vmax){
  return Vmax*(GLCx/Km-GLCp/Km)/(1+GLCx/Km+GLCp/Km);
}
