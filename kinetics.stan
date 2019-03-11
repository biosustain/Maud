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
