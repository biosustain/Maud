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

