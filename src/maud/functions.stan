vector unz_1d(vector[] mnsd, vector z){
  /* 
     Recover a real-valued vector from a 1xn vector z of z scores and a 2xn
     array mnsd of vectors consisting of the mean and standard deviation of
     each element. 
  */
  return mnsd[1] + mnsd[2] .* z;
}
vector unz_log_1d(vector[] mnsd, vector z){
  /* 
     Recover a positive-constrained vector from a 1xn vector z of z scores and
     a 2xn array mnsd of vectors consisting of the lognormal mean and standard
     deviation of each element. 
  */
  return exp(log(mnsd[1]) + mnsd[2] .* z);
}
vector[] unz_2d(vector[,] mnsd, vector[] z){
  /*
    Recover a mxn array of real-valued vectors from an mxn array z of vectors
    of z scores and a 2xmxn array mnsd of vectors consisting of mean and
    standard deviation arrays.
  */
  array[size(z)] vector [rows(z[1])] out;
  for (ex in 1:size(z)){
    out[ex] = unz_1d(mnsd[:,ex], z[ex]);
  }
  return out;
}
vector[] unz_log_2d(vector[,] mnsd, vector[] z){
  /*
    Recover a mxn array of positive-constrained vectors from an mxn array z of
    vectors of z scores and a 2xmxn array mnsd of vectors consisting of
    lognormal mean and standard deviation arrays.
  */
  array[size(z)] vector [rows(z[1])] out;
  for (ex in 1:size(z)){
    out[ex] = unz_log_1d(mnsd[:,ex], z[ex]);
  }
  return out;
}
vector get_keq(matrix S, vector dgf, int[] mic_to_met, vector water_stoichiometry){
  /*
      Calculate keqs from metabolite formation energies, assuming water's
      formation energy is known exactly.
  */
  real minus_RT = -0.008314 * 298.15;
  real dgf_water = -157.6;  // From http://equilibrator.weizmann.ac.il/metabolite?compoundId=C00001
  vector[cols(S)] delta_g = S' * dgf[mic_to_met] + water_stoichiometry * dgf_water;
  return exp(delta_g / minus_RT);
}
real get_Tr(vector metabolite,
            vector km,
            vector stoichiometry,
            real kcat,
            real keq){
  /* Tr coefficient in the modular rate law. */
  real plus_product = 1;
  real minus_product = 1;
  real k_minus = (kcat / keq);
  for (m in 1:rows(metabolite)){
    real multiplicand = (metabolite[m] / km[m]) ^ abs(stoichiometry[m]);
    k_minus *= km[m] ^ stoichiometry[m];
    if (stoichiometry[m] < 0)
      plus_product *= multiplicand;
    else
      minus_product *= multiplicand;
  }
  return kcat * plus_product - k_minus * minus_product;
}
real get_Tr_irreversible(vector metabolite,
                         vector km,
                         vector stoichiometry,
                         real kcat){
  /* Tr coefficient in the irreversible modular rate law. */
  real plus_product = 1;
  for (m in 1:rows(metabolite)){
    if (stoichiometry[m] < 0){
      real multiplicand = (metabolite[m] / km[m]) ^ abs(stoichiometry[m]);
      plus_product *= multiplicand;
    }
  }
  return kcat * plus_product;
}
real get_Dr_common_rate_law(vector metabolite, vector km, vector stoichiometry){
  /* Dr coefficient in the modular rate law. */
  real psi_plus = 1;
  real psi_minus = 1;
  for (m in 1:rows(metabolite)){
    real multiplicand = (1 + metabolite[m] / km[m]) ^ abs(stoichiometry[m]);
    if (stoichiometry[m] < 0)
      psi_plus *= multiplicand;
    else
      psi_minus *= multiplicand;
  }
  return psi_plus + psi_minus - 1;
}
real get_Dr_common_rate_law_irreversible(vector metabolite, vector km, vector stoichiometry){
  /* Dr coefficient in the modular rate law negating products. */
  real psi_plus = 1;
  real psi_minus = 1;
  for (m in 1:rows(metabolite)){
    if (stoichiometry[m] < 0){
      real multiplicand = (1 + metabolite[m] / km[m]) ^ abs(stoichiometry[m]);
      psi_plus *= multiplicand;
    }
  }
  return psi_plus;
}
int get_n_mic_for_edge(matrix S, int j, int edge_type){
  /*
    Get the number of active metabolites-in-compartment for a reaction or drain j,
    given stoichiometric matrix S.
  */
  int out = 0;
  for (s in S[,j]){
    if (edge_type != 3){ // not irreversible
      if (s != 0){
        out += 1;
      }
    }
    else {
      if (s < 0){
        out += 1;
      }
    }
  }
  return out;
}
int[] get_mics_for_edge(matrix S, int j, int edge_type){
  /*
    Get an array of active metabolites-in-compartment for a reaction or drain j, given
    stoichiometric matrix S.
  */
  int N_mic = get_n_mic_for_edge(S, j, edge_type);
  int out[N_mic];
  int pos = 1;
  for (i in 1:rows(S)){
    if (edge_type != 3){ // not irreversible
      if (S[i, j] != 0){
        out[pos] = i;
        pos += 1;
      }
    }
    else{
      if (S[i, j] < 0){
        out[pos] = i;
        pos += 1;
      }
    }
  }
  return out;
}

int check_steady_state(vector[] conc_balanced,
                       int e,
                       vector flux,
                       vector conc_init,
                       real[] timepoints,
                       vector conc_unbalanced,
                       vector conc_enzyme_experiment,
                       vector km,
                       vector drain,
                       vector kcat,
                       vector keq,
                       vector ki,
                       vector diss_t,
                       vector diss_r,
                       vector transfer_constant,
                       vector kcat_phos,
                       vector conc_phos_experiment){
  if ((max(fabs(conc_balanced[1]-conc_balanced[2])./conc_balanced[2]) > 0.001)){
    print("");
    print("Non-steady state in experiment ", e, ".");
    print("Balanced metabolite concentration at ", timepoints[1], " seconds: ", conc_balanced[1]);
    print("Balanced metabolite concentration at ", timepoints[2], " seconds: ", conc_balanced[2]);
    print("flux: ", flux);
    print("conc_init: ", conc_init);
    print("conc_unbalanced: ", conc_unbalanced);
    print("conc_enzyme_experiment: ", conc_enzyme_experiment);
    print("km: ", km);
    print("drain: ", drain);
    print("kcat: ", kcat);
    print("keq: ", keq);
    print("ki: ", ki);
    print("diss_t: ", diss_t);
    print("diss_r: ", diss_r);
    print("transfer_constant: ", transfer_constant);
    print("kcat_phos: ", kcat_phos);
    print("conc_phos_experiment: ", conc_phos_experiment);
    return 0;
  }
  else {
    return 1;
  }
}
