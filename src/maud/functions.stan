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
int get_n_mic_for_edge(matrix S, int j){
  /*
    Get the number of metabolites-in-compartment for a reaction or drain j,
    given stoichiometric matrix S.
  */
  int out = 0;
  for (s in S[,j]){
    if (s != 0){
      out += 1;
    }
  }
  return out;
}
int[] get_mics_for_edge(matrix S, int j){
  /*
    Get an array of metabolites-in-compartment for a reaction or drain j, given
    stoichiometric matrix S.
  */
  int N_mic = get_n_mic_for_edge(S, j);
  int out[N_mic];
  int pos = 1;
  for (i in 1:rows(S)){
    if (S[i, j] != 0){
      out[pos] = i;
      pos += 1;
    }
  }
  return out;
}
