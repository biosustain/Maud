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
vector get_dgrs(matrix S, vector dgf, int[] mic_to_met, vector water_stoichiometry){
  /*
      Calculate dgrs from metabolite formation energies, assuming water's
      formation energy is known exactly.
  */
  real minus_RT = -0.008314 * 298.15;
  real dgf_water = -157.6;  // From http://equilibrator.weizmann.ac.il/metabolite?compoundId=C00001
  vector[cols(S)] delta_g = S' * dgf[mic_to_met] + water_stoichiometry * dgf_water;
  return delta_g;
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

real get_met_km_ratio(real metabolite,
                      real km){
  /* Returns the ratio of metabolite concentration to km */
  return metabolite / km;
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

int get_n_sub_for_edge(matrix S, int j){
  /*
    Get the number of active metabolites-in-compartment for a reaction or drain j,
    given stoichiometric matrix S.
  */
  int out = 0;
  for (s in S[,j]){
    if (s < 0){
      out += 1;
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

int[] get_substrate_for_edge(matrix S, int j){
  /*
    Get an array of active metabolites-in-compartment for a reaction or drain j, given
    stoichiometric matrix S.
  */
  int N_mic = get_n_sub_for_edge(S, j);
  int out[N_mic];
  int pos = 1;
  for (i in 1:rows(S)){
    if (S[i, j] < 0){
      out[pos] = i;
      pos += 1;
    }
  }
  return out;
}

real get_reversibility(real dgrs, real reaction_quotient){
  /*
    Get a reversibility for single reaction given a standard gibbs energy
    and a reaction quotion.
  */
  real RT = 0.008314 * 298.15;
  real dgr = dgrs + RT*reaction_quotient;
  return 1-exp(dgr/RT);
}
real substrate_km_product(
  /*
    gets the numerate substate term by taking the product of
    concentration over the km values.
  */
  vector substrate_concs,
  vector substrate_kms
){
  return prod(substrate_concs./substrate_kms);
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
                       vector dgrs,
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
    print("dgrs: ", dgrs);
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

vector get_reaction_quotient(matrix S, vector conc){
  return S' * log(conc);
}

vector get_flux(vector conc,
                vector enz,
                vector km,
                vector drain,
                int[,] km_lookup,
                matrix S,
                int[] edge_type,
                int[] edge_to_drain,
                int[] edge_to_enzyme,
                vector kcat,
                vector dgrs,
                int[] ix_ci,
                int[] ix_ai,
                int[] ix_aa,
                int[] ix_pa,
                int[] ix_pi,
                int[] n_ci,
                int[] n_ai,
                int[] n_aa,
                int[] n_pa,
                int[] n_pi,
                vector ki,
                vector diss_t,
                vector diss_r,
                vector transfer_constant,
                int[] subunits,
                vector phos_kcat,
                vector phos_conc){
  vector[cols(S)] out;
  int pos_ci = 1;
  int pos_ai = 1;
  int pos_aa = 1;
  int pos_tc = 1;
  int pos_pa = 1;
  int pos_pi = 1;
  vector[cols(S)] reaction_quotient = get_reaction_quotient(S, conc);
  for (j in 1:cols(S)){
    int n_mic_j = get_n_mic_for_edge(S, j, edge_type[j]);
    int n_sub_j = get_n_sub_for_edge(S, j);
    int mics_j[n_mic_j] = get_mics_for_edge(S, j, edge_type[j]);
    int sub_j[n_sub_j] = get_substrate_for_edge(S, j);
    if (edge_type[j] == 1){  // reversible enzyme...
      vector[n_mic_j] km_j = km[km_lookup[mics_j, j]];
      vector[n_sub_j] km_j_substrate = km[km_lookup[sub_j, j]];
      real kcat_j = kcat[edge_to_enzyme[j]];
      real free_enzyme_ratio_denom = get_Dr_common_rate_law(conc[mics_j], km_j, S[mics_j, j]);
      if (n_ci[j] > 0){  // competitive inhibition
        int comp_inhs_j[n_ci[j]] = segment(ix_ci, pos_ci, n_ci[j]);
        vector[n_ci[j]] ki_j = segment(ki, pos_ci, n_ci[j]);
        free_enzyme_ratio_denom += sum(conc[comp_inhs_j] ./ ki_j);
        pos_ci += n_ci[j];
      }
      real free_enzyme_ratio = inv(free_enzyme_ratio_denom);
      real saturation_term = exp(log(substrate_km_product(conc[sub_j], km_j_substrate)) - log(free_enzyme_ratio_denom));
      real reversible_term = get_reversibility(dgrs[j], reaction_quotient[j]);
      out[j] = enz[edge_to_enzyme[j]] * kcat[edge_to_enzyme[j]] * saturation_term * reversible_term;
      if ((n_ai[j] > 0) || (n_aa[j] > 0)){  // allosteric regulation
        real Q_num = 1;
        real Q_denom = 1;
        if (n_ai[j] > 0){
          int allo_inhs_j[n_ai[j]] = segment(ix_ai, pos_ai, n_ai[j]);
          vector[n_ai[j]] diss_t_j = segment(diss_t, pos_ai, n_ai[j]);
          Q_num += sum(conc[allo_inhs_j] ./ diss_t_j);
          pos_ai += n_ai[j];
        }
        if (n_aa[j] > 0){
          int allo_acts_j[n_aa[j]] = segment(ix_aa, pos_aa, n_aa[j]);
          vector[n_aa[j]] diss_r_j = segment(diss_r, pos_aa, n_aa[j]);
          Q_denom += sum(conc[allo_acts_j] ./ diss_r_j);
          pos_aa += n_aa[j];
        }
        out[j] *= inv(1 + transfer_constant[pos_tc] * (free_enzyme_ratio * Q_num / Q_denom) ^ subunits[j]);
        pos_tc += 1;
      }
      if ((n_pi[j] > 0) || (n_pa[j] > 0)){  // phosphorylation
        real alpha = 0;
        real beta = 0;
        if (n_pa[j] > 0){
          int phos_acts_j[n_pa[j]] = segment(ix_pa, pos_pa, n_pa[j]);
          beta = sum(phos_kcat[phos_acts_j] .* phos_conc[phos_acts_j]);
          pos_pa += n_pa[j];
        }
        if (n_pi[j] > 0){
          int phos_inhs_j[n_pi[j]] = segment(ix_pi, pos_pi, n_pi[j]);
          alpha = sum(phos_kcat[phos_inhs_j] .* phos_conc[phos_inhs_j]);
          pos_pi += n_pi[j];
        }
        out[j] *= 1 / (1 + (alpha / beta) ^ subunits[j]);  // TODO: what if beta is zero and alpha is non-zero?
      }
    }
    else if (edge_type[j] == 2){  // drain...
      out[j] = drain[edge_to_drain[j]] * prod(conc[mics_j] ./ (conc[mics_j] + 1e-6));
    }
    else if (edge_type[j] == 3){  // irreversible modular rate law...
      vector[n_mic_j] km_j = km[km_lookup[mics_j, j]];
      vector[n_sub_j] km_j_substrate = km[km_lookup[sub_j, j]];
      real kcat_j = kcat[edge_to_enzyme[j]];
      real free_enzyme_ratio_denom = get_Dr_common_rate_law_irreversible(conc[mics_j], km_j, S[mics_j, j]);
      
      if (n_ci[j] > 0){  // competitive inhibition
        int comp_inhs_j[n_ci[j]] = segment(ix_ci, pos_ci, n_ci[j]);
        vector[n_ci[j]] ki_j = segment(ki, pos_ci, n_ci[j]);
        free_enzyme_ratio_denom += sum(conc[comp_inhs_j] ./ ki_j);
        pos_ci += n_ci[j];
      }
      real free_enzyme_ratio = inv(free_enzyme_ratio_denom);
      real saturation_term = exp(log(substrate_km_product(conc[sub_j], km_j_substrate)) - log(free_enzyme_ratio_denom));
      out[j] = enz[edge_to_enzyme[j]] * kcat[edge_to_enzyme[j]] * saturation_term;
      if ((n_ai[j] > 0) || (n_aa[j] > 0)){  // allosteric regulation
        real Q_num = 1;
        real Q_denom = 1;
        if (n_ai[j] > 0){
          int allo_inhs_j[n_ai[j]] = segment(ix_ai, pos_ai, n_ai[j]);
          vector[n_ai[j]] diss_t_j = segment(diss_t, pos_ai, n_ai[j]);
          Q_num += sum(conc[allo_inhs_j] ./ diss_t_j);
          pos_ai += n_ai[j];
        }
        if (n_aa[j] > 0){
          int allo_acts_j[n_aa[j]] = segment(ix_aa, pos_aa, n_aa[j]);
          vector[n_aa[j]] diss_r_j = segment(diss_r, pos_aa, n_aa[j]);
          Q_denom += sum(conc[allo_acts_j] ./ diss_r_j);
          pos_aa += n_aa[j];
        }
        out[j] *= inv(1 + transfer_constant[pos_tc] * (free_enzyme_ratio * Q_num / Q_denom) ^ subunits[j]);
        pos_tc += 1;
      }
      if ((n_pi[j] > 0) || (n_pa[j] > 0)){  // phosphorylation
        real alpha = 0;
        real beta = 0;
        if (n_pa[j] > 0){
          int phos_acts_j[n_pa[j]] = segment(ix_pa, pos_pa, n_pa[j]);
          beta = sum(phos_kcat[phos_acts_j] .* phos_conc[phos_acts_j]);
          pos_pa += n_pa[j];
        }
        if (n_pi[j] > 0){
          int phos_inhs_j[n_pi[j]] = segment(ix_pi, pos_pi, n_pi[j]);
          alpha = sum(phos_kcat[phos_inhs_j] .* phos_conc[phos_inhs_j]);
          pos_pi += n_pi[j];
        }
        out[j] *= 1 / (1 + (alpha / beta) ^ subunits[j]);  // TODO: what if beta is zero and alpha is non-zero?
      }
    }
    else reject("Unknown edge type ", edge_type[j]);
  }
  return out;
}
vector dbalanced_dt(real time,
                    vector current_balanced,
                    vector unbalanced,
                    int[] balanced_ix,
                    int[] unbalanced_ix,
                    vector enz,
                    vector km,
                    vector drain,
                    int[,] km_lookup,
                    matrix S,
                    int[] edge_type,
                    int[] edge_to_drain,
                    int[] edge_to_enzyme,
                    vector kcat,
                    vector dgrs,
                    int[] ix_ci,
                    int[] ix_ai,
                    int[] ix_aa,
                    int[] ix_pa,
                    int[] ix_pi,
                    int[] n_ci,
                    int[] n_ai,
                    int[] n_aa,
                    int[] n_pa,
                    int[] n_pi,
                    vector ki,
                    vector diss_t,
                    vector diss_r,
                    vector transfer_constant,
                    int[] subunits,
                    vector kcat_phos,
                    vector conc_phos){
  vector[rows(current_balanced)+rows(unbalanced)] current_concentration;
  current_concentration[balanced_ix] = current_balanced;
  current_concentration[unbalanced_ix] = unbalanced;
  vector[rows(S)] flux = get_flux(current_concentration,
                                  enz, km, drain, km_lookup, S, edge_type, edge_to_drain, edge_to_enzyme, kcat, dgrs,
                                  ix_ci, ix_ai, ix_aa, ix_pa, ix_pi, n_ci, n_ai, n_aa, n_pa, n_pi,
                                  ki, diss_t, diss_r, transfer_constant, subunits, kcat_phos, conc_phos);
  return (S * flux)[balanced_ix];
}
