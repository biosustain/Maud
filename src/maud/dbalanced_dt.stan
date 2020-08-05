int get_N_rxn_mics(matrix S, int i_rxn){
  int out = 0;
  for (s in S[,i_rxn]){
    if (s != 0){
      out += 1;
    }
  }
  return out;
}

int[] get_rxn_mics(matrix S, int i_rxn){
  int N_rxn_mics = get_N_rxn_mics(S, i_rxn);
  int out[N_rxn_mics];
  int ticker_pos = 1;
  for (i_met in 1:rows(S)){
    if (S[i_met, i_rxn] != 0){
      out[ticker_pos] = i_met;
      ticker_pos += 1;
    }
  }
  return out;
}

vector get_flux(vector conc_mic,
		vector conc_enz,
		vector km,
		int[,] km_lookup,
		matrix S,
		vector kcat,
		vector keq,
		int[] ci_ix,
		int[] ai_ix,
		int[] aa_ix,
		int[] n_ci,
		int[] n_ai,
		int[] n_aa,
		vector ki,
                vector dissociation_constant_t,
                vector dissociation_constant_r,
                vector transfer_constant){
  vector[cols(S)] flux;
  int pos_ci = 1;
  int pos_ai = 1;
  int pos_aa = 1;
  int pos_tc = 1;
  for (i_rxn in 1:cols(S)){
    int N_rxn_mics = get_N_rxn_mics(S, i_rxn);
    int rxn_mics[N_rxn_mics] = get_rxn_mics(S, i_rxn);
    vector[N_rxn_mics] km_rxn = km[km_lookup[rxn_mics, i_rxn]];
    int nci_rxn = n_ci[i_rxn];
    int nai_rxn = n_ai[i_rxn];
    int naa_rxn = n_aa[i_rxn];
    vector[nci_rxn] conc_ci;
    vector[nai_rxn] conc_ai;
    vector[naa_rxn] conc_aa;
    vector[nci_rxn] ki_rxn;
    vector[nai_rxn] diss_t_rxn;
    vector[naa_rxn] diss_r_rxn;
    real catalysis_factor;
    real allostery_factor = 1;
    real free_enzyme_ratio;
    int is_allosteric = ((nai_rxn > 0) || (naa_rxn > 0));
    if (nci_rxn != 0){
      conc_ci = conc_mic[segment(ci_ix, pos_ci, nci_rxn)];
      ki_rxn = segment(ki, pos_ci, nci_rxn);
    }
    if (nai_rxn != 0){
      conc_ai = conc_mic[segment(ai_ix, pos_ai, nai_rxn)];
      diss_t_rxn = segment(dissociation_constant_t, pos_ai, nai_rxn);
    }
    if (naa_rxn != 0){
      conc_aa = conc_mic[segment(aa_ix, pos_aa, naa_rxn)];
      diss_r_rxn = segment(dissociation_constant_r, pos_aa, naa_rxn);
    }
    free_enzyme_ratio = get_free_enzyme_ratio(conc_mic[rxn_mics],
                                              km_rxn,
                                              S[rxn_mics, i_rxn],
                                              conc_ci,
                                              ki_rxn);
    catalysis_factor = modular_rate_law(conc_mic[rxn_mics],
                                        km_rxn,
                                        S[rxn_mics, i_rxn],
                                        kcat[i_rxn],
                                        keq[i_rxn],
                                        conc_enz[i_rxn],
                                        conc_ci,
                                        ki_rxn);
    if (is_allosteric == 1){
      allostery_factor = get_allostery(conc_aa,
                                       conc_ai,
                                       free_enzyme_ratio,
                                       diss_r_rxn,
                                       diss_t_rxn,
                                       transfer_constant[pos_tc]);
    }
    flux[i_rxn] = catalysis_factor * allostery_factor;
    pos_ci += n_ci[i_rxn];
    pos_ai += n_ai[i_rxn];
    pos_aa += n_aa[i_rxn];
    pos_tc += is_allosteric;
  }
  return flux;
}

vector dbalanced_dt(real time,
		    vector current_balanced,
		    vector unbalanced,
		    int[] balanced_ix,
		    int[] unbalanced_ix,
		    vector enzyme_concentration,
		    vector km,
		    int[,] km_lookup,
		    matrix S,
		    vector kcat,
		    vector keq,
		    int[] ci_ix,  // index of competitive inhibitors (long-form ragged)
		    int[] ai_ix,  // index of allosteric inhibitors (long-form ragged)
		    int[] aa_ix,  // index of allosteric activators (long-form ragged)
		    int[] n_ci,   // number of competitive inhibitors for each reaction
		    int[] n_ai,   // number of allosteric inhibitors for each reaction
		    int[] n_aa,   // number of allosteric activators for each reaction
		    vector ki,    // vector of competitive inhibition constants
                    vector dissociation_constant_t,
                    vector dissociation_constant_r,
                    vector transfer_constant){
  vector[rows(current_balanced)+rows(unbalanced)] current_concentration;
  current_concentration[balanced_ix] = current_balanced;
  current_concentration[unbalanced_ix] = unbalanced;
  return S[balanced_ix] * get_flux(current_concentration,
				   enzyme_concentration,
				   km,
				   km_lookup,
				   S,
				   kcat,
				   keq,
				   ci_ix,
				   ai_ix,
				   aa_ix,
				   n_ci,
				   n_ai,
				   n_aa,
				   ki,
                                   dissociation_constant_t,
                                   dissociation_constant_r,
                                   transfer_constant);
}
