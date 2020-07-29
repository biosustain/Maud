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
		int[] n_ci,
		vector ki){
  vector[cols(S)] flux;
  int pos_ci = 1;
  for (i_rxn in 1:cols(S)){
    int N_rxn_mics = get_N_rxn_mics(S, i_rxn);
    int rxn_mics[N_rxn_mics] = get_rxn_mics(S, i_rxn);
    vector[N_rxn_mics] km_rxn = km[km_lookup[rxn_mics, i_rxn]];
    int nci_rxn = n_ci[i_rxn];
    vector[nci_rxn] conc_ci;
    vector[nci_rxn] ki_rxn;
    if (nci_rxn != 0){
      conc_ci = conc_mic[segment(ci_ix, pos_ci, nci_rxn)];
      ki_rxn = segment(ki, pos_ci, nci_rxn);
    }
    flux[i_rxn] = modular_rate_law(conc_mic[rxn_mics],
				   km_rxn,
				   S[rxn_mics, i_rxn],
				   kcat[i_rxn],
				   keq[i_rxn],
				   conc_enz[i_rxn],
				   conc_ci,
				   ki_rxn);
    pos_ci += n_ci[i_rxn];
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
		    int[] n_ci,   // number of competitive inhibitors for each reaction
		    vector ki){   // vector of competitive inhibition constants
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
				   n_ci,
				   ki);
}
