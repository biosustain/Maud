int get_N_enz_mics(matrix S, int i_enz){
  int out = 0;
  for (s in S[,i_enz]){
    if (s != 0){
      out += 1;
    }
  }
  return out;
}

int[] get_enz_mics(matrix S, int i_enz){
  int N_enz_mics = get_N_enz_mics(S, i_enz);
  int out[N_enz_mics];
  int ticker_pos = 1;
  for (i_met in 1:rows(S)){
    if (S[i_met, i_enz] != 0){
      out[ticker_pos] = i_met;
      ticker_pos += 1;
    }
  }
  return out;
}

vector get_flux_enz(vector conc_mic,
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
                vector transfer_constant,
                int[] subunits){
  vector[cols(S)] flux_enz;
  int pos_ci = 1;
  int pos_ai = 1;
  int pos_aa = 1;
  int pos_tc = 1;
  for (i_enz in 1:cols(S)){
    int N_enz_mics = get_N_enz_mics(S, i_enz);
    int enz_mics[N_enz_mics] = get_enz_mics(S, i_enz);
    vector[N_enz_mics] km_enz = km[km_lookup[enz_mics, i_enz]];
    int nci_enz = n_ci[i_enz];
    int nai_enz = n_ai[i_enz];
    int naa_enz = n_aa[i_enz];
    vector[nci_enz] conc_ci;
    vector[nai_enz] conc_ai;
    vector[naa_enz] conc_aa;
    vector[nci_enz] ki_enz;
    vector[nai_enz] diss_t_enz;
    vector[naa_enz] diss_r_enz;
    real catalysis_factor;
    real allostery_factor = 1;
    real free_enzyme_ratio;
    int is_allosteric = ((nai_enz > 0) || (naa_enz > 0));
    if (nci_enz != 0){
      conc_ci = conc_mic[segment(ci_ix, pos_ci, nci_enz)];
      ki_enz = segment(ki, pos_ci, nci_enz);
    }
    if (nai_enz != 0){
      conc_ai = conc_mic[segment(ai_ix, pos_ai, nai_enz)];
      diss_t_enz = segment(dissociation_constant_t, pos_ai, nai_enz);
    }
    if (naa_enz != 0){
      conc_aa = conc_mic[segment(aa_ix, pos_aa, naa_enz)];
      diss_r_enz = segment(dissociation_constant_r, pos_aa, naa_enz);
    }
    free_enzyme_ratio = get_free_enzyme_ratio(conc_mic[enz_mics],
                                              km_enz,
                                              S[enz_mics, i_enz],
                                              conc_ci,
                                              ki_enz);
    catalysis_factor = modular_rate_law(conc_mic[enz_mics],
                                        km_enz,
                                        S[enz_mics, i_enz],
                                        kcat[i_enz],
                                        keq[i_enz],
                                        conc_enz[i_enz],
                                        conc_ci,
                                        ki_enz);
    if (is_allosteric == 1){
      allostery_factor = get_allostery(conc_aa,
                                       conc_ai,
                                       free_enzyme_ratio,
                                       diss_r_enz,
                                       diss_t_enz,
                                       transfer_constant[pos_tc],
                                       subunits[i_enz]);
    }
    flux_enz[i_enz] = catalysis_factor * allostery_factor;
    pos_ci += n_ci[i_enz];
    pos_ai += n_ai[i_enz];
    pos_aa += n_aa[i_enz];
    pos_tc += is_allosteric;
  }
  return flux_enz;
}

vector get_flux_drain(matrix S_drain,
                      vector conc_mic,
                      vector drain){
  vector[cols(S_drain)] flux_drain;
  for (i_drain in 1:cols(S_drain)){
    int N_drain_mics = get_N_enz_mics(S_drain, i_drain);
    int drain_mics[N_drain_mics] = get_enz_mics(S_drain, i_drain);
    real drain_rate = drain_reaction(conc_mic[drain_mics], drain[i_drain]);
    flux_drain[i_drain] = drain_rate;
  }
  return flux_drain;
}

vector dbalanced_dt(real time,
                    vector current_balanced,
                    vector unbalanced,
                    int[] balanced_ix,
                    int[] unbalanced_ix,
                    vector enzyme_concentration,
                    vector km,
                    int[,] km_lookup,
                    matrix S_enz,
                    matrix S_drain,
                    matrix S_full,
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
                    vector transfer_constant,
                    int[] subunits,
                    vector drain){
  vector[rows(current_balanced)+rows(unbalanced)] current_concentration;
  vector[cols(S_enz)] flux_enz;
  vector[cols(S_drain)] flux_drain;

  current_concentration[balanced_ix] = current_balanced;
  current_concentration[unbalanced_ix] = unbalanced;

  flux_enz = get_flux_enz(current_concentration,
                          enzyme_concentration,
                          km,
                          km_lookup,
                          S_enz,
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
                          transfer_constant,
                          subunits);
  flux_drain = get_flux_drain(S_drain,
                              current_concentration,
                              drain);

  return S_full[balanced_ix] * append_row(flux_enz, flux_drain);
}
