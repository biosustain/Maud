functions {
  vector unfix(vector free, vector fixed, array[] int ix_free, array[] int ix_fixed){
    int N_out = rows(free) + rows(fixed);
    vector[N_out] out;
    out[ix_free] = free;
    out[ix_fixed] = fixed;
    return out;
  }
  vector unz_1d(array[] vector mnsd, vector z) {
    /*
      Recover a real-valued vector from a 1xn vector z of z scores and a 2xn
      array mnsd of vectors consisting of the mean and standard deviation of
      each element.
    */
    return mnsd[1] + mnsd[2] .* z;
  }

  vector unz_log_1d(array[] vector mnsd, vector z) {
    /*
      Recover a positive-constrained vector from a 1xn vector z of z scores and
      a 2xn array mnsd of vectors consisting of the lognormal mean and standard
      deviation of each element.
    */
    return exp(mnsd[1] + mnsd[2] .* z);
  }

  array[] vector unz_2d(array[,] vector mnsd, array[] vector z) {
    /*
      Recover a mxn array of real-valued vectors from an mxn array z of vectors
      of z scores and a 2xmxn array mnsd of vectors consisting of mean and
      standard deviation arrays.
    */
    array[size(z)] vector[rows(z[1])] out;
    for (ex in 1 : size(z)) {
      out[ex] = unz_1d(mnsd[ : , ex], z[ex]);
    }
    return out;
  }

  array[] vector unz_log_2d(array[,] vector mnsd, array[] vector z) {
    /*
      Recover a mxn array of positive-constrained vectors from an mxn array z of
      vectors of z scores and a 2xmxn array mnsd of vectors consisting of
      lognormal mean and standard deviation arrays.
    */
    array[size(z)] vector[rows(z[1])] out;
    for (ex in 1 : size(z)) {
      out[ex] = unz_log_1d(mnsd[ : , ex], z[ex]);
    }
    return out;
  }

  vector get_dgr(matrix S, vector dgf, real temperature,
                 array[] int mic_to_met, vector water_stoichiometry,
                 vector trans_charge, real psi) {
    /*
        Calculate dgr standard from metabolite formation energies, assuming water's
        formation energy is known exactly.
    */
    real minus_RT = -0.008314 * temperature;
    real dgf_water = -150.9; // From http://equilibrator.weizmann.ac.il/metabolite?compoundId=C00001
    real F = 96.5; // Faraday constant kJ/mol/V
    vector[cols(S)] dgrs = S' * dgf[mic_to_met]
                           + water_stoichiometry * dgf_water
                           + trans_charge * psi * F;
    return dgrs;
  }

  vector get_keq(matrix S, vector dgf, real temperature,
                 array[] int mic_to_met, vector water_stoichiometry,
                 vector trans_charge, real psi) {
    /*
        Calculate keqs from metabolite formation energies, assuming water's
        formation energy is known exactly.
    */
    real minus_RT = -0.008314 * temperature;
    vector[cols(S)] dgrs = get_dgr(S, dgf, temperature, mic_to_met,
                                   water_stoichiometry, trans_charge, psi);
    return exp(dgrs / minus_RT);
  }

  int check_steady_state(vector Sv, vector conc, real abs_thresh,
                         real rel_thresh) {
    /* Relative and absolute check for steady state. */
    vector[rows(conc)] rel_thresh_per_conc = conc * rel_thresh;
    int relative_check_failed = max(abs(Sv) - rel_thresh_per_conc) > 0;
    int absolute_check_failed = max(abs(Sv)) > abs_thresh;
    if (relative_check_failed) {
      print("Sv ", Sv, " not within ", rel_thresh_per_conc, " of zero.");
    }
    if (absolute_check_failed) {
      print("Sv ", Sv, " not within ", abs_thresh, " of zero.");
    }
    return (relative_check_failed || absolute_check_failed) ? 0 : 1;
  }

  int measure_ragged(array[,] int bounds, int i) {
    return bounds[i, 2] - bounds[i, 1] + 1;
  }

  array[] int extract_ragged(array[] int ix_long, array[,] int bounds, int i) {
    /*
      Extract the ith element of a ragged array stored in 1d array long.

      Make sure that members bounds[i, 1] to bounds[i, 2] of long form the
      required element!

     */
    return ix_long[bounds[i, 1] : bounds[i, 2]];
  }

  vector get_saturation(vector conc, vector km, vector free_enzyme_ratio,
                        array[] int sub_km_ix_by_edge_long,
                        array[,] int sub_km_ix_by_edge_bounds,
                        array[] int sub_by_edge_long,
                        array[,] int sub_by_edge_bounds,
                        array[] int edge_type) {
    int N_edge = size(sub_by_edge_bounds);
    vector[N_edge] prod_conc_over_km;
    for (f in 1 : N_edge) {
      if (edge_type[f] == 3) {
        prod_conc_over_km[f] = 1;
        continue;
      }
      int N_sub = measure_ragged(sub_by_edge_bounds, f);
      array[N_sub] int sub_ix = extract_ragged(sub_by_edge_long,
                                               sub_by_edge_bounds, f);
      array[N_sub] int sub_km_ix = extract_ragged(sub_km_ix_by_edge_long,
                                                  sub_km_ix_by_edge_bounds,
                                                  f);
      prod_conc_over_km[f] = prod(conc[sub_ix] ./ km[sub_km_ix]);
    }
    return prod_conc_over_km .* free_enzyme_ratio;
  }

  vector get_free_enzyme_ratio(vector conc, matrix S, vector km, vector ki,
                               array[] int edge_type, array[] int ci_mic_ix,
                               array[] int sub_km_ix_by_edge_long,
                               array[,] int sub_km_ix_by_edge_bounds,
                               array[] int prod_km_ix_by_edge_long,
                               array[,] int prod_km_ix_by_edge_bounds,
                               array[] int sub_by_edge_long,
                               array[,] int sub_by_edge_bounds,
                               array[] int prod_by_edge_long,
                               array[,] int prod_by_edge_bounds,
                               array[] int ci_ix_long,
                               array[,] int ci_ix_bounds) {
    /* Find the proportion of enzyme that is free, for each edge. */
    int N_edge = cols(S);
    vector[N_edge] denom;
    for (f in 1 : N_edge) {
      if (edge_type[f] == 3) {
        // drain
        denom[f] = 1;
        continue;
      }
      int N_sub = measure_ragged(sub_by_edge_bounds, f);
      int N_prod = measure_ragged(prod_by_edge_bounds, f);
      int N_ci = measure_ragged(ci_ix_bounds, f);
      array[N_sub] int sub_ix = extract_ragged(sub_by_edge_long,
                                               sub_by_edge_bounds, f);
      array[N_sub] int sub_km_ix = extract_ragged(sub_km_ix_by_edge_long,
                                                  sub_km_ix_by_edge_bounds,
                                                  f);
      array[N_prod] int prod_ix = extract_ragged(prod_by_edge_long,
                                                 prod_by_edge_bounds, f);
      vector[N_sub] sub_over_km = conc[sub_ix] ./ km[sub_km_ix];
      denom[f] = prod((rep_vector(1, N_sub) + sub_over_km)
                      ^ abs(S[sub_ix, f]));
      if (edge_type[f] == 1) {
        // reversible michaelis menten
        array[N_prod] int prod_km_ix = extract_ragged(prod_km_ix_by_edge_long,
                                                      prod_km_ix_by_edge_bounds,
                                                      f);
        vector[N_prod] prod_over_km = conc[prod_ix] ./ km[prod_km_ix];
        denom[f] += prod((rep_vector(1, N_prod) + prod_over_km)
                         ^ abs(S[prod_ix, f]))
                    - 1;
      }
      if (N_ci > 0) {
        array[N_ci] int ci_ix = extract_ragged(ci_ix_long, ci_ix_bounds, f);
        denom[f] += sum(conc[ci_mic_ix[ci_ix]] ./ ki[ci_ix]);
      }
    }
    return inv(denom);
  }

  vector get_reversibility(vector dgr, real temperature, matrix S,
                           vector conc, array[] int edge_type) {
    real RT = 0.008314 * temperature;
    int N_edge = cols(S);
    vector[N_edge] reaction_quotient = S' * log(conc);
    vector[N_edge] out;
    for (f in 1 : N_edge) {
      if (edge_type[f] == 1) {
        // reversible michaelis menten
        out[f] = 1 - exp((dgr[f] + RT * reaction_quotient[f]) / RT);
      } else {
        out[f] = 1;
      }
    }
    return out;
  }

  vector get_allostery(vector conc, // one per mic
                       vector free_enzyme_ratio, // one per edge
                       vector tc, // one per allosteric enzyme
                       vector dc, // one per allostery
                       vector subunits,
                       // one per edge
                       array[] int allostery_ix_long,
                       // - long and bounds encode a ragged
                       array[,] int allostery_ix_bounds,
                       //   array with one entry per edge
                       array[] int allostery_type,
                       // one per allostery
                       array[] int allostery_mic,
                       // one per allostery
                       array[] int edge_to_tc // one per edge
                       ) {
    int N_edge = size(allostery_ix_bounds);
    vector[N_edge] out = rep_vector(1, N_edge);
    for (f in 1 : N_edge) {
      int N_allostery = measure_ragged(allostery_ix_bounds, f);
      if (N_allostery == 0) {
        continue;
      }
      real Q_num = 1;
      real Q_denom = 1;
      real tc_edge = tc[edge_to_tc[f]];
      for (allostery in extract_ragged(allostery_ix_long,
                                       allostery_ix_bounds, f)) {
        real conc_over_dc = conc[allostery_mic[allostery]] / dc[allostery];
        if (allostery_type[allostery] == 1) {
          // activation
          Q_denom += conc_over_dc;
        } else {
          // inhibition
          Q_num += conc_over_dc;
        }
      }
      out[f] = inv(1
                   + tc_edge
                     * (free_enzyme_ratio[f] * Q_num / Q_denom) ^ subunits[f]);
    }
    return out;
  }

  vector get_phosphorylation(vector kcat_pme, vector conc_pme,
                             array[] int phos_ix_long,
                             array[,] int phos_ix_bounds,
                             array[] int phos_type, array[] int phos_pme,
                             vector subunits) {
    int N_edge = size(phos_ix_bounds);
    vector[N_edge] out = rep_vector(1, N_edge);
    for (f in 1 : N_edge) {
      int N_phos = measure_ragged(phos_ix_bounds, f);
      if (N_phos == 0) {
        continue;
      }
      real alpha = 0;
      real beta = 0;
      for (phos in extract_ragged(phos_ix_long, phos_ix_bounds, f)) {
        real kcat_times_conc = kcat_pme[phos_pme[phos]]
                               * conc_pme[phos_pme[phos]];
        if (phos_type[phos] == 2) {
          // inhibition
          alpha += kcat_times_conc;
        } else {
          beta += kcat_times_conc;
        }
      }
      out[f] = (beta / (alpha + beta)) ^ subunits[f];
    }
    return out;
  }

  vector get_drain_by_edge(vector drain, vector conc,
                           array[] int edge_to_drain,
                           array[] int sub_by_edge_long,
                           array[,] int sub_by_edge_bounds,
                           array[] int edge_type,
                           real drain_small_conc_corrector) {
    int N_edge = size(edge_type);
    vector[N_edge] out = rep_vector(1, N_edge);
    for (f in 1 : N_edge) {
      if (edge_type[f] == 3) {
        int N_sub = measure_ragged(sub_by_edge_bounds, f);
        array[N_sub] int subs = extract_ragged(sub_by_edge_long,
                                               sub_by_edge_bounds, f);
        out[f] = drain[edge_to_drain[f]]
                 * prod(conc[subs]
                        ./ (conc[subs] + drain_small_conc_corrector));
      }
    }
    return out;
  }

  vector get_vmax_by_edge(vector enzyme, vector kcat,
                          array[] int edge_to_enzyme, array[] int edge_type) {
    int N_edge = size(edge_to_enzyme);
    vector[N_edge] out = rep_vector(1, N_edge);
    for (f in 1 : N_edge) {
      if (edge_type[f] != 3) {
        out[f] = enzyme[edge_to_enzyme[f]] * kcat[edge_to_enzyme[f]];
      }
    }
    return out;
  }

  vector get_edge_flux(vector conc, vector enzyme, vector dgr, vector kcat,
                       vector km, vector ki, vector tc, vector dc,
                       vector kcat_pme, vector conc_pme, vector drain,
                       real temperature, real drain_small_conc_corrector,
                       matrix S, vector subunits, array[] int edge_type,
                       array[] int edge_to_enzyme, array[] int edge_to_drain,
                       array[] int ci_mic_ix,
                       array[] int sub_km_ix_by_edge_long,
                       array[,] int sub_km_ix_by_edge_bounds,
                       array[] int prod_km_ix_by_edge_long,
                       array[,] int prod_km_ix_by_edge_bounds,
                       array[] int sub_by_edge_long,
                       array[,] int sub_by_edge_bounds,
                       array[] int prod_by_edge_long,
                       array[,] int prod_by_edge_bounds,
                       array[] int ci_ix_long, array[,] int ci_ix_bounds,
                       array[] int allostery_ix_long,
                       array[,] int allostery_ix_bounds,
                       array[] int allostery_type, array[] int allostery_mic,
                       array[] int edge_to_tc, array[] int phos_ix_long,
                       array[,] int phos_ix_bounds,
                       array[] int phosphorylation_type,
                       array[] int phosphorylation_pme) {
    int N_edge = cols(S);
    vector[N_edge] vmax = get_vmax_by_edge(enzyme, kcat, edge_to_enzyme,
                                           edge_type);
    vector[N_edge] reversibility = get_reversibility(dgr, temperature, S,
                                                     conc, edge_type);
    vector[N_edge] free_enzyme_ratio = get_free_enzyme_ratio(conc, S, km, ki,
                                                             edge_type,
                                                             ci_mic_ix,
                                                             sub_km_ix_by_edge_long,
                                                             sub_km_ix_by_edge_bounds,
                                                             prod_km_ix_by_edge_long,
                                                             prod_km_ix_by_edge_bounds,
                                                             sub_by_edge_long,
                                                             sub_by_edge_bounds,
                                                             prod_by_edge_long,
                                                             prod_by_edge_bounds,
                                                             ci_ix_long,
                                                             ci_ix_bounds);
    vector[N_edge] saturation = get_saturation(conc, km, free_enzyme_ratio,
                                               sub_km_ix_by_edge_long,
                                               sub_km_ix_by_edge_bounds,
                                               sub_by_edge_long,
                                               sub_by_edge_bounds, edge_type);
    vector[N_edge] allostery = get_allostery(conc, free_enzyme_ratio, tc, dc,
                                             subunits, allostery_ix_long,
                                             allostery_ix_bounds,
                                             allostery_type, allostery_mic,
                                             edge_to_tc);
    vector[N_edge] phosphorylation = get_phosphorylation(kcat_pme, conc_pme,
                                                         phos_ix_long,
                                                         phos_ix_bounds,
                                                         phosphorylation_type,
                                                         phosphorylation_pme,
                                                         subunits);
    vector[N_edge] drain_by_edge = get_drain_by_edge(drain, conc,
                                                     edge_to_drain,
                                                     sub_by_edge_long,
                                                     sub_by_edge_bounds,
                                                     edge_type,
                                                     drain_small_conc_corrector);
    return vmax .* saturation .* reversibility .* allostery
           .* phosphorylation .* drain_by_edge;
  }

  vector dbalanced_dt(real time, vector current_balanced, vector unbalanced,
                      array[] int balanced_ix, array[] int unbalanced_ix,
                      vector enzyme, vector dgr, vector kcat, vector km,
                      vector ki, vector tc, vector dc, vector kcat_pme,
                      vector conc_pme, vector drain, real temperature,
                      real drain_small_conc_corrector, matrix S,
                      vector subunits, array[] int edge_type,
                      array[] int edge_to_enzyme, array[] int edge_to_drain,
                      array[] int ci_mic_ix,
                      array[] int sub_km_ix_by_edge_long,
                      array[,] int sub_km_ix_by_edge_bounds,
                      array[] int prod_km_ix_by_edge_long,
                      array[,] int prod_km_ix_by_edge_bounds,
                      array[] int sub_by_edge_long,
                      array[,] int sub_by_edge_bounds,
                      array[] int prod_by_edge_long,
                      array[,] int prod_by_edge_bounds,
                      array[] int ci_ix_long, array[,] int ci_ix_bounds,
                      array[] int allostery_ix_long,
                      array[,] int allostery_ix_bounds,
                      array[] int allostery_type, array[] int allostery_mic,
                      array[] int edge_to_tc,
                      array[] int phosphorylation_ix_long,
                      array[,] int phosphorylation_ix_bounds,
                      array[] int phosphorylation_type,
                      array[] int phosphorylation_pme) {
    vector[rows(current_balanced) + rows(unbalanced)] current_concentration;
    current_concentration[balanced_ix] = current_balanced;
    current_concentration[unbalanced_ix] = unbalanced;
    vector[cols(S)] edge_flux = get_edge_flux(current_concentration, enzyme,
                                              dgr, kcat, km, ki, tc, dc,
                                              kcat_pme, conc_pme, drain,
                                              temperature,
                                              drain_small_conc_corrector, S,
                                              subunits, edge_type,
                                              edge_to_enzyme, edge_to_drain,
                                              ci_mic_ix,
                                              sub_km_ix_by_edge_long,
                                              sub_km_ix_by_edge_bounds,
                                              prod_km_ix_by_edge_long,
                                              prod_km_ix_by_edge_bounds,
                                              sub_by_edge_long,
                                              sub_by_edge_bounds,
                                              prod_by_edge_long,
                                              prod_by_edge_bounds,
                                              ci_ix_long, ci_ix_bounds,
                                              allostery_ix_long,
                                              allostery_ix_bounds,
                                              allostery_type, allostery_mic,
                                              edge_to_tc,
                                              phosphorylation_ix_long,
                                              phosphorylation_ix_bounds,
                                              phosphorylation_type,
                                              phosphorylation_pme);
    return (S * edge_flux)[balanced_ix];
  }
  complex_vector get_complex_edge_flux_enzyme(vector conc,
                                              complex_vector enzyme,
                                              vector dgr, vector kcat,
                                              vector km, vector ki,
                                              vector tc, vector dc,
                                              vector kcat_pme,
                                              vector conc_pme, vector drain,
                                              real temperature,
                                              real drain_small_conc_corrector,
                                              matrix S, vector subunits,
                                              array[] int edge_type,
                                              array[] int edge_to_enzyme,
                                              array[] int edge_to_drain,
                                              array[] int ci_mic_ix,
                                              array[] int sub_km_ix_by_edge_long,
                                              array[,] int sub_km_ix_by_edge_bounds,
                                              array[] int prod_km_ix_by_edge_long,
                                              array[,] int prod_km_ix_by_edge_bounds,
                                              array[] int sub_by_edge_long,
                                              array[,] int sub_by_edge_bounds,
                                              array[] int prod_by_edge_long,
                                              array[,] int prod_by_edge_bounds,
                                              array[] int ci_ix_long,
                                              array[,] int ci_ix_bounds,
                                              array[] int allostery_ix_long,
                                              array[,] int allostery_ix_bounds,
                                              array[] int allostery_type,
                                              array[] int allostery_mic,
                                              array[] int edge_to_tc,
                                              array[] int phos_ix_long,
                                              array[,] int phos_ix_bounds,
                                              array[] int phosphorylation_type,
                                              array[] int phosphorylation_pme) {
    int N_edge = cols(S);
    complex_vector[N_edge] vmax = get_complex_vmax_by_edge(enzyme, kcat,
                                                           edge_to_enzyme,
                                                           edge_type);
    vector[N_edge] reversibility = get_reversibility(dgr, temperature, S,
                                                     conc, edge_type);
    vector[N_edge] free_enzyme_ratio = get_free_enzyme_ratio(conc, S, km, ki,
                                                             edge_type,
                                                             ci_mic_ix,
                                                             sub_km_ix_by_edge_long,
                                                             sub_km_ix_by_edge_bounds,
                                                             prod_km_ix_by_edge_long,
                                                             prod_km_ix_by_edge_bounds,
                                                             sub_by_edge_long,
                                                             sub_by_edge_bounds,
                                                             prod_by_edge_long,
                                                             prod_by_edge_bounds,
                                                             ci_ix_long,
                                                             ci_ix_bounds);
    vector[N_edge] saturation = get_saturation(conc, km, free_enzyme_ratio,
                                               sub_km_ix_by_edge_long,
                                               sub_km_ix_by_edge_bounds,
                                               sub_by_edge_long,
                                               sub_by_edge_bounds, edge_type);
    vector[N_edge] allostery = get_allostery(conc, free_enzyme_ratio, tc, dc,
                                             subunits, allostery_ix_long,
                                             allostery_ix_bounds,
                                             allostery_type, allostery_mic,
                                             edge_to_tc);
    vector[N_edge] phosphorylation = get_phosphorylation(kcat_pme, conc_pme,
                                                         phos_ix_long,
                                                         phos_ix_bounds,
                                                         phosphorylation_type,
                                                         phosphorylation_pme,
                                                         subunits);
    vector[N_edge] drain_by_edge = get_drain_by_edge(drain, conc,
                                                     edge_to_drain,
                                                     sub_by_edge_long,
                                                     sub_by_edge_bounds,
                                                     edge_type,
                                                     drain_small_conc_corrector);
    return vmax .* saturation .* reversibility .* allostery
           .* phosphorylation .* drain_by_edge;
  }

  complex_vector get_complex_vmax_by_edge(complex_vector enzyme, vector kcat,
                                          array[] int edge_to_enzyme,
                                          array[] int edge_type) {
    int N_edge = size(edge_to_enzyme);
    complex_vector[N_edge] out = rep_vector(1, N_edge);
    for (f in 1 : N_edge) {
      if (edge_type[f] != 3) {
        out[f] = enzyme[edge_to_enzyme[f]] * kcat[edge_to_enzyme[f]];
      }
    }
    return out;
  }
  vector get_flux_jacobian(complex_vector edge_flux,
                           complex complex_enzyme_step) {
    return get_imag(edge_flux) ./ get_imag(complex_enzyme_step);
  }
  matrix get_concentration_control_matrix(matrix S, matrix elasticity,
                                          array[] int balanced_mic_ix,
                                          int N_edge, int N_balanced) {
    matrix[N_balanced, N_edge] balanced_S = S[balanced_mic_ix,  : ];
    return -generalized_inverse(balanced_S * elasticity) * balanced_S;
  }
  matrix get_flux_control_matrix(matrix S, matrix elasticity,
                                 array[] int balanced_mic_ix, int N_edge,
                                 int N_balanced) {
    matrix[N_balanced, N_edge] balanced_S = S[balanced_mic_ix,  : ];
    matrix[N_edge, N_edge] I = identity_matrix(N_edge);
    return I
           - elasticity * generalized_inverse(balanced_S * elasticity)
             * balanced_S;
  }
  complex_vector get_complex_edge_flux_metabolite(complex_vector conc,
                                                  vector enzyme, vector dgr,
                                                  vector kcat, vector km,
                                                  vector ki, vector tc,
                                                  vector dc, vector kcat_pme,
                                                  vector conc_pme,
                                                  vector drain,
                                                  real temperature,
                                                  real drain_small_conc_corrector,
                                                  matrix S, vector subunits,
                                                  array[] int edge_type,
                                                  array[] int edge_to_enzyme,
                                                  array[] int edge_to_drain,
                                                  array[] int ci_mic_ix,
                                                  array[] int sub_km_ix_by_edge_long,
                                                  array[,] int sub_km_ix_by_edge_bounds,
                                                  array[] int prod_km_ix_by_edge_long,
                                                  array[,] int prod_km_ix_by_edge_bounds,
                                                  array[] int sub_by_edge_long,
                                                  array[,] int sub_by_edge_bounds,
                                                  array[] int prod_by_edge_long,
                                                  array[,] int prod_by_edge_bounds,
                                                  array[] int ci_ix_long,
                                                  array[,] int ci_ix_bounds,
                                                  array[] int allostery_ix_long,
                                                  array[,] int allostery_ix_bounds,
                                                  array[] int allostery_type,
                                                  array[] int allostery_mic,
                                                  array[] int edge_to_tc,
                                                  array[] int phos_ix_long,
                                                  array[,] int phos_ix_bounds,
                                                  array[] int phosphorylation_type,
                                                  array[] int phosphorylation_pme) {
    int N_edge = cols(S);
    vector[N_edge] vmax = get_vmax_by_edge(enzyme, kcat, edge_to_enzyme,
                                           edge_type);
    complex_vector[N_edge] reversibility = get_reversibility_complex(dgr,
                                                                    temperature,
                                                                    S, conc,
                                                                    edge_type);
    complex_vector[N_edge] free_enzyme_ratio = get_free_enzyme_ratio_complex(conc,
                                                                    S, km,
                                                                    ki,
                                                                    edge_type,
                                                                    ci_mic_ix,
                                                                    sub_km_ix_by_edge_long,
                                                                    sub_km_ix_by_edge_bounds,
                                                                    prod_km_ix_by_edge_long,
                                                                    prod_km_ix_by_edge_bounds,
                                                                    sub_by_edge_long,
                                                                    sub_by_edge_bounds,
                                                                    prod_by_edge_long,
                                                                    prod_by_edge_bounds,
                                                                    ci_ix_long,
                                                                    ci_ix_bounds);
    complex_vector[N_edge] saturation = get_saturation_complex(conc, km,
                                                               free_enzyme_ratio,
                                                               sub_km_ix_by_edge_long,
                                                               sub_km_ix_by_edge_bounds,
                                                               sub_by_edge_long,
                                                               sub_by_edge_bounds,
                                                               edge_type);
    complex_vector[N_edge] allostery = get_allostery_complex(conc,
                                                             free_enzyme_ratio,
                                                             tc, dc,
                                                             subunits,
                                                             allostery_ix_long,
                                                             allostery_ix_bounds,
                                                             allostery_type,
                                                             allostery_mic,
                                                             edge_to_tc);
    complex_vector[N_edge] phosphorylation = get_phosphorylation(kcat_pme,
                                                                 conc_pme,
                                                                 phos_ix_long,
                                                                 phos_ix_bounds,
                                                                 phosphorylation_type,
                                                                 phosphorylation_pme,
                                                                 subunits);
    complex_vector[N_edge] drain_by_edge = get_drain_by_edge_complex(drain,
                                                                    conc,
                                                                    edge_to_drain,
                                                                    sub_by_edge_long,
                                                                    sub_by_edge_bounds,
                                                                    edge_type,
                                                                    drain_small_conc_corrector);
    return vmax .* saturation .* reversibility .* allostery
           .* phosphorylation .* drain_by_edge;
  }
  complex_vector get_allostery_complex(complex_vector conc,
                                       // one per mic
                                       complex_vector free_enzyme_ratio,
                                       // one per edge
                                       vector tc,
                                       // one per allosteric enzyme
                                       vector dc,
                                       // one per allostery
                                       vector subunits,
                                       // one per edge
                                       array[] int allostery_ix_long,
                                       // - long and bounds encode a ragged
                                       array[,] int allostery_ix_bounds,
                                       //   array with one entry per edge
                                       array[] int allostery_type,
                                       // one per allostery
                                       array[] int allostery_mic,
                                       // one per allostery
                                       array[] int edge_to_tc // one per edge
                                       ) {
    int N_edge = size(allostery_ix_bounds);
    complex_vector[N_edge] out = rep_vector(1, N_edge);
    for (f in 1 : N_edge) {
      int N_allostery = measure_ragged(allostery_ix_bounds, f);
      if (N_allostery == 0) {
        continue;
      }
      complex Q_num = 1;
      complex Q_denom = 1;
      real tc_edge = tc[edge_to_tc[f]];
      for (allostery in extract_ragged(allostery_ix_long,
                                       allostery_ix_bounds, f)) {
        complex conc_over_dc = conc[allostery_mic[allostery]] / dc[allostery];
        if (allostery_type[allostery] == 1) {
          // activation
          Q_denom += conc_over_dc;
        } else {
          // inhibition
          Q_num += conc_over_dc;
        }
      }
      out[f] = 1
               / (1
                  + tc_edge
                    * (free_enzyme_ratio[f] * Q_num / Q_denom) ^ subunits[f]);
    }
    return out;
  }

  complex_vector get_reversibility_complex(vector dgr, real temperature,
                                           matrix S, complex_vector conc,
                                           array[] int edge_type) {
    real RT = 0.008314 * temperature;
    int N_edge = cols(S);
    int N_mets = rows(S);
    complex_vector[N_mets] log_conc;
    complex_vector[N_edge] reaction_quotient;
    complex_vector[N_edge] out;
    for (met in 1 : N_mets) {
      log_conc[met] = log(conc[met]);
    }
    reaction_quotient = S' * log_conc;
    for (f in 1 : N_edge) {
      if (edge_type[f] == 1) {
        // reversible michaelis menten
        out[f] = 1 - exp((dgr[f] + RT * reaction_quotient[f]) / RT);
      } else {
        out[f] = 1;
      }
    }
    return out;
  }

  complex_vector get_free_enzyme_ratio_complex(complex_vector conc, matrix S,
                                               vector km, vector ki,
                                               array[] int edge_type,
                                               array[] int ci_mic_ix,
                                               array[] int sub_km_ix_by_edge_long,
                                               array[,] int sub_km_ix_by_edge_bounds,
                                               array[] int prod_km_ix_by_edge_long,
                                               array[,] int prod_km_ix_by_edge_bounds,
                                               array[] int sub_by_edge_long,
                                               array[,] int sub_by_edge_bounds,
                                               array[] int prod_by_edge_long,
                                               array[,] int prod_by_edge_bounds,
                                               array[] int ci_ix_long,
                                               array[,] int ci_ix_bounds) {
    /* Find the proportion of enzyme that is free, for each edge. */
    int N_edge = cols(S);
    complex_vector[N_edge] denom;
    for (f in 1 : N_edge) {
      if (edge_type[f] == 3) {
        // drain
        denom[f] = 1;
        continue;
      }
      int N_sub = measure_ragged(sub_by_edge_bounds, f);
      int N_prod = measure_ragged(prod_by_edge_bounds, f);
      int N_ci = measure_ragged(ci_ix_bounds, f);
      array[N_sub] int sub_ix = extract_ragged(sub_by_edge_long,
                                               sub_by_edge_bounds, f);
      array[N_sub] int sub_km_ix = extract_ragged(sub_km_ix_by_edge_long,
                                                  sub_km_ix_by_edge_bounds,
                                                  f);
      array[N_prod] int prod_ix = extract_ragged(prod_by_edge_long,
                                                 prod_by_edge_bounds, f);
      complex_vector[N_sub] sub_over_km = conc[sub_ix] ./ km[sub_km_ix];
      denom[f] = prod((rep_vector(1, N_sub) + sub_over_km)
                      ^ abs(S[sub_ix, f]));
      if (edge_type[f] == 1) {
        // reversible michaelis menten
        array[N_prod] int prod_km_ix = extract_ragged(prod_km_ix_by_edge_long,
                                                      prod_km_ix_by_edge_bounds,
                                                      f);
        complex_vector[N_prod] prod_over_km = conc[prod_ix] ./ km[prod_km_ix];
        denom[f] += prod((rep_vector(1, N_prod) + prod_over_km)
                         ^ abs(S[prod_ix, f]))
                    - 1;
      }
      if (N_ci > 0) {
        array[N_ci] int ci_ix = extract_ragged(ci_ix_long, ci_ix_bounds, f);
        denom[f] += sum(conc[ci_mic_ix[ci_ix]] ./ ki[ci_ix]);
      }
    }
    return 1 / denom;
  }

  complex_vector get_saturation_complex(complex_vector conc, vector km,
                                        complex_vector free_enzyme_ratio,
                                        array[] int sub_km_ix_by_edge_long,
                                        array[,] int sub_km_ix_by_edge_bounds,
                                        array[] int sub_by_edge_long,
                                        array[,] int sub_by_edge_bounds,
                                        array[] int edge_type) {
    int N_edge = size(sub_by_edge_bounds);
    complex_vector[N_edge] prod_conc_over_km;
    for (f in 1 : N_edge) {
      if (edge_type[f] == 3) {
        prod_conc_over_km[f] = 1;
        continue;
      }
      int N_sub = measure_ragged(sub_by_edge_bounds, f);
      array[N_sub] int sub_ix = extract_ragged(sub_by_edge_long,
                                               sub_by_edge_bounds, f);
      array[N_sub] int sub_km_ix = extract_ragged(sub_km_ix_by_edge_long,
                                                  sub_km_ix_by_edge_bounds,
                                                  f);
      prod_conc_over_km[f] = prod(conc[sub_ix] ./ km[sub_km_ix]);
    }
    return prod_conc_over_km .* free_enzyme_ratio;
  }

  complex_vector get_drain_by_edge_complex(vector drain, complex_vector conc,
                                           array[] int edge_to_drain,
                                           array[] int sub_by_edge_long,
                                           array[,] int sub_by_edge_bounds,
                                           array[] int edge_type,
                                           real drain_small_conc_corrector) {
    int N_edge = size(edge_type);
    complex_vector[N_edge] out = rep_vector(1, N_edge);
    for (f in 1 : N_edge) {
      if (edge_type[f] == 3) {
        int N_sub = measure_ragged(sub_by_edge_bounds, f);
        array[N_sub] int subs = extract_ragged(sub_by_edge_long,
                                               sub_by_edge_bounds, f);
        out[f] = drain[edge_to_drain[f]]
                 * prod(conc[subs]
                        ./ (conc[subs] + drain_small_conc_corrector));
      }
    }
    return out;
  }

vector maud_ae_system(vector conc,
                      data real rel_tol_ode, data real abs_tol_ode, data int max_num_steps_ode,
                      vector conc_unbalanced,
                      array[] int balanced_mic_ix, array[] int unbalanced_mic_ix,
                      vector conc_enzyme, vector dgr, vector kcat, vector km,
                      vector ki, vector tc, vector dc, vector kcat_pme,
                      vector conc_pme, vector drain, real temperature,
                      real drain_small_conc_corrector, matrix S,
                      vector subunits, array[] int edge_type,
                      array[] int edge_to_enzyme, array[] int edge_to_drain,
                      array[] int ci_mic_ix,
                      array[] int sub_km_ix_by_edge_long,
                      array[,] int sub_km_ix_by_edge_bounds,
                      array[] int prod_km_ix_by_edge_long,
                      array[,] int prod_km_ix_by_edge_bounds,
                      array[] int sub_by_edge_long,
                      array[,] int sub_by_edge_bounds,
                      array[] int prod_by_edge_long,
                      array[,] int prod_by_edge_bounds,
                      array[] int ci_ix_long, array[,] int ci_ix_bounds,
                      array[] int allostery_ix_long,
                      array[,] int allostery_ix_bounds,
                      array[] int allostery_type, array[] int allostery_mic,
                      array[] int edge_to_tc,
                      array[] int phosphorylation_ix_long,
                      array[,] int phosphorylation_ix_bounds,
                      array[] int phosphorylation_type,
                      array[] int phosphorylation_pme){
  return conc - to_vector(ode_bdf_tol(dbalanced_dt, conc, 0, {1},
                          rel_tol_ode, abs_tol_ode, max_num_steps_ode,
                          conc_unbalanced,
                          balanced_mic_ix,
                          unbalanced_mic_ix,
                          conc_enzyme,
                          dgr, kcat, km, ki,
                          tc,
                          dc, kcat_pme,
                          conc_pme,
                          drain,
                          temperature,
                          drain_small_conc_corrector, S,
                          subunits, edge_type,
                          edge_to_enzyme, edge_to_drain,
                          ci_mic_ix, sub_km_ix_by_edge_long,
                          sub_km_ix_by_edge_bounds,
                          prod_km_ix_by_edge_long,
                          prod_km_ix_by_edge_bounds,
                          sub_by_edge_long,
                          sub_by_edge_bounds,
                          prod_by_edge_long,
                          prod_by_edge_bounds, ci_ix_long,
                          ci_ix_bounds, allostery_ix_long,
                          allostery_ix_bounds,
                          allostery_type, allostery_mic,
                          edge_to_tc,
                          phosphorylation_ix_long,
                          phosphorylation_ix_bounds,
                          phosphorylation_type,
                          phosphorylation_pme)[1]);
}

}
