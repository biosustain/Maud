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
      Calculate dgr standard from metabolite formation energies, assuming water's
      formation energy is known exactly.
  */
  real minus_RT = -0.008314 * 298.15;
  real dgf_water = -157.6;  // From http://equilibrator.weizmann.ac.il/metabolite?compoundId=C00001
  vector[cols(S)] dgrs = S' * dgf[mic_to_met] + water_stoichiometry * dgf_water;
  return dgrs;
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

int measure_ragged(int[] bounds, int i){
  return bounds[i, 2] - bounds[i, 1];
}

int extract_ragged(int[] long, int[] bounds, int i){
  /*
    Extract the ith element of a ragged array stored in 1d array long.

    Make sure that members bounds[i, 1] to bounds[i, 2] of long form the
    required element!

   */
  return long[bounds[i, 1]:bounds[i, 2]];
}

vector get_saturation(vector conc,
                      vector km,
                      vector free_enzyme_ratio,
                      int[] km_lookup,
                      int[] sub_by_flux_long,
                      int[] sub_by_flux_bounds
                      int[] flux_type){
  vector[N_flux] prod_conc_over_km;
  for (f in 1:N_flux){
    if (flux_type[f] == 2){
      prod_conc_over_km[f] = 1;
      continue;
    }
    int N_sub = measure_ragged(sub_by_flux_bounds, f);
    int[N_sub] sub_ix = extract_ragged(sub_by_flux_long, sub_by_flux_bounds, f);
    prod_conc_over_km[f] = prod(conc[sub_ix] ./ km[km_lookup[sub_ix, f]]);
  }
  return prod_conc_over_km * free_enzyme_ratio;
}

vector get_free_enzyme_ratio(vector conc,
                             matrix S,
                             vector km,
                             vector ki,
                             int[] flux_type;
                             int[] km_lookup,
                             int[] ki_lookup,
                             int[] sub_by_flux_long,
                             int[] sub_by_flux_bounds,
                             int[] prod_by_flux_long,
                             int[] prod_by_flux_bounds,
                             int[] ci_by_flux_long,
                             int[] ci_by_flux_bounds){
  /* Find the proportion of enzyme that is free, for each flux. */
  int N_km = rows(km);
  int N_flux = cols(S);
  vector[N_flux] out;
  for (f in 1:N_flux){
    if (flux_type[f] == 2){
      out[f] = 1;
      continue;
    }
    int N_sub = measure_ragged(sub_by_flux_bounds, f);
    int N_prod = measure_ragged(prod_by_flux_bounds, f);
    int N_ci = measure_ragged(ci_by_flux_bounds, f);
    int[N_sub] sub_ix = extract_ragged(sub_by_flux_long, sub_by_flux_bounds, f);
    int[N_prod] prod_ix = extract_ragged(prod_by_flux_long, prod_by_flux_bounds, f);
    real sub_over_km = conc[sub_ix] ./ km[km_lookup[sub_ix, f]];
    real prod_over_km = conc[prod_ix] ./ km[km_lookup[prod_ix, f]];
    out[f] = prod((rep_vector(1, N_sub) + sub_over_km) ^ fabs(S[sub_ix, f]));
    if (flux_type[f] == 1){
      out[f] += prod((rep_vector(1, N_prod) + prod_over_km) ^ fabs(S[prod_ix, f])) - 1;
    }
    if (N_ci > 0){
      int[N_ci] ci_ix = extract_ragged(ci_by_flux_long, ci_by_flux_bounds, f);
      out[f] += sum(conc[ci_ix] / ki[ki_lookup[ci_ix, f]]);
    }
  }
  return out;
}

vector get_reversibility(vector dgr, matrix S, vector conc, int[] flux_type){
  real RT = 0.008314 * 298.15;
  int N_flux = cols(S);
  vector[N_flux] reaction_quotient = S' * conc;
  vector[N_flux] out;
  for (f in 1:N_flux){
    if (flux_type[f] == 1)
      out[j] = 1 - exp(dgr[f] + RT * reaction_quotient[f]);
    else
      out[j] = 1;
  }
  return out;
}

vector get_allostery(vector conc,
                     vector free_enzyme_ratio,
                     vector tc,
                     vector dt,
                     vector dr,
                     vector subunits,
                     int[] dt_lookup,
                     int[] dr_lookup,
                     int[] flux_to_tc,
                     int[] ai_ix_long,
                     int[] ai_ix_bounds,
                     int[] aa_ix_long,
                     int[] aa_ix_bounds
){
  int N_flux = size(aa_ix_bounds);
  vector[N_flux] out = rep_vector(1, N_flux);
  for (f in 1:N_flux){
    int N_ai = measure_ragged(ai_bounds, f);
    int N_aa = measure_ragged(aa_bounds, f);
    real Q_num = 1;
    real Q_denom = 1;
    if (N_ai > 0) Q_num = 1 + sum(conc[ais] ./ dt[dt_lookup[ais, f]]);
    if (N_aa > 0) Q_denom = 1 + sum(conc[aas] ./ dr[dr_lookup[aas, f]]);
    if ((N_ai > 0) || (N_aa > 0)){
      real tc_f = flux_to_tc[f];
      out[j] = inv(1 + tc_f * free_enzyme_ratio[f] * Q_num / Q_denom) ^ subunits[f];
    }
  }
  return out;
}

vector get_phosphorylation(vector kcat_phos,
                           vector conc_phos,
                           int[] pa_ix_long,
                           int[] pa_ix_bounds,
                           int[] pi_ix_long,
                           int[] pi_ix_bounds,
                           vector subunits){
  int N_flux = size(pa_ix_bounds);
  vector[N_flux] out = rep_vector(1, N_flux);
  for (f in 1:N_flux){
    real alpha = 0;
    real beta = 0;
    int N_pa = measure_ragged(pa_ix_bounds, f);
    int N_pi = measure_ragged(pi_ix_bounds, f);
    if (N_pa > 0){
      int[N_pa] pas = extract_ragged(pa_ix_long, pa_ix_bounds, f);
      beta = sum(kcat_phos[pas] .* conc_phos[pas]);
    }
    if (N_pi > 0){
      int[N_pi] pis = extract_ragged(pi_ix_long, pi_ix_bounds, f);
      alpha = sum(kcat_phos[pis] .* conc_phos[pis]);
    }
    if ((N_pi > 0) || (N_pa > 0)){
      out[j] = inv(1 + (alpha / beta) ^ subunits[f]);  // TODO: what if beta is zero and alpha is non-zero?
    }
  }
  return out;
}

vector get_drain_flux(vector drain,
                      vector conc,
                      int[] flux_to_drain,
                      int[] sub_by_flux_long,
                      int[] sub_by_flux_bounds,
                      int[] flux_type){
  int N_flux = size(flux_type);
  vector[N_flux] out = rep_vector(1, N_flux);
  for (f in 1:N_flux){
    if (flux_type[f] == 2){
      int N_sub = measure_ragged(sub_by_flux_bounds, f);
      int subs[N_sub] = extract_ragged(sub_by_flux_long, sub_by_flux_bounds, f);
      out[f] = drain[flux_to_drain[j]] * prod(conc[subs] ./ (conc[subs] + 1e-6));
    }
  }
  return out;
}

vector get_flux(vector conc,
                vector enzyme,
                vector dgr,
                vector kcat,
                vector km,
                vector tc,
                vector dt,
                vector dr,
                vector kcat_phos,
                vector conc_phos,
                vector drain,
                matrix S,
                vector subunits,
                int[] flux_type,
                int[] flux_to_enzyme,
                int[] flux_to_tc,
                int[] flux_to_drain,
                int[] km_lookup,
                int[] ki_lookup,
                int[] dt_lookup,
                int[] dr_lookup,
                int[] sub_by_flux_long,
                int[] sub_by_flux_bounds,
                int[] prod_by_flux_long,
                int[] prod_by_flux_bounds,
                int[] ci_by_flux_long,
                int[] ci_by_flux_bounds,
                int[] ai_ix_long,
                int[] ai_ix_bounds,
                int[] aa_ix_long,
                int[] aa_ix_bounds);
                int[] pa_ix_long,
                int[] pa_ix_bounds,
                int[] pi_ix_long,
                int[] pi_ix_bounds){
  int N_flux = cols(S);
  vector[N_flux] vmax = get_vmax_by_flux(enzyme, kcat, flux_to_enzyme, flux_type);
  vector[N_flux] reversibility = get_reversibility(dgr, S, conc, flux_type);
  vector[N_flux] free_enzyme_ratio = get_free_enzyme_ratio(conc,
                                                           S,
                                                           km,
                                                           ki,
                                                           flux_type,
                                                           km_lookup,
                                                           ki_lookup,
                                                           sub_by_flux_long,
                                                           sub_by_flux_bounds,
                                                           prod_by_flux_long,
                                                           prod_by_flux_bounds,
                                                           ci_by_flux_long,
                                                           ci_by_flux_bounds);
  vector[N_flux] saturation = get_saturation(conc,
                                             km,
                                             free_enzyme_ratio,
                                             km_lookup,
                                             sub_by_flux_long,
                                             sub_by_flux_bounds,
                                             flux_type);
  vector[N_flux] allostery = get_allostery(conc,
                                           free_enzyme_ratio,
                                           tc,
                                           dt,
                                           dr,
                                           subunits,
                                           dt_lookup,
                                           dr_lookup,
                                           flux_to_tc,
                                           ai_ix_long,
                                           ai_ix_bounds,
                                           aa_ix_long,
                                           aa_ix_bounds);
  vector[N_flux] phosphorylation = get_phosphorylation(kcat_phos,
                                                       conc_phos,
                                                       pa_ix_long,
                                                       pa_ix_bounds,
                                                       pi_ix_long,
                                                       pi_ix_bounds,
                                                       subunits);
  vector[N_flux] drain_flux = get_drain_flux(drain,
                                             conc,
                                             flux_to_drain,
                                             sub_by_flux_long,
                                             sub_by_flux_bounds,
                                             flux_type);
  return vmax .* saturation .* reversibility .* allostery .* phosphorylation .* drain_flux;
}

vector dbalanced_dt(real time,
                    vector current_balanced,
                    vector unbalanced,
                    int[] balanced_ix,
                    int[] unbalanced_ix,
                    vector enzyme,
                    vector dgr,
                    vector kcat,
                    vector km,
                    vector tc,
                    vector dt,
                    vector dr,
                    vector kcat_phos,
                    vector conc_phos,
                    vector drain,
                    matrix S,
                    vector subunits,
                    int[] flux_type,
                    int[] flux_to_enzyme,
                    int[] flux_to_tc,
                    int[] flux_to_drain,
                    int[] km_lookup,
                    int[] ki_lookup,
                    int[] dt_lookup,
                    int[] dr_lookup,
                    int[] sub_by_flux_long,
                    int[] sub_by_flux_bounds,
                    int[] prod_by_flux_long,
                    int[] prod_by_flux_bounds,
                    int[] ci_by_flux_long,
                    int[] ci_by_flux_bounds,
                    int[] ai_ix_long,
                    int[] ai_ix_bounds,
                    int[] aa_ix_long,
                    int[] aa_ix_bounds);
                    int[] pa_ix_long,
                    int[] pa_ix_bounds,
                    int[] pi_ix_long,
                    int[] pi_ix_bounds
  vector[rows(current_balanced)+rows(unbalanced)] current_concentration;
  current_concentration[balanced_ix] = current_balanced;
  current_concentration[unbalanced_ix] = unbalanced;
  vector[rows(S)] flux = get_flux(current_concentration,
                                  enzyme,
                                  dgr,
                                  kcat,
                                  km,
                                  tc,
                                  dt,
                                  dr,
                                  kcat_phos,
                                  conc_phos,
                                  drain,
                                  S,
                                  subunits,
                                  flux_type,
                                  flux_to_enzyme,
                                  flux_to_tc,
                                  flux_to_drain,
                                  km_lookup,
                                  ki_lookup,
                                  dt_lookup,
                                  dr_lookup,
                                  sub_by_flux_long,
                                  sub_by_flux_bounds,
                                  prod_by_flux_long,
                                  prod_by_flux_bounds,
                                  ci_by_flux_long,
                                  ci_by_flux_bounds,
                                  ai_ix_long,
                                  ai_ix_bounds,
                                  aa_ix_long,
                                  aa_ix_bounds);
                                  pa_ix_long,
                                  pa_ix_bounds,
                                  pi_ix_long,
                                  pi_ix_bounds)
  return (S * flux)[balanced_ix];
}
