vector get_keq(matrix S,
               vector formation_energy,
               int[] mic_to_met,
               vector water_stoichiometry){
  /*
     Calculate keqs from metabolite formation energies, assuming water's
     formation energy is known exactly.
  */
  real minus_RT = -0.008314 * 298.15;
  real dgf_water = -157.6;  // From http://equilibrator.weizmann.ac.il/metabolite?compoundId=C00001
  vector[rows(S)] delta_g =
    S' * formation_energy[mic_to_met]
    + water_stoichiometry * dgf_water;
  return exp(delta_g / minus_RT);
}
