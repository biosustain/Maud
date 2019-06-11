real[] get_fluxes(real[] ode_metabolites,
                  real[] params,
                  real[] known_reals){
  // unpacking known values
  real S1 = known_reals[1];
  real S2 = known_reals[2];
  real C1 = known_reals[3];
  real C2 = known_reals[4];
  real temperature = known_reals[5];
  real gas_constant = known_reals [6];

  // defining metabolites
  real A = ode_metabolites[1];
  real B = ode_metabolites[2];
  real C = ode_metabolites[3];
  real D = ode_metabolites[4];
  real E = ode_metabolites[5];

  // thermodynamic parameters
  real RXN1_delta_g = params[1];
  real RXN2_delta_g = params[2];
  real RXN3_delta_g = params[3];
  real RXN4_delta_g = params[4];

  // kinetic parameters
  // RXN1
  real RXN1_Kcat1 = params[5];
  real RXN1_Kcat2 = params[6];
  real RXN1_Ka = params[7];

  // RXN2
  real RXN2_Kcat1 = params[8];
  real RXN2_Kcat2 = params[9];
  real RXN2_Ka = params[10];
  real RXN2_Kp = params[11];
  real RXN2_Kq = params[12];
  real RXN2_Kia = params[13];

  // RXN3
  real RXN3_Kcat1 = params[14];
  real RXN3_Kcat2 = params[15];
  real RXN3_Ka = params[16];

  // RXN4
  real RXN4_Kcat1 = params[17];
  real RXN4_Kcat2 = params[18];
  real RXN4_Ka = params[19];
  real RXN4_Kb = params[20];
  real RXN4_Kp = params[21];
  real RXN4_Kq = params[22];
  real RXN4_Kib = params[23];
  real RXN4_Kiq = params[24];

  // Influx
  real Influx_Kcat1 = params[25];
  real Influx_Kcat2 = params[26];
  real Influx_Ka = params[27];
  real Influx_Kp = params[28];

  // Outflux
  real Outflux_Kcat1 = params[29];
  real Outflux_Kcat2 = params[30];
  real Outflux_Ka = params[31];
  real Outflux_Kp = params[32];

  // Enzyme Concentrations
  real Influx_E0 = params[33];
  real RXN1_E0 = params[34];
  real RXN2_E0 = params[35];
  real RXN3_E0 = params[36];
  real RXN4_E0 = params[37];
  real Outflux_E0 = params[38];

  // equilibrium constants
  real Influx_Keq = (Influx_Kcat1/Influx_Kcat2)*(Influx_Kp/Influx_Ka);
  real RXN1_Keq = get_Keq(RXN1_delta_g, temperature, gas_constant);
  real RXN2_Keq = get_Keq(RXN2_delta_g, temperature, gas_constant);
  real RXN3_Keq = get_Keq(RXN3_delta_g, temperature, gas_constant);
  real RXN4_Keq = get_Keq(RXN4_delta_g, temperature, gas_constant);
  real Outflux_Keq = (Outflux_Kcat1/Outflux_Kcat2)*(Outflux_Kp/Outflux_Ka);

  // haldane_relationships
  real RXN2_Kip = get_Kip_ordered_unibi(RXN2_Keq, RXN2_Kia, RXN2_Kq, RXN2_Kcat1, RXN2_Kcat2);
  real RXN2_Kiq = get_Kiq_ordered_unibi(RXN2_Keq, RXN2_Ka, RXN2_Kp, RXN2_Kcat1, RXN2_Kcat2);

  real RXN4_Kip = get_Kip_ordered_bibi(RXN4_Keq, RXN4_Ka, RXN4_Kq, RXN4_Kib, RXN4_Kcat1, RXN4_Kcat2);
  real RXN4_Kia = get_Kip_ordered_bibi(RXN4_Keq, RXN4_Kb, RXN4_Kp, RXN4_Kiq, RXN4_Kcat1, RXN4_Kcat2);

  // fluxes
  real flux_Influx = uniuni(C1, A, Influx_E0 * Influx_Kcat1, Influx_E0 * Influx_Kcat2, Influx_Ka, Influx_Keq);
  real flux_RXN1 = uniuni(A, B, RXN1_E0 * RXN1_Kcat1, RXN1_E0 * RXN1_Kcat2, RXN1_Ka, RXN1_Keq);
  real flux_RXN2 = ordered_unibi(B, C, D, RXN2_E0 * RXN2_Kcat1, RXN2_E0 * RXN2_Kcat2, RXN2_Ka, RXN2_Kp, RXN2_Kq, RXN2_Kia, RXN2_Kip, RXN2_Kiq, RXN2_Keq);
  real flux_RXN3 = uniuni(C, D, RXN3_E0 * RXN3_Kcat1, RXN3_E0 * RXN3_Kcat2, RXN3_Ka, RXN3_Keq);
  real flux_RXN4 = ordered_bibi(D, S1, S2, E, RXN4_E0 * RXN4_Kcat1, RXN4_E0 * RXN4_Kcat2, RXN4_Ka, RXN4_Kb, RXN4_Kp, RXN4_Kq, RXN4_Kia, RXN4_Kib, RXN4_Kip, RXN4_Kiq, RXN4_Keq);
  real flux_Outflux = uniuni(E, C2, Outflux_E0 * Outflux_Kcat1, Outflux_E0 * Outflux_Kcat2, Outflux_Ka, Outflux_Keq);

  return {flux_Influx, flux_RXN1, flux_RXN2, flux_RXN3, flux_RXN4, flux_Outflux};
}

real[] get_odes(real[] fluxes){
  real influx = fluxes[1];
  real RXN1 = fluxes[2];
  real RXN2 = fluxes[3];
  real RXN3 = fluxes[4];
  real RXN4 = fluxes[5];
  real outflux = fluxes[6];

  return {
    influx - RXN1,   // A
    RXN1 - RXN2,   // B
    RXN2 - RXN3,   // C
    RXN2 + RXN3 - RXN4,   // D
    RXN4 - outflux   // E
  };
}

real[] steady_state_equation(real t,
                             real[] ode_metabolites,
                             real[] params,
                             real[] known_reals,
                             int[] known_ints){
   for (m in 1:size(ode_metabolites)){
     if (ode_metabolites[m] < 0){
       reject("Metabolite ", m, " is ", ode_metabolites[m], " but should be greater than zero");
     }
   }
  return get_odes(get_fluxes(ode_metabolites, params, known_reals));
}
