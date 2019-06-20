real[] get_fluxes(real[] metabolites,
                  real[] params,
                  real[] known_reals){
  // known reals
  real concentration_FBA = known_reals[1];
  real concentration_TDH = known_reals[2];
  real concentration_TPI = known_reals[3];
  real temperature = known_reals[4];
  real gas_constant = known_reals[5];
  // metabolites
  real M1 = metabolites[1];
  real M2 = metabolites[2];
  real M3 = metabolites[3];
  real M4 = metabolites[4];
  // thermodynamic parameters
  real delta_g_FBA = params[3];
  real delta_g_TDH = params[4];
  real delta_g_TPI = params[5];
  // kinetic parameters
  real influx_M1_alpha = params[7];
  real FBA_Kcat1 = params[8];
  real FBA_Kcat2 = params[9];
  real FBA_Ka = params[10];
  real FBA_Kp = params[11];
  real FBA_Kq = params[12];
  real FBA_Kia = params[13];
  real TDH_Kcat1 = params[14];
  real TDH_Kcat2 = params[15];
  real TDH_Ka = params[16];
  real TPI_Kcat1 = params[17];
  real TPI_Kcat2 = params[18];
  real TPI_Ka = params[19];
  real influx_M2_alpha = params[20];
  real influx_M4_alpha = params[21];
  // haldane_relationships
  real FBA_Keq = get_Keq(delta_g_FBA, temperature, gas_constant);
  real TDH_Keq = get_Keq(delta_g_TDH, temperature, gas_constant);
  real TPI_Keq = get_Keq(delta_g_TPI, temperature, gas_constant);
  real FBA_Kip = get_Kip_ordered_unibi(FBA_Keq, FBA_Kia, FBA_Kq, FBA_Kcat1, FBA_Kcat2);
  real FBA_Kiq = get_Kiq_ordered_unibi(FBA_Keq, FBA_Ka, FBA_Kp, FBA_Kcat1, FBA_Kcat2);
  // fluxes
  real influx_M1 = irr_mass_action(M1, influx_M1_alpha);
  real influx_M2 = irr_ass_action(M2, influx_M2_alpha);
  real FBA = ordered_unibi(M1, M2, M3,
                           concentration_FBA * FBA_Kcat1, concentration_FBA * FBA_Kcat2,
                           FBA_Ka, FBA_Kp, FBA_Kq,
                           FBA_Kia, FBA_Kip, FBA_Kiq,
                           FBA_Keq);
  real TDH = uniuni(M3, M4, concentration_TDH * TDH_Kcat1, concentration_TDH * TDH_Kcat2, TDH_Ka, TDH_Keq);
  real TPI = uniuni(DHAP, GAP, concentration_TPI * TPI_Kcat1, concentration_TPI * TPI_Kcat2, TPI_Ka, TPI_Keq);
  real influx_M4 = irr_mass_action(M4, influx_M4_alpha);
  return {influx_M1, influx_M2, FBA, TDH, TPI, influx_M4};
}
  
real[] get_odes(real[] fluxes){
  real influx_M1;
  real influx_M2;
  real FBA;
  real TDH;
  real TPI;
  real influx_M4;
  return {
    trasport_M1 - flux_FBA,
    influx_M2 + flux_FBA - flux_TPI,
    flux_FBA + flux_TPI - flux_TDH,
    influx_M4 + flux_TDH
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
