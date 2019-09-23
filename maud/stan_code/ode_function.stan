real[] ode_func(real t, real[] m, real[] p, real[] xr, int[] xi){
  vector[{{N_flux}}] fluxes = get_fluxes(m, p);
  return {
    {%- for ode in ode_stoic %}
    {{ode}}{{"," if not loop.last}}
    {%- endfor %}
  };
}
