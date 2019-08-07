
real[] get_odes(real[] fluxes){
  return {
      {% for ode in ode_stoic %}
        {{ode}}{{"," if not loop.last}}
      {%- endfor %}
  };
}
