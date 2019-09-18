vector get_fluxes(real[] m, real[] p, real[] xr){
  {%- for haldane in haldanes %}
  {{haldane}}
  {%- endfor %}
  {%- for fe in free_enzyme_ratio %} 
  {{fe}}
  {%- endfor %}
  return [
    {%- for flux in fluxes %} 
    {{flux}}{{"," if not loop.last-}} 
    {% endfor %}
  ]';
}
