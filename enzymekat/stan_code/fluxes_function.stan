vector get_fluxes(real[] m, real[] p){
  real empty_array[0];
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