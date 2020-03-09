functions{
#include big_k_rate_equations.stan
#include haldane_relationships.stan
#include allostery.stan
  {{fluxes_function}}
  {{ode_function}}
}
