vector steady_state_function(vector balanced, vector theta, data real[] xr, data int[] xi){
  int N_balanced = {{N_balanced}};
  real initial_time = 0;
  real time_step = {{time_step}};
  real balanced_new[{{N_balanced}}];
  balanced_new = integrate_ode_bdf(
    ode_func,
    to_array_1d(balanced),
    initial_time,
    rep_array(time_step, 1),
    to_array_1d(theta),
    xr,
    rep_array(0, 1),
    1e-8, 1e-12, 1e5
  )[1, ]; 
  return to_vector(balanced_new) - balanced;
}
