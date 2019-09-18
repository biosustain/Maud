vector steady_state_function(vector balanced, vector theta, real[] xr, int[] xi){
  int N_unbalanced = {{N_unbalanced}};
  int N_balanced = {{N_balanced}};
  real initial_time = 0;
  real time_step = {{time_step}};
  real conc[{{N_balanced + N_unbalanced}}];
  real balanced_new[{{N_balanced}}];
  conc[{ {{-balanced_codes|join(',')-}} }] = to_array_1d(balanced);
  conc[{ {{-unbalanced_codes|join(',')-}} }] = to_array_1d(theta[1:{{N_unbalanced}}]);
  balanced_new = integrate_ode_bdf(
    ode_func,
    conc,
    initial_time,
    rep_array(time_step, 1),
    to_array_1d(theta[N_unbalanced+1:]),
    xr,
    rep_array(0, 1),
    1e-5, 1e-5, 1e3
  )[1, { {{-balanced_codes|join(',')-}} }]; 
  return to_vector(balanced_new) - balanced;
}
