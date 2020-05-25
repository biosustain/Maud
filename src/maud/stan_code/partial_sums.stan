real partial_sum_conc(real[] y,
                      int start,
                      int end,
                      vector[] conc,
                      int[] experiment,
                      int[] mic,
                      vector sigma){
  real out = 0;
  for (c in start:end){
    out += lognormal_lpdf(y[c] | log(conc[experiment[c], mic[c]]), sigma[c]);
  }
  return out;
}
real partial_sum_enz(real[] y,
                     int start,
                     int end,
                     matrix enz_conc,
                     int[] experiment,
                     int[] enzyme,
                     vector sigma){
  real out = 0;
  for (c in start:end){
    out += lognormal_lpdf(y[c] | log(enz_conc[experiment[c], enzyme[c]]), sigma[c]);
  }
  return out;
}
real partial_sum_flux(real[] y,
                      int start,
                      int end,
                      matrix flux,
                      int[] experiment,
                      int[] reaction,
                      vector sigma){
  real out = 0;
  for (c in start:end){
    out += normal_lpdf(y[c] | flux[experiment[c], reaction[c]], sigma[c]);
  }
  return out;
}
