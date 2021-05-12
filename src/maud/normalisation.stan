vector unz_1d(vector[] priors, vector z){
  return priors[1] + priors[2] .* z;
}
vector unz_log_1d(vector[] priors, vector z){
  return exp(log(priors[1]) + priors[2] .* z);
}
vector[] unz_2d(vector[,] priors, vector[] z){
  array[size(z)] vector [rows(z[1])] out;
  for (ex in 1:size(z)){
    out[ex] = unz_1d(priors[:,ex], z[ex]);
  }
  return out;
}
vector[] unz_log_2d(vector[,] priors, vector[] z){
  array[size(z)] vector [rows(z[1])] out;
  for (ex in 1:size(z)){
    out[ex] = unz_log_1d(priors[:,ex], z[ex]);
  }
  return out;
}
