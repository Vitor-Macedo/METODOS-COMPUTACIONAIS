data{
  int j;
  vector[j] ni;
  vector[j-1]mi;
}
transformed data{
  real r;  
  r = sum(ni) - sum(mi);
}
parameters{
  vector[j] p;
  real N;
}
model{
  N ~ normal(50,50);
  p ~ uniform(0,1);
  target += tgamma(N+1) - tgamma(N - r + 1) + sum(ni .* p  + (N-ni).*(1-p)) + normal_lpdf(N|50,sqrt(50));
}
