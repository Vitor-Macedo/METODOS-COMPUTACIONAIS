data{
  int<lower=0> M;
  int<lower=0> C;
  int<lower=0,upper=min(M,C)> R;
}
parameters{
  real<lower=(C - R + M)> N;
}
model{
  N ~ exponential(.02);
  R ~ binomial(C, M / N);
}
