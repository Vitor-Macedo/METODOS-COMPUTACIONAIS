library(tidyverse)
dados <- data.frame(dias=c(0,8,28,41,63,79,97,117,135,154),
                    besouros=c(2,47,192,256,768,896,1120,896,1184,1024))
ft <- function(param,data,N0){
  t <- data$dias
  k <- param[1]
  r <- param[2]
  output <- (k*N0) / (N0 + (k-N0)*exp(-r*t)) 
  return(output)
}
# Função Objetivo
Q <- function(param,data,N0){
  ti <- data[,1]
  Ni <- data[,2]
  frk <- ft(param,data,N0)
  output <- sum((Ni-frk)^2)
  return(output)
}
# Derivadas parciais de 1° e 2° ordem da função f
dfdr <- function(param,data,N0){
  t <- data$dias
  K <- param[1]
  r <- param[2]
  output <- (K*N0*( (K-N0)*exp(-r*t)*t ))/ ((N0 + (K-N0)*exp(-r*t))^2)
  return(output)
}

dfdK <- function(param,data,N0){
  t <- data$dias
  K <- param[1]
  r <- param[2]
  output <- (N0*(-exp(-r*t)*N0+N0)) / (N0+exp(-r*t)*(K-N0))^2
  return(output)
}

d2fdr2 <- function(param,data,N0){
  t <- data$dias
  K <- param[1]
  r <- param[2]
  output <- (K*exp(-2*r*t)*N0*(t^2)*(-exp(r*t)*N0-N0+K)*(K-N0))/  ((N0+exp(-r*t)*(K-N0))^3)
  return(output)
}

d2fdK2 <- function(param,data,N0){
  t <- data$dias
  K <- param[1]
  r <- param[2]
  output <- (2*exp(-r*t)*N0*(-exp(-r*t)*N0+N0))/((N0+exp(-r*t)*(K-N0))^3)
  return(output)
}

d2fdKdr <- function(param,data,N0){
  t <- data$dias
  K <- param[1]
  r <- param[2]
  output <- (exp(-r*t)*N0*N0*t*(exp(-r*t)*N0-N0+2*K-K*exp(-r*t)))/((N0+exp(-r*t)*(K-N0))^3)
  return(output)
}

# Veto gradriente da função objetivo
dQdr <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dQdr=-2*besouros*dfdr(param = param,data=data,N0 = N0) +
             2*aux*dfdr(param = param,data=data,N0 = N0))
  output <- sum(df$dQdr)
  return(output)
}

dQdK <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dQdk=-2*besouros*dfdK(param = param,data=data,N0 = N0) +
             2*aux*dfdK(param = param,data=data,N0 = N0))
  output <- sum(df$dQdk)
  return(output)
}

gradQ <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  gradR <- dQdr(param,data,N0)
  gradK <- dQdK(param,data,N0)
  gradiente <- c(gradK,gradR)
  return(gradiente)
}

# Matriz Hessiana da função objetivo
d2Qdr2 <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           d2Qdr2=-2*besouros*d2fdr2(param,data,N0)+2*dfdr(param,data,N0)^2 +
             2*ft(param,data,N0)*d2fdr2(param,data,N0))
  output <- sum(df$d2Qdr2)
  return(output)
}

d2Qdk2 <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           d2QdK2=-2*besouros*d2fdK2(param,data,N0)+2*dfdK(param,data,N0)^2 +
             2*ft(param,data,N0)*d2fdK2(param,data,N0))
  output <- sum(df$d2QdK2)
  return(output)
}

d2QdKdr <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           d2QdKdr=-2*besouros*d2fdKdr(param,data,N0)+
             2*dfdK(param,data,N0)*dfdr(param,data,N0)+
             2*ft(param,data,N0)*d2fdKdr(param,data,N0))
  output <- sum(df$d2QdKdr)
  return(output)
}

hessianQ <- function(param,data,N0){
  A11 <- d2Qdk2(param,data,N0)
  A22 <- d2Qdr2(param,data,N0)
  A12 <- A21 <- d2QdKdr(param,data,N0)
  hessiana <- matrix(data = c(A11,A21,A12,A22),2,2)
  return(hessiana)
}

# Algorítimo de Newton Raphson
NR <- function(theta_k,data,N0,epsilon,criterium){
  iter=0
  while (criterium>epsilon) {
    iter <- iter + 1
    theta_k_plus_1 <- theta_k - solve(hessianQ(param=theta_k,data=data,N0=2))%*%gradQ(param=theta_k,data=data,N0=2)
    criterium <- abs(Q(param = theta_k_plus_1,data = data,N0=2)/Q(param = theta_k,data = data,N0=2) - 1)
    theta_k <- theta_k_plus_1
  }
  output <- list(iter=iter,theta=theta_k,Q_func=Q(param=theta_k,data,N0))
  return(output)
}

NR(theta_k = c(900,.1),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
NR(theta_k = c(500,.1),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
NR(theta_k = c(1024,.15),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
NR(theta_k = c(2000,.15),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
NR(theta_k = c(1024,.5),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
# Algorítimo de line search
ft_ls <- function(gama,param,data,N0){
  t <- data$dias
  k <- param[1]
  r <- param[2]
  p <- solve(hessianQ(param=param,data=data,N0=2))%*%gradQ(param=param,data=data,N0=2)
  pk <- p[1];pr=p[2]
  output <- ((k-gama*pk)*N0) / (N0 + (k-gama*pk-N0)*exp(-(r-gama*pr)*t)) 
  return(output)
}
Q_ls <- function(gama,param,data,N0){
  ti <- data$dias
  Ni <- data$besouros
  frk <- ft_ls(gama,param,data,N0)
  output <- sum((Ni-frk)^2)
  if(gama>0.5){
    return(output)
  }else{
    return(Inf)
  }
}

ls <- function(theta_k,data,N0,epsilon,criterium){
  iter=0
  while (criterium>epsilon) {
    iter <- iter + 1
    gama <- optimize(f = Q_ls,interval = c(0,10),param=theta_k,data=data,N0=N0)$minimum
    theta_k_plus_1 <- theta_k - gama*solve(hessianQ(param=theta_k,data=data,N0=2))%*%gradQ(param=theta_k,data=data,N0=2)
    criterium <- abs(Q(param = theta_k_plus_1,data = data,N0=2)/Q(param = theta_k,data = data,N0=2) - 1)
    theta_k <- theta_k_plus_1
  }
  output <- list(iter=iter,theta=theta_k,Q_func=Q(param=theta_k,data,N0))
  return(output)
}
ls(theta_k = c(900,.1),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
ls(theta_k = c(500,.1),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
ls(theta_k = c(1024,.15),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
ls(theta_k = c(2000,.15),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
ls(theta_k = c(1024,1),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)
# Escore de fisher
llk <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  ti <- data[,1]
  Ni <- data[,2]
  frk <- ft(param,data,N0)
  output <- sum(-.5*log(2*pi*sigma2) -(1/(2*sigma2))*(log(Ni)-log(frk))^2 )
  if(sigma2>0 & r>0 & K>0){
    return(-output)
  }else{
    return(Inf)
  }
} 

dldK <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dfdK=dfdK(param,data,N0),
           dldK=(-1/(2*sigma2))* (-2*log(besouros)*(1/aux)*dfdK + 2*log(aux)*(1/aux)*dfdK ) )
  output <- sum(df$dldK)
  return(output)
}

dldr <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dfdr=dfdr(param,data,N0),
           dldr=(-1/(2*sigma2))* (-2*log(besouros)*(1/aux)*dfdr + 2*log(aux)*(1/aux)*dfdr ) )
  output <- sum(df$dldr)
  return(output)
}

dldsigma2 <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  ti <- data[,1]
  Ni <- data[,2]
  frk <- ft(param,data,N0)
  output <- sum(-(1/(2*sigma2)) + (1/(2*sigma2^2))*(log(Ni)-log(frk))^2 )
  return(output)
}

gradllk <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  gradR <- dldr(param,data,N0)
  gradK <- dldK(param,data,N0)
  gradsigma2 <- dldsigma2(param,data,N0)
  gradiente <- c(gradK,gradR,gradsigma2)
  return(gradiente)
}

Ed2ldr2 <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dfdr=dfdr(param,data,N0),
           d2fdr2=d2fdr2(param,data,N0),
           d2ldr2=(-1/(2*sigma2))* (-2*log(aux)*((1/aux)*d2fdr2-(aux^(-2))*(dfdr^2)) +
                  2*(((aux^(-2))*dfdr-(aux^(-2))*dfdr*log(aux))*dfdr + log(aux)*(1/aux)*d2fdr2)) )
  output <- sum(df$d2ldr2)
  return(output)
}

Ed2ldK2 <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dfdK=dfdK(param,data,N0),
           d2fdK2=d2fdK2(param,data,N0),
           d2ldK2=(-1/(2*sigma2))* (-2*log(aux)*((1/aux)*d2fdK2-(aux^(-2))*(dfdK^2)) +
                                      2*(((aux^(-2))*dfdK-(aux^(-2))*dfdK*log(aux))*dfdK + log(aux)*(1/aux)*d2fdK2)  ) )
  output <- sum(df$d2ldK2)
  return(output)
}

Ed2ldsigma22 <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  ti <- data[,1]
  Ni <- data[,2]
  frk <- ft(param,data,N0)
  output <- sum((1/(2*sigma2^2)) - (1/(sigma2^3))*(sigma2))
  return(output)
}

Ed2ldKdr <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dfdK=dfdK(param,data,N0),
           dfdr=dfdr(param,data,N0),
           d2fdKdr=d2fdKdr(param,data,N0),
           d2ldKdr=(-1/(2*sigma2))* (-2*log(aux)*((1/aux)*d2fdKdr-(aux^(-2))*(dfdK*dfdr)) +
                                      2*(((aux^(-2))*dfdK-(aux^(-2))*dfdK*log(aux))*dfdK + log(aux)*(1/aux)*d2fdKdr)))
  output <- sum(df$d2ldKdr)
  return(output)
}

Ed2ldrdsigma2 <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dfdr=dfdr(param,data,N0),
           d2ldrdsigma2=(1/(2*sigma2^2))* (-2*log(aux)*(1/aux)*dfdr + 2*log(aux)*(1/aux)*dfdr ) )
  output <- sum(df$d2ldrdsigma2)
  return(output)
}

Ed2ldKdsigma2 <- function(param,data,N0){
  K <- param[1]
  r <- param[2]
  sigma2 <- param[3]
  df <- data%>%
    mutate(aux=ft(param = param,data=data,N0 = N0),
           dfdk=dfdr(param,data,N0),
           d2ldkdsigma2=(1/(2*sigma2^2))* (-2*log(aux)*(1/aux)*dfdk + 2*log(aux)*(1/aux)*dfdk))
  output <- sum(df$d2ldkdsigma2)
  return(output)
}

fisher_inf <- function(param,data,N0){
  A11 <- Ed2ldK2(param,data,N0)
  A22 <- Ed2ldr2(param,data,N0)
  A33 <- Ed2ldsigma22(param,data,N0)
  A12 <- A21 <- Ed2ldKdr(param,data,N0)
  A13 <- A31 <- Ed2ldKdsigma2(param,data,N0)
  A23 <- A32 <- Ed2ldrdsigma2(param,data,N0)
  FI <- matrix(c(A11,A21,A31,A12,A22,A32,A13,A23,A33),3,3)
  return(FI)
}

theta_k=c(1000,.1,1)
data=dados[-1,]
N0=2
EF <- function(theta_k,data,N0,epsilon,criterium){
  iter=0
  while (criterium>epsilon) {
    iter <- iter + 1
    theta_k_plus_1 <- theta_k + solve(-fisher_inf(param=theta_k,data=data,N0=2))%*%gradllk(param=theta_k,data=data,N0=2)
    criterium <- abs(llk(param = theta_k_plus_1,data = data,N0=2)/llk(param = theta_k,data = data,N0=2) - 1)
    theta_k <- theta_k_plus_1
  }
  output <- list(iter=iter,theta=theta_k,llk=llk(param=theta_k,data,N0))
  return(output)
}
EF(theta_k = c(1000,.1,1),data = dados[-1,],N0 = 2,epsilon = .00001,criterium = 1)

# Teste das derivadas do gradiente Q e da hessiana em pontos arbitrários
# numDeriv::grad(Q,x=c(1000,0.5),data=dados[-1,],N0=2)
# gradQ(param=c(1000,0.5),data=dados[-1,],N0=2)
# numDeriv::hessian(Q,x=c(1000,0.5),data=dados[-1,],N0=2)
# hessianQ(param=c(1000,.5),data=dados[-1,],N0=2)
# ft_vec <- function(param,t,N0){
#   k <- param[1]
#   r <- param[2]
#   output <- (k*N0) / (N0 + (k-N0)*exp(-r*t))
#   return(output)
# }
# ggplot(dados,aes(x=dias,y=besouros))+
#   geom_point()+
#   geom_function(fun = ft_vec,args=list(param=c(1033.5623,0.1179361),N0=2))+
#   theme_bw()
# ggsave(filename = "D:/UNICAMP/Disciplinas/COMPUTACIONAL/Atividade2/resolucao/imgs/NR.pdf",width = 8,height = 4)
# Q <- function(param,data,N0){
#   ti <- data[,1]
#   Ni <- data[,2]
#   frk <- ft(param,data,N0)
#   output <- sum((Ni-frk)^2)
#   return(output)
# }
grid <- expand.grid(r=seq(-.5,2,.1),k=seq(1,1500,10))
df <- grid%>%
  rowwise()%>%
  mutate(valueQ=llk(param = c(k,r,1),data = dados[-1,],N0=2))
df%>%
  ggplot(aes(x=k,y=r,z=valueQ))+
  geom_contour_filled()+theme_bw()
# ggsave(filename = "D:/UNICAMP/Disciplinas/COMPUTACIONAL/Atividade2/resolucao/imgs/space.pdf",width = 8,height = 4)
theta_k=c(1000,.5,1)
data=dados[-1,]
N0=2
iter=0
iter <- iter + 1
theta_k_plus_1 <- theta_k - solve(fisher_inf(param=theta_k,data=data,N0=2))%*%gradllk(param=theta_k,data=data,N0=2)
theta_k_plus_1%>%round(4)
criterium <- abs(llk(param = theta_k_plus_1,data = data,N0=2)/llk(param = theta_k,data = data,N0=2) - 1)
theta_k <- theta_k_plus_1
