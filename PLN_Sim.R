# Simulation of the estimation by composite likelihood of Poisson-log-Normale parameters

#Simulating Poisson Log-Normale Count.

PLN_dist <- function(mu,Sigma){
  Donnees=rep(0,length(mu))
  Z=mvrnorm(n = 1, rep(0,length(mu)), Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  for (k in 1:length(mu)){
   Donnees[k]=rpois(1,exp(mu[k]+Z[k]))}
  return(Donnees)
}

#Exemple :
Sigma=matrix(c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,3,0,0,0,0,0,0,3), nrow=6,ncol=6)
mu=c(0.1,3.4,9,3.02,2.2,1.45)
Obs=PLN_dist(mu,Sigma)