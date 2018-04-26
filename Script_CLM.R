#New script to define the function 

library(MASS)
library(PLNmodels)
library(poilog)
library(nloptr)


##################################################################################################
#Notation

#On note Y la matrice des observations. Y est une matrice de taille n*p où n est le nombre d'observations (ie le nombre de répétitions) et p le nombre de variables observées
#On note X la matrice des covariables. X est une matrice de taille n*d ou d représente le nombre de covariables observées.
#On note param le vecteur contenant les paramètres de la loi Poisson log normale.
#   Ceux-ci sont disposés dans l'ordre suivant : les p*d premiers paramètres correspondent aux coefficients de corrélation avec les covariables. Il sont rangés par colonne : les d premières variables sont mu_11 jusqu'à mu_d1 (ie le coeff de corrélation pour chaque covariable par rapport à la première variable observée).
#                                                Les p paramètres suivant correspondent aux facteurs de variance de chacune des variables observées
#                                                Enfin tous les paramètres restants correspondent aux paramètres de covariation, rangé par ligne : Sigma_12...Sigma_1p,Sigma_23...Sigma_2p...Sigma_p-1,p


#################################################################################################
#Definition of the composite-likelihood function :
CL_f<-function(param,Y,X){
  n=nrow(Y)
  p=ncol(Y)
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covariance
  Xmu=X%*%Mu
  CL<-0
  for(i in 1:n){
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        CL<-CL+log(dbipoilog(Y[i,j],Y[i,k],mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k]))
      }
    }
  }
return(CL)
}

neg_CL<-function(param,Y,X){
  return(-CL_f(param,Y,X))
}

#Definition of the gradient function of the composite likelihood :
grad_CL_f<-function(param,Y,X){
  n=nrow(Y)
  p=ncol(Y)
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)
  Xmu=X%*%Mu
  grad_CL=rep(0,length(param))
  for (j in 1:p){#On va dériver par rapport à chacun des paramètres de moyenne
    for (k in 1:p){
      if (k != j){
        for (i in 1:n){
          pjki=dbipoilog(Y[i,j],Y[i,k],mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
          pjki1=dbipoilog(Y[i,j]+1,Y[i,k],mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
          grad_CL[((j-1)*d+1):(j*d)]=grad_CL[((j-1)*d+1):(j*d)]+(X[i,]/pjki)*(Y[i,j]*pjki - (Y[i,j]+1)*pjki1 )
        }
      }
    }
  }
  for (j in 1:p){#On dérive par rapport à chacun des éléments de variance
    for(k in 1:p){
      if (k != j){
        for (i in 1:n){
          pjki=dbipoilog(Y[i,j],Y[i,k],mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
          pjk1i=dbipoilog(Y[i,j],Y[i,k]+1,mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
          pjk2i=dbipoilog(Y[i,j],Y[i,k]+2,mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
          grad_CL[d*p+j]<-grad_CL[d*p+j]+ ((1/pjki)*(Y[i,k]*Y[i,k]*pjki-((Y[i,k]+1)*(Y[i,k]+1)+Y[i,k]*(Y[i,k]+1))*pjk1i+(Y[i,k]+1)*(Y[i,k]+2)*pjk2i))
        }
      }
    }
  }
  compteur=1
  for(j in 1:(p-1)){ #On dérive ici par rapport au terme de covariance
    for(k in (j+1):p){
      for (i in 1:n){
        pjki=dbipoilog(Y[i,j],Y[i,k],mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
        pj1ki=dbipoilog(Y[i,j]+1,Y[i,k],mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
        pjk1i=dbipoilog(Y[i,j],Y[i,k]+1,mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
        pj1k1i=dbipoilog(Y[i,j]+1,Y[i,k]+1,mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
        grad_CL[d*p+p+compteur]<-grad_CL[d*p+p+compteur]+ (Y[i,k]*Y[i,j]*pjki-Y[i,k]*(Y[i,j]+1)*pj1ki-(Y[i,k]+1)*Y[i,j]*pjk1i+(Y[i,k]+1)*(Y[i,j]+1)*pj1k1i)
      }
      compteur=compteur+1
    }
  }
  return(grad_CL)
}

neg_grad_CL<-function(param,Y,X){
  return(-grad_CL_f(param,Y,X))
}

########################################################################################################################################################################
### Optimisation des paramètres
ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
opts <- list( "algorithm" = paste("NLOPT_LD", ctrl$method, sep="_"),
              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
              "print_level" = max(0,ctrl$trace-1))
param_optimaux<-nloptr(x0=Beta.vem.tmp, eval_f=neg_CL, eval_grad_f=neg_grad_CL, 
                       # lb = c(rep(-Inf, d), 0), ub = rep(Inf, (d+1)), 
                       opts=opts, Y=Y, X=X)

#################################################################################################################################################################
##Simulation de donnees
#Pour les paramètres additionnels on notera ici n le nombre de données à simulées 
Observations_simulees<-function(n,p,X,param){
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covariance
  Xmu=X%*%Mu
  Obs=c()
  for(i in 1:n){
    Y=rep(0,p)
    Z<-mvrnorm(n = 1, mu=rep(0,p), Sigma=Rho, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    for (k in 1:p){
      Y[k]=rpois(1,exp(Xmu[i,k]+Z[k]))
    }
    Obs=rbind(Obs,Y)
  }
  return(Obs)
}



Observation_rpoilog<-function(n,p,X,param){
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covariance
  Xmu=X%*%Mu
  Obs=c()
  for(i in 1:n){
    Y=rpoilog(S=p,Xmu[i,],Rho)
    Obs=rbind(Obs,Y)
  }
  return(Obs)
}


########################################################################################################################
##Approximation

#model 
ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
opts <- list( "algorithm" = paste("NLOPT_LD", ctrl$method, sep="_"),
              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
              "print_level" = max(0,ctrl$trace-1))

param_0<-function(Y,X){ #Fonction pour trouver les paramètres initiaux à partir de la vem
  p=ncol(Y)
  Model<-PLN(Y~-1+X)
  Sigma=Model$model_par$Sigma
  param=c(Model$model_par$Theta,diag(Sigma))
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      a=Sigma[j,k]
      param=c(param,a)
    }
  }
  return(param)
}

#######################################################################################################################
## Application numérique

param=c(2,1.5,1,1,-0.5)
p=2
n=100
X=matrix(1,n,1)

Obs_1=Observations_simulees(n,p,X,param)
x_0=param_0(Obs_1,X)


ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
opts <- list( "algorithm" = paste("NLOPT_LD", ctrl$method, sep="_"),
              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
              "print_level" = max(0,ctrl$trace-1))
param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL, 
                       # lb = c(rep(-Inf, d), 0), ub = rep(Inf, (d+1)), 
                       opts=opts, Y=Obs_1, X=X)
print(param_optimaux$solution)
param_optimaux_MLE<-bipoilogMLE(Obs_1,startVals = x_0)
#results :
#param_optimaux$solution
#[1]  2.1385876  1.5187045  0.9962633  1.2455938 -0.1394972
#x_0
#[1]  2.1406709  1.5261414  1.0368881  1.2523913 -0.1396207


###Application numérique 2 :
param=c(2,1.5,1,1,1,1,-0.5,0,-0.5)
p=3
n=20
X=matrix(1,n,1)

Obs_2=Observations_simulees(n,p,X,param)
x_0=param_0(Obs_2,X)

ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
opts <- list( "algorithm" = paste("NLOPT_LD", ctrl$method, sep="_"),
              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
              "print_level" = max(0,ctrl$trace-1))
param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL, 
                       # lb = c(rep(-Inf, d), 0), ub = rep(Inf, (d+1)), 
                       opts=opts, Y=Obs_1, X=X,lb=c(0,0,0,1.0e-4,1.0e-4,1.0e-4,-1,-1,-1),ub=c(Inf,Inf,Inf,Inf,Inf,Inf,1,1,1))

#> param_optimaux$solution
#[1]  2.1676392  1.3031737  1.1499010  0.7337491  1.0157961  1.8099489 -0.4129029 -0.0132766
#[9] -0.6377163
#> x_0
#[1]  2.04036268  1.13056145  1.00790819  0.78550596  0.90326739  1.63775291 -0.41271376
#[8] -0.01361177 -0.63722910

##################################################################################################################################
##Application numérique de vérification :
#Perturber le résultats du PLN pour voir si dans ces conditions on bouge.
#Calculer si effectivement avec les paramètres finaux on trouve bien le maximum de la fonction de vraisemblance

param=c(2,1.5,1,1,-0.5)
p=2
n=20
X=matrix(1,n,1)

Obs_1=Observations_simulees(n,p,X,param)
Obs_1=Observation_rpoilog(n,p,X,param)
x_0=param_0(Obs_1,X)
ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
opts <- list( "algorithm" = "NLOPT_LD_MMA",
              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
              "print_level" = max(0,ctrl$trace-1))
param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL, 
                       lb = c(-Inf,-Inf,1.0e-4,1.0e-4,-0.999999), ub = c(Inf,Inf,Inf,Inf,0.9999999), 
                       opts=opts, Y=Obs_1, X=X)
print(grad_CL_f(param_optimaux$solution,Obs_1,X))
param_optim<-optim(par = x_0,fn = neg_CL,gr=neg_grad_CL,method="L-BFGS-B",lower =c(-Inf,-Inf,1.0e-7,1.0e-7,-0.999999),Y=Obs_1,X=X)
print(grad_CL_f(param_optim$par,Obs_1,X))
param_optim_MLE=bipoilogMLE(n1=Obs_1,startVals = x_0 )
print(grad_CL_f(param_optim_MLE$par,Obs_1,X))


#Ici on a regardé par rapport au critère du gradient. On peut aussi comparer à l'optimum réel de la fonction.

##Rédigeons un exemple dans le cas où il y a des covariables.
n=30
p=3
d=2
X1=runif(n,-1,1)
X2=runif(n,0.5,1.3)
X=cbind(X1,X2)
param=c(0.3,0.4,-1.2,0.1,1.4,-0.02,0.2,0.2,0.3,-0.1,0,-0.1)
Obs_3=Observations_simulees(n,p,X,param)
x_0=param_0(Obs_3,X)

ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
opts <- list( "algorithm" = "NLOPT_LD_MMA",
              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
              "print_level" = max(0,ctrl$trace-1))
param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL, 
                       lb = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,1.0e-7,1.0e-7,1.0e-7,-0.999999,-0.999999,-0.999999), ub = c(rep(Inf,9),0.99999,0.99999,0.99999), 
                       opts=opts, Y=Obs_3, X=X)
print(grad_CL_f(param_optimaux$solution,Obs_3,X))
param_optim<-optim(par = x_0,fn = neg_CL,gr=neg_grad_CL,method="L-BFGS-B",lower =c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,1.0e-7,1.0e-7,1.0e-7,-0.999999,-0.999999,-0.999999),upper = c(rep(Inf,9),0.99999,0.99999,0.99999),Y=Obs_3,X=X)
print(grad_CL_f(param_optim$par,Obs_3,X))

