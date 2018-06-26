######################################################################################################
## Script Simu François
#######################################################################################################
library(MASS)
library(PLNmodels)
library(poilog)
library(nloptr)
library(reshape2)
library(ggplot2)

## Function composite likelihood
CL_f<-function(param,Y,O,X){
  n=nrow(Y)
  p=ncol(Y)
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covariance
  Xmu=O+X%*%Mu
  Corr=cov2cor(Rho)
  CL<-0
  for(i in 1:n){
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        CL<-CL+log(dbipoilog(Y[i,j],Y[i,k],mu1=Xmu[i,j],mu2=Xmu[i,k],sig1 = sqrt(Sigma[j]),sig2 = sqrt(Sigma[k]),rho = Corr[j,k]))
      }
    }
  }
  return(CL)
}

CL_f_uni<-function(param,Y,O,X,d,p,Xmu,Sigma,Rho){
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covariance
  Xmu=O+X%*%Mu
  Corr=cov2cor(Rho)
  CL<-0
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      CL<-CL+log(dbipoilog(Y[j],Y[k],mu1=Xmu[j],mu2=Xmu[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Corr[j,k]))
    }
  }
  return(CL)
}

neg_CL<-function(param,Y,O,X){
  return(-CL_f(param,Y,O,X))
}


Observations_simulees_bis<-function(n,p,X,O,param){
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covariance
  Xmu=X%*%Mu
  expOXMu=exp(O+Xmu)
  Z=mvrnorm(n,rep(0,p),Rho)
  means_mat=expOXMu*exp(Z)
  Y = matrix(rpois(n*p,means_mat ), n, p)
  return(list(Y=Y,Z=Z,mean=means_mat))
}

##Détermination des paramètres initiaux
param_0<-function(Y,O,X){ #Fonction pour trouver les paramètres initiaux à partir de la vem
  p=ncol(Y)
  Model<-PLN(Y~-1+X+offset(O))
  Sigma=Model$model_par$Sigma #Matrice de variance-covariance
  param=c(Model$model_par$Theta,diag(Sigma))
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      a=Sigma[j,k]
      param=c(param,a)
    }
  }
  return(param)
}

#Fonction qui renvoie la matrice de variance-covariance à partir des distances, alpha et sigma
var_cov_spatial<-function(points_x,points_y,alpha,sigma){
  p=length(points_x)
  var=matrix(0,p,p)
  for(i in 1:(p-1)){
    for (j in (i+1):p){
      var[i,j]<-exp(-alpha*sqrt((points_x[i]-points_x[j])^2 + (points_y[i]-points_y[j])^2))
    }
  }
  var<-var+t(var)
  diag(var)<-rep(1,p)
  return(sigma*var)
}




##Optimisation en utilisant d'abord nloptr qui est rapide, puis optim pour avoir une estimation de la hessienne

#Simulation des paramètres
p=5
d=2
n=50
nb_simu=250
O=matrix(1,n,p)
X=matrix(runif(n,-1,1),n,1)
if(d>1){
  for (j in 1:(d-1)){
    X<-cbind(X,runif(n,-1,1))
  }
}
mu<-runif(p*d,0,1)
#On simule les distances
points_x=runif(p,-1,1)
points_y=runif(p,-1,1)

alpha=1.0
sigma=3.0
Var_Cov=var_cov_spatial(points_x,points_y,alpha,sigma)
param<-c(mu,diag(Var_Cov))
for (j in 1:(p-1)){
  for (k in (j+1):p){
    a=Var_Cov[j,k]
    param=c(param,a)
  }
}


param_estim=c()
param_estim_norm=c()
diff_hessian=c()
##Simulation des observations
for (k in 1:nb_simu){
  Obs_3=Observations_simulees_bis(n,p,X,O,param)$Y
  x_0=param_0(Obs_3,O,X)
  x_init=c(x_0[1:(p*d)]) #PLN ne nous sort pas le modèle spatialisé. L'idée est donc d'approcher sigma et alpha. 
  #Pour sigma on fait la moyenne des variances (cf ci-dessous)
  #Pour alpha on fait :
  var_es<-c()
  compteur=1
  for (i in 1:(p-1)){
    for(j in (i+1):p){
      var_es<-c(var_es,-log(abs(x_0[p*d+p+compteur]/mean(x_0[(p*d+1):(p*d+p)])))*(1/ sqrt((points_x[i]-points_x[j])^2+(points_y[i]-points_y[j])^2)) )
      compteur= compteur+1
    }
  }
  x_init=c(x_init,mean(var_es),mean(x_0[(p*d+1):(p*d+p)]))
  ##paramètres de optim
  opts <- list("xtol_rel"=1*10^(-6), "ftol_rel"=0.0,"ftol_abs"=0.0,"maxeval"=1000,"check.derivatives"=FALSE,"algorithm"="NLOPT_GN_DIRECT_L")
  fonction_a_optimiser<-function(param){
    mu_trans=param[1:(p*d)]
    alpha_trans=param[p*d+1]
    sigma_trans=param[p*d+2]
    Var_Cov=var_cov_spatial(points_x,points_y,alpha_trans,sigma_trans)
    param_CL<-c(mu_trans,diag(Var_Cov))
    for (j in 1:(p-1)){
      for (i in (j+1):p){
        a=Var_Cov[j,i]
        param_CL=c(param_CL,a)
      }
    }
    return(neg_CL(param_CL,Obs_3,O,X))
  }
  ##Optimisation avec nloptr 
  lower=c(rep(-Inf,p*d),10^(-5),10^(-5))
  upper=c(rep(Inf,p*d),10^9,10^9)
  param_optimaux_nl<-nloptr(x_init,eval_f=fonction_a_optimiser,lb=lower,ub=upper,eval_grad_f=NULL,opts = opts)
  ##Optimisation avec optim
  param_optimaux<-optim(param_optimaux_nl$solution,fn=fonction_a_optimiser,eval_grad_f=NULL,method = "SANN",hessian = TRUE,opts = opts)
  param_estim=rbind(param_estim,param_optimaux$solution)
  Vp_inf=ginv(param_optimaux$hessian)
  grad=param_optimaux$gradient%*%t(param_optimaux$gradient)
  Godambe=Vp_inf%*%grad%*%Vp_inf
  param_estim_norm<-rbind(param_estim_norm,(param_optimaux$solution-param)/sqrt(diag(Godambe)))
}
E=data.frame(param_estim_norm)
E_prime=melt(E)
pltt<-ggplot(data=E_prime,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
  theme_bw()+ggtitle("Optimisation wrt spatial parameters")
