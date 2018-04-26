#Installation des packages 
library(MASS)
library(PLNmodels)
library(nloptr)
library(poilog)
#Simulating Poisson Log-Normale Count.

PLN_dist <- function(mu,Sigma){
  Donnees<-rep(0,length(mu))
  Z<-mvrnorm(n = 1, rep(0,length(mu)), Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  for (k in 1:length(mu)){
    Donnees[k]=rpois(1,exp(mu[k]+Z[k]))}
  return(Donnees)
}



Observation<-function(mu,Sigma,n){Obs=PLN_dist(mu,Sigma)
  for (k in 1:n){
    Obs=rbind(Obs, PLN_dist(mu,Sigma))}
  return(Obs)}

Observation_bis<-function(mu,Sigma,n){Obs=rpoilog(length(mu),mu,Sigma)
  for (k in 1:n){
    Obs=rbind(Obs, rpoilog(length(mu),mu,Sigma))}
  return(Obs)}

##En utilisant le package poilog :
vdbipoilog<-Vectorize(dbipoilog,vectorize.args =c("n1","n2") )


  CL_etape1<-function(param,Data){
    a=c()
    for (k in 1:nrow(Data)){
      a=c(a,dbipoilog(n1=Data[k,1],n2=Data[k,2],mu1=param[1],mu2=param[2],sig1 =max(param[3],1.0e-4),sig2 = max(param[4],1.0e-4), rho = param[5]))
    }
    a
  }#Cette fonction renvoie un vecteur avec la valeur de la densité pour chaque observation.
  eval_function<-function(param,Data){
    -sum(log(CL_etape1(param,Data)))
  }
  grad_function<-function(param,Data){
    Data_1=cbind(Data[,1]+1,Data[,2])
    Data_11=cbind(Data[,1]+2,Data[,2])
    Data_2=cbind(Data[,1],Data[,2]+1)
    Data_22=cbind(Data[,1],Data[,2]+2)
    Data_12=cbind(Data_1[,1],Data_2[,2])
    a<-CL_etape1(param,Data)
    a_k<-CL_etape1(param,Data_1)
    a_kk<-CL_etape1(param,Data_11)
    a_j<-CL_etape1(param,Data_2)
    a_jj<-CL_etape1(param,Data_22)
    a_jk<-CL_etape1(param,Data_12)
    mu_k=Data[,1]-Data_1[,1]* a_k/a #derivee par rapport à mu_k
    mu_j=Data[,2]-Data_2[,2]* a_j/a
    sigma_kk=Data[,2]*Data[,2]-(Data_2[,2]*Data_2[,2] + Data[,2]*Data_2[,2])*a_j/a + Data_2[,2]*Data_22[,2]*a_jj/a
    sigma_jj=Data[,1]*Data[,1]-(Data_1[,1]*Data_1[,1]+Data[,1]*Data_1[,1])*a_k/a+ Data_1[,1]*Data_11[,1]*a_kk/a
    sigma_jk=Data[,1]*Data[,2]-Data_1[,1]*Data[,2]*a_k/a-Data[,1]*Data_12[,2]*a_j/a+Data_1[,1]*Data_2[,2]*a_jk/a
    c(sum(-mu_k),sum(-mu_j),sum(-sigma_kk),sum(-sigma_jj),sum(-sigma_jk))  
  }
  opts <- list("algorithm"="NLOPT_LD_TNEWTON")
  param_optimaux=nloptr(param_O,
                        eval_function,
                        grad_function,
                        lb=c(-Inf,-Inf,0,1.0e-4,-1),
                        ub=c(+Inf,+Inf,+Inf,+Inf,1),
                        Data=Data_test,
                        opts=opts)  
  
  ##Test 1 :
  Sigma=matrix(c(1,-0.5,-0.5,1),2,2)
  mu=c(2,1.5)
  Obs_1=Observation(mu,Sigma,100)
  model=PLN(Obs_1)
  mu_O=model$model_par$Theta[1:2]
  Sigma_O=model$model_par$Sigma[1:2,1:2]
  param_O=c(mu_O,Sigma_O[1,1],Sigma_O[2,2],Sigma_O[1,2])
  opts <- list("algorithm"="NLOPT_LD_TNEWTON_RESTART",
               "xtol_rel"=1.0e-4)
  param_optimaux=nloptr(param_O,
                        eval_function,
                        grad_function,
                        lb=c(-Inf,-Inf,0,1.0e-4,-1),
                        ub=c(+Inf,+Inf,+Inf,+Inf,1),
                        Data=Obs_1,
                        opts=opts) 
  #Resultats : mu_O  [1] 1.922906 1.561254
  #Sigma_O :          [,1]      [,2]
 #          [1,]  1.171024 -0.684009
 #          [2,] -0.684009  1.303741
 # Minimization using NLopt version 2.4.2 
  
#  NLopt solver status: -1 ( NLOPT_FAILURE: Generic failure code. )
  
#  Number of Iterations....: 37 
# Termination conditions:  xtol_rel: 1e-04 
#  Number of inequality constraints:  0 
# Number of equality constraints:    0 
#  Current value of objective function:  663.043101411342 
 # Current value of controls: 2.134131 1.306078 1.085642 1.22029 -0.5285077
  
##Faisons plusieurs simulations pour voir les écarts :
  nombre_donnees=c(1:5)*10
  Sigma=matrix(c(1,-0.5,-0.5,1),2,2)
  mu=c(2,1.5)
  opts <- list("algorithm"="NLOPT_LD_TNEWTON_RESTART",
               "xtol_rel"=1.0e-4)
  Result_CL=list(c(),c(),c(),c(),c())
  Result_PLN=list(c(),c(),c(),c(),c())
  for(k in nombre_donnees){
    for (j in 1:20){ 
      Obs_1=Observation(mu,Sigma,k)
      model=PLN(Obs_1)
      mu_O=model$model_par$Theta[1:2]
      Sigma_O=model$model_par$Sigma[1:2,1:2]
      param_O=c(mu_O,Sigma_O[1,1],Sigma_O[2,2],Sigma_O[1,2])
      param_optimaux=nloptr(param_O,
                            eval_function,
                            grad_function,
                            lb=c(-Inf,-Inf,0,1.0e-4,-1),
                            ub=c(+Inf,+Inf,+Inf,+Inf,1),
                            Data=Obs_1,
                            opts=opts) 
      Result_PLN[[k/10]]=rbind(Result_PLN[[k/10]],param_O)
      Result_CL[[k/10]]=rbind(Result_CL[[k/10]],param_optimaux$solution)
    }
  }
  
 #Fonction MLE du package poi-log 
  bipoilogMLE(Obs_1,n2 = NULL,
              startVals = param_O,
              nboot = 0, zTrunc = TRUE, file = NULL,
              method = "BFGS", control = list(maxit=1000))
  
  #En utilisant optim et non nloptr
  param_optim<- optim(param_O,eval_function,
                      gr=grad_function,
                      method="L-BFGS-B",
                      lower=c(0,0,0,1.0e-4,-1),
                      upper=c(Inf,Inf,Inf,Inf,1),
                      Data=Obs_1)
  
  Sigma=matrix(c(4,-0.5,-0.5,4),2,2)
  mu=c(4,5)
  Obs_1=Observation(mu,Sigma,k)
  
  nombre_donnees=c(1:5)*10
  Sigma=matrix(c(1,-0.5,-0.5,1),2,2)
  mu=c(2,1.5)
  opts <- list("algorithm"="NLOPT_LD_TNEWTON_RESTART",
               "xtol_rel"=1.0e-4)
  Result_CL=list(c(),c(),c(),c(),c())
  Result_PLN=list(c(),c(),c(),c(),c())
  for(k in nombre_donnees){
    for (j in 1:20){ 
      Obs_1=Observation(mu,Sigma,k)
      model=PLN(Obs_1)
      mu_O=model$model_par$Theta[1:2]
      Sigma_O=model$model_par$Sigma[1:2,1:2]
      param_O=c(mu_O,Sigma_O[1,1],Sigma_O[2,2],Sigma_O[1,2])
      param_optim<- optim(param_O,eval_function,
                          gr=grad_function,
                          method="L-BFGS-B",
                          lower=c(0,0,0,1.0e-4,-1),
                          upper=c(Inf,Inf,Inf,Inf,1),
                          Data=Obs_1)
      Result_PLN[[k/10]]=rbind(Result_PLN[[k/10]],param_O)
      Result_CL[[k/10]]=rbind(Result_CL[[k/10]],param_optim$par)
    }
  }
  
#les resultats sauvés en image le sont pour 13 données
  mpg_1<-cbind((Result_CL[[1]][1:13,1]-mu[1])/mu[1],(Result_CL[[2]][1:13,1]-mu[1])/mu[1],(Result_CL[[3]][1:13,1]-mu[1])/mu[1],(Result_CL[[4]][1:13,1]-mu[1])/mu[1],(Result_CL[[5]][1:13,1]-mu[1])/mu[1]) 
  boxplot(mpg_1[,1],mpg_1[,2],mpg_1[,3],mpg_1[,4],mpg_1[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à mu1 par CL",
          +         xlab="Nombre de données", ylab="Ecart relatif")
  
  mpg_2<-cbind((Result_CL[[1]][1:13,2]-mu[2])/mu[2],(Result_CL[[2]][1:13,2]-mu[2])/mu[2],(Result_CL[[3]][1:13,2]-mu[2])/mu[2],(Result_CL[[4]][1:13,2]-mu[2])/mu[2],(Result_CL[[5]][1:13,2]-mu[2])/mu[2]) 
  boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à mu2 par CL",
          +         xlab="Nombre de données", ylab="Ecart relatif")
  
  mpg<-cbind((Result_CL[[1]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_CL[[2]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_CL[[3]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_CL[[4]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_CL[[5]][1:13,3]-Sigma[1,1])/Sigma[1,1]) 
   boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à Sigma11 par CL",
            +         xlab="Nombre de données", ylab="Ecart relatif")
   
   mpg<-cbind((Result_CL[[1]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_CL[[2]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_CL[[3]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_CL[[4]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_CL[[5]][1:13,4]-Sigma[2,2])/Sigma[2,2]) 
    boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à Sigma22 par CL",
             +         xlab="Nombre de données", ylab="Ecart relatif")
    
    mpg<-cbind((Result_CL[[1]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_CL[[2]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_CL[[3]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_CL[[4]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_CL[[5]][1:13,5]-Sigma[1,2])/Sigma[1,2]) 
    boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à Sigma12 par CL",
            +         xlab="Nombre de données", ylab="Ecart relatif")
    
     mpg<-cbind((Result_PLN[[1]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_PLN[[2]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_PLN[[3]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_PLN[[4]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_PLN[[5]][1:13,5]-Sigma[1,2])/Sigma[1,2]) 
     boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à Sigma12 par PLN", xlab="Nombre de données", ylab="Ecart relatif")
     mpg<-cbind((Result_PLN[[1]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_PLN[[2]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_PLN[[3]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_PLN[[4]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_PLN[[5]][1:13,4]-Sigma[2,2])/Sigma[2,2]) 
     boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à Sigma22 par PLN", xlab="Nombre de données", ylab="Ecart relatif")
     mpg<-cbind((Result_PLN[[1]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_PLN[[2]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_PLN[[3]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_PLN[[4]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_PLN[[5]][1:13,3]-Sigma[1,1])/Sigma[1,1])
     boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à Sigma11 par PLN", xlab="Nombre de données", ylab="Ecart relatif")
     mpg<-cbind((Result_PLN[[1]][1:13,2]-mu[2])/mu[2],(Result_PLN[[2]][1:13,2]-mu[2])/mu[2],(Result_PLN[[3]][1:13,2]-mu[2])/mu[2],(Result_PLN[[4]][1:13,2]-mu[2])/mu[2],(Result_PLN[[5]][1:13,2]-mu[2])/mu[2]) 
     boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à mu2 par PLN", xlab="Nombre de données", ylab="Ecart relatif")
    mpg<-cbind((Result_PLN[[1]][1:13,1]-mu[1])/mu[1],(Result_PLN[[2]][1:13,1]-mu[1])/mu[1],(Result_PLN[[3]][1:13,1]-mu[1])/mu[1],(Result_PLN[[4]][1:13,1]-mu[1])/mu[1],(Result_PLN[[5]][1:13,1]-mu[1])/mu[1])
    boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="bisque",main="Ecart relatif à mu1 par PLN", xlab="Nombre de données", ylab="Ecart relatif")
    
    
    nombre_donnees=c(1:5)*10
    Sigma=matrix(c(1,-0.5,-0.5,1),2,2)
    mu=c(2,1.5)
    opts <- list("algorithm"="NLOPT_LD_TNEWTON_RESTART",
                 "xtol_rel"=1.0e-4)
    Result_MLE=list(c(),c(),c(),c(),c())
    Result_PLN=list(c(),c(),c(),c(),c())
    for(k in nombre_donnees){
      for (j in 1:15){ 
        Obs_1=Observation(mu,Sigma,k)
        model=PLN(Obs_1)
        mu_O=model$model_par$Theta[1:2]
        Sigma_O=model$model_par$Sigma[1:2,1:2]
        param_O=c(mu_O,Sigma_O[1,1],Sigma_O[2,2],Sigma_O[1,2])
        param_optim<- bipoilogMLE(Obs_1,startVals=param_O,
                                nboot = 0, zTrunc = TRUE, file = NULL, 
                            method = "BFGS", control = list(maxit=1000))
        Result_PLN[[k/10]]=rbind(Result_PLN[[k/10]],param_O)
        Result_MLE[[k/10]]=rbind(Result_MLE[[k/10]],param_optim$par)
      }
    }
    
    mpg<-cbind((Result_MLE[[1]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_MLE[[2]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_MLE[[3]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_MLE[[4]][1:13,5]-Sigma[1,2])/Sigma[1,2],(Result_MLE[[5]][1:13,5]-Sigma[1,2])/Sigma[1,2]) 
    boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="light blue",main="Ecart relatif à Sigma12 par MLE", xlab="Nombre de données", ylab="Ecart relatif")
    mpg<-cbind((Result_MLE[[1]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_MLE[[2]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_MLE[[3]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_MLE[[4]][1:13,4]-Sigma[2,2])/Sigma[2,2],(Result_MLE[[5]][1:13,4]-Sigma[2,2])/Sigma[2,2]) 
    boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="light blue",main="Ecart relatif à Sigma22 par MLE", xlab="Nombre de données", ylab="Ecart relatif")
    mpg<-cbind((Result_MLE[[1]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_MLE[[2]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_MLE[[3]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_MLE[[4]][1:13,3]-Sigma[1,1])/Sigma[1,1],(Result_MLE[[5]][1:13,3]-Sigma[1,1])/Sigma[1,1])
    boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="light blue",main="Ecart relatif à Sigma11 par MLE", xlab="Nombre de données", ylab="Ecart relatif")
    mpg<-cbind((Result_PLN[[1]][1:13,2]-mu[2])/mu[2],(Result_MLE[[2]][1:13,2]-mu[2])/mu[2],(Result_MLE[[3]][1:13,2]-mu[2])/mu[2],(Result_MLE[[4]][1:13,2]-mu[2])/mu[2],(Result_MLE[[5]][1:13,2]-mu[2])/mu[2]) 
    boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="light blue",main="Ecart relatif à mu2 par MLE", xlab="Nombre de données", ylab="Ecart relatif")
    mpg<-cbind((Result_MLE[[1]][1:13,1]-mu[1])/mu[1],(Result_MLE[[2]][1:13,1]-mu[1])/mu[1],(Result_MLE[[3]][1:13,1]-mu[1])/mu[1],(Result_MLE[[4]][1:13,1]-mu[1])/mu[1],(Result_MLE[[5]][1:13,1]-mu[1])/mu[1])
    boxplot(mpg[,1],mpg[,2],mpg[,3],mpg[,4],mpg[,5],names = c("10","20","30","40","50"),col="light blue",main="Ecart relatif à mu1 par MLE", xlab="Nombre de données", ylab="Ecart relatif")