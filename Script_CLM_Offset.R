#New script to define the function 

library(MASS)
library(PLNmodels)
library(poilog)
library(nloptr)
library(reshape2)
library(ggplot2)

##################################################################################################
#Notation

#On note Y la matrice des observations. Y est une matrice de taille n*p où n est le nombre d'observations (ie le nombre de répétitions) et p le nombre de variables observées
#On note X la matrice des covariables. X est une matrice de taille n*d ou d représente le nombre de covariables observées.
#On note param le vecteur contenant les paramètres de la loi Poisson log normale.
#   Ceux-ci sont disposés dans l'ordre suivant : les p*d premiers paramètres correspondent aux coefficients de corrélation avec les covariables. Ils sont rangés par colonne : les d premières variables sont mu_11 jusqu'à mu_d1 (ie le coeff de corrélation pour chaque covariable par rapport à la première variable observée).
#                                                Les p paramètres suivant correspondent aux facteurs de variance de chacune des variables observées
#                                                Enfin tous les paramètres restants correspondent aux paramètres de covariation, rangé par ligne : Sigma_12...Sigma_1p,Sigma_23...Sigma_2p...Sigma_p-1,p


#################################################################################################
#Definition of the composite-likelihood function :

#La matrice O est une matrice n*p qui correspond à l'intercept (?)
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

#Definition of the gradient function of the composite likelihood :
grad_CL_f<-function(param,Y,O,X){
  n=nrow(Y)
  p=ncol(Y)
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=Rho+t(Rho)
  Xmu=O+X%*%Mu
  grad_CL=rep(0,length(param))
  #On calcule une fois pour toute les termes en bipoilog dont on a besoin
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

grad_CL_uni<-function(Y,X,O,d,p,Xmu,Sigma,Rho){
  #Fonction qui calcule le gradient pour une observation.
  ##Y est le vecteur des observations
  ##X est le vecteur des covariables
  ##O est le vecteur des offsets
  ##Xmu est la combinaison des covariables et facteurs de régression
  ## Sigma est le vecteur des variances
  ## Rho est la matrice des variances convariances
  #On commence par calculer tous les termes dont on a besoin 
  terms<-matrix(0,3*p,3*p)
  Corr=cov2cor(Rho)
  for(j in 2:p){
    for (k in 1:(j-1)){
      for (l in 1:3){
        for (m in 1:3){
          if (l+m<5){
            terms[(3*(k-1)+l),(3*(j-1)+m)]<- dbipoilog(Y[k]+l-1,Y[j]+m-1,mu1=Xmu[k],mu2= Xmu[j],sig1 =  sqrt(Sigma[k]),sig2=sqrt(Sigma[j]),rho = Corr[k,j])
          }
        }
      }
    }
  }
  terms<-terms+t(terms) 
  gradCL_calc<-rep(0,d*p+p+0.5*p*(p-1))
  for (j in 1:p){#dérivées par rapport aux termes de moyenne
    for (k in 1:p){
      if (j != k){
        gradCL_calc[((j-1)*d+1):(d*j)]<-gradCL_calc[((j-1)*d+1):(d*j)]+(X/terms[3*(j-1)+1,3*(k-1)+1])*(Y[j]*terms[3*(j-1)+1,3*(k-1)+1]-(Y[j]+1)*terms[3*(j-1)+2,3*(k-1)+1])
        
      }
    }
  }
  for (j in 1:p){#Dérivées par rapport aux termes de variance
    for (k in 1:p){
      if (j != k) {
        gradCL_calc[d*p+j]<-gradCL_calc[d*p+j]+ (Y[k]^2-(Y[k]+1)*(2*Y[k]+1)*(terms[3*(j-1)+1,3*(k-1)+2]/terms[3*(j-1)+1,3*(k-1)+1])+(Y[k]+1)*(Y[k]+2)*(terms[3*(j-1)+1,3*(k-1)+3]/terms[3*(j-1)+1,3*(k-1)+1]))
      }
    }
  }
  compteur=1
  liste_indices = list()
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      liste_indices[[compteur]]<-c(j,k)
      compteur=compteur+1
    }
  }
  for (l in 1:(0.5*p*(p-1))){
    j=liste_indices[[l]][1]
    k=liste_indices[[l]][2]
    gradCL_calc[d*p+p+l]<-Y[k]*Y[j]-(Y[k]+1)*Y[j]*(terms[3*(j-1)+1,3*(k-1)+2]/terms[3*(j-1)+1,3*(k-1)+1])-Y[k]*(Y[j]+1)*(terms[3*(j-1)+2,3*(k-1)+1]/terms[3*(j-1)+1,3*(k-1)+1]) + (Y[j]+1)*(Y[k]+1)*terms[3*(j-1)+2,3*(k-1)+2]/terms[3*(j-1)+1,3*(k-1)+1]
    
  }
  return(gradCL_calc)
}


grad_CL_opt<-function(param,Y,O,X){
  #Le vecteur param contient tous les paramètres dans l'ordre suivant : coefficents de régression, variances, covrariances
  n=nrow(Y)
  p=ncol(Y)
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=Rho+t(Rho)
  diag(Rho)<-Sigma #Rho est ici une matrice de variance-covariance
  Xmu=O+X%*%Mu
  grad_CL=rep(0,length(param))
  for (i in 1:n){#On répète sur les observations
    grad_CL<-grad_CL+grad_CL_uni(Y[i,],X[i,],O[i,],d,p,Xmu[i,],Sigma,Rho)
  }
  return(grad_CL)
}


neg_grad_CL<-function(param,Y,O,X){
  return(-grad_CL_opt(param,Y,O,X))
}

neg_grad_CL_n<-function(param,Y,O,X){
  return(-grad_CL_f(param,Y,O,X))
}
########################################################################################################################################################################################################################################################################
##Simulation de donnees
# Générer une matrice de variance-covariance
var_cov<-function(p){
  points_x=runif(p,-1,1)
  points_y=runif(p,-1,1)
  var=matrix(0,p,p)
  for(i in 1:(p-1)){
    for (j in (i+1):p){
      var[i,j]<-exp(-sqrt((points_x[i]-points_x[j])^2 + (points_y[i]-points_y[j])^2))
    }
  }
  var<-var+t(var)
  diag(var)<-rep(1,p)
  return(var)
}

#Pour les paramètres additionnels on notera ici n le nombre de données à simulées 

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
######################################################################################################################################
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


# ####################################################################################################################################
# ##Simulation
param=c(2,1.5,1,1,-0.5)
p=2
n=20
X=matrix(1,n,1)
O=matrix(1,n,p)

Obs_1=Observations_simulees_bis(n,p,X,O,param)
# x_0=param_0(Obs_1,O,X)
# ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
#              xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
# opts <- list( "algorithm" = "NLOPT_LD_MMA",
#               "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
#               "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
#               "print_level" = max(0,ctrl$trace-1))
# param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL, 
#                        lb = c(-Inf,-Inf,1.0e-4,1.0e-4,-0.999999), ub = c(Inf,Inf,Inf,Inf,0.9999999), 
#                        opts=opts, Y=Obs_1, X=X,O=O)
# print(grad_CL_f(param_optimaux$solution,Obs_1,O,X))
# param_optim<-optim(par = x_0,fn = neg_CL,gr=neg_grad_CL,method="L-BFGS-B",lower =c(-Inf,-Inf,1.0e-7,1.0e-7,-0.999999),Y=Obs_1,X=X,O=O)
# print(grad_CL_f(param_optim$par,Obs_1,O,X))
# 
# 
# 
# 
# n=30
# p=3
# d=2
# X1=runif(n,-1,1)
# X2=runif(n,0.5,1.3)
# X=cbind(X1,X2)
# param=c(0.3,0.4,-1.2,0.1,1.4,-0.02,0.2,0.2,0.3,-0.1,0,-0.1)
# O=matrix(1,n,p)
# Obs_3=Observations_simulees_bis(n,p,X,O,param)
# x_0=param_0(Obs_3,O,X)
# 
# ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
#              xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
# opts <- list( "algorithm" = "NLOPT_LD_MMA",
#               "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
#               "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
#               "print_level" = max(0,ctrl$trace-1))
# param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL, 
#                        lb = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,1.0e-7,1.0e-7,1.0e-7,-0.999999,-0.999999,-0.999999), ub = c(rep(Inf,9),0.99999,0.99999,0.99999), 
#                        opts=opts, Y=Obs_3, X=X,O=O)
# print(grad_CL_f(param_optimaux$solution,Obs_3,O,X))
# param_optim<-optim(par = x_0,fn = neg_CL,gr=neg_grad_CL,method="L-BFGS-B",lower =c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,1.0e-7,1.0e-7,1.0e-7,-0.999999,-0.999999,-0.999999),upper = c(rep(Inf,9),0.99999,0.99999,0.99999),Y=Obs_3,X=X,O)
# print(grad_CL_f(param_optim$par,Obs_3,O,X))

# #######################################################################################################
# ##Simulation de données pour voir la distribution de l'estimateur.
# 
# nb_simu = 250
# n=30
# p=3
# d=2
# X1=runif(n,-1,1)
# X2=runif(n,0.5,1.3)
# X=cbind(X1,X2)
# param=c(0.3,0.4,-1.2,0.1,1.4,-0.02,0.2,0.2,0.3,-0.1,0,-0.1)
# O=matrix(1,n,p)
# ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
#              xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
# opts <- list( "algorithm" = "NLOPT_LD_MMA",
#               "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
#               "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
#               "print_level" = max(0,ctrl$trace-1))
# 
# 
# param_estim=c()
# for (k in 1:nb_simu){
#   Obs_3=Observations_simulees_bis(n,p,X,O,param)
#   x_0=param_0(Obs_3,O,X)
#   param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL,
#                          lb = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,1.0e-7,1.0e-7,1.0e-7,-0.999999,-0.999999,-0.999999), ub = c(rep(Inf,9),0.99999,0.99999,0.99999),
#                          opts=opts, Y=Obs_3, X=X,O=O)
#   param_estim=rbind(param_estim,param_optimaux$solution-param)
# }
# for (j in 1:ncol(param_estim)){
#   param_estim[,j]<-param_estim[,j]/sqrt(var(param_estim[,j]))
# }
# write.table(param_estim, file = "/home/teuliere/PLN-Cl/Estim_250", append = FALSE, quote = TRUE, sep = " ",
#             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#             col.names = TRUE, qmethod = c("escape", "double"))
# 
# setwd(dir="/home/teuliere/PLN-Cl")
# Er=read.table("Estim_250",sep=" ",header = TRUE)
# Er=t(Er)
# for (k in 1:nrow(Er)){
#   Er[k,]<-Er[k,]*sqrt(var(Er[k,]))
# }
# for (j in 1:ncol(Er)){
#   Er[,j]<-Er[,j]/sqrt(var(Er[,j]))
# }
# Er2=melt(Er)
# ggplot(data=Er2,aes(x=value))+geom_histogram()+facet_wrap(~Var2,scales="free")+
# theme_bw()
# 
# #Estim_10 a été simulé avec les paramètres suivants : nb_simu=10 , n=30, p=3, d=2, X1=runif(n,-1,1), X2=runif(n,0.5,1.3),X=cbind(X1,X2), param=c(0.3,0.4,-1.2,0.1,1.4,-0.02,0.2,0.2,0.3,-0.1,0,-0.1), O=matrix(1,n,p)
# #ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
# #             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
# #opts <- list( "algorithm" = "NLOPT_LD_MMA",
# #              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
# #              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
# #              "print_level" = max(0,ctrl$trace-1))
