#################################################################################################
## FONCTION FOR SIMULATIONS OF THE CL MAXIMUM ESTIMATOR OF PLN

# ELsa TEULIERE Juin 2018

#Notations
#On note Y la matrice des observations. Y est une matrice de taille n*p où n est le nombre d'observations (ie le nombre de répétitions) et p le nombre de variables observées
#On note X la matrice des covariables. X est une matrice de taille n*d ou d représente le nombre de covariables observées.
#On note param le vecteur contenant les paramètres de la loi Poisson log normale.
#   Ceux-ci sont disposés dans l'ordre suivant : les p*d premiers paramètres correspondent aux coefficients de corrélation avec les covariables. Ils sont rangés par colonne : les d premières variables sont mu_11 jusqu'à mu_d1 (ie le coeff de corrélation pour chaque covariable par rapport à la première variable observée).
#                                                Les p paramètres suivant correspondent aux facteurs de variance de chacune des variables observées
#                                                Enfin tous les paramètres restants correspondent aux paramètres de covariation, rangé par ligne : Sigma_12...Sigma_1p,Sigma_23...Sigma_2p...Sigma_p-1,p


#################################################################################################

#Library

library(MASS)
library(PLNmodels)
library(poilog)
library(nloptr)
library(reshape2)
library(ggplot2)
################################################################################################
## Calcul of the composite likelihood



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
}

neg_CL<-function(param,Y,O,X){
  return(-CL_f(param,Y,O,X))
}

##########################################################################################################################
##Gradient

##Pour une observation
grad_CL_uni_st<-function(Y,X,O,d,p,Xmu,Sigma,Rho){
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
    gradCL_calc[d*p+p+l]<-2*(Y[k]*Y[j]-(Y[k]+1)*Y[j]*(terms[3*(j-1)+1,3*(k-1)+2]/terms[3*(j-1)+1,3*(k-1)+1])-Y[k]*(Y[j]+1)*(terms[3*(j-1)+2,3*(k-1)+1]/terms[3*(j-1)+1,3*(k-1)+1]) + (Y[j]+1)*(Y[k]+1)*terms[3*(j-1)+2,3*(k-1)+2]/terms[3*(j-1)+1,3*(k-1)+1])
    
  }
  return(gradCL_calc)
}


grad_CL_opt_st<-function(param,Y,O,X){
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
    grad_CL<-grad_CL+grad_CL_uni_st(Y[i,],X[i,],O[i,],d,p,Xmu[i,],Sigma,Rho)
  }
  return(grad_CL)
}


neg_grad_CL_st<-function(param,Y,O,X){
  return(-grad_CL_opt_st(param,Y,O,X))
}
################################################################################################################
##Matrice de Gadambe


## Fonction calculant pour une observation la valeur de la matrice hessienne et du gradient en les parametres param

Hessian_corr_st<-function(param,Y,X,O,d,p){ #Ici on rentre les paramètres d et p pour éviter ensuite dans la boucle sur les observations d'avoir à les recalculer
  #mise en forme des paramètre 
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=Rho+t(Rho)
  Xmu=O+X%*%Mu
  diag(Rho)<-Sigma
  Rho=cov2cor(Rho)
  #On commence par calculer tous les termes dont on aura besoin 
  terms<-matrix(0,5*p,5*p)
  for(j in 2:p){
    for (k in 1:(j-1)){
      for (l in 1:5){
        for (m in 1:5){
          if (l+m<7){
            terms[(5*(k-1)+l),(5*(j-1)+m)]<- dbipoilog(Y[k]+l-1,Y[j]+m-1,Xmu[k],Xmu[j],sqrt(Sigma[k]),sqrt(Sigma[j]),Rho[k,j])
          }
        }
      }
    }
  }
  terms<-terms+t(terms) #On peut faire cela car tous les blocs diagonaux sont nuls
  
  #Maintenant on rempli la matrice hessienne
  hess_mat<-matrix(0,d*p+p+0.5*p*(p-1),d*p+p+0.5*p*(p-1))
  #On commence par le bloc des ddérivées par rapport aux termes de moyenne
  #Les dérivées non croisées
  for (j in 1:p){ #on parcourt tous les paramètres mu
    bloc<-0
    for (k in 1:p){#puis tous les couples ou j intervient
      if (j != k){
        bloc<-bloc+ X%*%t(X) *(1/terms[5*(j-1)+1,5*(k-1)+1]^2)*((Y[j]+1)*(Y[j]+2)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+3,5*(k-1)+1]-(Y[j]+1)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+2,5*(k-1)+1]-(Y[j]+1)^2 * terms[5*(j-1)+2,5*(k-1)+1]^2)
      }
    }
    hess_mat[(d*(j-1)+1):(d*j),(d*(j-1)+1):(d*j)] <- bloc
  }
  #Les termes de dérivée croisée
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      hess_mat[(d*(j-1)+1):(d*j),(d*(k-1)+1):(d*k)]<- X%*%t(X) * (Y[j]+1)*(Y[k]+1)*(1/terms[5*(j-1)+1,5*(k-1)+1]^2)*(terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+2,5*(k-1)+2]-terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2])
    }
  }
  #On s'attaque maintenant au bloc des dérivées par rapport aux termes de variance
  #Tout d'abord les dérivées non croisées
  for (j in 1:p){
    for(k in 1:p){
      if (j != k){
        hess_mat[d*p+j,d*p+j]<-hess_mat[d*p+j,d*p+j]+(1/terms[5*(j-1)+1,5*(k-1)+1]^2)*(-(Y[k]+1)*(2*Y[k]+1)^2*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]+(4*Y[k]^2+12*Y[k]+7)*(Y[k]+1)*(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+3]-2*(Y[k]+1)*(Y[k]+2)*(Y[k]+3)*(2*Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+4]-(Y[k]+1)^2*(2*Y[k]+1)^2*terms[5*(j-1)+1,5*(k-1)+2]^2+2*(Y[k]+1)^2*(Y[k]+2)*(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+3]+(Y[k]+1)*(Y[k]+2)*(Y[k]+3)*(Y[k]+4)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+5]-(Y[k]+1)^2*(Y[k]+2)^2*terms[5*(j-1)+1,5*(k-1)+3]^2)
      }
    }
  }
  #Maintenant les termes de dérivée croisée avec les termes de moyenne
  #On commence par les dérivées mu_j sigma_jj
  for(j in 1:p){
    vect=0
    for (k in 1:p){
      if(j !=k){
        vect<- vect + (Y[j]+1)*(Y[k]+1)*(1/terms[5*(j-1)+1,5*(k-1)+1]^2)*((2*Y[k]+1)*(terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+2,5*(k-1)+2]-terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2])+(Y[k]+2)*(terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+3]-terms[5*(j-1)+2,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]))
      }
    }
    hess_mat[((j-1)*d+1):(j*d), d*p+j]<-X*vect
  }
  #On calcule maintenant les termes mu_j sigma_kk ou k et j diffèrent
  for (j in 1:p){#indice du sigma par rapport auquel on dérive
    for (k in 1:p){#indice du mu par rapport auquel on dérive
      if (k != j){
        hess_mat[((k-1)*d+1):(k*d),d*p+j]<- X * ((Y[k]+1)*(1/terms[5*(j-1)+1,5*(k-1)+1]^2)*(-(2*Y[k]+1) *terms[5*(j-1)+1,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1]+(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+1]*((2*Y[k]+3)* terms[5*(j-1)+1,5*(k-1)+3]-(Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+4])+(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]*((Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+3]-(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2])))
      }
    }
  }
  #Enfin on calcule les termes du type sigma_j sigma_k
  for (j in 2:p){
    for (k in 1:(j-1)){
      hess_mat[(d*p+k),d*p+j]<-(Y[j]+1)*(Y[k]+1)*(1/terms[5*(j-1)+1,5*(k-1)+1])^2*((2*Y[j]+1)*((2*Y[k]+1)*(terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1]-terms[5*(j-1)+1,5*(k-1)+2]*terms[5*(j-1)+2,5*(k-1)+1])+(Y[k]+2)*(terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+3]-terms[5*(j-1)+2,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]))+(Y[j]+2)*((2*Y[k]+1)*(terms[5*(j-1)+3,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]-terms[5*(j-1)+3,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1])+(Y[k]+2)*(terms[5*(j-1)+3,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]-terms[5*(j-1)+1,5*(k-1)+3]*terms[5*(j-1)+3,5*(k-1)+1])))
    }
  }
  #Reste maintenant les dérivées par rapport aux termes de covariance
  #On va créer une liste contenant les indices des termes de covariances
  compteur=1
  liste_indices = list()
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      liste_indices[[compteur]]<-c(j,k)
      compteur=compteur+1
    }
  }
  #On rempli d'abord la diagonale du bloc
  for (l in 1:(0.5*p*(p-1))){
    j=liste_indices[[l]][1]
    k=liste_indices[[l]][2]
    hess_mat[d*p+p+l,d*p+p+l]<- (4.0/terms[5*(j-1)+1,5*(k-1)+1]^2)*(-(Y[j]+1)*Y[k]^2*terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+1]+(Y[j]+1)*(Y[j]+2)*Y[k]^2*terms[5*(j-1)+3,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+1]+(Y[j]+1)*(Y[k]+1) * (2*Y[k]*Y[j]+2*Y[k]+2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1]-(Y[j]+1)*(Y[j]+2)*(Y[k]+1)*(2*Y[k]+1)*terms[5*(j-1)+3,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1]-(Y[j]+1)^2*Y[k]^2*terms[5*(j-1)+2,5*(k-1)+1]^2-2*(Y[j]+1)*(Y[k]+1)*Y[j]*Y[k]*terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]+2* (Y[j]+1)^2*(Y[k]+1)*Y[k]*terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+2,5*(k-1)+1]-Y[j]^2*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]+Y[j]^2*(Y[k]+1)*(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]-(Y[k]+1)*(Y[k]+2)*(Y[j]+1)*(2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]-Y[j]^2*(Y[k]+1)^2*terms[5*(j-1)+1,5*(k-1)+2]^2+2*(Y[k]+1)^2*Y[j]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+2]+(Y[j]+1)*(Y[k]+1)*(Y[j]+2)*(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+3,5*(k-1)+3]-(Y[j]+1)^2*(Y[k]+1)^2*terms[5*(j-1)+2,5*(k-1)+2]^2)
  }
  #On s'attaque aux termes croisés avec la moyenne
  for (l in 1 : (0.5*p*(p-1))){
    j=liste_indices[[l]][1]
    k=liste_indices[[l]][2]
    #On remplit tout d'abord les blocs mu_j sigma_jk
    hess_mat[((j-1)*d+1):(j*d),d*p+p+l]<-2*X*((Y[j]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]-Y[k]*terms[5*(j-1)+2,5*(k-1)+1]+(Y[j]+2)*Y[k]*terms[5*(j-1)+3,5*(k-1)+1]-(Y[j]+2)*(Y[k]+1)*terms[5*(j-1)+3,5*(k-1)+2]) +terms[5*(j-1)+2,5*(k-1)+1] *(-Y[j]*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]-(Y[j]+1)*Y[k]*terms[5*(j-1)+2,5*(k-1)+1]+(Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]) )
    #On remplit maintenant les mu_k sigma_jk
    hess_mat[((k-1)*d+1):(k*d),d*p+p+l] <- 2*X*((Y[j]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[k]+1)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]-Y[j]*terms[5*(j-1)+1,5*(k-1)+2]+(Y[k]+2)*Y[j]*terms[5*(j-1)+1,5*(k-1)+3]-(Y[k]+2)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+3]) +terms[5*(j-1)+1,5*(k-1)+2] *(-Y[k]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]-(Y[k]+1)*Y[j]*terms[5*(j-1)+1,5*(k-1)+2]+(Y[k]+1)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]) )
    #les termes sigma_jj sigma_jk
    hess_mat[d*p+j,(d+1)*p+l]<-(2*(Y[k]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[j]+1)*(2*Y[k]^2+3*Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]-3*(Y[j]+1)*(Y[k]+2)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+3]+(Y[j]+1)*(Y[k]+2)*(Y[k]+3)*terms[5*(j-1)+2,5*(k-1)+4]-Y[j]*(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]+Y[j]*(Y[k]+2)*(2*Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+3]- Y[j]*(Y[k]+2)*(Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+4])+(Y[k]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]+Y[j]*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]-(Y[k]+1)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2])*((Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+3]-(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]))
    #les termes sigma_kk sigma_jk
    hess_mat[d*p+k,(d+1)*p+l]<-(2*(Y[j]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[k]+1)*(2*Y[j]^2+3*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]-3*(Y[k]+1)*(Y[j]+2)*(Y[j]+1)*terms[5*(j-1)+3,5*(k-1)+2]+(Y[k]+1)*(Y[j]+2)*(Y[j]+3)*terms[5*(j-1)+4,5*(k-1)+2]-Y[k]*(2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]+Y[k]*(Y[j]+2)*(2*Y[j]+3)*terms[5*(j-1)+3,5*(k-1)+1]- Y[k]*(Y[j]+2)*(Y[j]+3)*terms[5*(j-1)+4,5*(k-1)+1])+(Y[j]*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]+Y[k]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]-(Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2])*((Y[j]+2)*terms[5*(j-1)+3,5*(k-1)+1]-(2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]))
    #Il nous reste les termes croisés de covariance qui sont tous nuls
  }
  gradCL_calc<-rep(0,d*p+p+0.5*p*(p-1))
  for (j in 1:p){#dérivées par rapport aux termes de moyenne
    for (k in 1:p){
      if (j != k){
        gradCL_calc[((j-1)*d+1):(d*j)]<-gradCL_calc[((j-1)*d+1):(j*d)]+X*(Y[j]-(Y[j]+1)*(terms[5*(j-1)+2,5*(k-1)+1]/terms[5*(j-1)+1,5*(k-1)+1]))
        
      }
    }
  }
  for (j in 1:p){#Dérivées par rapport aux termes de variance
    for (k in 1:p){
      if (j != k) {
        gradCL_calc[d*p+j]<-gradCL_calc[d*p+j]+ (Y[k]^2-(Y[k]+1)*(2*Y[k]+1)*(terms[5*(j-1)+1,5*(k-1)+2]/terms[5*(j-1)+1,5*(k-1)+1])+(Y[k]+1)*(Y[k]+2)*(terms[5*(j-1)+1,5*(k-1)+3]/terms[5*(j-1)+1,5*(k-1)+1]))
      }
    }
  }
  for (l in 1:(0.5*p*(p-1))){
    j=liste_indices[[l]][1]
    k=liste_indices[[l]][2]
    gradCL_calc[d*p+p+l]<-2*(Y[k]*Y[j]-(Y[k]+1)*Y[j]*(terms[5*(j-1)+1,5*(k-1)+2]/terms[5*(j-1)+1,5*(k-1)+1])-Y[k]*(Y[j]+1)*(terms[5*(j-1)+2,5*(k-1)+1]/terms[5*(j-1)+1,5*(k-1)+1]) + (Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]/terms[5*(j-1)+1,5*(k-1)+1])
    
  }
  return(list(hess_mat,gradCL_calc))
}

## Donne les deux termes de la matrice de godambe estimé sur toutes les observations
Estimateurs_esp_corr_st<-function(param,Yobs,Xobs,Oobs){#On va calculer la moyenne empirique sur toutes les observation de la matrice hessienne
  n=nrow(Yobs)
  p=ncol(Yobs)
  d=ncol(Xobs)
  Esp_hess=matrix(0,(p*d+p+0.5*p*(p-1)),(p*d+p+0.5*p*(p-1)))
  Esp_gradCL = matrix(0,(p*d+p+0.5*p*(p-1)),(p*d+p+0.5*p*(p-1)))
  for (k in 1:n){
    objets<-Hessian_corr_st(param,Yobs[k,],Xobs[k,],Oobs[k,],d,p)
    Esp_hess<-Esp_hess+objets[[1]]
    Esp_gradCL<-Esp_gradCL+objets[[2]]%*%t(objets[[2]])
  }
  return(list((1/n)*Esp_hess,(1/n)*Esp_gradCL))
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
