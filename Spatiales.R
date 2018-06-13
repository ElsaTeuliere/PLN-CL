library(MASS)
library(PLNmodels)
library(poilog)
library(nloptr)
library(reshape2)
library(ggplot2)

#####################################################################################################
##L'idée de ce script est juste d'utiliser ce que nous avons déjà fait mais au cas plus simple où nos paramètres sont sous forme 
##sigma_jk= r exp(- alpha d_jk).

#Dans ce cas là l'idée est d'optimiser en fonction des paramètres de régression et en fonction des r et alpha.
#On ajoute donc un argument à la fonction : la matrice D qui comporte les distances entre deux parcelles.

#Ici param sera de la forme c(mu11...mu1d,mu21.....mupd,r,alpha)
#La matrice des distances D aura la forme p p .
###################################################################################################################


#La fonction à optimiser qui prend directement sur l'ensemble des n observations
CL_f_sp<-function(param,D,Y,O,X){
  n=nrow(Y)
  p=ncol(Y)
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Rho=param[p*d+1]*exp(- param[p*d+1]*D) #matrice de variance/covariance
  Sigma=diag(Rho) 
  Xmu=O+X%*%Mu
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

#Composite Likelihood pour une observation
CL_f_uni_sp<-function(param,Y,O,X,d,p,Xmu,Sigma,Rho){
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covariance
  Xmu=O+X%*%Mu
  CL<-0
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      CL<-CL+log(dbipoilog(Y[j],Y[k],mu1=Xmu[j],mu2=Xmu[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k]))
    }
  }
  return(CL)
}

# negative likelihood 
neg_CL_sp<-function(param,D,Y,O,X){
  return(-CL_f_sp(param,D,Y,O,X))
}

#############################################################################################################################
# Gradient

#Pour une observation. On sommera ensuite sur toutes les observations

grad_CL_uni_sp<-function(Y,X,O,d,p,Xmu,Sigma,Rho){
  #On commence par calculer tous les termes dont on a besoin 
  terms<-matrix(0,3*p,3*p)
  for(j in 2:p){
    for (k in 1:(j-1)){
      for (l in 1:3){
        for (m in 1:3){
          if (l+m<5){
            terms[(3*(k-1)+l),(3*(j-1)+m)]<- dbipoilog(Y[k]+l-1,Y[j]+m-1,Xmu[k],Xmu[j],Sigma[k],Sigma[j],Rho[k,j])
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


grad_CL_opt_sp<-function(param,D,Y,O,X){
  n=nrow(Y)
  p=ncol(Y)
  d=ncol(X)
  Mu=matrix(param[1:(p*d)],d,p)
  Xmu=O+X%*%Mu
  Rho=param[p*d+1]*exp(- param[p*d+1]*D) #matrice de variance/covariance
  Sigma=diag(Rho) 
  grad_CL=rep(0,length(param))
  for (i in 1:n){#On répète sur les observations
    grad_CL<-grad_CL+grad_CL_uni(Y[i,],X[i,],O[i,],d,p,Xmu[i,],Sigma,Rho)
  }
  return(grad_CL)
}


neg_grad_CL_opt_sp<-function(param,D,Y,O,X){
  return(-grad_CL_opt_sp(param,D,Y,O,X))
}

#################################################################################################
##Matrice hessienne

Hessian_corr_sp<-function(param,D,Y,X,O,d,p){ #Ici on rentre les paramètres d et p pour éviter ensuite dans la boucle sur les observations d'avoir à les recalculer
  #mise en forme des paramètre 
  Mu=matrix(param[1:(p*d)],d,p)
  Xmu=O+X%*%Mu
  Rho=param[p*d+1]*exp(- param[p*d+1]*D) #matrice de variance/covariance
  Sigma=diag(Rho) 
  #On commence par calculer tous les termes dont on aura besoin 
  terms<-matrix(0,5*p,5*p)
  for(j in 2:p){
    for (k in 1:(j-1)){
      for (l in 1:5){
        for (m in 1:5){
          if (l+m<7){
            terms[(5*(k-1)+l),(5*(j-1)+m)]<- dbipoilog(Y[k]+l-1,Y[j]+m-1,Xmu[k],Xmu[j],Sigma[k],Sigma[j],Rho[k,j])
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
        hess_mat[d*p+j,d*p+j]<-hess_mat[d*p+j,d*p+j]+(1/terms[5*(j-1)+1,5*(k-1)+1]^2)*(-(Y[k]+1)*(2*Y[k]+1)^2*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]+2*(2*Y[k]+1)*(Y[k]+1)^2*(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+3]-2*(Y[k]+1)*(Y[k]+2)*(Y[k]+3)*(2*Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+4]-(Y[k]+1)^2*(2*Y[k]+1)^2*terms[5*(j-1)+1,5*(k-1)+2]^2+2*(Y[k]+1)^2*(Y[k]+2)*(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+3]+(Y[k]+1)*(Y[k]+2)*(Y[k]+3)*(Y[k]+4)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+5]-(Y[k]+1)^2*(Y[k]+2)^2*terms[5*(j-1)+1,5*(k-1)+3]^2)
      }
    }
  }
  #Maintenant les termes de dérivée croisée avec les termes de moyenne
  #On commence par les dérivées mu_j sigma_jj
  for(j in 1:p){
    vect=0
    for (k in 1:p){
      if(j !=k){
        vect<- vect + (Y[j]+1)*(Y[k]+2)*(1/terms[5*(j-1)+1,5*(k-1)+1]^2)*((2*Y[k]+1)*(terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+2,5*(k-1)+2]-terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2])+(Y[k]+2)*(terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+3]-terms[5*(j-1)+2,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]))
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
      hess_mat[(d*p+k),d*p+j]<-(Y[j]+1)*(Y[k]+1)*((2*Y[j]+1)*((2*Y[k]+1)*(terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1]-terms[5*(j-1)+1,5*(k-1)+2]*terms[5*(j-1)+2,5*(k-1)+1])+(Y[k]+2)*(terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+3]-terms[5*(j-1)+2,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]))+(Y[j]+2)*((2*Y[k]+1)*(terms[5*(j-1)+3,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]-terms[5*(j-1)+3,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1])+(Y[k]+2)*(terms[5*(j-1)+3,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]-terms[5*(j-1)+1,5*(k-1)+3]*terms[5*(j-1)+3,5*(k-1)+1])))
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
    hess_mat[d*p+p+l,d*p+p+l]<- (1/terms[5*(j-1)+1,5*(k-1)+1]^2)*(-(Y[j]+1)*Y[k]^2*terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+1]+(Y[j]+1)*(Y[j]+2)*Y[k]^2*terms[5*(j-1)+3,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+1]+(Y[j]+1)*(Y[k]+1) * (2*Y[k]*Y[j]+2*Y[k]+2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1]-(Y[j]+1)*(Y[j]+2)*(Y[k]+1)*(2*Y[k]+1)*terms[5*(j-1)+3,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1]-(Y[j]+1)^2*Y[k]^2*terms[5*(j-1)+2,5*(k-1)+1]^2-2*(Y[j]+1)*(Y[k]+1)*Y[j]*Y[k]*terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]+2* (Y[j]+1)^2*(Y[k]+1)*Y[k]*terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+2,5*(k-1)+1]-Y[j]^2*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]+Y[j]^2*(Y[k]+1)*(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]-(Y[k]+1)*(Y[k]+2)*(Y[j]+1)*(2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]-Y[j]^2*(Y[k]+1)^2*terms[5*(j-1)+1,5*(k-1)+2]^2+2*(Y[k]+1)^2*Y[j]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+2]+(Y[j]+1)*(Y[k]+1)*(Y[j]+2)*(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+3,5*(k-1)+3]-(Y[j]+1)^2*(Y[k]+1)^2*terms[5*(j-1)+2,5*(k-1)+2]^2)
  }
  #On s'attaque aux termes croisés avec la moyenne
  for (l in 1 : (0.5*p*(p-1))){
    j=liste_indices[[l]][1]
    k=liste_indices[[l]][2]
    #On remplit tout d'abord les blocs mu_j sigma_jk
    hess_mat[((j-1)*d+1):(j*d),d*p+p+l]<-X*((Y[j]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]-Y[k]*terms[5*(j-1)+2,5*(k-1)+1]+(Y[j]+2)*Y[k]*terms[5*(j-1)+3,5*(k-1)+1]-(Y[j]+2)*(Y[k]+1)*terms[5*(j-1)+3,5*(k-1)+2]) +terms[5*(j-1)+2,5*(k-1)+1] *(-Y[j]*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]-(Y[j]+1)*Y[k]*terms[5*(j-1)+2,5*(k-1)+1]+(Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]) )
    #On remplit maintenant les mu_k sigma_jk
    hess_mat[((k-1)*d+1):(k*d),d*p+p+l] <- X*((Y[j]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[k]+1)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]-Y[j]*terms[5*(j-1)+1,5*(k-1)+2]+(Y[k]+2)*Y[j]*terms[5*(j-1)+1,5*(k-1)+3]-(Y[k]+2)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+3]) +terms[5*(j-1)+1,5*(k-1)+2] *(-Y[k]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]-(Y[k]+1)*Y[j]*terms[5*(j-1)+1,5*(k-1)+2]+(Y[k]+1)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]) )
    #les termes sigma_jj sigma_jk
    hess_mat[d*p+j,(d+1)*p+l]<-((Y[k]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[j]+1)*(2*Y[k]^2+3*Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]-3*(Y[j]+1)*(Y[k]+2)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+3]+(Y[j]+1)*(Y[k]+2)*(Y[k]+3)*terms[5*(j-1)+2,5*(k-1)+4]-Y[j]*(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]+Y[j]*(Y[k]+2)*(2*Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+3]- Y[j]*(Y[k]+2)*(Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+4])+(Y[k]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]+Y[j]*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]-(Y[k]+1)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2])*((Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+3]-(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]))
    #les termes sigma_kk sigma_jk
    hess_mat[d*p+k,(d+1)*p+l]<-((Y[j]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[k]+1)*(2*Y[j]^2+3*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]-3*(Y[k]+1)*(Y[j]+2)*(Y[j]+1)*terms[5*(j-1)+3,5*(k-1)+2]+(Y[k]+1)*(Y[j]+2)*(Y[j]+3)*terms[5*(j-1)+4,5*(k-1)+2]-Y[k]*(2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]+Y[k]*(Y[j]+2)*(2*Y[j]+3)*terms[5*(j-1)+3,5*(k-1)+1]- Y[k]*(Y[j]+2)*(Y[j]+3)*terms[5*(j-1)+4,5*(k-1)+1])+(Y[j]*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]+Y[k]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]-(Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2])*((Y[j]+2)*terms[5*(j-1)+3,5*(k-1)+1]-(2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]))
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
    gradCL_calc[d*p+p+l]<-Y[k]*Y[j]-(Y[k]+1)*Y[j]*(terms[5*(j-1)+1,5*(k-1)+2]/terms[5*(j-1)+1,5*(k-1)+1])-Y[k]*(Y[j]+1)*(terms[5*(j-1)+2,5*(k-1)+1]/terms[5*(j-1)+1,5*(k-1)+1]) + (Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]/terms[5*(j-1)+1,5*(k-1)+1]
    
  }
  return(list(hess_mat,gradCL_calc))
}

Estimateurs_esp_corr<-function(param,Yobs,Xobs,Oobs){#On va calculer la moyenne empirique sur toutes les observation de la matrice hessienne
  n=nrow(Yobs)
  p=ncol(Yobs)
  d=ncol(Xobs)
  Esp_hess=matrix(0,(p*d+p+0.5*p*(p-1)),(p*d+p+0.5*p*(p-1)))
  Esp_gradCL = matrix(0,(p*d+p+0.5*p*(p-1)),(p*d+p+0.5*p*(p-1)))
  for (k in 1:n){
    objets<-Hessian_corr(param,Yobs[k,],Xobs[k,],Oobs[k,],d,p)
    Esp_hess<-Esp_hess+objets[[1]]
    Esp_gradCL<-Esp_gradCL+objets[[2]]%*%t(objets[[2]])
  }
  return(list((1/n)*Esp_hess,(1/n)*Esp_gradCL))
}