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


####################################################################################################################
## Composite likelihood du modèle spatial

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


var_cov_spatial_dist<-function(mat_dist,alpha,sigma){
  p=length(points_x)
  var=matrix(0,p,p)
  for(i in 1:(p-1)){
    for (j in (i+1):p){
      var[i,j]<-exp(-alpha*mat_dist)
    }
  }
  var<-var+t(var)
  diag(var)<-rep(1,p)
  return(sigma*var)
}


n=100
nb_simu=250
p=5
d=3
O=matrix(1,n,p)
X=matrix(runif(n,-1,1),n,1)
if(d>1){
  for (j in 1:(d-1)){
    X<-cbind(X,runif(n,-1,1))
  }
}
mu<-runif(p*d,0,1)
alpha=5
sigma=3
#On simule les distances
points_x=runif(p,-1,1)
points_y=runif(p,-1,1)

var=matrix(0,p,p)
Var_Cov=var_cov_spatial(points_x,points_y,alpha,sigma)
param<-c(mu,diag(Var_Cov))
for (j in 1:(p-1)){
  for (k in (j+1):p){
    a=Var_Cov[j,k]
    param=c(param,a)
  }
}

Obs_3_temp=Observations_simulees_bis(n,p,X,O,param)
Obs_3=Obs_3_temp$Y
Z=Obs_3_temp$Z
#On trace d'abord la composite likelihood par rapport au parametre sigma :
abscisse=seq(1,5,length=10)
ordonnees=c()
for (k in 1:length(abscisse)){
  vc=var_cov_spatial(points_x,points_y,alpha,abscisse[k])
  param_profil<-c(mu,diag(vc))
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      a=vc[j,k]
      param_profil=c(param_profil,a)
    }
  }
  ordonnees= c(ordonnees,CL_f(param_profil,Obs_3,O,X))
}
plot(abscisse,ordonnees)
title("CL wrt sigma")

#On trace de même pour alpha
abscisse=seq(0.5,10,length=20)
ordonnees=c()
for (k in 1:length(abscisse)){
  vc=var_cov_spatial(points_x,points_y,abscisse[k],sigma)
  print(vc)
  param_profil<-c(mu,diag(vc))
  for (j in 1:(p-1)){
    for (i in (j+1):p){
      a=vc[j,i]
      param_profil=c(param_profil,a)
    }
  }
  print(param_profil)
  ordonnees= c(ordonnees,CL_f(param_profil,Obs_3,O,X))
}
plot(abscisse,ordonnees)
title("CL wrt alpha (alpha=5)")

#####################################################################################################################################
## Optimisation sans gradient 

##################################################################################################################################
##Optimisation par une méthode sans gradients

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
  Obs_3=Observations_simulees_bis(n,p,X,O,param)
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
    mu=param[1:(p*d)]
    alpha=param[p*d+1]
    sigma=param[p*d+2]
    Var_Cov=var_cov_spatial(points_x,points_y,alpha,sigma)
    param_CL<-c(mu,diag(Var_Cov))
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
    param_optimaux<-nloptr(x_init,eval_f=fonction_a_optimiser,lb=lower,ub=upper,eval_grad_f=NULL,opts = opts)
    param_estim=rbind(param_estim,param_optimaux$solution)
}
E=data.frame(param_estim-param)
E_prime=melt(E)
pltt<-ggplot(data=E_prime,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
  theme_bw()

#Optimisation avec optim en récupérant le gradient et la hessienne
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
  Obs_3=Observations_simulees_bis(n,p,X,O,param)
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
    mu=param[1:(p*d)]
    alpha=param[p*d+1]
    sigma=param[p*d+2]
    Var_Cov=var_cov_spatial(points_x,points_y,alpha,sigma)
    param_CL<-c(mu,diag(Var_Cov))
    for (j in 1:(p-1)){
      for (i in (j+1):p){
        a=Var_Cov[j,i]
        param_CL=c(param_CL,a)
      }
    }
    return(neg_CL(param_CL,Obs_3,O,X))
  }
  ##Optimisation avec optim
  param_optimaux<-optim(x_init,fn=fonction_a_optimiser,eval_grad_f=NULL,method = "SANN",hessian = TRUE,opts = opts)
  param_estim=rbind(param_estim,param_optimaux$solution)
  Vp_inf=ginv(param_optimaux$hessian)
  grad=param_optimaux$gradient%*%t(param_optimaux$gradient)
  Godambe=Vp_inf%*%grad%*%Vp_inf
  param_estim_norm<-rbind(param_estim_norm,(param_optimaux$par-param)/sqrt(diag(Godambe)))
}
E=data.frame(param_estim_norm)
E_prime=melt(E)
pltt<-ggplot(data=E_prime,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
  theme_bw()

##Optimisation en utilisant d'abord nloptr qui est rapide, puis optim pour avoir une estimation de la hessienne

#Dans ce cas, on n'a pas besoin de contraintes car la matrice des variances, covariances spatiales sera toutjours défénie positive
#Simulation des paramètres
p=5
d=2
n=50
nb_simu=100
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
  param_optimaux<-optim(param_optimaux_nl$solution,fn=fonction_a_optimiser,method = "SANN",hessian = TRUE)
  param_estim=rbind(param_estim,param_optimaux$solution)
  Vp_inf=ginv(param_optimaux$hessian)
  grad=param_optimaux$gradient%*%t(param_optimaux$gradient)
  Godambe=Vp_inf%*%grad%*%Vp_inf
  param_estim_norm<-rbind(param_estim_norm,(param_optimaux$par-param)/sqrt(diag(Godambe)))
}
E=data.frame(param_estim_norm)
E_prime=melt(E)
pltt<-ggplot(data=E_prime,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
  theme_bw()+ggtitle("Optimisation wrt spatial parameters")
