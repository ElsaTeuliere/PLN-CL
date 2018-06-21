source('/home/teuliere/PLN-Cl/Script_CLM_Offset.R')

##Code des dérivées d'ordre 2 pour la matrice de Godambe

###############################################################################
#Calcul de la matrice pour une observation :
#Cette fonction va directement nous sortir le gradient et la matrice hessienne

Hessian<-function(param,Y,X,O,d,p){ #Ici on rentre les paramètres d et p pour éviter ensuite dans la boucle sur les observations d'avoir à les recalculer
  #mise en forme des paramètre 
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)
  diag(Rho)<-Sigma
  Rho=cov2cor(Rho)
  Xmu=O+X%*%Mu
  
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
          hess_mat[d*p+j,d*p+j]<-hess_mat[d*p+j,d*p+j]+(1/terms[5*(j-1)+1,5*(k-1)+1]^2)*(-(Y[k]+1)*(2*Y[k]+1)^2*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]+(4*Y[k]^2+12*Y[k]+7)*(Y[k]+1)*(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+3]-(Y[k]+3)*(4*Y[k]+6)*(Y[k]+1)*(Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+4]-(Y[k]+1)^2*(2*Y[k]+1)^2*terms[5*(j-1)+1,5*(k-1)+2]^2+2*(Y[k]+1)*(Y[k]+2)*(Y[k]+3)*(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+2]+(Y[k]+1)*(Y[k]+2)*(Y[k]+3)*(Y[k]+4)*terms[5*(j-1)+1,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+5]-(Y[k]+1)^2*(Y[k]+2)^2*terms[5*(j-1)+1,5*(k-1)+3]^2)
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
    hess_mat[(d*p+k),d*p+j]<-(Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+1]^(-2)*((2*Y[j]+1)*((2*Y[k]+1)*(terms[5*(j-1)+2,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1]-terms[5*(j-1)+1,5*(k-1)+2]*terms[5*(j-1)+2,5*(k-1)+1])+(Y[k]+2)*(terms[5*(j-1)+2,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+3]-terms[5*(j-1)+2,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]))+(Y[j]+2)*((2*Y[k]+1)*(terms[5*(j-1)+3,5*(k-1)+1]*terms[5*(j-1)+1,5*(k-1)+2]-terms[5*(j-1)+3,5*(k-1)+2]*terms[5*(j-1)+1,5*(k-1)+1])+(Y[k]+2)*(terms[5*(j-1)+3,5*(k-1)+3]*terms[5*(j-1)+1,5*(k-1)+1]-terms[5*(j-1)+1,5*(k-1)+3]*terms[5*(j-1)+3,5*(k-1)+1])))
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
  hess_mat[d*p+j,(d+1)*p+l]<-((Y[k]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[j]+1)*(2*Y[k]^2+3*Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2]-3*(Y[j]+1)*(Y[k]+2)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+3]+(Y[j]+1)*(Y[k]+2)*(Y[k]+3)*terms[5*(j-1)+2,5*(k-1)+4]+Y[j]*(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]+Y[j]*(Y[k]+2)*(2*Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+3]- Y[j]*(Y[k]+2)*(Y[k]+3)*terms[5*(j-1)+1,5*(k-1)+4])+(Y[k]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]+Y[j]*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]-(Y[k]+1)*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2])*((Y[k]+2)*terms[5*(j-1)+1,5*(k-1)+3]-(2*Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]))
  #les termes sigma_kk sigma_jk
  hess_mat[d*p+k,(d+1)*p+l]<-((Y[j]+1)/terms[5*(j-1)+1,5*(k-1)+1]^2) *(terms[5*(j-1)+1,5*(k-1)+1]*((Y[k]+1)*(2*Y[j]^2+3*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+2]-3*(Y[k]+1)*(Y[j]+2)*(Y[j]+1)*terms[5*(j-1)+3,5*(k-1)+2]+(Y[k]+1)*(Y[j]+2)*(Y[j]+3)*terms[5*(j-1)+4,5*(k-1)+2]+Y[k]*(2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]+Y[k]*(Y[j]+2)*(2*Y[j]+3)*terms[5*(j-1)+3,5*(k-1)+1]- Y[k]*(Y[j]+2)*(Y[j]+3)*terms[5*(j-1)+4,5*(k-1)+1])+(Y[j]*(Y[k]+1)*terms[5*(j-1)+1,5*(k-1)+2]+Y[k]*(Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]-(Y[j]+1)*(Y[k]+1)*terms[5*(j-1)+2,5*(k-1)+2])*((Y[j]+2)*terms[5*(j-1)+3,5*(k-1)+1]-(2*Y[j]+1)*terms[5*(j-1)+2,5*(k-1)+1]))
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


 Estimateurs_esp<-function(param,Yobs,Xobs,Oobs){#On va calculer la moyenne empirique sur toutes les observation de la matrice hessienne
   n=nrow(Yobs)
   p=ncol(Yobs)
   d=ncol(Xobs)
   Esp_hess=matrix(0,(p*d+p+0.5*p*(p-1)),(p*d+p+0.5*p*(p-1)))
   Esp_gradCL = matrix(0,(p*d+p+0.5*p*(p-1)),(p*d+p+0.5*p*(p-1)))
   for (k in 1:n){
     objets<-Hessian(param,Yobs[k,],Xobs[k,],Oobs[k,],d,p)
     Esp_hess<-Esp_hess+objets[[1]]
     Esp_gradCL<-Esp_gradCL+objets[[2]]%*%t(objets[[2]])
   }
   return(list((1/n)*Esp_hess,(1/n)*Esp_gradCL))
 }

 
 symetrisation<-function(M){
   for (k in 2:nrow(M)){
     for (j in 1:(k-1)){
       if (M[k,j]==0){
         M[k,j]=M[j,k]
       }
     }
   }
   return (M)
 }
 
 
 
 
 Hessian_corr<-function(param,Y,X,O,d,p){ #Ici on rentre les paramètres d et p pour éviter ensuite dans la boucle sur les observations d'avoir à les recalculer
   #mise en forme des paramètre 
   Mu=matrix(param[1:(p*d)],d,p)
   Sigma=param[(p*d+1):(p*d+p)]
   Rho=matrix(0,nrow=p,ncol=p)
   Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
   Rho=t(Rho)
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
   return(list((1/n)*E,(1/n)*Esp_gradCL))
 }
##Test
param=c(2,1.5,1,1,-0.5)
p=2
n=20
X=matrix(1,n,1)
O=matrix(1,n,p)

Obs_1=Observations_simulees_bis(n,p,X,O,param)
Hessian(param,Obs_1[1,],X[1,],O[1,],1,2)