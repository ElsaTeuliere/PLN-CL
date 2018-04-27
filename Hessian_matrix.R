##Code des dérivées d'ordre 2 pour la matrice de Godambe

###############################################################################
#Calcul de la matrice pour une observation :

Hessian<-function(param,Yi,X,O,d,p){ #Ici on rentre les paramètres d et p pour éviter ensuite dans la boucle sur les observations d'avoir à les recalculer
  #mise en forme des paramètre 
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=t(Rho)
  Xmu=O+X%*%Mu
  
  #On commence par calculer tous les termes dont on aura besoin 
  terms<-matrix(4*p,4*p)
  for(j in 2:p){
    for (k in 1:(j-1){
      for (l in 1:5){
        for (m in 1:5){
          if (l+m<6){
            terms[(5*(k-1)+l),(5*(j-1)+m)]<- dbipoilog(Y[k]+l-1,Y[j]+l-1,Xmu[k],Xmu[j],Sigma[k],Sigma[j],Rho[k,j])
          }
        }
      }
    }
  }
#Maintenant on rempli la matrice hessienne
#On commence par le bloc des ddérivées par rapport aux termes de moyenne
  

}