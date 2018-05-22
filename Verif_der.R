#########################################################################################
## Vérification des formules de dérivées
########################################################################################

source('/home/teuliere/PLN-Cl/Script_CLM_Offset.R')
source('/home/teuliere/PLN-Cl/Hessian_matrix.R')

approx_der<-function(param,Y,O,X,epsilon){
  #Cette fonction permet d'approcher la valeur du gradient  en utilisant un pas epsilon
  approx_grad<-c()
  for (j in 1:length(param)){
    param_per<-param
    param_per[j]<-param[j]+epsilon
    approx_grad<-c(approx_grad,(CL_f(param_per,Obs,O,X)-CL_f(param,Obs,O,X))/epsilon)
  }
  return(approx_grad)
}

approx_hes<-function(param,Y,O,X,epsilon){
  #Cette fonction permet d'approcher la valeur du gradient  en utilisant un pas epsilon
  approx_hess<-matrix(0,length(param),length(param))
  for (j in 1:length(param)){
      param_per<-param
      param_per[j]<-param[j]+epsilon
      approx_hess[,j]<-(grad_CL_opt(param_per,Obs,O,X)-grad_CL_opt(param,Obs,O,X))/epsilon
    }
  return(approx_hess)
}


norme_L2<-function(x){
  return(sqrt(sum(x^2)))
}
##############################################################################################
#Script de test

p=5
d=3
n=100
O=matrix(1,n,p)
X=matrix(runif(n,-1,1),n,1)
if(d>1){
  for (j in 1:(d-1)){
    X<-cbind(X,runif(n,-1,1))
  }
}
epsilon_set<-c(10^(-3),10^(-4),10^(-5),10^(-6))
liste_epsilon<-list()
liste_param<-list()
liste_indice_param<-list()
liste_Obs<-list()
liste_grad_approx<-list()
liste_grad<-list()
liste_delta_grad<-list()
c=1
for(epsilon in epsilon_set){
  for(k in 1:50){
    mu<-runif(p*d,0,3)
    Var_Cov<-var_cov(p)
    param<-c(mu,diag(Var_Cov))
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        a=Var_Cov[j,k]
        param=c(param,a)
      }
    }
    for (j in 1:10){
      Obs<-Observations_simulees_bis(n,p,X,O,param)
      grad_approx<-approx_der(param,Obs,O,X,epsilon)
      grad<-grad_CL_opt(param,Obs,O,X)
      liste_epsilon[[c]]<-epsilon
      liste_param[[c]]<-param
      liste_indice_param[[c]]<-k
      liste_Obs[[c]]<-Obs
      liste_grad_approx[[c]]<-grad_approx
      liste_grad[[c]]<-grad
      liste_delta_grad[[c]]<- norme_L2(grad-grad_approx)/norme_L2(grad)
      c<-c+1
    }
  }
  tableau<-cbind(liste_epsilon,liste_indice_param,liste_delta_grad)
  data<-data.frame(tableau)
  plt<-ggplot(data=data,aes(x=liste_indice_param,y=liste_delta_grad))+geom_boxplot()+facet_wrap(~liste_epsilon,scales="free")+
    theme_bw()
}

############################################################################################
#Tentons de trouver les erreurs dans la fonction gradient ----> RESOLU

#Testons d'abord la création de la matrice terms
param<-c( 0.1430776, 0.4540851, 1.8289825, 0.6347092, 0.8830140, 2.8492039, 1.8304776, 2.4976766, 0.3186815,0.3693200, 0.5114818, 1.4293980, 1.8165842, 2.0046715, 2.2155389, 1.0000000, 1.0000000, 1.0000000,1.0000000, 1.0000000, 0.2642435, 0.2950962, 0.5789222, 0.3901704, 0.6662627, 0.2586593, 0.6323566,0.2452045, 0.7436866, 0.3291081)
n=100
p=5
d=3
Mu=matrix(param[1:(p*d)],d,p)
Sigma=param[(p*d+1):(p*d+p)]
Rho=matrix(0,nrow=p,ncol=p)
Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
Rho=Rho+t(Rho)
Xmu=O+X%*%Mu

Obs=Observations_simulees_bis(n,p,X,O,param)

Y=Obs[1,]
Xmu1=Xmu[1,]

terms<-matrix(0,3*p,3*p)
for(j in 2:p){
  for (k in 1:(j-1)){
    for (l in 1:3){
      for (m in 1:3){
        if (l+m<5){
          terms[(3*(k-1)+l),(3*(j-1)+m)]<- dbipoilog(Y[k]+l-1,Y[j]+m-1,Xmu1[k],Xmu1[j],Sigma[k],Sigma[j],Rho[k,j])
        }
      }
    }
  }
}
terms<-terms+t(terms) 

j=1
k=2
dbipoilog(Y[j],Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
dbipoilog(Y[j]+1,Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
dbipoilog(Y[j],Y[k]+1,mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
dbipoilog(Y[j]+1,Y[k]+1,mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])

j=2
k=5
dbipoilog(Y[j],Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
dbipoilog(Y[j]+1,Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
dbipoilog(Y[j],Y[k]+1,mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
dbipoilog(Y[j]+1,Y[k]+1,mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])

##La matrice terms est bonne

##Calcul d'une étape du gradient
j=1
k=5
pjki=dbipoilog(Y[j],Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
pjki1=dbipoilog(Y[j]+1,Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
(X[1,]/pjki)*(Y[j]*pjki - (Y[j]+1)*pjki1 )
(X[1,]/terms[3*(j-1)+1,3*(k-1)+1])*(Y[j]*terms[3*(j-1)+1,3*(k-1)+1]-(Y[j]+1)*terms[3*(j-1)+2,3*(k-1)+1])

##OK

## Calcul de la somme sur k
j=4
gradCL_calc<-rep(0,d*p+p+0.5*p*(p-1))
for (k in 1:p){
  if (j != k){
    gradCL_calc[((j-1)*d+1):(d*j)]<-gradCL_calc[((j-1)*d+1):(d*j)]+(X[1,]/terms[3*(j-1)+1,3*(k-1)+1])*(Y[j]*terms[3*(j-1)+1,3*(k-1)+1]-(Y[j]+1)*terms[3*(j-1)+2,3*(k-1)+1])
  } 
}

grad_CL=rep(0,length(param))
for (k in 1:p){
  if (k != j){ 
    pjki=dbipoilog(Y[j],Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
    pjki1=dbipoilog(Y[j]+1,Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
    grad_CL[((j-1)*d+1):(j*d)]=grad_CL[((j-1)*d+1):(j*d)]+(X[1,]/pjki)*(Y[j]*pjki - (Y[j]+1)*pjki1 )
  }
}

gradCL_calc
grad_CL

##OK

# Boucle sur j
gradCL_calc<-rep(0,d*p+p+0.5*p*(p-1))
for (j in 1:p){
  for (k in 1:p){
    if (j != k){
      gradCL_calc[((j-1)*d+1):(d*j)]<-gradCL_calc[((j-1)*d+1):(d*j)]+(X[1,]/terms[3*(j-1)+1,3*(k-1)+1])*(Y[j]*terms[3*(j-1)+1,3*(k-1)+1]-(Y[j]+1)*terms[3*(j-1)+2,3*(k-1)+1])
    }
  } 
}

grad_CL=rep(0,length(param))
for (j in 1:p){
  for (k in 1:p){
    if (k != j){ 
      pjki=dbipoilog(Y[j],Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
      pjki1=dbipoilog(Y[j]+1,Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
      grad_CL[((j-1)*d+1):(j*d)]=grad_CL[((j-1)*d+1):(j*d)]+(X[1,]/pjki)*(Y[j]*pjki - (Y[j]+1)*pjki1 )
    }
  }
}

gradCL_calc
grad_CL

##OK

##Boucle sur n
gradCL_calc<-rep(0,d*p+p+0.5*p*(p-1))
for (i in 1:n){
  Y=Obs[i,]
  Xmu1=Xmu[i,]
  gradCL_calc<-gradCL_calc+grad_CL_uni(Y,X[i,],O[i],d,p,Xmu1,Sigma,Rho)
}

grad_CL=rep(0,length(param))
for (j in 1:p){
  for (k in 1:p){
    if (k != j){ 
      for (i in 1:n){
        Y=Obs[i,]
        Xmu1=Xmu[i,]
        pjki=dbipoilog(Y[j],Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
        pjki1=dbipoilog(Y[j]+1,Y[k],mu1=Xmu1[j],mu2=Xmu1[k],sig1 = Sigma[j],sig2 = Sigma[k],rho = Rho[j,k])
        grad_CL[((j-1)*d+1):(j*d)]=grad_CL[((j-1)*d+1):(j*d)]+(X[i,]/pjki)*(Y[j]*pjki - (Y[j]+1)*pjki1 )
      }
    }
  }
}


gradCL_calc
grad_CL
####
grad_CLcalc=rep(0,length(param))
for (i in 1:n){#On répète sur les observations
  grad_CLcalc<-grad_CLcalc+grad_CL_uni(Obs[i,],X[i,],O[i,],d,p,Xmu[i,],Sigma,Rho)
}

grad_CLcalc
grad_CL


###
grad_CLcalc<-grad_CL_opt(param,Obs,O,X)

grad_CLcalc
grad_CL


#################################################################################################################
##Vérification des dérivées secondes
p=5
d=3
n=100
O=matrix(1,n,p)
X=matrix(runif(n,-1,1),n,1)
if(d>1){
  for (j in 1:(d-1)){
    X<-cbind(X,runif(n,-1,1))
  }
}
param<-c( 0.1430776, 0.4540851, 1.8289825, 0.6347092, 0.8830140, 2.8492039, 1.8304776, 2.4976766, 0.3186815,0.3693200, 0.5114818, 1.4293980, 1.8165842, 2.0046715, 2.2155389, 1.0000000, 1.0000000, 1.0000000,1.0000000, 1.0000000, 0.2642435, 0.2950962, 0.5789222, 0.3901704, 0.6662627, 0.2586593, 0.6323566,0.2452045, 0.7436866, 0.3291081)
hes_approx<-approx_hes(param,Obs,O,X,epsilon)
hes<-n*Estimateurs_esp(param,Obs,X,O)[[1]]
epsilon=10^(-6)

##Les dériées secondes ont l'air OK too


