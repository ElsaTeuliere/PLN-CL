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


approx_hes_uni<-function(param,Y,O,X,epsilon){
  #Cette fonction permet d'approcher la valeur du gradient  en utilisant un pas epsilon
  approx_hess<-matrix(0,length(param),length(param))
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=Rho+t(Rho)
  Xmu=O+X%*%Mu
  for (j in 1:length(param)){
    param_per<-param
    param_per[j]<-param[j]+epsilon
    Mu_per=matrix(param_per[1:(p*d)],d,p)
    Sigma_per=param_per[(p*d+1):(p*d+p)]
    Rho_per=matrix(0,nrow=p,ncol=p)
    Rho_per[lower.tri(Rho_per,diag=F)]<-param_per[(p*d+p+1):length(param_per)]
    Rho_per=Rho_per+t(Rho_per)
    Xmu_per=O+X%*%Mu_per
    approx_hess[,j]<-(grad_CL_uni(Y,X,O,d,p,Xmu_per,Sigma_per,Rho_per)-grad_CL_uni(Y,X,O,d,p,Xmu_per,Sigma_per,Rho_per))/epsilon
  }
  return(approx_hess)
}


approx_uni_bis<-function(param,Y,O,X,epsilon){
  d=length(X)
  p=length(O)
  approx_hess<-matrix(0,length(param),length(param))
  Mu=matrix(param[1:(p*d)],d,p)
  Sigma=param[(p*d+1):(p*d+p)]
  Rho=matrix(0,nrow=p,ncol=p)
  Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
  Rho=Rho+t(Rho)
  Xmu=O+X%*%Mu
  #On va remplir la matrice comme on fait pour construire la matrice hessienne
  approx_hess<-matrix(0,length(param),length(param))
  #On commence par les termes de moyenne
  termes_derivees<-grad_CL_uni(Y,X,O,d,p,Xmu,Sigma,Rho)
  terme_original<-CL_f_uni(param,Y,O,X,d,p)
  for (j in 1:(d*p)){ #on fixe le paramètre de moyenne
    param_per_1<-param
    param_per_1[j]<-param[j]+epsilon
    approx_hess[j,j]<-(2/epsilon^2)*(CL_f_uni(param_per_1,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[j])
    if(j<(d*p)){
      for (k in (j+1):(d*p)){
        param_per_2<-param_per_1
        param_per_2[k]<-param_per_1[k]+epsilon
        approx_hess[j,k]<-(1/epsilon^2)*(CL_f_uni(param_per_2,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[j]-epsilon*termes_derivees[k])
      }
    }
  }
  #On rempli maintenant les termes de variances
  for (j in 1:p){ #on fixe le paramètre de variance que l'on étudie
    indice=d*p+j
    param_per_1<-param
    param_per_1[indice]<-param[indice]+epsilon
    #dérivée non croisée
    approx_hess[indice,indice]<-(2/epsilon^2)*(CL_f_uni(param_per_1,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[indice])
    #termes de dérivée croisée avec la moyenne
    for (k in 1:(d*p)){
      param_per_2<-param_per_1
      param_per_2[k]<-param_per_1[k]+epsilon
      approx_hess[k,indice]<-(1/epsilon^2)*(CL_f_uni(param_per_2,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[indice]-epsilon*termes_derivees[k])
    }
    #termes de dérivée croisés avec la variance
    if(j<p){
      for (k in (j+1):p){
        indice_prime=d*p+k
        param_per_2<-param_per_1
        param_per_2[indice_prime]<-param_per_1[indice_prime]+epsilon
        approx_hess[indice,indice_prime]<-(1/epsilon^2)*(CL_f_uni(param_per_2,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[indice]-epsilon*termes_derivees[indice_prime])
      }
    }
  }
  #On dérive maintenant par rapport aux termes de covariance
  for (l in 1:(0.5*p*(p-1))){
    indice=d*p+p+l
    param_per_1<-param
    param_per_1[indice]<-param[indice]+epsilon
    approx_hess[indice,indice]<-(2/epsilon^2)*(CL_f_uni(param_per_1,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[indice])
    #termes croisés avec la moyenne
    for (k in 1:(d*p)){
      param_per_2<-param_per_1
      param_per_2[k]<-param_per_1[k]+epsilon
      approx_hess[k,indice]<-(1/epsilon^2)*(CL_f_uni(param_per_2,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[indice]-epsilon*termes_derivees[k])
    }
    for (k in 1:p){
      indice_prime=d*p+k
      param_per_2<-param_per_1
      param_per_2[indice_prime]<-param_per_1[indice_prime]+epsilon
      approx_hess[indice_prime,indice]<-(1/epsilon^2)*(CL_f_uni(param_per_2,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[indice]-epsilon*termes_derivees[indice_prime])
    }
    if(l<(0.5*p*(p-1))){
      for (k in (l+1):(0.5*p*(p-1))){
        indice_prime=d*p+p+k
        param_per_2<-param_per_1
        param_per_2[indice_prime]<-param_per_1[indice_prime]+epsilon
        approx_hess[indice,indice_prime]<-(1/epsilon^2)*(CL_f_uni(param_per_2,Y,O,X,d,p)-terme_original-epsilon*termes_derivees[indice]-epsilon*termes_derivees[indice_prime])
      }
    }
  }
  return(approx_hess)
}

hes_approx_tot<-function(param,Y,O,X,epsilon){
  d=ncol(X)
  p=ncol(Y)
  n=nrow(Y)
  hes<-0
  for (i in 1:n) {
    hes<- hes+approx_uni_bis(param,Y[i,],O[i,],X[i,],epsilon)
  }
  return((1/n)*hes)
}

norme_L2<-function(x){
  return(sqrt(sum(x^2)))
}

#######################################################################################################################
##En calculant Godambe avec les bonnes observations

param_estim=c()
param_estim_norm<-c()
c=0
tps<-0
grad_sortie_PLN<-c()
grad_sortie_nloptr<-c()
for (k in 1:nb_simu){
  Obs_3=Observations_simulees_bis(n,p,X,O,param)
  b_time<-Sys.time()
  x_0=param_0(Obs_3,O,X)
  lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-0.999,0.5*p*(p-1)))
  ub = c(rep(Inf,p*d+p),rep(0.999,0.5*p*(p-1)))
  if (all(x_0<=ub) & all(x_0>=lb)){
    grad_sortie_PLN<-c(grad_sortie_PLN,norme_L2(grad_CL_opt(x_0,Obs_3,O,X)))
    param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL,
                           lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-0.999,0.5*p*(p-1))), ub = c(rep(Inf,p*d+p),rep(0.999,0.5*p*(p-1))),
                           opts=opts, Y=Obs_3, X=X,O=O)
    param_estim=rbind(param_estim,param_optimaux$solution)
    grad_sortie_nloptr<-c(grad_sortie_nloptr,norme_L2(grad_CL_opt(param_optimaux$solution,Obs_3,O,X)))
    c<-c+1
    nombre_iterations[length(nombre_iterations)+1]<-param_optimaux$iter
    V_inf=Estimateurs_esp(param_optimaux$solution,Obs_3,X,O)
    Vp_inf=ginv(symetrisation(V_inf[[1]]))
    Godambe=Vp_inf%*%V_inf[[2]]%*%Vp_inf
    param_estim_norm<-rbind(param_estim_norm,(param_optimaux$solution-param)/sqrt(diag(Godambe)))
    tps<-tps+(Sys.time()-b_time)
  }
}
temps_moyen[[length(temps_moyen)+1]]<-tps/c
name_table=paste(c('values',as.character(p),as.character(d),"_new"), collapse = "")
name_table_norm=paste(c('values',as.character(p),as.character(d),"_norm_new"), collapse = "")
name_graph=paste(c('graph',as.character(p),as.character(d),"_new"), collapse = "")
write.table(param_estim_norm, file = paste(c("/home/teuliere/PLN-Cl/",name_table_norm),collapse=""), append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))
write.table(param_estim, file = paste(c("/home/teuliere/PLN-Cl/",name_table),collapse=""), append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))
E_20=read.table(name_table_norm,sep="\t",header = TRUE)
E_202=melt(E_20)
plt<-ggplot(data=E_202,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
  theme_bw()
png(paste(name_graph,".png",sep=""),width=1000, height=800,res=110)
print(plt)
dev.off()

E=read.table(name_table,sep="\t",header=TRUE)
for (i in 1:nb_simu){
  E[j,]<-E[j,]-param
}
E_prime=melt(E)
pltt<-ggplot(data=E_prime,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
  theme_bw()
print(pltt)


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
epsilon=10^(-6)
hes_approx<-approx_uni_bis(param,Y,O,X,epsilon)
hes<-Hessian(param,Y,X,O,d,p)[[1]]
delta_hes<-matrix(NA,nrow(hes),ncol(hes))
for (i in 1:nrow(hes)){
  for (j in 1:ncol(hes)){
    if (hes[i,j]!=0){
      delta_hes[i,j]<-abs(hes_approx[i,j]-hes[i,j])/hes[i,j]
    }
  }
}


###########
##Peut-être que l'erreur vient du fait qu'on ne rerempli pas la matrice dont on a rempli en gros que le triangle sup.



##################################################################################################################################
##Vérification grâce à l'approximation numérique fournie par optim

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
mu<-runif(p*d,0,2)
Var_Cov<-var_cov(p)
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

  ##Optimisation avec optim
  lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-0.999,0.5*p*(p-1)))
  ub = c(rep(Inf,p*d+p),rep(0.999,0.5*p*(p-1)))
  if (all(x_0<=ub) & all(x_0>=lb)){
    param_optimaux<-optim(x_0,neg_CL,neg_grad_CL,method = "L-BFGS-B",lower=lb,upper=ub,hessian = TRUE, Y=Obs_3,O=O,X=X)
    param_estim=rbind(param_estim,param_optimaux$par)
    V_inf=Estimateurs_esp(param_optimaux$par,Obs_3,X,O)
    Vp_inf=ginv(param_optimaux$hessian)
    Godambe=Vp_inf%*%V_inf[[2]]%*%Vp_inf
    param_estim_norm<-rbind(param_estim_norm,(param_optimaux$par-param)/sqrt(diag(Godambe)))
    diff_hessian=c(diff_hessian, max(abs(symetrisation(V_inf[[1]])-param_optimaux$hessian)))
  }
}
E_prime=melt(data.frame(param_estim_norm))
pltt<-ggplot(data=E_prime,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
  theme_bw()
