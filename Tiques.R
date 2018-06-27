library(sp)
library(stringr)
setwd(dir="/home/teuliere/PLN-Cl")
load("Elsa_data_tranche.Rdata")


## Essai numero 1, en prenant la structure des sites, et en considérant les tranches comme des répétitions
#On a donc 220 variables d'interêt et 16 répétitions.
################################################################################################
##Mise en forme du jeu de données 
################################################################################################
tranche2=as.data.frame(tranche)


#Creation matrice distance

mat_dist <- as.matrix(dist(coordinates(tranche)))
 for (i in 1:(ncol(mat_dist)-1)){
   for(j in (i+1):ncol(mat_dist)){
     if (identical(tranche$N6[i],tranche$N6[j])){}else{
       mat_dist[i,j]<-10^9
       mat_dist[j,i]<-10^9
     }
   }
 }


#Creation de la matrice des observations
Obs_nymphe=tranche$nb_nymphe
#Creation de la matrice des covariables
X=cbind(as.numeric(tranche$cl_sat_def_parcelle),as.numeric(tranche$cl_deer_ibw_ab_parcelle_corr),as.numeric(tranche$cl_boar_ibw_ab_parcelle_corr),as.numeric(tranche$cl_dist_chip_intro_parcelle),as.numeric(tranche$vegetation),as.numeric(tranche$forest_cover),as.numeric(tranche$superficial_soil),as.numeric(tranche$order_sampling))
#Creation de la matrice des Offsets
O=rep(1,3504)

#Paramétrisation p et d
p=nrow(tranche)
d=ncol(X)
#################################################################################################
##Analyse des données
################################################################################################
#On utilise le modèle spatial donc l'optimisation sans gradient

x_0=param_0(Obs_nymphe,O,X)
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
    Var_Cov=var_cov_spatial_dist(mat_dist,alpha_trans,sigma_trans)
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
  param_estim=param_optimaux$solution
  Vp_inf=ginv(param_optimaux$hessian)
  grad=param_optimaux$gradient%*%t(param_optimaux$gradient)
  Godambe=Vp_inf%*%grad%*%Vp_inf
  param_estim_norm<-(param_optimaux$par-param)/sqrt(diag(Godambe))

