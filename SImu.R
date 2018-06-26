#Script pour effectuer des simulations avec la variance infinie calculées avec notre fonction Hessian

source('/home/teuliere/PLN-Cl/Script_CLM_Offset.R')
source('/home/teuliere/PLN-Cl/Hessian_matrix.R')

# #######################################################################################################
# ##Simulation de données pour voir la distribution de l'estimateur.
# 
# nb_simu = 20
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
#                          lb = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,1.0e-4,1.0e-4,1.0e-4,-0.999,-0.999,-0.999), ub = c(rep(Inf,9),0.999,0.999,0.999),
#                          opts=opts, Y=Obs_3, X=X,O=O)
#   param_estim=rbind(param_estim,param_optimaux$solution)
# }
# param_estim_norm<-param_estim
# for (j in 1:nrow(param_estim)){#Ici on calcule la matrice de Godambe pour chaque jeu de paramètres estimés
#   V_inf=Estimateurs_esp(param_estim[j,],Obs_3,X,O)
#   Vp_inf=ginv(V_inf[[1]])
#   Godambe=Vp_inf%*%V_inf[[2]]%*%Vp_inf
#   param_estim_norm[j,]<-(param_estim[j,]-param)/diag(Godambe)
# }
# write.table(param_estim_norm, file = "/home/teuliere/PLN-Cl/Estim_T_20", append = FALSE, quote = TRUE, sep = "\t",
#             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
#             col.names = TRUE, qmethod = c("escape", "double"))
# 
# setwd(dir="/home/teuliere/PLN-Cl")
# E_20=read.table("Estim_T_20",sep="\t",header = TRUE)
# E_202=melt(E_20)
# ggplot(data=E_202,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
# theme_bw()

##################################################################################################################################
##Simulation de grande ampleur
setwd(dir="/home/teuliere/PLN-Cl")
n=100
nb_simu=250
p_set=c(5,7,10,15,30,50)
d_set=1:3
temps_moyen<-list()
nombre_iterations<-list()
for (p in p_set){
  O=matrix(1,n,p)
  for (d in d_set){
    X=matrix(runif(n,-1,1),n,1)
    if(d>1){
      for (j in 1:(d-1)){
        X<-cbind(X,runif(n,-1,1))
      }
    }
    mu<-runif(p*d,0,1)
    Var_Cov<-var_cov(p)
    param<-c(mu,diag(Var_Cov))
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        a=Var_Cov[j,k]
        param=c(param,a)
      }
    }
    ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
                 xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
    opts <- list( "algorithm" = "NLOPT_LD_MMA",
                  "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
                  "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
                  "print_level" = max(0,ctrl$trace-1))
    param_estim=c()
    c=0
    tps<-0
    for (k in 1:nb_simu){
      Obs_3=Observations_simulees_bis(n,p,X,O,param)
      b_time<-Sys.time()
      x_0=param_0(Obs_3,O,X)
      lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-0.999,0.5*p*(p-1)))
      ub = c(rep(Inf,p*d+p),rep(0.999,0.5*p*(p-1)))
      if (all(x_0<=ub) & all(x_0>=lb)){
        param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL,
                              lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-0.999,0.5*p*(p-1))), ub = c(rep(Inf,p*d+p),rep(0.999,0.5*p*(p-1))),
                              opts=opts, Y=Obs_3, X=X,O=O)
        param_estim=rbind(param_estim,param_optimaux$solution)
        c<-c+1
        tps<-tps+(Sys.time()-b_time)
        nombre_iterations[length(nombre_iterations)+1]<-param_optimaux$iter
      }
    }
    temps_moyen[[length(temps_moyen)+1]]<-tps/c
    param_estim_norm<-param_estim
    for (j in 1:nrow(param_estim)){#Ici on calcule la matrice de Godambe pour chaque jeu de paramètres estimés
      V_inf=Estimateurs_esp(param_estim[j,],Obs_3,X,O)
      Vp_inf=ginv(symetrisation(V_inf[[1]]))
      Godambe=Vp_inf%*%V_inf[[2]]%*%Vp_inf
      param_estim_norm[j,]<-(param_estim[j,]-param)/sqrt(diag(Godambe))
    }
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
  }
}

E_bis=read.table(name_table,sep="\t",header=TRUE)
E_bis_norm=E_bis
for (j in 1:nb_simu){
  param_estim=E_bis[j,]
  V_inf_approx=symetrisation(hes_approx_tot(param,Obs_3,O,X,epsilon=10^(-6)))
  Vp_inf_approx=ginv(V_inf_approx)
  V_inf=Estimateurs_esp(param_estim,Obs_3,X,O)
  Godambe=Vp_inf_approx%*%V_inf[[2]]%*%Vp_inf_approx
  E_bis_norm[j,]<- (E_bis[j,]-param)/sqrt(diag(Godambe))
}


#####################################################################################################
##En calculant Godambe avec les bonnes observations
######################################################################################################

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
Var_Cov<-var_cov(p)
param<-c(mu,diag(Var_Cov))
for (j in 1:(p-1)){
  for (k in (j+1):p){
    a=Var_Cov[j,k]
    param=c(param,a)
  }
}
temps_moyen<-list()
nombre_iterations<-list()
param_estim=c()
param_estim_norm<-c()
c=0
tps<-0
grad_sortie_PLN<-c()
grad_sortie_nloptr<-c()
for (k in 1:nb_simu){
  Obs_3=Observations_simulees_bis(n,p,X,O,param)$Y
  b_time<-Sys.time()
  x_0=param_0(Obs_3,O,X)
  lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-Inf,0.5*p*(p-1)))
  ub = c(rep(Inf,p*d+p),rep(+Inf,0.5*p*(p-1)))
    param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL,
                           lb = lb, ub = ub,
                           opts=list("xtol_rel"=1*10^(-6),"algorithm"="NLOPT_LD_TNEWTON"),
                           Y=Obs_3, X=X,O=O)
    param_estim=rbind(param_estim,param_optimaux$par)
    c<-c+1
    nombre_iterations[length(nombre_iterations)+1]<-param_optimaux$iter
    V_inf=Estimateurs_esp_corr(param_optimaux$par,Obs_3,X,O)
    hess=symetrisation(V_inf[[1]])
    A=diag(hess)
    hess_ok=hess+t(hess)-diag(A)
    Vp_inf=ginv(hess_ok)
    Godambe=Vp_inf%*%V_inf[[2]]%*%Vp_inf
    param_estim_norm<-rbind(param_estim_norm,(sqrt(n)*(param_optimaux$par-param))/sqrt(diag(Godambe)))
    tps<-tps+(Sys.time()-b_time)
}
temps_moyen[[length(temps_moyen)+1]]<-tps/c
name_table=paste(c('values',as.character(p),as.character(d),"_new"), collapse = "")
name_table_norm=paste(c('values',as.character(p),as.character(d),"_norm_corr"), collapse = "")
name_graph=paste(c('graph',as.character(p),as.character(d),"_corr"), collapse = "")
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







###########################################################################################""
#Mesurer le temps d'optimisation
##########################################################################################

p=10
d=3
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
ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
opts <- list( "algorithm" = "NLOPT_LD_MMA",
              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
              "print_level" = max(0,ctrl$trace-1))
param_estim=c()

Obs_3=Observations_simulees_bis(n,p,X,O,param)
x_0=param_0(Obs_3,O,X)

#Avec la fonction gradient optimisée
bo_time<-Sys.time()
param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL,
                       lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-0.999,0.5*p*(p-1))), ub = c(rep(Inf,p*d+p),rep(0.999,0.5*p*(p-1))),
                       opts=opts, Y=Obs_3, X=X,O=O)
param_estim=rbind(param_estim,param_optimaux$solution)
eo_time<-Sys.time()

#Avec la fonction gradient non optimisée
b_time<-Sys.time()
param_optimaux_n<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL_n,
                       lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-0.999,0.5*p*(p-1))), ub = c(rep(Inf,p*d+p),rep(0.999,0.5*p*(p-1))),
                       opts=opts, Y=Obs_3, X=X,O=O)
param_estim=rbind(param_estim,param_optimaux_n$solution)
e_time<-Sys.time()

print("fonction optimisée :")
eo_time-bo_time

print("fonction non optimisée")
e_time-b_time
