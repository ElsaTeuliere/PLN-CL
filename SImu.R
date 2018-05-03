#Script pour effectuer des simulations avec la variance infinie calculées avec notre fonction Hessian

source('/home/teuliere/PLN-Cl/Script_CLM_Offset.R')
source('/home/teuliere/PLN-Cl/Hessian_matrix.R')

#######################################################################################################
##Simulation de données pour voir la distribution de l'estimateur.

nb_simu = 20
n=30
p=3
d=2
X1=runif(n,-1,1)
X2=runif(n,0.5,1.3)
X=cbind(X1,X2)
param=c(0.3,0.4,-1.2,0.1,1.4,-0.02,0.2,0.2,0.3,-0.1,0,-0.1)
O=matrix(1,n,p)
ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
             xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
opts <- list( "algorithm" = "NLOPT_LD_MMA",
              "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
              "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
              "print_level" = max(0,ctrl$trace-1))


param_estim=c()
for (k in 1:nb_simu){
  Obs_3=Observations_simulees_bis(n,p,X,O,param)
  x_0=param_0(Obs_3,O,X)
  param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL,
                         lb = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,1.0e-4,1.0e-4,1.0e-4,-0.999,-0.999,-0.999), ub = c(rep(Inf,9),0.999,0.999,0.999),
                         opts=opts, Y=Obs_3, X=X,O=O)
  param_estim=rbind(param_estim,param_optimaux$solution)
}
param_estim_norm<-param_estim
for (j in 1:nrow(param_estim)){#Ici on calcule la matrice de Godambe pour chaque jeu de paramètres estimés
  V_inf=Estimateurs_esp(param_estim[j,],Obs_3,X,O)
  Vp_inf=ginv(V_inf[[1]])
  Godambe=Vp_inf%*%V_inf[[2]]%*%Vp_inf
  param_estim_norm[j,]<-(param_estim[j,]-param)/diag(Godambe)
}
write.table(param_estim_norm, file = "/home/teuliere/PLN-Cl/Estim_T_20", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))

setwd(dir="/home/teuliere/PLN-Cl")
E_20=read.table("Estim_T_20",sep="\t",header = TRUE)
E_202=melt(E_20)
ggplot(data=E_202,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
theme_bw()

##################################################################################################################################
##Simulation de grande ampleur
setwd(dir="/home/teuliere/PLN-Cl")
n=100
nb_simu=250
p_set=1:10
d_set=1:3
for (p in p_set){
  O=matrix(1,n,p)
  for (d in d_set){
    X=runif(n,-1,1)
    for (j in 2:d){
      X<-cbind(X,runif(n,-1,1))
    }
    param<-c(round(runif(p*d,0,3)),runif(p,0,1),runif(0.5*p*(p-1),-1,1))
    Sigma=param[(p*d+1):(p*d+p)]
    Rho=matrix(0,nrow=p,ncol=p)
    Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
    Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covariance
    while (min(eigen(Rho,symmetric = TRUE)$values) <0){
      param<-c(round(runif(p*d,0,3)),runif(p,0,1),runif(0.5*p*(p-1),-1,1))
      Sigma=param[(p*d+1):(p*d+p)]
      Rho=matrix(0,nrow=p,ncol=p)
      Rho[lower.tri(Rho,diag=F)]<-param[(p*d+p+1):length(param)]
      Rho=t(Rho)+Rho+diag(Sigma,p,p) #matrice de variance/covarianc
    }
    ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
                 xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
    opts <- list( "algorithm" = "NLOPT_LD_MMA",
                  "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
                  "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
                  "print_level" = max(0,ctrl$trace-1))
    param_estim=c()
    for (k in 1:nb_simu){
      Obs_3=Observations_simulees_bis(n,p,X,O,param)
      x_0=param_0(Obs_3,O,X)
      param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL, eval_grad_f=neg_grad_CL,
                             lb = c(rep(-Inf,p*d),rep(1.0e-4,p),rep(-0.999,0.5*p*(p-1))), ub = c(rep(Inf,p*d+p),rep(0.999,0.5*p*(p-1))),
                             opts=opts, Y=Obs_3, X=X,O=O)
      param_estim=rbind(param_estim,param_optimaux$solution)
    }
    param_estim_norm<-param_estim
    for (j in 1:nrow(param_estim)){#Ici on calcule la matrice de Godambe pour chaque jeu de paramètres estimés
      V_inf=Estimateurs_esp(param_estim[j,],Obs_3,X,O)
      Vp_inf=ginv(V_inf[[1]])
      Godambe=Vp_inf%*%V_inf[[2]]%*%Vp_inf
      param_estim_norm[j,]<-(param_estim[j,]-param)/diag(Godambe)
    }
    name_table=paste(c('values',as.character(p),as.character(d)), collapse = "")
    name_graph=paste(c('graph',as.character(p),as.character(d)), collapse = "")
    write.table(param_estim_norm, file = paste(c("/home/teuliere/PLN-Cl/",name_table),collapse=""), append = FALSE, quote = TRUE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"))
    E_20=read.table(name_table,sep="\t",header = TRUE)
    E_202=melt(E_20)
    plt<-ggplot(data=E_202,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales="free")+
      theme_bw()
    png(paste(name_graph,".png",sep=""),width=1000, height=800,res=110)
    print(plt)
    dev.off()
  }
}
