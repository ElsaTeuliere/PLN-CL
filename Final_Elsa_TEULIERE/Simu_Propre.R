######################################################################################################
## Script pour les simulations classiques
######################################################################################################

source('./Fonction_CL_PLN.R')




n=100
nb_simu=250
p=5
d=3
O=matrix(1,n,p)
X=matrix(runif(n,-5,5),n,1)
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
  param_optimaux<-nloptr(x0=x_0, eval_f=neg_CL,
                         eval_g_ineq = constraint, 
                         lb = lb, ub = ub,
                         opts = list("algorithm" = "NLOPT_GN_DIRECT_L",
                                     "print_level" = 2,
                                     "xtol_rel"=1.0e-7,
                                     "xtol_abs"=1.0e-7,
                                     "maxeval"=1000), Y=Obs_3, X=X,O=O)
  param_estim=rbind(param_estim,param_optimaux$solution)
  c<-c+1
  nombre_iterations[length(nombre_iterations)+1]<-param_optimaux$iter
  V_inf=Estimateurs_esp_corr_st(param_optimaux$solution,Obs_3,X,O)
  hess=symetrisation(V_inf[[1]])
  A=diag(hess)
  hess_ok=hess+t(hess)-diag(A)
  Vp_inf=ginv(hess_ok)
  Godambe=Godambe + Vp_inf%*%V_inf[[2]]%*%Vp_inf
}
matrix_param=matrix(rep(param),nrow=250,ncol=length(param),byrow=TRUE)
param_estim_norm<-sqrt(nb_simu)(param_estim-matrix_param)/diag((1/nb_simu) * Godambe)
temps_moyen[[length(temps_moyen)+1]]<-tps/c
name_table=paste(c('values',as.character(p),as.character(d),"_new"), collapse = "")
name_table_norm=paste(c('values',as.character(p),as.character(d),"_norm_st_DIRECT_L"), collapse = "")
name_graph=paste(c('graph',as.character(p),as.character(d),"_st_DIRETC_L_tol"), collapse = "")
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