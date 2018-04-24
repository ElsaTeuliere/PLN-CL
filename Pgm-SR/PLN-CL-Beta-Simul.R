# Composite likelihood for the Poisson log normal

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('PLN-CL-Beta-Functions.R')
library(PLNmodels)
library(nloptr)
library(mvtnorm)
library(poilog)

# Parms list
p.list = c(3, 10, 10); n.list = c(50, 50, 100); d.list = c(2, 2, 2); sim.nb = length(p.list)
rho = .5
B = 100; 

for (s in 1:sim.nb){
# for (s in sim.nb:sim.nb){
   p = p.list[s]; n = n.list[s]; d = d.list[s]; 
   P = p*(p-1)/2
   SimName = paste0('PLN-CL-Simul-p', p, '-n', n, '-d', d)
   
   # Model parms
   Sigma = matrix(rho, p, p); diag(Sigma) = 1; 
   sigma = sqrt(diag(Sigma)); Rho = cov2cor(Sigma)
   O = matrix(1, n, p); 
   W = matrix(1, p, p)
   
   Beta.vem = Beta.vem.sd = Beta.vem.stat = Beta.cl.pval = matrix(0, B, d*p)
   Beta.cl = Beta.cl.sd = Beta.cl.stat = Beta.vem.pval = matrix(0, B, d*p)
   x.list = seq(-3, +3, length.out=101); phi.list = dnorm(x.list)
   Iter.vem = Iter.cl = Time.vem = Time.vem.sd = Time.cl = Time.cl.sd = rep(0, B)
   for (b in 1:B){
      # b = 1
      cat('\nSimul ', b, ':')
      # Data simulation
      X = matrix(rnorm(d*n), n, d)
      Beta = matrix(rnorm(d*p), d, p)
      expOXbeta = exp(O + X %*% Beta)
      Z = rmvnorm(n, sigma=Sigma)
      Y = matrix(rpois(n*p, expOXbeta*exp(Z)), n, p)
      # plot(log(1+expOXbeta), log(1+Y)); abline(0, 1)
      Beta.true = as.vector(Beta)
      CL.true = CL.Beta(Beta.true, Y, X, O, sigma, Rho)
      d1.CL.true = d1.CL.Beta(Beta.true, Y, X, O, sigma, Rho)
      d1.CL.true
      d1.CL.Beta.old(Beta.true, Y, X, O, sigma, Rho); 
      
      # Fit PLN model
      Time.vem[b] = system.time(PLN <- PLN(Y ~ -1 + X + offset(O)))[[1]]
      Iter.vem[b] = PLN$optim_par$iterations
      Beta.vem.tmp = as.vector(t(PLN$model_par$Theta))
      CL.init = CL.Beta(Beta.vem.tmp, Y, X, O, sigma, Rho)
      d1.CL.init = d1.CL.Beta(Beta.vem.tmp, Y, X, O, sigma, Rho)
      cat('\n', Beta.vem.tmp, '/', max(abs(d1.CL.init)), '/', CL.init)
      
      # Approximate asymptotic covariance matrix
      d1.CL = rep(0, d*p)
      J = H = G = matrix(0, d*p, d*p)
      Time.vem.sd[b] = system.time(
         invisible(sapply(1:n, function(i){
            d1.CLi <- d1.CLi.Beta(Beta.vem.tmp, Y[i, ], X[i, ], O[i, ], sigma, Rho)
            d1.CL <<- d1.CL + d1.CLi 
            J <<- J + d1.CLi%o%d1.CLi / n
         })))[[1]]
      H <- - d2.CL.Beta(Beta.vem.tmp, Y, X, O, sigma, Rho) / n
      H.inv <- solve(H)
      G <- (H %*% solve(J) %*% H)
      Var.vem <- (H.inv %*% J %*% H.inv) / n
      
      # Results & test
      Sd.vem = sqrt(diag(Var.vem))
      T.vem = (Beta.vem.tmp - Beta.true)/Sd.vem
      # Pval.vem = 2*pnorm(abs(T.vem), lower.tail=F)
      Pval.vem = pnorm(abs(T.vem))
      print(rbind(Beta.true, Beta.vem.tmp, Sd.vem, T.vem, Pval.vem))
      print(Pval.vem)
      
      # Max CL: nloptr
      ctrl <- list(ftol_rel = ifelse(n < 1.5*p, 1e-6, 1e-8), ftol_abs = 0,
                   xtol_rel = 1e-4, xtol_abs = 1e-4, maxeval = 10000, method = "MMA")
      opts <- list( "algorithm" = paste("NLOPT_LD", ctrl$method, sep="_"),
                    "maxeval" = ctrl$maxeval, "ftol_rel" = ctrl$ftol_rel,
                    "ftol_abs" = ctrl$ftol_abs, "xtol_rel" = ctrl$xtol_rel,
                    "print_level" = max(0,ctrl$trace-1))
      Time.cl[b] = system.time(MCL <- nloptr(x0=Beta.vem.tmp, eval_f=negCL.Beta, eval_grad_f=d1.negCL.Beta, 
                  # lb = c(rep(-Inf, d), 0), ub = rep(Inf, (d+1)), 
                  opts=opts, Y=Y, X=X, O=O, sigma=sigma, Rho=Rho, W=W))[[1]]
      Iter.cl[b] = MCL$iterations
      Beta.cl.tmp = MCL$solution
      CL.cl = CL.Beta(Beta.cl.tmp, Y, X, O, sigma, Rho); 
      d1.CL.cl = d1.CL.Beta(Beta.cl.tmp, Y, X, O, sigma, Rho); 
      # print(cbind(rbind(Beta.true, Beta.vem.tmp, Beta.cl.tmp),
      # rbind(d1.CL.true, d1.CL.init, d1.CL.cl),
      # c(CL.true, CL.init, CL.cl)))
      cat('\n', Beta.cl.tmp, '/', max(abs(d1.CL.cl)), '/', CL.cl, '\n')
      cat(CL.true, CL.init, CL.cl, '\n')
      
      # Approximate asymptotic covariance matrix
      d1.CL = rep(0, d*p)
      J = H = G =  matrix(0, d*p, d*p)
      Time.cl.sd[b] = system.time(
         invisible(sapply(1:n, function(i){
            d1.CLi = d1.CLi.Beta(Beta.cl.tmp, Y[i, ], X[i, ], O[i, ], sigma, Rho)
            d1.CL <<- d1.CL + d1.CLi 
            J <<- J + d1.CLi%o%d1.CLi / n
         }))
      )[[1]]
      H = - d2.CL.Beta(Beta.cl.tmp, Y, X, O, sigma, Rho) / n
      H.inv = solve(H)
      G = (H %*% solve(J) %*% H)
      Var.cl = (H.inv %*% J %*% H.inv) / n
      
      # Results & test
      Sd.cl = sqrt(diag(Var.cl))
      T.cl = (Beta.cl.tmp - Beta.true)/Sd.cl
      # Pval.cl = 2*pnorm(abs(T.cl), lower.tail=F)
      Pval.cl = pnorm(abs(T.cl))
      print(rbind(Beta.true, Beta.cl.tmp, Sd.cl, T.cl, Pval.cl))
      print(Pval.cl)
      
      Beta.vem[b, ] = Beta.vem.tmp; Beta.cl[b, ] = Beta.cl.tmp
      Beta.vem.sd[b, ] = Sd.vem; Beta.vem.stat[b, ] = T.vem; Beta.vem.pval[b, ] = Pval.vem
      Beta.cl.sd[b, ] = Sd.cl; Beta.cl.stat[b, ] = T.cl; Beta.cl.pval[b, ] = Pval.cl
      
      # Plot
      if (b%%round(sqrt(B))==0){
         par(mfrow=c(2, 3))
         H = hist(as.vector(Beta.vem.stat[1:b, ]), breaks=sqrt(d*p*b))
         lines(x.list, b*d*p*mean(diff(H$breaks))*phi.list, col=2, lwd=2)
         H = hist(as.vector(Beta.vem.pval[1:b, ]), breaks=sqrt(d*p*b))
         lines(c(0, 1), rep(2*b*d*p*mean(diff(H$breaks)), 2), col=2, lwd=2)
         qqnorm(as.vector(Beta.vem.stat[1:b, ])); abline(0, 1)
         
         H = hist(as.vector(Beta.cl.stat[1:b, ]), breaks=sqrt(d*p*b))
         lines(x.list, b*d*p*mean(diff(H$breaks))*phi.list, col=2, lwd=2)
         H = hist(as.vector(Beta.cl.pval[1:b, ]), breaks=sqrt(d*p*b))
         lines(c(0, 1), rep(2*b*d*p*mean(diff(H$breaks)), 2), col=2, lwd=2)
         qqnorm(as.vector(Beta.cl.stat[1:b, ])); abline(0, 1)
      }
      print(c(Time.vem[b], Time.vem.sd[b], Time.cl[b], Time.cl.sd[b], Iter.vem[b], Iter.cl[b]))
   }
   save(Sigma, Beta, X, O, Beta.vem, Beta.cl, 
        Beta.vem.sd, Beta.vem.stat, Beta.vem.pval, 
        Beta.cl.sd, Beta.cl.stat, Beta.cl.pval, 
        Time.vem, Time.vem.sd, Time.cl, Time.cl.sd,
        Iter.vem, Iter.cl, 
        file=paste0(SimName, '.Rdata'))
   print(c(mean(Time.vem), mean(Time.vem.sd), mean(Time.cl), mean(Time.cl.sd)))
}
