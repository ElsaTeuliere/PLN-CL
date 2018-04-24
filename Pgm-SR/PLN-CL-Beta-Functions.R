###############################################################################
# CL with known Sigma = (sigma2, Rho)
CL.Beta <- function(Beta.vec, Y, X, O, sigma, Rho, W=matrix(1, ncol(Y), ncol(Y))){
   p = ncol(Y); d = ncol(X); 
   Mu = O + X%*%matrix(Beta.vec, d, p)
   CL = 0; 
   sapply(1:(p-1), function(j){sapply((j+1):p, function(k){
      sapply(1:nrow(Y), function(i){
         CL <<- CL + W[j, k]*log(dbipoilog(n1=Y[i, j], n2=Y[i, k], mu1=Mu[i, j], mu2=Mu[i, k], 
                                   sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k]))
      })
   })})
   return(CL)
}
negCL.Beta <- function(Beta.vec, Y, X, O, sigma, Rho, W=matrix(1, ncol(Y), ncol(Y))){
   return(-CL.Beta(Beta.vec, Y, X, O, sigma, Rho, W))
}

###############################################################################
# derivative of CLi with known Sigma = (sigma2, Rho)
d1.CLi.Betaj <- function(Beta.vec, j, Yi, Xi, Oi, sigma, Rho, W=matrix(1, length(Yi), length(Yi))){
   # Beta.vec = as.vector(Beta0)
   p = length(Yi); d = length(Xi); 
   mu = Oi + Xi%*%matrix(Beta.vec, d, p)
   d1CLi.Betaj = matrix(0, d, p); 
   sapply(1:p, function(k){
      if (k != j){
         pjk = dbipoilog(n1=Yi[j], n2=Yi[k], mu1=mu[j], mu2=mu[k], 
                         sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
         pj1k = dbipoilog(n1=Yi[j]+1, n2=Yi[k], mu1=mu[j], mu2=mu[k], 
                          sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
         d1CLi.Betaj[, j] <<- d1CLi.Betaj[, j] + W[j, k] * (Yi[j]*pjk - (Yi[j]+1)*pj1k) * Xi / pjk
         }})
   return(as.vector(d1CLi.Betaj))
}

###############################################################################
# derivative of CLi with known Sigma = (sigma2, Rho)
d1.CLi.Beta <- function(Beta.vec, Yi, Xi, Oi, sigma, Rho, W=matrix(1, length(Yi), length(Yi))){
   # Beta.vec = as.vector(Beta0)
   d1CLi.Beta = rep(0, length(Beta.vec)); 
   sapply(1:length(Yi), function(j){
      d1CLi.Beta <<- d1CLi.Beta + d1.CLi.Betaj(Beta.vec, j, Yi, Xi, Oi, sigma, Rho, W)
      # print(d1CLi.Beta)
   })
   return(d1CLi.Beta)
}

###############################################################################
# derivative of CLi with known Sigma = (sigma2, Rho)
d1.CLi.Beta.old <- function(Beta.vec, Yi, Xi, Oi, sigma, Rho, W=matrix(1, length(Yi), length(Yi))){
   # Beta.vec = as.vector(Beta0)
   p = length(Yi); d = length(Xi); 
   mu = Oi + Xi%*%matrix(Beta.vec, d, p)
   d1CLi.Beta = matrix(0, d, p); 
   sapply(1:(p-1), function(j){sapply((j+1):p, function(k){
      pjk = dbipoilog(n1=Yi[j], n2=Yi[k], mu1=mu[j], mu2=mu[k], 
                      sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
      pj1k = dbipoilog(n1=Yi[j]+1, n2=Yi[k], mu1=mu[j], mu2=mu[k], 
                       sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
      pjk1 = dbipoilog(n1=Yi[j], n2=Yi[k]+1, mu1=mu[j], mu2=mu[k], 
                       sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
      d1CLi.Beta[, j] <<- d1CLi.Beta[, j] + W[j, k] * (Yi[j]*pjk - (Yi[j]+1)*pj1k) * Xi / pjk
      d1CLi.Beta[, k] <<- d1CLi.Beta[, k] + W[j, k] * (Yi[k]*pjk - (Yi[k]+1)*pjk1) * Xi / pjk
   })})
   return(as.vector(d1CLi.Beta))
}

###############################################################################
# derivative of CL with known Sigma = (sigma2, Rho)
d1.CL.Beta <- function(Beta.vec, Y, X, O, sigma, Rho, W=matrix(1, ncol(Y), ncol(Y))){
   d1CL.Beta = rep(0, length(Beta.vec))
   sapply(1:nrow(Y), function(i){d1CL.Beta <<- 
      d1CL.Beta + d1.CLi.Beta(Beta.vec, Y[i, ], X[i, ], O[i, ], sigma, Rho, W)})
   return(d1CL.Beta)
}
d1.CL.Beta.old <- function(Beta.vec, Y, X, O, sigma, Rho, W=matrix(1, ncol(Y), ncol(Y))){
   d1CL.Beta = rep(0, length(Beta.vec))
   sapply(1:nrow(Y), function(i){d1CL.Beta <<- 
      d1CL.Beta + d1.CLi.Beta.old(Beta.vec, Y[i, ], X[i, ], O[i, ], sigma, Rho, W)})
   return(d1CL.Beta)
}
d1.negCL.Beta <- function(Beta.vec, Y, X, O, sigma, Rho, W=matrix(1, ncol(Y), ncol(Y))){
   return(-d1.CL.Beta(Beta.vec, Y, X, O, sigma, Rho, W))
}

###############################################################################
# - second derivative of CLi with known Sigma = (sigma2, Rho)
d2.CLi.Beta <- function(Beta.vec, Yi, Xi, Oi, sigma, Rho, W=matrix(1, length(Yi), length(Yi))){
   # Beta.vec = as.vector(Beta0)
   mu = Oi + Xi%*%matrix(Beta.vec, d, p)
   d2CLi.Beta = matrix(0, p*d, p*d); 
   # image(d2CLi.Beta);
   sapply(1:(p-1), function(j){
      posj = ((j-1)*d+1):(j*d)
      sapply((j+1):p, function(k){
         posk = ((k-1)*d+1):(k*d)
         pjk = dbipoilog(n1=Yi[j], n2=Yi[k], mu1=mu[j], mu2=mu[k], 
                         sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
         pj1k = dbipoilog(n1=Yi[j]+1, n2=Yi[k], mu1=mu[j], mu2=mu[k], 
                          sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
         pj2k = dbipoilog(n1=Yi[j]+2, n2=Yi[k], mu1=mu[j], mu2=mu[k], 
                          sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
         pjk1 = dbipoilog(n1=Yi[j], n2=Yi[k]+1, mu1=mu[j], mu2=mu[k], 
                          sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
         pjk2 = dbipoilog(n1=Yi[j], n2=Yi[k]+2, mu1=mu[j], mu2=mu[k], 
                          sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
         pj1k1 = dbipoilog(n1=Yi[j]+1, n2=Yi[k]+1, mu1=mu[j], mu2=mu[k], 
                           sig1=sigma[j], sig2=sigma[k], rho=Rho[j, k])
         Ci = matrix(0, 2, 2)
         Ci[1, 1] = Yi[j]^2*pjk - (2*Yi[j]+1)*(Yi[j]+1)*pj1k + (Yi[j]+1)*(Yi[j]+2)*pj2k
         Ci[2, 2] = Yi[k]^2*pjk - (2*Yi[k]+1)*(Yi[k]+1)*pjk1 + (Yi[k]+1)*(Yi[k]+2)*pjk2
         Ci[1, 2] = Yi[j]*Yi[k]*pjk - (Yi[j]+1)*Yi[k]*pj1k - Yi[j]*(Yi[k]+1)*pjk1 + 
            (Yi[j]+1)*(Yi[k]+1)*pj1k1
         Ci[2, 1] = Ci[1, 2]
         Bi = c( Yi[j]*pjk - (Yi[j]+1)*pj1k, Yi[k]*pjk - (Yi[k]+1)*pjk1 )
         d2CLi.Beta[c(posj, posk), c(posj, posk)] <<-
            d2CLi.Beta[c(posj, posk), c(posj, posk)] +
            W[j, k] * kronecker((Ci/pjk - (Bi%o%Bi)/pjk^2), (Xi%o%Xi))
         # image(1:(d*p), 1:(d*p), sign(d2CLi.Beta)*log(abs(d2CLi.Beta)), xlab=paste(j, k));
         # print(Bi)
      })})
   return(d2CLi.Beta)
}

###############################################################################
# - second derivative of CL with known Sigma = (sigma2, Rho)
d2.CL.Beta <- function(Beta.vec, Y, X, O, sigma, Rho, W=matrix(1, ncol(Y), ncol(Y))){
   # print(Beta.vec); print(Y); print(X); print(O); print(sigma); print(Rho)
   d2CL.Beta = matrix(0, length(Beta.vec), length(Beta.vec))
   sapply(1:nrow(Y), function(i){d2CL.Beta <<- d2CL.Beta + 
      d2.CLi.Beta(Beta.vec, Y[i, ], X[i, ], O[i, ], sigma, Rho, W)})
   return(d2CL.Beta)
}
d2.negCL.Beta <- function(Beta.vec, Y, X, O, sigma, Rho, W=matrix(1, ncol(Y), ncol(Y))){
   return(-d2.CL.Beta(Beta.vec, Y, X, O, sigma, Rho, W))
}
