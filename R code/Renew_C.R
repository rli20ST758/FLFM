
Renew_C <- function(Yi, beta, G, sigma2, knots=35, p=3, m=2, Nsubj, L, J){
  ############################################################ 
  #Arguments:
  #Yi -- Yi with (K,J) dimension for each subject
  #beta -- the beta estimated in the previous step
  #G -- the estimated self covariance matrices in the previous step
  #sigma2 -- the sigma2 estimated in the previous step
  #knots -- list of two vectors of knots or number of equidistant knots for all dimensions; defaults to 35
  #p -- degrees of B-splines; defaults to 3; details see fbps
  #m -- order of differencing penalty; defaults to 2; details see fbps
  #Nsubj -- number of subjects
  #L -- number of reduced dimensions
  #J -- number of variables
  ############################################################ 
  
  t <- ncol(G[[1]])
  C <- list()
  beta <- as.matrix(beta)
  
  for (l in 1:L) {
    
    # Sample covariance
    YY <- lapply(1:Nsubj, function(x){ 
      Y_star <- Yi[[x]] %*% beta[,l]
      Y_star %*% t(Y_star) })
    Cov_Y <- Reduce('+', YY) / Nsubj
    # Delete variance caused by self part
    beta_G <- lapply(1:J, function(x){ beta[x,l]^2 * G[[x]] })  
    self <- Reduce('+', beta_G) 
    # Delete noise for diagnal elements
    noise <- sum(beta[,l]^2 * sigma2) * diag(rep(1,t))
    
    # Raw covariance matrix for shared part
    Cl <- Cov_Y - self - noise
    # Smoothing raw Cl  by the fast bivariate P-spline method
    BP_Cl <- fbps.cov(as.matrix(Cl), diag.remove=FALSE, knots=knots, p=p, m=m)
    C[[l]] <- BP_Cl$cov
    }
  
  return(C)
}


############################################################
# Renew_C2() for common G for all functional outcomes
############################################################
Renew_C2 <- function(Yi, beta, G, gamma, sigma2, knots=35, p=3, m=2, Nsubj, L, J){
  t <- ncol(G)
  C <- list()
  beta <- as.matrix(beta)
  
  for (l in 1:L) {
    
    # Sample covariance
    YY <- lapply(1:Nsubj, function(x){ 
      Y_star <- Yi[[x]] %*% beta[,l]
      Y_star %*% t(Y_star) })
    Cov_Y <- Reduce('+', YY) / Nsubj
    # Delete variance caused by self part
    beta_G <- lapply(1:J, function(x){ beta[x,l]^2 * gamma[x]^2 * G })  
    self <- Reduce('+', beta_G) 
    # Delete noise for diagnal elements
    noise <- sum(beta[,l]^2 * sigma2) * diag(rep(1,t))
    
    # Raw covariance matrix for shared part
    Cl <- Cov_Y - self - noise
    # Smoothing raw Cl  by the fast bivariate P-spline method
    BP_Cl <- fbps.cov(as.matrix(Cl), diag.remove=FALSE, knots=knots, p=p, m=m)
    C[[l]] <- BP_Cl$cov
  }
  
  return(C)
}
