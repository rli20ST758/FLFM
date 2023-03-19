
Renew_G_sigma <- function(Y_resid, C, beta, sigma2, knots=35, p=3, m=2, Nsubj, L, J){
  ############################################################ 
  #Arguments:
  #Y_resid -- (observed Y - mean function)
  #C -- the C estimated in the previous step
  #beta -- the beta estimated in the previous step
  #sigma2 <- sigma2_old
  #knots -- list of two vectors of knots or number of equidistant knots for all dimensions; defaults to 35
  #p -- degrees of B-splines; defaults to 3; details see fbps
  #m -- order of differencing penalty; defaults to 2; details see fbps
  #Nsubj -- number of subjects
  #L -- number of reduced dimensions
  #J -- number of variables
  ############################################################ 
  
  t <- ncol(C[[1]])
  G <- list()
  beta <- as.matrix(beta)

  for (j in 1:J) {
    
    ##############################################
    # Estimate Gj
    ##############################################
    Cov_Y <- t(Y_resid[[j]]) %*% Y_resid[[j]] / Nsubj
    # Variance caused by shared part
    share <- Reduce('+', lapply(1:L, function(x){beta[j,x]^2 * C[[x]] })) 
    # Noise for diagnal elements
    noise <- diag(rep(sigma2[j],t))
    
    # Raw covariance matrix for self part
    Gj <- Cov_Y - share - noise
    # Smoothing raw Gj  by the fast bivariate P-spline method
    BP_Gj <- fbps.cov(as.matrix(Gj), diag.remove=FALSE, knots=knots, p=p, m=m)
    G[[j]] <- BP_Gj$cov

    ##############################################
    # Estimate sigma
    ##############################################
    Cov_Y <- sum(Y_resid[[j]]^2) / Nsubj
    # Variance caused by shared part
    share <- Reduce('+', lapply(1:L, function(x){beta[j,x]^2 * sum(diag(C[[x]]))})) 
    # Variance caused by self part
    self <- sum(diag(G[[j]]))
    # Estimate sigma
    sigma2[j] <- (Cov_Y - share - self) / t

  }
  
  return(list(G=G, sigma2=sigma2))
}


