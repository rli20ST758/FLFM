
########################################################################################
# This function is for common individual covariance functions for outcome-specific terms
#########################################################################################

Renew_G_gamma_sigma <- function(Y_resid, C, beta, gamma, sigma2, knots=35, p=3, m=2, Nsubj, L, J){
  ############################################################ 
  #Arguments:
  #Y_resid -- (observed Y - mean function)
  #C -- the C estimated in the previous step
  #beta -- the beta estimated in the previous step
  #gamma <- gamma_old
  #sigma2 <- sigma2_old
  #knots -- list of two vectors of knots or number of equidistant knots for all dimensions; defaults to 35
  #p -- degrees of B-splines; defaults to 3; details see fbps
  #m -- order of differencing penalty; defaults to 2; details see fbps
  #Nsubj -- number of subjects
  #L -- number of reduced dimensions
  #J -- number of variables
  ############################################################ 
  
  t <- ncol(C[[1]])
  beta <- as.matrix(beta)
  
  ##############################################
  # Estimate G
  ##############################################
  Ej <- lapply(1:J, function(j){
    Cov_Y <- t(Y_resid[[j]]) %*% Y_resid[[j]] / Nsubj
    # Variance caused by shared part
    share <- Reduce('+', lapply(1:L, function(x){beta[j,x]^2 * C[[x]] })) 
    # Noise for diagnal elements
    noise <- diag(rep(sigma2[j],t))
    # Ej for self part
    Cov_Y - share - noise
  })
  G <- Reduce("+", lapply(1:J, function(j){gamma[j]^2 * Ej[[j]]}))/sum(gamma^4)
  # Smoothing raw G  by the fast bivariate P-spline method
  BP_G <- fbps.cov(as.matrix(G), diag.remove=FALSE, knots=knots, p=p, m=m)
  G <- BP_G$cov
  
  for (j in 1:J) {
    ##############################################
    # Estimate gamma
    ##############################################
    gamma[j] <- sqrt(sum(diag(t(G)%*%Ej[[j]]))/sum(diag(t(G)%*%G)))
    
    ##############################################
    # Estimate sigma
    ##############################################
    Cov_Y <- sum(Y_resid[[j]]^2) / Nsubj
    # Variance caused by shared part
    share <- Reduce('+', lapply(1:L, function(x){beta[j,x]^2 * sum(diag(C[[x]]))})) 
    # Variance caused by self part
    self <- gamma[j]^2 * sum(diag(G))
    # Estimate sigma
    sigma2[j] <- (Cov_Y - share - self) / t
  }
  
  gamma <- gamma / gamma[1]
  return(list(G=G, gamma=gamma, sigma2=sigma2))
}


