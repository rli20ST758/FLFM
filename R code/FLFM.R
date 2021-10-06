
FLFM <- function(Y, Nsubj, argvals, L, sparse=FALSE, rho=0, pve1=0.95, pve2=0.95,
                initial=list(C=NULL,beta=NULL,G=NULL,sigma2=NULL),
                knots=35, p=3, m=2, tol=1e-4, Nmax=200, fit=FALSE) {
  ############################################################ 
  #Arguments:
  #argvals -- the argument values of the function evaluations 
  #Y -- observed data. A list of length J containing densely observed multivariate (J-dimensional) functional data. 
  #Y[[j]] is an (Nsubj x length(argvals)) matrix of functional data for Nsubj subjects observed on argvals.
  #Nsubj -- number of subjects
  #L -- number of reduced dimensions
  #sparse -- whether to update beta by sparse pca, default to FALSE. 
  #rho -- the penalty parameter for sparse pca, only useful when "sparse=TRUE". 
  #pve1 -- the proportion of variance explained for joint terms, defuault to 0.95. 
  #pve2 -- the proportion of variance explained for outcome-specific terms, defuault to 0.95. 
  #initial --  the initial values of parameters for the iteration. It can be from cross-validation for L. 
  #            When "sparse=TRUE", it highly recommends to set initial values as the estimation when "sparse=FALSE".
  #knots -- the number of equidistant knots for all dimensions; defaults to 35
  #p -- degrees of B-splines; defaults to 3; details see fbps
  #m -- order of differencing penalty; defaults to 2; details see fbps
  #tol -- tolerance for iteration, default to 10e-4
  #Nmax -- the number of maximal iteration, default to 200
  #fit -- logical with default FALSE. If TRUE, get the fitted values and scores
  ############################################################ 
  require(Rcpp)
  require(RSpectra)
  require(refund)
  require(face)
  require(splines)
  require(MASS)
  require(fda)
  require(Matrix)
  require(fields)
  require(rARPACK)
  require(quadprog)
  require(elasticnet)
  require(pracma)
  require(mgcv)
  
  # source("./R code/fbps.cov.R") # for smooth covariance matrices
  # source("./R code/spca_seqADMM.R") # for sparse beta
  # sourceCpp("./R code/eigendecomp.cpp") # for sparse beta
  # sourceCpp("./R code/MatrixMtp.cpp") # for sparse beta
  
  t <- length(argvals)
  J <- length(Y)
  ####################################
  ####step 1: Demean by fpca.face
  ####################################
  mu <- matrix(NA, nr=J, nc=t)
  sigma0 <- rep(NA, J)
  Y_resid <- list()
  for(j in 1:J){
    face_fit <- fpca.face(Y=Y[[j]], argvals=argvals, var=T)
    mu[j,] <- face_fit$mu
    # Demean for each variate
    Y_resid[[j]] <- Y[[j]] - matrix(rep(mu[j,],Nsubj), byrow=TRUE, nr=Nsubj)
    sigma0[j] <- face_fit$sigma2
  }
  
  #Yi for each subject 
  Yi <- lapply(1:Nsubj, function(x){ 
    do.call(cbind, lapply(1:J, function(j){ Y_resid[[j]][x,] }) )
  }) 
  
  
  ############################################
  ####step 2: Estimate parameters by iteration
  ############################################
  # initialize the values of parameters
  diff_C <- 1
  diff_beta <- 1
  diff_G <- 1
  diff_sigma2 <- 1
  if(!is.null(initial$C)) {
    C_new <- initial$C
  } else{
    C_new <- lapply(1:L, function(x){matrix(rnorm(t^2), nr=t, nc=t)})
  }
  if(!is.null(initial$beta)) {
    beta_new <- initial$beta
  } else{
    beta_new <- beta_init_func(Y=Y_resid, L=L)
  }
  if(!is.null(initial$G)) {
    G_new <- initial$G
  } else{
    G_new <- lapply(1:J, function(x){matrix(0, nr=t, nc=t)})
  }
  if(!is.null(initial$sigma2)) {
    sigma2_new <- initial$sigma2
  } else{
    sigma2_new <- sigma0
  }
  k = 0
  
  
  while( (diff_C>tol | diff_beta>tol | diff_G>tol | diff_sigma2>tol) & k <= Nmax){
    
    # Estimate each Cl
    C_old <- C_new
    C_new <- Renew_C(Yi=Yi, beta=beta_new, G=G_new, sigma2=sigma2_new, knots=knots, m=m, p=p, Nsubj=Nsubj, L=L, J=J)
    num_C <- unlist(lapply(1:L, function(x){norm(C_new[[x]]-C_old[[x]],"F")^2}))
    den_C <- unlist(lapply(1:L, function(x){norm(C_old[[x]],"F")^2}))
    diff_C <- sqrt(sum(num_C)/sum(den_C))
   
    # Estimate beta
    beta_old <- beta_new
    beta_new <- Renew_beta(Yi=Yi, C=C_new, beta=beta_old, G=G_new, sigma2=sigma2_new, Nsubj=Nsubj, L=L, J=J, sparse=sparse, rho=rho)
    diff_beta <- norm(beta_old-beta_new,"F")/norm(beta_old,"F")
    
    # Estimate each Gj and sigma2
    G_old <- G_new
    sigma2_old <- sigma2_new
    G_sigma <- Renew_G_sigma(Y_resid=Y_resid, C=C_new, beta=beta_new, sigma2=sigma2_old, knots=knots, m=m, p=p, Nsubj=Nsubj, L=L, J=J)
    G_new <- G_sigma$G
    sigma2_new <- G_sigma$sigma2
    num_G <- unlist(lapply(1:J, function(x){norm(G_new[[x]]-G_old[[x]],"F")^2}))
    den_G <- unlist(lapply(1:J, function(x){norm(G_old[[x]],"F")^2}))
    diff_G <- sqrt(sum(num_G)/sum(den_G))
    diff_sigma2 <- norm(sigma2_old-sigma2_new,"F")/norm(sigma2_old,"F")
    
    k = k + 1 
    c(diff_C,diff_beta,diff_G,diff_sigma2)
    } #end iteration 
    
  
  #####################################################################################
  ####step 3: Get estimated results and eigen components for covariance functions
  #####################################################################################
  # Estimated Cl
  C_temp<- C_new 
  order <- order(unlist(lapply(1:L, function(x){sum(diag(C_temp[[x]]))})), decreasing=TRUE)
  C_hat0 <- lapply(1:L, function(x){C_temp[[order[x]]]})
  # Estimated beta
  beta_hat <- as.matrix(beta_new)
  beta_hat[which(abs(beta_hat)<1e-5)] <- 0
  beta_hat <- as.matrix(beta_hat[,order])
  # Estimated Gj
  G_hat0 <- G_new
  # Estimated sigma2
  sigma2_hat <- sigma2_new
  
  
  # Get eigen components of C by pve1
  eigens_C <- lapply(C_hat0, function(x){eigen(x)})
  lambda_C0 <- lapply(eigens_C, function(x){(x$value/(t-1))[which(x$value>0)]})
  order_lambda_C <- sort(unlist(lambda_C0), decreasing=TRUE) #order all lambda_C
  lambda_C_pve <- order_lambda_C/sum(order_lambda_C) #the proportion of variance explained by each lambda_C
  lambda_C_star <- order_lambda_C[which(cumsum(lambda_C_pve)>pve1)[1]]
  lambda_C <- lapply(lambda_C0,function(x){ #eigenvalues
    num <- sum(x>=lambda_C_star)
    if(num==0){return(NULL)}
    return(x[1:num])})
  phi <- lapply(1:L, function(x){ #eigenvectors
    eigenvec <- (eigens_C[[x]])$vectors * sqrt(t-1)
    if(length(lambda_C[[x]])==0) NULL
    else as.matrix(eigenvec[,1:length(lambda_C[[x]])]) })
  C_hat <- lapply(1:L, function(x){
    num <- length(lambda_C[[x]])
    if(num==0){return(NULL)}
    if(num==1){return(lambda_C[[x]]*phi[[x]]%*%t(phi[[x]]))}
    return(Reduce("+", lapply(1:num, function(y){lambda_C[[x]][y]*phi[[x]][,y]%*%t(phi[[x]][,y])})))
  })

  # Get eigen components of G by pve2
  eigens_G <- lapply(G_hat0, function(x){eigen(x)})
  lambda_G0 <- lapply(eigens_G, function(x){(x$value/(t-1))[which(x$value>0)]})
  order_lambda_G <- sort(unlist(lambda_G0), decreasing=TRUE) #order all lambda_G
  lambda_G_pve <- order_lambda_G/sum(order_lambda_G)
  lambda_G_star <- order_lambda_G[which(cumsum(lambda_G_pve)>pve2)[1]]
  lambda_G <- lapply(lambda_G0,function(x){ #eigenvalues
    num <- sum(x>=lambda_G_star)
    if(num==0){return(NULL)}
    return(x[1:num])})
  psi <- lapply(1:J, function(x){ #eigenvectors
    eigenvec <- (eigens_G[[x]])$vectors * sqrt(t-1)
    if(length(lambda_G[[x]])==0) NULL
    else as.matrix(eigenvec[,1:length(lambda_G[[x]])]) })
  G_hat <- lapply(1:J, function(x){
    num <- length(lambda_G[[x]])
    if(num==0){return(NULL)}
    if(num==1){return(lambda_G[[x]]*psi[[x]]%*%t(psi[[x]]))}
    return(Reduce("+", lapply(1:num, function(y){lambda_G[[x]][y]*psi[[x]][,y]%*%t(psi[[x]][,y])})))
  })

  
  ######################################
  ####step 4: Estimate scores
  ######################################
  Yhat <- NULL
  xi <- NULL
  eta <- NULL

  # set NULL to 0
  lambda_C_temp = lambda_C
  lambda_G_temp = lambda_G
  phi_temp = phi
  psi_temp = psi
  for (x in 1:L) {
    if(length(lambda_C[[x]])==0){
      lambda_C_temp[[x]] = 1e-6
      phi_temp[[x]] = matrix(0, nc=1, nr=t)
    }
  }
  for (x in 1:J) {
    if(length(lambda_G[[x]])==0){
      lambda_G_temp[[x]] = 1e-6
      psi_temp[[x]] = matrix(0, nc=1, nr=t)
    }
  }
  
  if(fit){
    Mat1 <- do.call(cbind, lapply(1:L,function(x){kronecker(beta_hat[,x], phi_temp[[x]])}))
    Mat2 <- do.call(cbind, lapply(1:L,function(x){kronecker(beta_hat[,x]/sigma2_hat, phi_temp[[x]])}))
    Mat3 <- lapply(1:J, function(x){t(psi_temp[[x]])%*%psi_temp[[x]]/sigma2_hat[x]})
    Mat4 <- lapply(1:J, function(x){t(psi_temp[[x]])/sigma2_hat[x]})
    invGamma1 <- bdiag(lapply(lambda_C_temp,function(x){if(length(x)>1) diag(1/x) else 1/x}))
    invGamma2 <- bdiag(lapply(lambda_G_temp,function(x){if(length(x)>1) diag(1/x) else 1/x}))
    A <- t(Mat1) %*% Mat2 + invGamma1
    C <- bdiag(Mat4) %*% Mat1
    D <- bdiag(Mat3) + invGamma2
    invD <- solve(D)
    MatE <- solve(A-t(C)%*%invD%*%C)
    MatF <- - MatE%*%t(C)%*%invD
    MatG <- - invD%*%C%*%MatE
    MatH <- invD - invD%*%C%*%MatF
    part_xi <- MatE%*%t(Mat2) + MatF%*%bdiag(Mat4)
    part_eta <- MatG%*%t(Mat2) + MatH%*%bdiag(Mat4)
    y_vec <- do.call(cbind, lapply(Yi, function(x){as.vector(x)}))
    xi <- part_xi %*% y_vec
    eta <- part_eta %*% y_vec
    
    part1 <- Mat1%*%xi
    part2 <- bdiag(psi_temp) %*% eta
    yhat <- matrix(part1+part2, nc=t, byrow=T) + kronecker(rep(1,Nsubj),mu)
    idx <- rep(1:J, Nsubj)
    Yhat <- lapply(1:J, function(x){yhat[idx==x,]})
  }
  
  # delete scores corresponding to lambda=NULL
  idxC <- which(unlist(lambda_C_temp)==1e-6)
  if(length(idxC) > 0)  xi=xi[-idxC,]
  idxG <- which(unlist(lambda_G_temp)==1e-6)
  if(length(idxG) > 0)  eta=eta[-idxG,]
  
  
  return(list(Yhat=Yhat, meanfunctions=mu, xi=xi, eta=eta, C=C_hat, beta=beta_hat, G=G_hat, 
              sigma2=sigma2_hat, lambda_C=lambda_C, phi=phi, lambda_G=lambda_G, psi=psi, k=k))
}




beta_init_func <- function(Y, L){
  #   ############################################################ 
  #   #Arguments:
  #   #Y <- Y_residual (observed Y - mean function)
  #   #argval - number of points at which functions are observed
  #   #L -- number of reduced dimensions
  #   ############################################################  
  J <- length(Y)
  Nsubj <- nrow(Y[[1]])
  Y_vec <- do.call('cbind', lapply(1:J, function(x){colSums(Y[[x]])}))
  beta_init <- eigs(cov(Y_vec),k=L)$vectors
}



Renew_C <- function(Yi, beta, G, sigma2, knots=35, p=3, m=2, Nsubj, L, J){
  ############################################################ 
  #Arguments:
  #Yi -- Yi with (K,J) dimension for each subject
  #beta -- the beta estimated in the previous step
  #G -- the estimated self covariance matrices in the previous step
  #sigma2 -- the sigma2 estimated in the previous step
  #knots -- the number of equidistant knots for all dimensions; defaults to 35
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



Renew_beta <- function(Yi, C, beta, G, sigma2, Nsubj, L, J, sparse, rho){
  ############################################################ 
  #Arguments:
  #Yi -- Yi with (K,J) dimension for each subject
  #C -- the C estimated in the previous step
  #beta -- beta_old
  #G -- the estimated self covariance matrices in the previous step
  #sigma2 -- the sigma2 estimated in the previous step
  #Nsubj -- number of subjects
  #L -- number of reduced dimensions
  #J -- number of variables
  #sparse -- whether to update beta by sparse pca, default to FALSE. 
  #rho -- the penalty parameter for sparse pca, only useful when "sparse=TRUE". 
  ############################################################ 
  
  if(L==1){
    
    # Compute the sample covairance matrix for updating beta
    Yi_C_Yi <- lapply(1:Nsubj, function(x) {
      temp <- t(Yi[[x]]) %*% C[[1]]
      temp %*% Yi[[x]]})
    Y_C_Y <- Reduce('+', Yi_C_Yi)
    # Variance caused by self part
    C_Gj <- unlist( lapply(1:J, function(x){sum(C[[1]]*G[[x]])}) )
    self <- diag(C_Gj)
    # Delete noise for diagnal elements
    noise  <- sum(diag(C[[1]])) * diag(sigma2)
    
    # Update beta 
    A <- self + noise - Y_C_Y/Nsubj 
    if(sparse==TRUE){
      S <- t(A)%*%A / mean(diag(t(A)%*%A))
      H <- LA_seqadmm(S=S,eta=100,rho=rho,t=10,K=20,etastep=2) #sparse pca
      beta_temp <- getEigenDecomp(H)[[2]][,ncol(H)]
    } else{
      beta_temp <- eigs_sym(A, 1, which="SA")$vectors
    }
    beta <- sign(beta_temp[1]) * beta_temp
    
  } else {
    
    # Compute tr( t(Cl)%*%Cl')
    ClCl <- matrix(0, nrow=L, ncol=L)
    for (l1 in 1:L) {
      for (l2 in l1:L) { 
        ClCl[l1,l2] <- sum(C[[l1]]*C[[l2]])
      }
    }
    ClCl <- ClCl + t(ClCl) - diag(diag(ClCl))
    
    ##############################################
    # Update beta column by column
    ##############################################
    for(l in 1:L){ 
      
      # Compute the intersect part with beta_l'
      Cl_beta0 <- lapply((1:L)[-l], function(x){ClCl[l,x] * beta[,x] %*% t(beta[,x])}) 
      Cl_beta <- Reduce('+', Cl_beta0)
      # Compute the sample covairance matrix for updating beta
      Yi_Cl_Yi <- lapply(1:Nsubj, function(x) {
        temp <- t(Yi[[x]]) %*% C[[l]]
        temp %*% Yi[[x]]})
      Y_Cl_Y <- Reduce('+', Yi_Cl_Yi)
      # Variance caused by self part
      Cl_Gj <- unlist( lapply(1:J, function(x){sum(C[[l]]*G[[x]])}) )
      self <- diag(Cl_Gj)
      # Delete noise for diagnal elements
      noise <- sum(diag(C[[l]])) * diag(sigma2)
      
      # Update beta 
      A <- Cl_beta + self + noise - Y_Cl_Y/Nsubj 
      if(sparse==TRUE){
        S <- -A / mean(abs(diag(A)))
        if(l==1) { # sparse beta
          H <- LA_seqadmm(S=S,eta=100,rho=rho,t=10,K=20,etastep=2)
          beta_temp <- getEigenDecomp(H)[[2]][,ncol(H)]
          beta[,l] <- sign(beta_temp[1]+1e-6) * beta_temp
        } else {
          PrevPi <- beta[,1:(l-1)]%*%t(beta[,1:(l-1)])
          PrevPi_Eig <- getEigenDecomp(diag(1,J)-PrevPi)
          PrevPi_Eig[[1]] <- rev(PrevPi_Eig[[1]])
          PrevPi_Eig[[2]] <- PrevPi_Eig[[2]][,ncol(PrevPi_Eig[[2]]):1]
          H <- LA_seqadmm(S=S,PrevPi_Eig=PrevPi_Eig,PrevPi_d=l-1,eta=100,rho=rho,t=10,K=20,etastep=2)
          beta_temp <- getEigenDecomp(H)[[2]][,ncol(H)]
          beta[,l] <- sign(beta_temp[1]+1e-6) * beta_temp
        }
      } else{
        if(l!=1){ # general beta
          P <- diag(J) - beta[,1:(l-1)] %*% t(beta[,1:(l-1)])
          A <- P %*% A %*% P
        }
        beta_temp <- eigs_sym(A, 1, which="SA")$vectors
        beta[,l] <- sign(beta_temp[1]+1e-6) * beta_temp
      } # end if sparse
      
    } # end for loop
  } # end if else
  return(beta)
}



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


