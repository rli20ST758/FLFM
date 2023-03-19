#' The R function FLFM() is for Functional Latent Factor Model
#' 
#' @param Y observed data. A list of length J containing densely observed multivariate (J-dimensional) functional data.  
#'        Y[[j]] is an (Nsubj x length(argvals)) matrix of functional data for Nsubj subjects observed on argvals
#' @param Nsubj number of subjects
#' @param argvals the argument values of the function evaluations 
#' @param L number of reduced dimensions. L=0 is only allowed for common_G = FALSE
#' @param sparse Logical with default FALSE. If TRUE, beta are sparse
#' @param rho the penalty parameter for sparse pca, only useful when sparse is TRUE
#' @param pve1 the proportion of variance explained for joint terms, default to 0.95 
#' @param pve2 the proportion of variance explained for outcome-specific terms, default to 0.95
#' @param common_G Logical with default FALSE. If TRUE, use common individual covariance functions for outcome-specific terms
#' @param initial the initial values of parameters for the iteration. It can be from cross-validation for L
#'                For sparse=TRUE, it highly recommends to set initial values as the estimation for sparse=FALSE
#' @param knots list of two vectors of knots or number of equidistant knots for all dimensions; defaults to 35
#' @param p degrees of B-splines; defaults to 3; details see fbps
#' @param m order of differencing penalty; defaults to 2; details see fbps
#' @param tol tolerance for iteration, default to 10e-4
#' @param Nmax the number of maximal iteration, default to 200
#' @param fit Logical with default FALSE. If TRUE, get the fitted values and scores
#' 
#' @return 
#' @name FLFM
#' @author Ruonan Li and Luo Xiao
#' @export
#' 

FLFM <- function(Y, Nsubj, argvals, L, sparse = FALSE, rho = 0, pve1 = 0.95, pve2 = 0.95, common_G = FALSE,
                 initial = list(C=NULL, beta=NULL, G=NULL, gamma=NULL, sigma2=NULL),
                 knots = 35, p = 3, m = 2, tol = 1e-4, Nmax = 200, fit = FALSE) {
  
  if(common_G == FALSE) {
    LFM(Y = Y, Nsubj = Nsubj, argvals = argvals, L = L, sparse = sparse, 
        rho = rho, pve1 = pve1, pve2 = pve2, initial = initial,
        knots = knots, p = p, m = m, tol = tol, Nmax = Nmax, fit = fit)
  } else{
    LFM2(Y = Y, Nsubj = Nsubj, argvals = argvals, L = L, sparse = sparse, 
         rho = rho, pve1 = pve1, pve2 = pve2, initial = initial,
         knots = knots, p = p, m = m, tol = tol, Nmax = Nmax, fit = fit)
  }
  
}



####################################################################################### 
# FLFM using different individual covariance functions for outcome-specific terms
#######################################################################################
LFM <- function(Y, Nsubj, argvals, L, sparse = FALSE, rho = 0, pve1 = 0.95, pve2 = 0.95, 
                initial = list(C=NULL, beta=NULL, G=NULL, gamma=NULL, sigma2=NULL),
                knots = 35, p = 3, m = 2, tol = 1e-4, Nmax = 200, fit = FALSE) {
  
  t <- length(argvals)
  J <- length(Y)
  mu <- matrix(NA, nr=J, nc=t)
  sigma0 <- rep(NA, J)
  Y_resid <- list()
  
  ############################################
  ####if L is 0
  ############################################
  if(L==0){
    lambda_G <- list()
    G_hat <- list()
    psi <- list()
    zeta <- list()
    Yhat <- list()
    for(j in 1:J){
      face_fit <- fpca.face(Y=Y[[j]], argvals=argvals, pve=pve1, var=T)
      mu[j,] <- face_fit$mu
      # Demean for each variate
      Y_resid[[j]] <- Y[[j]] - matrix(rep(mu[j,],Nsubj), byrow=TRUE, nr=Nsubj)
      sigma0[j] <- face_fit$sigma2
      lambda_G[[j]] <- face_fit$evalues/t
      psi[[j]] <- face_fit$efunctions * sqrt(t)
      zeta[[j]] <- t(face_fit$scores/sqrt(t))
      Yhat[[j]] <- face_fit$Yhat
      G_hat[[j]] <- Reduce("+", lapply(1:length(lambda_G[[j]]), function(l){
        lambda_G[[j]][l]*psi[[j]][,l]%*%t(psi[[j]][,l])}))
    }
    sigma2_hat <- sigma0
    zeta <- Reduce(rbind, zeta)
    k <- 1
    beta_hat <- NULL
    lambda_C <- NULL
    phi <- NULL
    xi <- NULL
    C_hat <- NULL
  } 
  
  
  if(L > 0){
    ####################################
    ####step 1: Demean by fpca.face
    ####################################
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
    # Initilize the values of parameters
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
    lambda_G <- lapply(1:J,function(j){ #eigenvalues
      idx <- which(cumsum(lambda_G0[[j]]/sum(lambda_G0[[j]]))>pve2)[1]
      return(lambda_G0[[j]][1:idx]) })
    psi <- lapply(1:J, function(x){ #eigenvectors
      eigenvec <- (eigens_G[[x]])$vectors * sqrt(t-1)
      as.matrix(eigenvec[,1:length(lambda_G[[x]])])
    })
    G_hat <- lapply(1:J, function(x){
      num <- length(lambda_G[[x]])
      return(Reduce("+", lapply(1:num, function(y){lambda_G[[x]][y]*psi[[x]][,y]%*%t(psi[[x]][,y])})))
    })
    
    
    ######################################
    ####step 4: Estimate scores
    ######################################
    Yhat <- NULL
    xi <- NULL
    zeta <- NULL
    
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
      part_zeta <- MatG%*%t(Mat2) + MatH%*%bdiag(Mat4)
      y_vec <- do.call(cbind, lapply(Yi, function(x){as.vector(x)}))
      xi <- part_xi %*% y_vec
      zeta <- part_zeta %*% y_vec
      
      part1 <- Mat1%*%xi
      part2 <- bdiag(psi_temp) %*% zeta
      yhat <- matrix(part1+part2, nc=t, byrow=T) + kronecker(rep(1,Nsubj),mu)
      idx <- rep(1:J, Nsubj)
      Yhat <- lapply(1:J, function(x){yhat[idx==x,]})
    }
    
    # delete scores corresponding to lambda=NULL
    idxC <- which(unlist(lambda_C_temp)==1e-6)
    if(length(idxC) > 0)  xi=xi[-idxC,]
    idxG <- which(unlist(lambda_G_temp)==1e-6)
    if(length(idxG) > 0)  zeta=zeta[-idxG,]
  }
  
  
  return(list(Yhat=Yhat, meanfunctions=mu, xi=xi, zeta=zeta, C=C_hat, beta=beta_hat, G=G_hat,
              common_G = FALSE, sigma2=sigma2_hat, lambda_C=lambda_C,
              phi=phi, lambda_G=lambda_G, psi=psi, k=k))
}



####################################################################################### 
# FLFM using common individual covariance functions for outcome-specific terms
#######################################################################################
LFM2 <- function(Y, Nsubj, argvals, L, sparse = FALSE, rho = 0, pve1 = 0.95, pve2 = 0.95, 
                 initial = list(C=NULL, beta=NULL, G=NULL, gamma=NULL, sigma2=NULL),
                 knots = 35, p = 3, m = 2, tol = 1e-4, Nmax = 200, fit = FALSE) {
  
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
  # Initilize the values of parameters
  diff_C <- 1
  diff_beta <- 1
  diff_G <- 1
  diff_gamma <- 1
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
    G_new <- matrix(0, nr=t, nc=t)
  }
  if(!is.null(initial$gamma)) {
    gamma_new <- initial$gamma
  } else{
    gamma_new <- rep(1,J)
  }
  if(!is.null(initial$sigma2)) {
    sigma2_new <- initial$sigma2
  } else{
    sigma2_new <- sigma0
  }
  k = 0
  
  
  while( (diff_C>tol | diff_beta>tol | diff_G>tol | diff_gamma>tol | diff_sigma2>tol) & k <= Nmax){
    
    # Estimate each Cl
    C_old <- C_new
    C_new <- Renew_C2(Yi=Yi, beta=beta_new, G=G_new, gamma=gamma_new, sigma2=sigma2_new, knots=knots, m=m, p=p, Nsubj=Nsubj, L=L, J=J)
    num_C <- unlist(lapply(1:L, function(x){norm(C_new[[x]]-C_old[[x]],"F")^2}))
    den_C <- unlist(lapply(1:L, function(x){norm(C_old[[x]],"F")^2}))
    diff_C <- sqrt(sum(num_C)/sum(den_C))
    
    # Estimate beta
    beta_old <- beta_new
    beta_new <- Renew_beta2(Yi=Yi, C=C_new, beta=beta_old, G=G_new, gamma=gamma_new, sigma2=sigma2_new, Nsubj=Nsubj, L=L, J=J, sparse=sparse, rho=rho)
    diff_beta <- norm(beta_old-beta_new,"F")/norm(beta_old,"F")
    
    # Estimate each Gj, gamma and sigma2
    G_old <- G_new
    gamma_old <- gamma_new
    sigma2_old <- sigma2_new
    G_gamma_sigma <- Renew_G_gamma_sigma(Y_resid=Y_resid, C=C_new, beta=beta_new, gamma=gamma_old, sigma2=sigma2_old, knots=knots, m=m, p=p, Nsubj=Nsubj, L=L, J=J)
    G_new <- G_gamma_sigma$G
    gamma_new <- G_gamma_sigma$gamma
    sigma2_new <- G_gamma_sigma$sigma2
    diff_G <- norm(G_new-G_old,"F")/norm(G_old,"F")
    diff_gamma <- norm(gamma_old-gamma_new,"F")/norm(gamma_old,"F")
    diff_sigma2 <- norm(sigma2_old-sigma2_new,"F")/norm(sigma2_old,"F")
    
    k = k + 1 
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
  # Estimated gamma
  gamma_hat <- gamma_new
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
  eigens_G <- eigen(G_hat0)
  lambda_G0 <- (eigens_G$value/(t-1))[which(eigens_G$value>0)]
  idx <- which(cumsum(lambda_G0/sum(lambda_G0))>pve2)[1]
  lambda_G <- lambda_G0[1:idx]
  eigenvec <- (eigens_G)$vectors * sqrt(t-1)
  psi <- as.matrix(eigenvec[,1:length(lambda_G)])
  G_hat <- Reduce("+", lapply(1:length(lambda_G), function(y){lambda_G[y]*psi[,y]%*%t(psi[,y])}))
  ######################################
  ####step 4: Estimate scores
  ######################################
  Yhat <- NULL
  xi <- NULL
  zeta <- NULL
  
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
  
  if(fit){
    Mat1 <- do.call(cbind, lapply(1:L,function(x){kronecker(beta_hat[,x], phi_temp[[x]])}))
    Mat2 <- do.call(cbind, lapply(1:L,function(x){kronecker(beta_hat[,x]/sigma2_hat, phi_temp[[x]])}))
    Mat3 <- lapply(1:J, function(x){gamma_hat[x]^2*t(psi_temp)%*%psi_temp/sigma2_hat[x]})
    Mat4 <- lapply(1:J, function(x){gamma_hat[x]*t(psi_temp)/sigma2_hat[x]})
    invGamma1 <- bdiag(lapply(lambda_C_temp,function(x){if(length(x)>1) diag(1/x) else 1/x}))
    invGamma2 <- bdiag(lapply(1:J,function(x){if(length(lambda_G_temp)>1) diag(1/lambda_G_temp) else 1/lambda_G_temp}))
    A <- t(Mat1) %*% Mat2 + invGamma1
    C <- bdiag(Mat4) %*% Mat1
    D <- bdiag(Mat3) + invGamma2
    invD <- solve(D)
    MatE <- solve(A-t(C)%*%invD%*%C)
    MatF <- - MatE%*%t(C)%*%invD
    MatG <- - invD%*%C%*%MatE
    MatH <- invD - invD%*%C%*%MatF
    part_xi <- MatE%*%t(Mat2) + MatF%*%bdiag(Mat4)
    part_zeta <- MatG%*%t(Mat2) + MatH%*%bdiag(Mat4)
    y_vec <- do.call(cbind, lapply(Yi, function(x){as.vector(x)}))
    xi <- part_xi %*% y_vec
    zeta <- part_zeta %*% y_vec
    
    part1 <- Mat1%*%xi
    part2 <- bdiag(lapply(1:J, function(x){gamma_hat[x]*psi_temp})) %*% zeta
    yhat <- matrix(part1+part2, nc=t, byrow=T) + kronecker(rep(1,Nsubj),mu)
    idx <- rep(1:J, Nsubj)
    Yhat <- lapply(1:J, function(x){yhat[idx==x,]})
  }
  
  # delete scores corresponding to lambda=NULL
  idxC <- which(unlist(lambda_C_temp)==1e-6)
  if(length(idxC) > 0)  xi=xi[-idxC,]
  
  
  return(list(Yhat=Yhat, meanfunctions=mu, xi=xi, zeta=zeta, C=C_hat, beta=beta_hat,
              G=G_hat, gamma=gamma_hat, common_G = TRUE, sigma2=sigma2_hat, 
              lambda_C=lambda_C, phi=phi, lambda_G=lambda_G, psi=psi, k=k))
}





