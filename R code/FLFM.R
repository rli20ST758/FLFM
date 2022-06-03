
LFM <- function(Y, Nsubj, argvals, L, sparse=FALSE, rho=0, pve1=0.95, pve2=0.95,
                initial=list(C=NULL,beta=NULL,G=NULL,sigma2=NULL),
                knots=35, p=3, m=2, tol=1e-4, Nmax=200, fit=FALSE) {
  ############################################################ 
  #Arguments:
  #argvals -- the argument values of the function evaluations 
  #Y -- observed data. A list of length J containing densely observed multivariate (J-dimensional) functional data. 
  #Y[[j]] is an (Nsubj x length(argvals)) matrix of functional data for Nsubj subjects observed on argvals.
  #Nsubj -- number of subjects
  #L -- number of reduced dimensions
  #sparse -- whether to update beta by sparse pca with default FALSE. 
  #rho -- the penalty parameter for sparse pca, only useful when sparse is TRUE. 
  #pve1 -- the proportion of variance explained for joint terms, defu=ault to 0.95. 
  #pve2 -- the proportion of variance explained for outcome-specific terms, defu=ault to 0.95. 
  #initial --  the initial values of parameters for the iteration. It can be from cross-validation for L. 
  #            For sparse=TRUE, it highly recommends to set initial values as the estimation for sparse=FALSE.
  #knots -- list of two vectors of knots or number of equidistant knots for all dimensions; defaults to 35
  #p -- degrees of B-splines; defaults to 3; details see fbps
  #m -- order of differencing penalty; defaults to 2; details see fbps
  #tol -- tolerance for iteration, default to 10e-4
  #Nmax -- the number of maximal iteration, default to 200
  #fit -- Logical with default FALSE. If TRUE, get the fitted values and scores
  ############################################################ 
  
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
  # Get eigen components of G by pve2
  # eigens_G <- lapply(G_hat0, function(x){eigen(x)})
  # lambda_G0 <- lapply(eigens_G, function(x){(x$value/(t-1))[which(x$value>0)]})
  # order_lambda_G <- sort(unlist(lambda_G0), decreasing=TRUE) #order all lambda_G
  # lambda_G_pve <- order_lambda_G/sum(order_lambda_G)
  # lambda_G_star <- order_lambda_G[which(cumsum(lambda_G_pve)>pve2)[1]]
  # lambda_G <- lapply(lambda_G0,function(x){ #eigenvalues
  #   num <- sum(x>=lambda_G_star)
  #   if(num==0){return(NULL)}
  #   return(x[1:num])})
  # psi <- lapply(1:J, function(x){ #eigenvectors
  #   eigenvec <- (eigens_G[[x]])$vectors * sqrt(t-1)
  #   if(length(lambda_G[[x]])==0) NULL
  #   else as.matrix(eigenvec[,1:length(lambda_G[[x]])]) })
  # G_hat <- lapply(1:J, function(x){
  #   num <- length(lambda_G[[x]])
  #   if(num==0){return(NULL)}
  #   if(num==1){return(lambda_G[[x]]*psi[[x]]%*%t(psi[[x]]))}
  #   return(Reduce("+", lapply(1:num, function(y){lambda_G[[x]][y]*psi[[x]][,y]%*%t(psi[[x]][,y])})))
  # })

  
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
  
  
  return(list(Yhat=Yhat, meanfunctions=mu, xi=xi, zeta=zeta, C=C_hat, beta=beta_hat, G=G_hat, 
              sigma2=sigma2_hat, lambda_C=lambda_C, phi=phi, lambda_G=lambda_G, psi=psi, k=k))
}




# If do projection to get fitted values and scores
# if(fit){
#   # Assign global eigenfunctions
#   efuncs1 <- list()
#   index <- 0
#   for(l in 1:L){
#     #efuncs1[[l]] <- kronecker(beta_hat[,l], phi[[l]])
#     efuncs1[[l]] <- kronecker(rep(1,J), phi[[l]])
#     for(m in 1:ncol(efuncs1[[l]])){
#       index <- index + 1
#       assign(paste("Efunc1", index, sep=""), efuncs1[[l]][,m], envir=parent.frame()) #make it global variable
#     }
#   }
#   
#   max_n <- max(unlist(lapply(lambda_G, function(x)length(x))))
#   efuncs2 <- do.call(rbind, lapply(1:J, function(x){
#     num <- length(lambda_G[[x]])
#     if(num<max_n){return(cbind(psi[[x]],matrix(0,nr=t,nc=(max_n-num))))}
#     return(psi[[x]])
#   }))
#   for(r in 1:max_n){
#     assign(paste("Efunc2",r,sep=""), efuncs2[,r], envir = parent.frame()) #make it global variable
#   }
#   
#   
#   yhat <- list()
#   xi <- matrix(NA, nrow=Nsubj, ncol=index)
#   eta0 <- matrix(NA, nrow=Nsubj, ncol=max_n*J)
#   g_var <- rep(as.factor(paste(1:J,sep="")),each=t) #indicator for each subjects's variables
#   #phimodel <- paste("s(Efunc1",1:index,", bs='re') ", collapse="+", sep="")
#   phimodel <- paste("s(g_var, by=Efunc1",1:index,", bs='re') ", collapse="+", sep="")
#   psimodel <- paste("s(g_var, by=Efunc2",1:max_n,", bs='re') ", collapse="+", sep="")
#   flags <- rep(1:L, unlist(lapply(1:L,function(x){length(lambda_C[[x]])})))
#   for (i in 1:Nsubj) {
#     # Assign model
#     y <- as.vector(Yi[[i]])
#     formula <- paste("y ~ ", phimodel, "+", psimodel)
#     # Fit model and estimate scores
#     fitting <- bam(as.formula(formula), drop.intercept=T)
#     #xi[i,] <- fitting$coefficients[1:index]
#     #eta0[i,] <- fitting$coefficients[-(1:index)]
#     xi[i,] <- unlist(lapply(1:index, function(x){
#       temp <- fitting$coefficients[((x-1)*J+1):(x*J)]
#       sum(temp)/sum(beta_hat[,flags[x]])
#       }))
#     eta0[i,] <- fitting$coefficients[-(1:(index*J))]
#     
#     yhat[[i]] <- matrix(fitting$fitted.values, nc=t, byrow=T) + mu
#   }
#   
#   Yhat <- lapply(1:J, function(x){do.call("rbind",lapply(1:Nsubj, function(sub){yhat[[sub]][x,]}))})
#   eta <- lapply(1:J, function(x){
#     start <- max_n*(x-1)+1
#     end <- max_n*(x-1) + length(lambda_G[[x]])
#     eta0[,start:end]})
# }