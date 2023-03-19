
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
  #sparse -- Logical with default FALSE. If TRUE, beta are sparse
  #rho -- the penalty parameter for sparse pca, only useful when sparse is TRUE. 
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
        #S <- t(A)%*%A / mean(diag(t(A)%*%A))  
        S <- -A / mean(abs(diag(A)))
        if(l==1) { # sparse beta
          H <- LA_seqadmm(S=S,eta=100, rho=rho, t=10, K=20, etastep=2)
          beta_temp <- getEigenDecomp(H)[[2]][,ncol(H)]
          beta[,l] <- sign(beta_temp[1]+1e-6) * beta_temp
        } else {
          PrevPi <- beta[,1:(l-1)]%*%t(beta[,1:(l-1)])
          PrevPi_Eig <- getEigenDecomp(diag(1,J)-PrevPi)
          PrevPi_Eig[[1]] <- rev(PrevPi_Eig[[1]])
          PrevPi_Eig[[2]] <- PrevPi_Eig[[2]][,ncol(PrevPi_Eig[[2]]):1]
          H <- LA_seqadmm(S=S, PrevPi_Eig=PrevPi_Eig, PrevPi_d=l-1, eta=100, rho=rho, t=10, K=20, etastep=2)
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


##########################################################
# Renew_beta2() for common G for all functional outcomes
##########################################################
Renew_beta2 <- function(Yi, C, beta, G, gamma, sigma2, Nsubj, L, J, sparse, rho){
  
  if(L==1){
    
    # Compute the sample covairance matrix for updating beta
    Yi_C_Yi <- lapply(1:Nsubj, function(x) {
      temp <- t(Yi[[x]]) %*% C[[1]]
      temp %*% Yi[[x]]})
    Y_C_Y <- Reduce('+', Yi_C_Yi)
    # Variance caused by self part
    C_Gj <- gamma^2 * sum(C[[1]]*G)   
    self <- diag(C_Gj)
    # Delete noise for diagnal elements
    noise  <- sum(diag(C[[1]])) * diag(sigma2)
    
    # Update beta 
    A <- self + noise - Y_C_Y/Nsubj 
    if(sparse==TRUE){
      S <- t(A)%*%A / mean(diag(t(A)%*%A))
      H <- LA_seqadmm(S=S, eta=100, rho=rho, t=10, K=20, etastep=2) #sparse pca
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
      Cl_Gj <- gamma^2 * sum(C[[l]]*G)  
      self <- diag(Cl_Gj)
      # Delete noise for diagnal elements
      noise <- sum(diag(C[[l]])) * diag(sigma2)
      
      # Update beta 
      A <- Cl_beta + self + noise - Y_Cl_Y/Nsubj 
      if(sparse==TRUE){
        #S <- t(A)%*%A / mean(diag(t(A)%*%A))  
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
          H <- LA_seqadmm(S=S, PrevPi_Eig=PrevPi_Eig, PrevPi_d=l-1, eta=100, rho=rho, t=10, K=20, etastep=2)
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





