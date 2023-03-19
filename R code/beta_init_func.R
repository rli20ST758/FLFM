

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


# beta_init_func <- function(Y, argval, L){
#   ############################################################ 
#   #Arguments:
#   #Y <- Y_residual (observed Y - mean function)
#   #argval - number of points at which functions are observed
#   #L -- number of reduced dimensions
#   ############################################################ 

#   J <- length(Y)
#   Nsubj <- nrown(Y[[1]])
#   argval <- ncol(Y[[1]])
#
#   cov <- lapply(1:J, function(x){  #compute C_jj
#     crossprod(Y[[x]]) / Nsubj
#   })
#   C0 <- Reduce('+', cov) #C0 is the sum of C_jj
#   
#   if(L==1){
#     # Compute beta square
#     beta.sq <- unlist(lapply(1:J, function(x){sum(cov[[x]])/sum(C0)}))
#     # Sign
#     cov.1j <- lapply(1:J, function(x){
#       t(Y[[1]]) %*% Y[[x]] / Nsubj
#       })
#     # Decide the sign of beta
#     beta.sign <- sign(unlist(lapply(1:J, function(x){sum(cov.1j[[x]])/ sum(C0)})))
#     # Get initial point of beta
#     beta.init <- sqrt(beta.sq) * beta.sign
#   } else{
#     eigendec <- eigen(C0)  #get eigenvalues and eigenfunctions of estimated corvariance
#     npc <- which(cumsum(eigendec$values)/sum(eigendec$values) > 0.95)[1] #numnber of principle components
#     phi <- eigendec$vectors[,1:npc]
#     
#     # lambda_{jjk} = t(phi_k) %*% C_{jj} %*% phi_k
#     lambda.jjk <- unlist(lapply(1:J, function(x){colSums( (t(cov[[x]]) %*% phi) * phi / argval^2  )}))
#     lambda.jjk <- matrix(lambda.jjk, nc=npc, nr=J, byrow=T)
#     lambda.colsum <- colSums(lambda.jjk)
#     
#     # Find the best group for phi (random 50 times) and estimate beta.sq
#     beta.sq <-matrix(NA, nrow=J, ncol=L)
#     residual <- 1e5
#     # Get the best eigenfunction group assigment
#     for(i in 1:50){
#       # Randomly assign eigenfunction to latent level Cl
#       # Make sure each Cl has its eigenfunctions
#       group.phi0 <- NULL
#       while(sum(1:L %in% group.phi0) < L){
#         group.phi0 <- as.integer(runif(npc, min = 1, max = L+1))
#       }
#       
#       beta.sq0 <- matrix(NA, nrow=J, ncol=L)
#       residual0 <- 0
#       # Solve the quadratic problem to get estimate of beta.square
#       for (j in 1:L) {
#         Y.data <- rowSums(as.matrix(lambda.jjk[ ,group.phi0==j]))
#         X.data <- sum(lambda.colsum[group.phi0==j]) * diag(J)
#         beta.sq0[,j] <- solve.QP(Dmat=crossprod(X.data), dvec=t(X.data)%*%Y.data, 
#                                  Amat=matrix(1, nrow=J, ncol=1), bvec=1, meq=T, factorized=FALSE)$solution
#         residual0 <- norm(Y.data - X.data %*% beta.sq0[,j], "2")/norm(Y.data,"2") + residual0
#       }
#       # Save the assigned group of phi with the minimal residual
#       if (residual0 < residual){
#         residual <- residual0
#         group.phi <- group.phi0
#         beta.sq <- beta.sq0
#       }
#     } #end i
#     
#     # Choose j0
#     j0 <- which.max(apply(beta.sq, 1, min)) 
#     # cov_{j0j}
#     cov.j0j <- lapply(1:J, function(x){
#       t(Y[[j0]]) %*% Y[[x]] / Nsubj
#     })
#     # lambda_{j0jk} = t(phi_k) %*% C_{j0j} %*% phi_k
#     lambda.j0jk <- unlist(lapply(1:J, function(x){colSums( (t(cov.j0j[[x]]) %*% phi) * phi / argval^2  )}))
#     lambda.j0jk <- matrix(lambda.j0jk, nc=npc, nr=J, byrow=T)
#     # Decide the sign of beta
#     beta.sign <- do.call(cbind, lapply(1:L, function(x){ sign(rowSums(as.matrix(lambda.j0jk[ ,group.phi==x]))) } ) )
#     # Change the sign relative t0 beta_{1l}
#     beta.sign <- beta.sign %*% diag(sign(beta.sign[1,]))
#     # Get initial point of beta
#     beta.init0 <- sqrt(beta.sq) * beta.sign
#     # Project to a orthogonal space by Gram-Schmidt Process
#     beta.init <- gramSchmidt(beta.init0)$Q
#     
#   }
#   return(beta.init)
# }
