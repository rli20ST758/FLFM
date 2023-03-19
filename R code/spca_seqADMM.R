########################################
LA_seqadmm <- function(H20=NULL,S,PrevPi_Eig=NULL,PrevPi_d=NULL,eta=NULL,rho=NULL,t,K,etastep,eps=1e-2){
  
  norm1 <- rep(0,K-2)
  norm2 <- rep(0,K-2)
  for (i in 1:K){
    
    H <- seqADMM_perstage(H20=H20,S=S,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,eta=eta,rho=rho,t=t)    
    H1 <- H[[1]]
    H2 <- H[[2]]
    W <- H[[3]]
    eta <- eta*etastep
    
    if(i>2){
      # stopping criterion
      norm1[i-2] <- (norm(H1-H2,'F'))^2 #primal residual
      norm2[i-2] <- (norm(H2-H20,'F'))^2  #dual residual
      maxnorm <- max(norm1[i-2],norm2[i-2])
      if(maxnorm<eps^2){
        break
      }
    }
    # reset parameter
    H20 <- H2 #store primal variable
  }
  
  if(i==K){
    return(-1) #not converge
  } else{
    return(H2)
  }
}



########################################
########################################
seqADMM_perstage <- function(H20=NULL,S,PrevPi_Eig=NULL,PrevPi_d=NULL,eta=NULL,rho=NULL,t){
  
  p <- dim(S)[1]  #dimension of sample covariance matrix  
  #starting value
  niter=0
  #if no initial value of H20 is given, then the default value is zero.
  if(is.null(H20)==TRUE){
    H20 <- matrix(0,p,p)
  }
  H2 <- matrix(0,p,p) ## used to store  value for H2 at the current step
  W0 <- matrix(0,p,p) ##intial value for dual variable , used to store value of dual variable at the previous step
  W <- matrix(0,p,p)  ##used to store value of dual variable at the current step
  H1_sum <- matrix(0,p,p) ## used to store the sum of all H1 in each iteration of this stage
  H2_sum <- matrix(0,p,p) ## used to store the sum of all H2 in each iteration of this stage
  W_sum <- matrix(0,p,p) ## used to store the sum of all dual variable W in each iteration of this stage
  
  # iteration
  while(niter<t){
    # update primal variable
    # update H_1
    mat <- H20 - (1/eta)*(W0-S)
    mat <- (mat+t(mat))/2 #force symmetric
    H1 <- FantopeProj(mat=mat,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d)
    # update H_2
    H2 <- H2Prox(W0,H1,eta,rho)                                     
    # update dual variable W
    W <- W0 + eta*(H1-H2)
    
    # reset H20,W0
    H20 <- H2
    W0 <- W
    # update sum
    H1_sum <- H1_sum + H1
    H2_sum <- H2_sum + H2
    W_sum <- W_sum + W
    niter <- niter+1
  }
  return(list(H1_sum/t,H2_sum/t,W_sum/t))
}


########################################
########################################
FantopeProj <- function(mat,PrevPi_Eig=NULL,PrevPi_d=NULL){
  
  # PrevPi==NULL, then do fantope projection, otherwise, do deflated fantope projection
  if(is.null(PrevPi_Eig)==FALSE){

    U <- PrevPi_Eig[[2]]
    p <- dim(U)[1]
    U <- U[,(1:(p-PrevPi_d))] #U: the orthogonal complement basis of PrevPi
    mat <- eigenMapMatMult3(t(U),mat,U)
    mat <- (t(mat)+mat)/2 #force to be symmetric
  }
  
  U_eig <- eigs_sym(A=mat,k=10)
  U_value <- U_eig$values
  U_vec <- U_eig$vectors
  
  theta <- GetTheta(U_value) #get theta by solving corresponding piecewise linear equation
  
  # compute the "eigenvalues" of the projection
  new.values <- unlist(sapply(U_value,FUN=function(u){
    min(max(u-theta,0),1)
  }))
  newmat <- eigenMapMatMult3(U_vec,diag(new.values),t(U_vec))
  
  if(is.null(PrevPi_Eig)==FALSE){
    newmat <- eigenMapMatMult3(U,newmat,t(U))
  }
  
  return(newmat)
}

########################################
########################################
H2Prox <- function(W,H1,eta,rho){
  
  p <- dim(W)[1] 
  M <- (1/eta)*W + H1
  A <- matrix(0,p,p)
  
  # first step:elementwise soft thresholding
  for(i in 1:p){
    for(j in 1:p){
      A[i,j] <- SoftThreshold(M[i,j],rho/eta)
    }
  }
  return(A)
}





##############################################################
# Using interpolation to solve the piecewise linear equation
##############################################################
GetTheta <- function(r){
  # concatenation vector r-1 and r
  knots <- as.numeric(as.character((as.data.frame(table(c(r-1,r)))[,1])))
  index <- find_first(r,knots,simplex_sum) #return the index of the right endpoint 
  ia <- knots[index] #right endpoint
  ib <- knots[index-1]
  fa <- simplex_sum(r,ia)
  fb <- simplex_sum(r,ib)
  theta <- ia + (ib-ia)*(1-fa)/(fb-fa) #interpolation
  return(theta)
}

################################################################
# Find the right endpoint of the interval which contains theta
################################################################
find_first <- function(r,knots,f){
  nk <- length(knots)
  for (i in 1:nk){
    if(f(r,knots[i])<1){
      return(i)
      break
    }
  }
}

########################################
########################################
simplex_sum <- function(r,theta){
  vec <- sapply(r,FUN=function(u){
    min(max(u-theta,0),1)
  })
  return(sum(vec))
}

########################################
########################################
SoftThreshold <- function(x,lambda){
  value <- sign(x)*max(abs(x)-lambda,0)
  return(value)
}


