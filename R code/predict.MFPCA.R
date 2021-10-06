
predict.MFPCA <- function(object, data, newdata, argvals, pve){
  ############################################################ 
  #Arguments:
  #object -- an object returned by LFM
  #score vectors. It is used for calculating predictions
  #data -- observed data for fit model
  #newdata -- a data list consisting for which predicted values are desired. 
  #The list of length J contains densely observed multivariate (J-dimensional) functional data.
  #argvals -- the argument values of the function evaluations
  #pve <- the proportion of variance explained, defu=ault to 0.95
  ############################################################ 

  J <- length(data)
  N_new <- nrow(newdata[[1]])
  vectors <- object$vectors
  eigenfucntions <- object$functions
  meanfucntions <- object$meanFunction

  
  temp0 <- lapply(1:J, function(x){funData(argvals=argvals, X = data[[x]])})
  sim <- multiFunData(temp0)
  temp1 <- lapply(1:J, function(x){funData(argvals=argvals, X = newdata[[x]])})
  sim_new <- multiFunData(temp1)

  # Get predicted scores for newdata
  scores <- c()
  for (j in 1:J) {
    u_fit <- PACE(sim@.Data[[j]], predData=sim_new@.Data[[j]], pve=pve)
    scores <- cbind(scores, u_fit$scores)
  }
  scores_MFPCA <- scores %*% vectors

  
  # Get fitted values for newdata
  Ypred <- list()
  for (j in 1:J) {
    Ypred[[j]] <- scores_MFPCA %*% eigenfucntions[[j]]@X + kronecker(rep(1,N_new), meanfucntions[[j]]@X)
  }
  
  return(list(Ypred=Ypred, scores=scores_MFPCA))
  
}












                   