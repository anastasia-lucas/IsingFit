#' Run glmnet and capture summary features
#'
#' @description This is an internal function to run glmnet and extracts relevant model 
#' features—intercepts, betas, and lambdas
#'
#' @param x x as in glmnet model
#' @param y y as in glmnet model
#' @param family family as in glmnet model
#' @return list
#' @noRd

run_glmnet <- function(x, y, family){
  l <- list()
  a <- glmnet::glmnet(x = x, y = y, family = family)
  l$intercepts <- a$a0
  l$betas <- a$beta
  l$lambdas <- a$lambda
  l$nlambdas <- length(l$lambdas)
  return(l)
}

#' Get counts of beta != 0 to create J
#'
#' @description This is an internal function to create matrix J
#'
#' @param x dgCMatrix
#' @param length length of list (= maxlambdas)
#' @return vector
#' @noRd

get_counts <- function(x, length){
  j <- matrix(0, length, 1)
  j[1:ncol(x), 1] <- colSums(x!=0)
  return(j)
}

#' Create column of P matrix
#'
#' @description This is an internal function to create the P matrix
#'
#' @param betas betas
#' @param intercepts intercepts
#' @param N nrow(x)
#' @param x x
#' @maxlambdas maxlambdas
#' @return vector
#' @noRd

construct_P <- function(betas, intercepts, N, x, maxlambdas){
    betas <- as.matrix(betas)
    y <- matrix(0, nrow=N , ncol= ncol(betas))
    xi <- x[,-i]
    NB <- nrow(betas) # number of rows in beta
    for (bb in 1:NB){   # bb <- 1
      y <- y + betas[rep(bb,N),] * xi[,bb]
    }
    y <- matrix(intercepts, nrow=N, ncol=ncol(y) , byrow=TRUE ) + y
    # number of NAs
    n_NA <- maxlambdas-ncol(y)
    if (n_NA > 0 ){
      for ( vv in 1:n_NA){
        y <- cbind( y , NA )
      }
    }
    # calculate P matrix
    sum_logl_M <- colSums(log(exp(y*x[,i])/(1+exp(y))))
    
    # need to return the column sums & then can combine
    return(sum_logl_M)
}