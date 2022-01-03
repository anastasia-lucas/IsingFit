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
  j <- matrix(0, maxlambdas, 1)
  j[1:ncol(x), 1] <- colSums(x!=0)
  return(j)
}
