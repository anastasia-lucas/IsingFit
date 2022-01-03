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
