ComparableNetworks <- function(x1, x2, strong.restr = FALSE, AND=TRUE, gamma=.25,
                               plot = TRUE, progressbar = TRUE){
  ifelse(nrow(x1) < nrow(x2), {ss <- x1 ; ll <- x2}, 
{ss <- x2 ; ll <- x1})

p <- ncol(ss)
n <- nrow(ss)
res.ss <- IsingFit(ss, lowerbound.lambda = NA, AND = AND, gamma = gamma,
                   plot = plot, progressbar = progressbar)
lowerbound.lambda <- min(res.ss$lambda.values)
if(strong.restr == TRUE) {
  lowerbound.lambda <- sqrt(log(p)/n)
  res.ssstr <- IsingFit(ss,lowerbound.lambda=lowerbound.lambda,AND=AND,gamma=gamma,
                        plot=plot,progressbar=progressbar)
}
if(strong.restr == FALSE) res.ssstr = NULL
res.ll <- IsingFit(ll,lowerbound.lambda=lowerbound.lambda,AND=AND,gamma=gamma,
                   plot=plot,progressbar=progressbar)
res <- list(res.smallest = res.ss, res.smalleststr = res.ssstr,
                                      res.largest = res.ll)
return(res)
}