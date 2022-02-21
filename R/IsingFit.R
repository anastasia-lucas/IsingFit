#' Run Ising model
#'
#' @description This network estimation procedure eLasso, which is based on the
#' Ising model, combines l1-regularized logistic regression with model selection
#' based on the Extended Bayesian Information Criterion (EBIC). EBIC is a fit
#' measure that identifies relevant relationships between variables. The
#' resulting network consists of variables as nodes and relevant relationships
#' as edges. Can deal with binary data.
#'
#' @usage IsingFit(x, family='binomial', AND = TRUE, gamma = 0.25, plot = TRUE,
#'                 ncores = 1, forcecores = FALSE, lowerbound.lambda = NA,...)
#'
#' @param x Input matrix. The dimension of the matrix is nobs x nvars; each row
#' is a vector of observations of the variables. Must be cross-sectional data.
#' @param family The default is 'binomial', treating the data as binary.
#' Currently, this procedure is only supported for binary data.
#' @param AND Logical. Can be TRUE of FALSE to indicate whether the AND-rule or
#' the OR-rule should be used to define the edges in the network. Defaults to TRUE.
#' @param gamma A value of hyperparameter gamma in the extended BIC. Can be
#' anything between 0 and 1. Defaults to .25.
#' @param plot Logical. Should the resulting network be plotted?
#' @param progressbar Logical. Should the pbar be plotted in order to see the
#' progress of the estimation procedure?
#' @param lowerbound.lambda The minimum value of tuning parameter lambda
#' (regularization parameter). Can be used to compare networks that are based on
#' different sample sizes. The lowerbound.lambda is based on the number of
#' observations in the smallest group n: sqrt(log(p)/n). p is the number of
#' variables, that should be the same in both groups. When both networks are
#' estimated with the same lowerbound for lambda (based on the smallest group),
#' the two networks can be directly compared.
#' @param ncores number of cores to use for distributed computing, default: 1
#' @param forcecores logical, allows IsingFit to use more than 75% of all
#' currently available cores, default: FALSE
#' @param \dots Arguments sent to \code{qgraph}.
#'
#' @return IsingFit returns (invisibly) a 'IsingFit' object that contains the
#' following items:
#' \item{weiadj }{The weighted adjacency matrix.}
#' \item{thresholds }{Thresholds of the variables.}
#' \item{q }{The object that is returned by qgraph (class 'qgraph').}
#' \item{gamma }{The value of hyperparameter gamma.}
#' \item{AND }{A logical indicating whether the AND-rule is used or not. If not,
#' the OR-rule is used.}
#' \item{time }{The time it took to estimate the network.}
#' \item{asymm.weights }{The (asymmetrical) weighted adjacency matrix before
#' applying the AND/OR rule.}
#' \item{lambda.values }{The values of the tuning parameter per node that
#' ensured the best fitting set of neighbors.}
#'
#' @references
#' \itemize{
#'  \item Chen, J., & Chen, Z. (2008). Extended bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.
#'  \item Foygel, R., & Drton, M. (2011). Bayesian model choice and information criteria in sparse generalized linear models. arXiv preprint arXiv:1112.5635.
#'  \item Ravikumar, P., Wainwright, M. J., & Lafferty, J. D. (2010). High-dimensional Ising model selection using l1-regularized logistic regression. The Annals of Statistics, 38, 1287 - 1319.
#'  \item van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for constructing networks from binary data. Scientific Reports 4, 5918; DOI:10.1038/srep05918.
#' }
#'
#' @examples
#' ### Simulate dataset ###
#' # Input:
#' N <- 6 # Number of nodes
#' nSample <- 1000 # Number of samples
#'
#' # Ising parameters:
#' set.seed(123)
#' Graph <- matrix(sample(0:1,N^2,TRUE,prob = c(0.8, 0.2)),N,N) * runif(N^2,0.5,2)
#' Graph <- pmax(Graph,t(Graph))
#' diag(Graph) <- 0
#' Thresh <- -rowSums(Graph) / 2
#'
#' # Simulate:
#' set.seed(123)
#' Data <- IsingSampler::IsingSampler(nSample, Graph, Thresh)
#'
#' ### Fit using IsingFit ###
#' Res <- IsingFit(Data, family='binomial', plot=FALSE)
#' @importFrom foreach %dopar%

IsingFit <-
  # TODO: add progressbar back in
  function(x, family='binomial', AND = TRUE, gamma = 0.25, plot = TRUE,
           progressbar = TRUE, ncores = 1, forcecores = FALSE,
           lowerbound.lambda = NA,...){
    t0 <- Sys.time()

    # set names
    xNames <- colnames(x)

    if(is.null(xNames)) {
      xNames <- rep(paste0("X", 1:ncol(x)))
    }

    if (family!='binomial')
      stop ("This procedure is currently only supported for binary (family='binomial') data")

    ###### Check to prevent error of lognet() in package glmnet ######
    # I think all of the checklognets() function can be replaced by the following
    NodesToAnalyze <- apply(x, 2, function (x) length(unique(x))) > 1 # require >= 2 unique values

    names(NodesToAnalyze) <- xNames
    if (!any(NodesToAnalyze)) stop("No variance in dataset")
    if (any(!NodesToAnalyze)) {
      warning(paste("Nodes with too little variance (not allowed):",paste(colnames(x)[!NodesToAnalyze],collapse = ", ")))
    }
    ######## End check ######

    x <- as.matrix(x)
    allthemeans <- colMeans(x)
    x <- x[,NodesToAnalyze,drop=FALSE]
    nvar <- ncol(x)
    # p <- nvar - 1
    # intercepts <- betas <- lambdas <- list(vector,nvar)
    # nlambdas <- rep(0,nvar)

    # Run glmnet on all of the columns
    if(isTRUE(forcecores)){
      cores = ncores
    } else {
      cores = min(ncores,
                  max(1, floor(parallel::detectCores()*0.75)))
      if(cores < ncores){
        print(paste("Warning: using", cores, "cores due to resource availability"))
      }
    }

    doParallel::registerDoParallel(cores = cores)

    glmnetRes <- foreach::foreach(i = 1:nvar, .packages = c('foreach',
                                                            'parallel',
                                                            'doParallel')) %dopar% {
      # TODO %do% is not being exported by foreach package here
      IsingFit:::run_glmnet(x[,-i], x[,i], family = family)
    }


    maxlambdas <- max(unlist(lapply(glmnetRes, function (x) x$nlambdas)),
                      na.rm = TRUE)

    if (progressbar==TRUE) pb <- txtProgressBar(max=nvar, style = 3)

    P <- logl <- sumlogl <- J <- matrix(0, maxlambdas, nvar)

    # TODO check to make sure this should be getting N != 0 and not sum(betas!=0)
    doParallel::registerDoParallel(cores = cores)

    J <- foreach::foreach(i = 1:nvar, .packages = c('foreach',
                                                    'parallel',
                                                    'doParallel',
                                                    'Matrix')) %dopar% {
                                                              # TODO %do% is not being exported by foreach package here
                                                              IsingFit:::get_counts(glmnetRes[[i]]$betas,
                                                                                    length = maxlambdas)
                                                            }

    J <- do.call(cbind, J)

    # construct p matrix
    doParallel::registerDoParallel(cores = cores)

    P_M <- foreach::foreach(i = 1:nvar, .packages = c('foreach',
                                                      'parallel',
                                                      'doParallel',
                                                      'Matrix')) %dopar% {
                                                      # TODO %do% is not being exported by foreach package here
                                                      IsingFit:::construct_P(betas = glmnetRes[[i]]$betas,
                                                                             intercepts = glmnetRes[[i]]$intercepts,
                                                                             N = nrow(x),
                                                                             x = x,
                                                                             i = i,
                                                                             maxlambdas = maxlambdas)
                                                    }

    sumlogl <- do.call(cbind, P_M)

    sumlogl[sumlogl==0]=NA
    penalty <- J * log(nrow(x)) + 2 * gamma * J * log(nvar - 1)
    EBIC <- -2 * sumlogl + penalty

    lambda.mat <- matrix(NA,nrow(EBIC),ncol(EBIC))
    for (i in 1:nvar){
      lambda.mat[,i] <- c(glmnetRes[[i]]$lambdas,rep(NA,nrow(EBIC)-length(glmnetRes[[i]]$lambdas)))
    }

    if(!is.na(lowerbound.lambda)){
      EBIC <- EBIC/(lambda.mat>=lowerbound.lambda)*1
    }

    lambda.opt <- apply(EBIC,2,which.min)
    lambda.val <- rep(NA,nvar)
    thresholds <- 0
    for(i in 1:length(lambda.opt)){
      lambda.val[i] <- lambda.mat[lambda.opt[i],i]
      thresholds[i] <- glmnetRes[[i]]$intercepts[lambda.opt[i]]
    }
    weights.opt <- matrix(,nvar,nvar)
    for (i in 1:nvar){
      weights.opt[i,-i] <- glmnetRes[[i]]$betas[,lambda.opt[i]]
    }
    asymm.weights <- weights.opt
    diag(asymm.weights)=0
    if (AND==TRUE) {
      adj <- weights.opt
      adj <- (adj!=0)*1
      EN.weights <- adj * t(adj)
      EN.weights <- EN.weights * weights.opt
      meanweights.opt <- (EN.weights+t(EN.weights))/2
      diag(meanweights.opt) <- 0
    } else {
      meanweights.opt <- (weights.opt+t(weights.opt))/2
      diag(meanweights.opt) <- 0
    }
    graphNew <- matrix(0,length(NodesToAnalyze),length(NodesToAnalyze))
    graphNew[NodesToAnalyze,NodesToAnalyze] <- meanweights.opt
    colnames(graphNew) <- rownames(graphNew) <- xNames
    threshNew <- ifelse(allthemeans > 0.5, -Inf, Inf)
    threshNew[NodesToAnalyze] <- thresholds
    if (plot==TRUE) notplot=FALSE else notplot=TRUE
    q <- qgraph(graphNew,layout='spring',labels=names(NodesToAnalyze),DoNotPlot=notplot,...)
    Res <- list(weiadj = graphNew, thresholds = threshNew, q = q, gamma = gamma,
                AND = AND, time = Sys.time() - t0, asymm.weights = asymm.weights,
                lambda.values = lambda.val)
    class(Res) <- "IsingFit"
    return(Res)
  }

## Methods:
plot.IsingFit <- function(object,...) qgraph(object$q,DoNotPlot = FALSE, ...)

print.IsingFit <- function(x)
{
  cat("Estimated network:\n")

  print(round(x$weiadj,2))

  cat("\n\nEstimated Thresholds:\n")

  print(x$thresholds)
}

summary.IsingFit <- function(object)
{
  cat("\tNetwork Density:\t\t", round(mean(object$weiadj[upper.tri(object$weiadj)]!=0),2),"\n",
      "Gamma:\t\t\t",round(object$gamma,2),"\n",
      "Rule used:\t\t",ifelse(object$AND,"And-rule","Or-rule"),"\n",
      "Analysis took:\t\t",format(object$time,format="%s"),"\n"
  )
}



#
# print(fit)
# plot(fit)
# summary(fit)
# export(fit)

