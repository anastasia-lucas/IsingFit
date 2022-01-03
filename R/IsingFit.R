IsingFit <-
  function(x, family='binomial', AND = TRUE, gamma = 0.25, plot = TRUE, 
           progressbar = TRUE, ncores = 1, lowerbound.lambda = NA,...){
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
    cores = min(ncores, 
                max(1, floor(parallel::detectCores()*0.75)))
    if(cores < ncores){ 
      print(paste("Warning: using", cores, "cores due to resource availability"))
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
    
    # for (i in 1:nvar)
    # {
    #   # TODO requires Matrix library to be loaded why ???
    #   J[1:ncol(betas[[i]]),i] <- colSums(betas[[i]]!=0)
    # }
    
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
    
    # TODO continue
    
    logl_M <- P_M <- array(0, dim=c(nrow(x),maxlambdas, nvar) )
    N <- nrow(x)
    for (i in 1:nvar){  # i <- 1
      betas.ii <- as.matrix( betas[[i]] )
      int.ii <- intercepts[[i]]
      y <- matrix( 0 , nrow=N , ncol= ncol(betas.ii) ) 
      xi <- x[,-i]
      NB <- nrow( betas.ii) # number of rows in beta
      for (bb in 1:NB){   # bb <- 1
        y <- y + betas.ii[rep(bb,N),] * xi[,bb]
      }
      y <- matrix( int.ii , nrow=N , ncol=ncol(y) , byrow=TRUE ) + y
      # number of NAs
      n_NA <- maxlambdas-ncol(y)
      if (n_NA > 0 ){ 
        for ( vv in 1:n_NA){ 
          y <- cbind( y , NA ) 
        } 
      }
      # calculate P matrix
      P_M[,,i] <- exp(y*x[,i])/(1+exp(y))
      logl_M[,,i] <- log(P_M[,,i])  
      if (progressbar==TRUE) setTxtProgressBar(pb, i)
    }
    
    logl_Msum <- colSums( logl_M , 1, na.rm=FALSE )
    if (progressbar==TRUE) close(pb)
    sumlogl <- logl_Msum 
    sumlogl[sumlogl==0]=NA
    penalty <- J * log(nrow(x)) + 2 * gamma * J * log(p)
    EBIC <- -2 * sumlogl + penalty
    
    lambda.mat <- matrix(NA,nrow(EBIC),ncol(EBIC))
    for (i in 1:nvar){
      lambda.mat[,i] <- c(lambdas[[i]],rep(NA,nrow(EBIC)-length(lambdas[[i]])))
    }
    
    if(!is.na(lowerbound.lambda)){
      EBIC <- EBIC/(lambda.mat>=lowerbound.lambda)*1
    }
    
    lambda.opt <- apply(EBIC,2,which.min)
    lambda.val <- rep(NA,nvar)
    thresholds <- 0
    for(i in 1:length(lambda.opt)){
      lambda.val[i] <- lambda.mat[lambda.opt[i],i]
      thresholds[i] <- intercepts[[i]][lambda.opt[i]]
    }
    weights.opt <- matrix(,nvar,nvar)
    for (i in 1:nvar){
      weights.opt[i,-i] <- betas[[i]][,lambda.opt[i]]
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

