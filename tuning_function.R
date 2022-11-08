#' Check loss function
#' @param t, the evaluate value 
#' @param tau, the tau-th quantile 
#' @return the check loss
rho.func <- function(t, tau){
    tau * ifelse(t>0, t, 0) - (1 - tau) * ifelse(t<0, t, 0)
}

#' This function is used for choosing tuning parameter 
#' 1. only smoothing is imposed, eta.n
#' 2. only functional MCP penalty is imposed, lambda 1
#' 3. impose the hierarchy penalty,that is, both lambda1 and lambda 2 are not equal to zero
#' @param method, specific the method to tuning, "OLS", “QR”
#' @param tuning.type, specific the tuning type, "validation"
#' @param dat, a list including the training and validation and testing set
#' @param Jmat, Jmat matrix
#' @param Psi, the augmented design matrix, n by (Ls*p+p)
#' @param y, the response, n by 1
#' @param tau, the tau-th quantile 
#' @param eta.n.vec, the eta.n for tuning, smoothing spline parameter 
#' @param lambda.1.vec, the lambda.1 for tuning
#' @param lambda.2.vec, the lambda.2 for tuning
#' @param beta.basis, the basis for beta
#' @param p, one plus the number of scalar covariates
#' @param V, the smoothness penality weight matrix
#' @param W, the proposed penality weight matrix
#' @param s, the regularized parameter for MCP
#' @param tol, the threshold for convergence
#' @param cutoff, the cutoff for small estimated theta's
#' @param max.iter, the max iteration times
#' @param epsilon, small value to control the MM precision
#' 
#' @return tuning results, a list including 
#' score: error score, 
#' eta.n: the selected best eta.n, 
#' lambda.1, the selected best lambda.1, 
#' lambda.2, the selected best lambda.1, 
#' min.score, the minmum error score

QR.OLS.tune <- function(method = "QR", tuning.type = "validation", 
                        dat, Jmat, Psi, y, tau, eta.n.vec, lambda.1.vec, lambda.2.vec, 
                        beta.basis, p, V, W, s, tol, cutoff, max.iter, epsilon)
{
    if((method %in% c("QR", "OLS")) == FALSE) print("undefined loss function type")
    if((tuning.type %in% c("validation")) == FALSE) print("undefined tuning type")
    
    #' manipulate the lambda and eta.n
    Ls    <- beta.basis$nbasis
    n.lam <- max(length(lambda.1.vec), length(lambda.2.vec))
    n.eta <- length(eta.n.vec)
    
    if(n.lam == 0)
    {
        lambda.1.vec <- 0
        lambda.2.vec <- 0
        n.lam        <- 1
    }
    
    #' calculate Psi and y on the validation data
    if(tuning.type == "validation")
    {
        X.valid <- crossprod(dat$valid$x$coefs, Jmat)
        Z.valid <- dat$valid$Z
        n.valid <- nrow(X.valid)
        
        U.valid <- matrix(0, nr = n.valid, nc = (p-1)*Ls)
        for(i in 1:n.valid) U.valid[i,] = kronecker(Z.valid[i,-1], X.valid[i,])
        
        Psi.valid <- cbind(X.valid, U.valid, Z.valid)
        y.valid   <- dat$valid$y
    }
    
    #' calculate the error score for all lambda and eta.n
    err.score = matrix(0, nr = n.lam, nc = n.eta)
    
    for(r in 1:n.lam)
    {
        for(u in 1:n.eta)
        {
            if(method == "QR")
            {
                fit.tuning <- QR.alg(Psi = Psi, y = y, tau = tau, eta.n = eta.n.vec[u], 
                                     lambda.1 = lambda.1.vec[r], lambda.2 = lambda.2.vec[r], 
                                     beta.basis = beta.basis, p = p, V = V, W = W, s = s, 
                                     tol = tol, cutoff = cutoff, max.iter = max.iter, epsilon = epsilon)
            }
            
            if(method == "OLS")
            {
                fit.tuning <- OLS.alg(Psi = Psi, y = y, eta.n = eta.n.vec[u], 
                                      lambda.1 = lambda.1.vec[r], lambda.2 = lambda.2.vec[r],
                                      beta.basis = beta.basis, p = p, V = V, W = W, s = s, 
                                      tol = tol, cutoff = cutoff, max.iter = max.iter)
            }
            
            # for validation method
            if(tuning.type == "validation")
            {
                beta.coef   <- c(fit.tuning$bVecTilde, fit.tuning$gammaVecTilde)
                y.valid.hat <- Psi.valid %*% beta.coef
                resid.valid <- y.valid - y.valid.hat
                
                #' OLS residual square 
                if(method == "OLS") err = sum(resid.valid^2)
                
                #' quantile loss
                if(method == "QR")
                {
                    rss.temp = vector(length = n.valid)
                    for(i in 1:n.valid) rss.temp[i] = rho.func(t = resid.valid[i], tau = tau)
                    err = sum(rss.temp)
                }
            }
            
            err.score[r,u] = err
        }
    }
    #' minimum error score index
    h = which(err.score == min(err.score), arr.ind = T) 
    
    #' return list
    tuning <- list(score = err.score,
                   eta.n = eta.n.vec[h[nrow(h),2]],
                   lambda.1 = lambda.1.vec[h[nrow(h), 1]], 
                   lambda.2 = lambda.2.vec[h[nrow(h), 1]],
                   min.score = min(err.score))
    
    return(tuning)
}
