#' This function is used for calculate the quantile + smoothing spline method, without interaction penalty
#' @param thetaVecHat, the initial value for theta
#' @param tau, tau-th quantile
#' @param V.star, a diagonal matrix
#' @param Psi, the design matrix,c(X, U, Z)
#' @param y, the response, n by 1 
#' @param eta.n, tuning parameter for smooothing spline
#' @param tol, the threshold for convergence
#' @param max.iter, the max iteration times
#' @param epsilon, small value to control the MM precision
#' 
#' @return thetaVecHat, the estimate for thetaVecHat, 
#' delta, the vec norm difference between last two estimate
#' it, number of iterations

QR.smoothing.spline <- function(thetaVecHat, tau, V.star, Psi, y, eta.n, tol = tol, max.iter = 300, epsilon = epsilon){
    n <- nrow(Psi)
    delta <- 1
    it <- 0
    pPsi <- nrow(thetaVecHat)
    
    while(it < max.iter & delta > tol){
        it <- it + 1
        thetaVecHat0 <- thetaVecHat
        
        R <- y-Psi %*% thetaVecHat0
        A <- 1/(epsilon + abs(R))
        Vrho <- 1-2*tau-R*A
        #updata the thetaVec
        DeltaNum <- crossprod(Psi, Vrho) + 4*n*eta.n * V.star %*% thetaVecHat0
        DeltaDen <- crossprod(Psi, diag(A[,1])) %*% Psi + 4*n*eta.n*V.star
        
        thetaVecHat <- thetaVecHat0 - MASS::ginv(DeltaDen) %*% DeltaNum
        delta <- vec.norm.func(thetaVecHat-thetaVecHat0)/pPsi
    }
    est.list <- list("thetaVecHat" = thetaVecHat, "delta" = delta, "it" = it)
}

#' This function is used for calculate the quantile + smoothing spline + interaction penalty
#' @param thetaVecHat, the initial value for theta
#' @param tau, the tau-th quantile
#' @param V.star, quasi diagonal matrix
#' @param Psi, the design matrix,c(X, U, Z)
#' @param y, the response, n by 1 
#' @param beta.basis, the B-spline basis for beta
#' @param eta.n, tuning parameter for smooothing spline
#' @param lambda.1, the paramter for first penalty
#' @param lambda.2, the paramter for second penalty 
#' @param n, the length of y
#' @param p, the number of scalar covariates
#' @param W, W matrix
#' @param s, the regularization parameter for MCP penalty
#' @param tol, the threshold for convergence 
#' @param max.iter, the max iteration times
#' @param epsilon, small value to control the MM precision
#'
#' @return a list: thetaVecHat, the estimate for thetaVecHat. 
#' it, the number of iteration when convergence 

QR.interaction.lqa <- function(thetaVecHat, tau, V.star, Psi, y, beta.basis, eta.n, 
                                    lambda.1, lambda.2, p, W, s, tol, max.iter = 300, epsilon = epsilon)
{
    #' calculate parameters
    Ls     <- beta.basis$nbasis
    Mn     <- length(beta.basis$params) + 1
    ds     <- Ls-Mn
    sMoT   <- sqrt(Mn/diff(getbasisrange(beta.basis)))
    pPsi   <- nrow(thetaVecHat)
    n      <- nrow(Psi)
    
    #' initialization
    thetaNormj <- matrix(0, Mn, p)
    bZeroMat   <- matrix(FALSE, 1, pPsi)
    thetaNorm  <- rep(Inf, p)
    
    delta <- 1
    it <- 0
    
    #' loop
    while(it < max.iter & delta > tol)
    {
        it <- it + 1
        thetaVecHat0 <- thetaVecHat
        
        #' caluclate the norm for whole estimator
        thetaNormOld <- thetaNorm
        for(j in 1:p) thetaNorm[j] <- vec.norm.func(thetaVecHat0[(1:Ls)+(j-1)*Ls,])
        delta <- max((thetaNormOld - thetaNorm)^2)
        
        #' calculate the thetaNormj
        for(j in 1:p)
        {
            for(k in 1:Mn)
            {
                index <- c(k:(k+ds))
                b <- thetaVecHat0[index+(j-1)*Ls]
                thetaNormj[k,j] <- sqrt(t(b) %*% W[,,k] %*% b)
            }
        }
        
        #' calculate W.2.tilde
        thetaNormjSum <- sqrt(rowSums(thetaNormj^2))
        W.Hat <- matrix(0, Ls, Ls)
        
        if(lambda.2 > 0)
        {
            for(k in 1:Mn)
            {
                index <- c(k:(k+ds))
                c2k <- DMCP.func((thetaNormjSum[k] * sMoT), s, lambda = lambda.2)
                
                if(c2k != 0){
                    if(thetaNormjSum[k] < tol) 
                    {
                        index.temp <- prod.func(index,c(1:(p-1)),Ls)
                        bZeroMat[index.temp] <- TRUE
                    }else{ 
                        W.Hat[index, index] <- W.Hat[index, index] + c2k/2 * (sMoT/thetaNormjSum[k]) * W[,,k]
                    }
                }
            }
            W.tilde.2 <- V.star.func(p, W.Hat)
        }
        
        bZeroVec <- bZeroMat
        bNonZeroVec <- !bZeroVec
        
        #' calculate W.1.tilde
        W.tilde.1 <- matrix(0, nr = pPsi, nc = pPsi)
        
        for(j in 1:p){
            if(lambda.2 > 0 && j == 1){
                next
            }else{
                WHatj1 <- matrix(0, Ls, Ls)
                
                for(k in 1:Mn){
                    index <- c(k:(k+ds))
                    c1kj <- DMCP.func((thetaNormj[k,j]*sMoT), s, lambda = lambda.1)
                    
                    if(c1kj != 0)
                    {
                        if(thetaNormj[k,j] < tol){
                            bZeroVec[Ls*(j-1) + index] <- TRUE
                        }else{
                            WHatj1[index,index] <- WHatj1[index, index] + c1kj/2*(sMoT/thetaNormj[k,j])*W[,,k]
                        }
                    }
                }
                W.tilde.1[(1:Ls)+(j-1)*Ls, (1:Ls)+(j-1)*Ls] <- WHatj1
            }
        }
        bNonZeroVec <- !bZeroVec
        Psi.temp <- Psi[, bNonZeroVec]
        thetaVecHat0.temp <- thetaVecHat0[bNonZeroVec,]
        
        if(lambda.2 > 0){
            W.tilde.temp <- W.tilde.1[bNonZeroVec, bNonZeroVec] + W.tilde.2[bNonZeroVec, bNonZeroVec]
        }else{
            W.tilde.temp <- W.tilde.1[bNonZeroVec, bNonZeroVec]
        }
        
        V.star.temp <- V.star[bNonZeroVec, bNonZeroVec]
        
        #' updata the estimator 
        R <- y - Psi.temp %*% thetaVecHat0.temp
        A <- 1/(epsilon + abs(R))
        Vrho <- 1 - 2 * tau - R * A
        
        #' updata the thetaVec
        DeltaNum <- crossprod(Psi.temp, Vrho) + 4 * n * W.tilde.temp %*% thetaVecHat0.temp +
            4 * n * eta.n * V.star.temp %*% thetaVecHat0.temp
        
        DeltaDen <- crossprod(Psi.temp, diag(A[,1])) %*% Psi.temp + 4 * n * W.tilde.temp +
            4 * n * eta.n *V.star.temp
        thetaVecHat <- matrix(0, pPsi, 1)
        thetaVecHat[bNonZeroVec,1] <- thetaVecHat0.temp - MASS::ginv(DeltaDen) %*% DeltaNum
    }
    est.list <- list("thetaVecHat" = thetaVecHat, "delta" = delta, "it" = it)
}


#' This function is used for calculate the the Alt.4, Alt.5, and proposed methods in our main paper
#' @param Psi, the augmented design matrix, n by (Ls*p+p)
#' @param y, the response, n by 1
#' @param tau, the tau-th quantile
#' @param eta.n, tuning parameter for smooothing spline
#' @param lambda.1, the parameter for first penalty
#' @param lambda.2, the parameter for second penalty
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
#' @return est, a list including 
#' bVecTilde: the estimate for b, 
#' gammaVecTilde: the estimated gamma, 
#' beta.est.list, the estimated beta function,
#' fitted.value, the fitted y

QR.alg <- function(Psi, y, tau, eta.n, lambda.1, lambda.2, beta.basis, p, V, W, s, tol, cutoff, max.iter, epsilon)
{
    #' calculate parameters
    Ls     <- beta.basis$nbasis
    n      <- nrow(Psi)
    V.star <- V.star.func(p = p, V = V)
    
    #' caluculate the initialization theta
    thetaVecHat <- MASS::ginv(crossprod(Psi, Psi) + n*eta.n*V.star) %*% crossprod(Psi,y)
    thetaTilde  <- thetaVecHat
    
    #' caluculate the OLS smoothing estimator
    if(lambda.1 == 0 & lambda.2 == 0){
        est.QR.SS<- QR.smoothing.spline(thetaVecHat, tau, V.star, Psi, y, eta.n, 
                                                  tol = tol, max.iter = 300, epsilon = epsilon)
        thetaTilde <- est.QR.SS$thetaVecHat
    }
    
    #if sparse,perform NISF 
    if(lambda.1 != 0 || lambda.2 != 0){
        est.QR.ilqa <- QR.interaction.lqa(thetaVecHat, tau, V.star, Psi, y, beta.basis, eta.n, 
                                                    lambda.1, lambda.2, p, W, s, tol, max.iter = 300, epsilon = epsilon)
        thetaTilde <- est.QR.ilqa$thetaVecHat
    }
    
    bZero <- (abs(thetaTilde) < cutoff)
    thetaTilde[bZero] <- 0
    bNonZero <- !bZero
    
    #' the estimated b and gamma
    bVecTilde <- thetaTilde[1:(p*Ls)]
    gammaVecTilde <- thetaTilde[(p*Ls+1):(p*Ls+p)]
    
    #' the estimated beta function
    beta.list <- list()
    for(j in 1:p) beta.list[[j]] <- fd(coef = bVecTilde[(1+(j-1)*Ls):(Ls+(j-1)*Ls)], basisobj = beta.basis)
    
    #' the fitted response
    fitted.values <- Psi %*% thetaTilde
    
    #' return list
    est <- list( "bVecTilde" = bVecTilde, "gammaVecTilde"= gammaVecTilde, 
                 "beta.est.list" = beta.list, "fitted.values" = fitted.values)
    return(est)
}

