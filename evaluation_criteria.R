#' This function return estimated beta function given the coefficients and beta.basis
#' @param beta.basis, the basis for beta
#' @param beta.coef, the coefficients for beta
#' @return z, the beta_t functions
beta.hat.func <- function(t, beta.basis, beta.coef)
{
    beta.Hat <- eval.basis(t,beta.basis) %*% beta.coef
    z <- t(beta.Hat[,1])
    return(z)
}


PMSE <- function(dat, Jmat, est, method, tau = tau)
{
    Ls     <- ncol(Jmat)
    X.test <- crossprod(dat$test$x$coefs, Jmat)
    Z.test <- dat$test$Z
    y.test <- dat$test$y
    
    n.test <- length(y.test)
    U.test <- matrix(0, nr = n.test, nc = (p-1)*Ls)
    for(i in 1:n.test){
        U.test[i,] <- kronecker(Z.test[i,-1],X.test[i,])
    }
    y.hat.test <- X.test %*% est$bVecTilde[1:Ls] + U.test %*% est$bVecTilde[(Ls+1):(p*Ls)] + Z.test %*% est$gammaVecTilde
    if(method == "OLS"){
        PMSE <- 1/n.test * sum((y.test - y.hat.test)^2) 
    }else if(method == "QR"){
        PMSE <- 1/n.test * sum((rho.func((y.test - y.hat.test), tau = tau)))
    }
    return(PMSE)
}

#' This function calculate the integrated squared error on the given region
#' @param beta.func, the true beta function
#' @param region.list, the evaluated regions
#' @param beta.basis, the B-spline basis for beta
#' @param beta.coef, the estimate beta coefficient, b
#' 
#' @return ISE, 

ISE.func <- function(beta.func, region.list, beta.basis, beta.coef)
{
    LL <- 0
    ISE_all <- 0
    LL <- vector()
    for(i in 1:length(region.list))
    {
        rg1 <- region.list[[i]][1]
        rg2 <- region.list[[i]][2]
        LL[i] <- rg2 - rg1
        ISE_all <- ISE_all + integrate(function(t) {(beta.func(t) - beta.hat.func(t, beta.basis, beta.coef))^2}, rg1, rg2)$value       
    }
    ISE <- ISE_all/sum(LL)
    return(ISE)
}

#' This function is used to calculate the correct identified proportions for given region
#' @param region, the evaluated region
#' @param bVecTilde, the estimated b vector
#' @param p, the number of scalar covariates
#' @param beta.basis, the B-spline basis for beta
#' @param region.type, indicator the type of region, including null regions and nonnull regions
#' 
#' @return corr.ratio, the corrected identified ratio.

region.corr.ration.func <- function(region, bVecTilde, p, beta.basis, region.type = "null.region")
{
    Ls <- beta.basis$nbasis
    if(region.type %in% c("null.region", "nonnull.region") == FALSE){
        print("unknown region type")
    }
    region.seq <- list()
    corr.ratio <- vector()
    
    for(i in 1:p){
        region.seq[[i]] = vector()
        for(j in 1:length(region[[i]])){
            region.seq[[i]] = c(region.seq[[i]], seq(region[[i]][[j]][1], region[[i]][[j]][2], by = 0.001))
        }
        
        temp.vec = lapply(region.seq[[i]], beta.hat.func, beta.basis = beta.basis,
                           beta.coef = bVecTilde[(1:Ls)+(i-1)*Ls]) %>% unlist
        if(region.type == "null.region"){
            corr.ratio[i] = sum(abs(temp.vec) <= 0.01)/length(region.seq[[i]])
        }else if(region.type == "nonnull.region"){
            corr.ratio[i] = sum(abs(temp.vec) > 0.01)/length(region.seq[[i]])
        }
    }
    return(corr.ratio)
}





