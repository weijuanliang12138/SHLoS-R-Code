################ Beta function for simulation 1 ################
#' Define the beta_0 functions
beta.func.10 <- function(t)
{
  t1 <- t[t<=0.3]
  t3 <- t[t>=0.7]
  y <- matrix(0,1,length(t))
  y[t<=0.3] <- 2*(1-t1)*sin(2*pi*(t1+0.2))
  y[t>=0.7] <- 2*t3*sin(2*pi*(t3-0.2))
  y
} 

#' Define the beta_1 function
beta.func.11 <- function(t)
{
  t1 <- t[t<=0.3]
  t3 <- t[t>=0.7]
  y <- matrix(0,1,length(t))
  y[t<=0.3] <- 0
  y[t>=0.7] <- 2*t3*sin(2*pi*(t3-0.2))
  y
} 

#' Define the beta_2 function
beta.func.12 <- function(t)
{
  t1 <- t[t<=0.3]
  t3 <- t[t>=0.7]
  y <- matrix(0,1,length(t))
  y[t<=0.3] <- 2*(1-t1)*sin(2*pi*(t1+0.2))
  y[t>=0.7] <- 0
  y
}


################ Beta function for simulation 2 ################
#' Define the beta_0 functions
a = 5
beta.func.20 <- function(t)
{
    t1 = t[t<=0.3 & t > 0.2]
    t2 = t[t<=0.6 & t > 0.5]
    t3 = t[t<=0.8 & t > 0.7]
    y  = matrix(0,1,length(t))
    y[t<=0.3 & t > 0.2] = a*sin(pi*(t1-0.2)/0.1)
    y[t<=0.6 & t > 0.5] = -a*0.6*sin(pi*(t2-0.5)/0.1)
    y[t<=0.8 & t > 0.7] = a*0.7*sin(pi*(t3-0.7)/0.1)
    y
} 

#' Define the beta_1 functions
beta.func.21 <- function(t)
{
    t1 = t[t<=0.3 & t > 0.2]
    t2 = t[t<=0.6 & t > 0.5]
    y  = matrix(0,1,length(t))
    y[t<=0.3 & t > 0.2] = 0.4*a*((t1-0.25)^2 - 0.05^2)/0.05^2
    y[t<=0.6 & t > 0.5] = a*sin(pi*(t2-0.5)/0.1)
    y
} 

#' Define the beta_2 functions
beta.func.22 <- function(t)
{
    t2 = t[t<=0.6 & t > 0.5]
    t3 = t[t<=0.8 & t > 0.7]
    y <- matrix(0,1,length(t))
    y[t<=0.6 & t > 0.5] = 0.5*a*sin(pi*(t2-0.5)/0.1)
    y[t<=0.8 & t > 0.7] = 0.8*a*((t3-0.75)^2 - 0.05^2)/0.05^2
    y
}


################ Beta function for simulation 3 ################
#' Define the beta_0 functions
beta.func.30 <- function(t)
{
    t1 = t[t<=0.15 & t>0.125] 
    t2 = t[t<=0.2 & t>0.175]
    t3 = t[t<=0.35  & t>0.325]
    t4 = t[t<=0.625 & t>0.6]
    t5 = t[t<=0.725 & t>0.7]
    t6 = t[t<=0.825 & t>0.8]
    t7 = t[t<=0.9 & t>0.875]
    
    y  = matrix(0,1,length(t))
    y[t<=0.15 & t>0.125] = 0.4*2*a*((t1-0.1375)^2 - 0.0125^2)/0.0125^2
    y[t<=0.2 & t>0.175] = 0.7*2*a*sin(pi*(t2-0.175)/0.025)
    y[t<=0.35  & t>0.325] = -0.6*2*a*sin(pi*(t3-0.325)/0.025)
    y[t<=0.625 & t>0.6] = 0.8*2*a*sin(pi*(t4-0.6)/0.025)
    y[t<=0.725 & t>0.7] = -a*2*sin(pi*(t5-0.7)/0.025)
    y[t<=0.825 & t>0.8] = 0.5*2*a*sin(pi*(t6-0.8)/0.025)
    y[t<=0.9 & t>0.875] = -0.7*2*a*sin(pi*(t7-0.875)/0.025)
    
    y
} 
#plot(beta.func.30, ylim = c(-a, a))

#' Define the beta_1 functions
beta.func.31 <- function(t)
{
    t1 = t[t<=0.15 & t>0.125] 
    t3 = t[t<=0.35  & t>0.325]
    t5 = t[t<=0.725 & t>0.7]
    t7 = t[t<=0.9 & t>0.875]
    
    y  = matrix(0,1,length(t))
    y[t<=0.15 & t>0.125] = 2*a*sin(pi*(t1-0.125)/0.025)
    y[t<=0.35  & t>0.325] = 0.6*2*a*sin(pi*(t3-0.325)/0.025)
    y[t<=0.725 & t>0.7] = 0.8*2*a*((t5-0.7125)^2 - 0.0125^2)/0.0125^2
    y[t<=0.9 & t>0.875] = 0.9*2*a*sin(pi*(t7-0.875)/0.025)
    
    y
} 
#plot(beta.func.31, ylim = c(-a, a))

#' Define the beta_2 functions
beta.func.32 <- function(t)
{
 
    t2 = t[t<=0.2 & t>0.175]
    t4 = t[t<=0.625 & t>0.6]
    t6 = t[t<=0.825 & t>0.8]


    y  = matrix(0,1,length(t))
    y[t<=0.2 & t>0.175] = 0.5*2*a*sin(pi*(t2-0.175)/0.025)
    y[t<=0.625 & t>0.6] = 2*a*((t4-0.6125)^2 - 0.0125^2)/0.0125^2
    y[t<=0.825 & t>0.8] = 0.7*2*a*sin(pi*(t6-0.8)/0.025)
    y
}



# \int_0^1 B_j(t) * beta(t)dt
#' This function is used to calculate the inner product of f(t) and bfun(t)
#' @param f, function
#' @param basis, the basis for bfun(t)
#' @param j, evaluate at j
#' 
#' @return y$value, the inner product of f(t) and bfun(t)
inner.prod <- function(f,basis,j)
{
  rng <- getbasisrange(basis)
  knots <- c(rng[1],basis$params,rng[2])
  nbasis <- basis$nbasis
  norder <- basis$nbasis - length(knots) + 2
  
  a <- rng[1]
  if(j-norder > 0) a <- knots[j-norder+1]
  
  b <- rng[2]
  if (j <= nbasis-norder) b <- knots[j+1]
  
  bfun <- function(t)
  {
    mat <- eval.basis(t,basis)
    z <- t(mat[,j])
  }
  
  y <- integrate(function(t) {f(t)*bfun(t)},a,b)
  y$value
}


preprocess.x <- function(X,D)
{
  if(is.matrix(X) || is.fd(X)){
    X <- X %>% list()
    K <- 1
  }else K <- length(X)
  
  xfd.list <- list()
  for(k in 1:K){
    x <- X[[k]]
    if(!is.fd(x)){
      if(is.matrix(D) || is.vector(D)){
        x <- Data2fd(argvals = D, y = x)
      }else{
        x <- Data2fd(argvals = D[[k]], y = x)
      }
    }
    if(!is.fd(x))
      stop("X must be a matrix, fd object or a list of matrices/fd objects")
    xfd.list[[k]] <- x
  }
  xfd.list
}

#' This function is use to calculate the weight matrix W
#' @param basis, the basis to evaluated
#' 
#' @return W, the weight matrix
compute.weights <- function(basis)
{
  L <- basis$nbasis
  rng <- getbasisrange(basis)
  breaks <- c(rng[1],basis$params,rng[2])
  M <- length(breaks) - 1
  norder <- L-M+1
  W <- array(0,dim=c(norder,norder,M))
  for (j in 1:M)
  {
    temp <- inprod(basis,basis,rng=c(breaks[j],breaks[j+1]))
    W[,,j] <- temp[j:(j+norder-1),j:(j+norder-1)]
  }
  W
}


#' This function is use to calculate the vector L2 norm
#' @param v, vector to evaluate.
#' 
#' @return the L2 norm of v
vec.norm.func <- function(v)
{
  sqrt(sum(v^2))
}

#' The derivative of MCP function
#' @param t, the value to evaluate
#' @param s, the regularization parameter
#' @param lambda, the tunning parameter
#' 
#' @return the derivative of MCP function evaluate at t, with regularization parameter s, and tuning parameter lambda
DMCP.func <- function(t, s, lambda){
  ifelse(t >= (lambda*s), 0, (lambda-t/s))
}


#' This function is used to calculate the main effect and its corresponding interaction index.
#' @param a, a index vector
#' @param b, the number of interactions plus the main effect
#' @param Ls, the number of B-spline basis
#' 
#' @return the index vector
prod.func <- function(a, b, Ls){
  L <- length(b)
  ret <- a
  for(l in 1:L){
    ret = cbind(ret,a + b[l]*Ls) 
  }
  ret %>% as.vector
}

#' This function is used to calculate the V.star function
#' @param p, the number of scalar covariates plus 1
#' @param V, orignal V matrix
#' 
#' @return V.star matrix
V.star.func <- function(p ,V){
    Ls <- nrow(V)
    ncoef <- p*Ls+p
    V.star = matrix(0, nr = ncoef, nc = ncoef)
    
    for(i in 1:p) V.star[((1:Ls)+Ls*(i-1)),((1:Ls)+Ls*(i-1))] = V
    return(V.star)
}



