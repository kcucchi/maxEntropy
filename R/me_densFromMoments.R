#' probability density function from distribution moments using maximum entropy
#' 
#' \code{me_densFromMoments} calculates the probability density function corresponding to moments using the principle of maximum entropy
#' 
#' Some more detailed description. Authors : Karina Cucchi, Changhong Wang 
#' 
#' @param mu a vector containing moments
#' @param x the vector where to evaluate the me probability density function
#' @param lambda (optional) a vector for starting values of Lagrangian parameters. If the value of lambda is not provided by the user, a uniform distribution is used.
#' @param eps a number corresponding to the tolerance for the stopping rule. 
#' The iterations stop when the relative error is less that eps.
#' 
#' @return 
#' a list containing 
#' \enumerate{
#' \item a vector of Lagrangian parameters,
#' \item a vector containing values of the pdf at locations indicated in input parameter \code{x}
#' \item a vector containing values for entropy calculated at each iteration
#' } 
#' 
#' @examples
#' x=seq(-6,6,by=0.1)
#' mu = c(0.0000,5.4300,2.6779,51.3181,50.0401)
#' output = me_densFromMoments(mu = mu,x = x,eps = 1e-30)
#' lambda=output[[1]]
#' print(paste0('lambda = [',paste(round(lambda,3),collapse = ' '),']'))
#' 
#' @export
me_densFromMoments <- function(mu,x,xEval=NULL,lambda=NULL,eps=1e-5,verbose=F){

  mu = c(1,mu); # lambda_0 is one
  lx = length(x);
  xmin = x[1]; xmax = x[lx]; dx = x[2]-x[1]
  
  # if lambda not defined by user, set to uniform distribution
  if(is.null(lambda)){
    lambda = numeric(length=length(mu))
    lambda[1] = log(xmax-xmin)
  }
  N = length(lambda)
  M = 2*N-1 
  
  # needed for calculation of vandermonde matrix
  if(!("matrixcalc" %in% rownames(installed.packages()))) {install.packages("matrixcalc")}
  require(matrixcalc)
  fin = matrixcalc::vandermonde.matrix(x,M) # from package matrixcalc
  
  # start iterations to calculate lambda
  entr = numeric(0)
  iter=0
  while(TRUE){
    iter=iter+1
    if(verbose){cat(paste0('\n----------------\niter=',iter))}
    
    p=exp(-(fin[,1:N]%*%lambda)); # calculate p(x)
    
    # calculate Gn moments
    G=numeric(length=M)
    for(n in 1:M){
      G[n]=dx*sum(fin[,n]*p) # elementwise multiplication
    }
    
    # maximum entropy
    entr[iter]=-lambda %*% G[1:N]; # Calculate the entropy value
    if(verbose){cat(paste0('\nEntropy = ',entr[iter]))}
    
    # Calculate gnk
    gnk = array(0,dim=c(N,N))
    for (i in 1:N){ # Matrix G is a Hankel matrix
      gnk[,i]=-G[i:(N+i-1)]
    }
    
    #
    v = mu-G[1:N]
    delta=solve(gnk,v)
    lambda = lambda + delta
    
    # Convergence criteria
    if(all(abs(delta/lambda)<eps)){break}
    if(iter>2){
      if(abs((entr[iter]-entr[iter-1])/entr[iter])<eps){
        if(verbose){cat('\n CONVERGENCE REACHED !!')}
        break
      }
    }
  }
  
  # evaluate density at xEval with final lambda coefficients
  if(is.null(xEval)){xEval=x}
  # check if xEval in range of x
  rangex = range(x)
  if(xEval>= rangex[1] && xEval<=rangex[2]){
    fin = matrixcalc::vandermonde.matrix(xEval,M) # from package matrixcalc
    pEval = exp(-(fin[,1:N]%*%lambda))
  }else{pEval=0} # if not set to zero to avoid infinity problems

  
  if(verbose){cat('\n--------------- END -------------------')}
  
  return(list(lambda=lambda,p=p,pEval=pEval,xEval=xEval,entr=entr))
  
}


# write.table(q,file='realizations.txt',col.names=F,row.names=F,quote=F)