#' probability density function evaluation from distribution samples using maximum entropy
#' 
#' \code{me_densFromSamples} calculates the probability density function corresponding to distribution samples using the principle of maximum entropy
#' 
#' Some more detailed description. Authors : Karina Cucchi, Changhong Wang 
#' 
#' @param q the vector containing samples from the distribution
#' @param nbMoments an integer corresponding to the number of moments to account for in the calculation of the me pdf
#' @param x is a vector containing values to evaluate the pdf and necessary to calculate the entropy.
#' If x is not defined, the default is seq(min(q)-sd(q),max(q)+sd(q),length=1000)
#' @param xEval vector containing where to evaluate the probability density function.
#' If xEval is not specified, it is set to x.
#' @param eps a number corresponding to the tolerance for the stopping rule.
#' @param verbose a boolean indicating whether to allow for printing outputs
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
#' q=c(3.132,3.132,3.132,-1.625,-1.821,1.026,-3.003,-2.332,
#' 0.827,-0.437,0.305,-3.433,1.865,-1.131,-2.828,3.921,-0.730)
#' output = me_densFromSamples(q=q,nbMoments=5,x = x,eps = 1e-30)
#' lambda=output[[1]]
#' print(paste0('lambda = [',paste(round(lambda,3),collapse = ' '),']'))
#' histData = hist(x = q,freq = F,xlim=range(x),col = 'lightblue')
#' lines(x,p,col='red')
#' 
#' @export
me_densFromSamples <- function(q,nbMoments,x=NULL,xEval=NULL,lambda=NULL,eps=1e-5,verbose=F){
  
  # center and normalize data
  mean_q = mean(q)
  sd_q = sd(q)
  q_center = (q - mean_q) / sd_q
  
  # calculate first nbMoments moments of samples contained in q
  q_org_mom = numeric(length = nbMoments)
  q_n = length(q_center)
  for(i in 1:nbMoments){
    q_org_mom[i]=sum(q_center^i)/q_n
  }
  if(verbose){cat(paste0('moments : ',paste(round(q_org_mom,5),collapse = ' ')))}
  
  # call function calculating density from moments
  # define x vector if not defined in arguments
  if(is.null(x)){x_center=seq(from=min(q_center)-sd(q_center),
                       to=max(q_center)+sd(q_center),length=1000)
  }else{x_center = (x - mean_q)/sd_q }
  # if evaluation points are not specified, defined to be x
  if(is.null(xEval)){
    xEval_center=x_center
  }else{xEval_center = (xEval - mean_q)/sd_q}
  output_me = me_densFromMoments(mu = q_org_mom,lambda=lambda,eps=eps,
                                 x = x_center,xEval=xEval_center,
                                 verbose=verbose)

  return(output_me)
  
}

