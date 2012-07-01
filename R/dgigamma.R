#' generalized inverse gamma density
#'
#' @param x vector of quantiles
#' @param alpha alpha parameter
#' @param beta beta parameter
#' @param gamma gamma parameter
#' @param mu mu parameter
#' @export
dgigamma <- function(x,alpha,beta,gamma=1,mu=0){
  r = rep(0,length(x))
  ok = x>mu
  bxmu = beta/(x[ok]-mu)
  d = exp(-(bxmu)^gamma)* gamma * bxmu^(1+alpha*gamma)
  r[ok]=d/(beta*gamma(alpha))
  r
}

.base_pgigamma <- function(q,alpha,beta,gamma=1,mu=0){
  force(q);force(alpha);force(beta);force(gamma);force(mu)
  f = function(x){dgigamma(x,alpha,beta,gamma,mu)}
  integrate(f,lower=0,upper=q)$value
}

#' distribution function for generalized inverse gamma
#'
#' @param q vector of quantiles
#' @param alpha alpha parameter
#' @param beta beta parameter
#' @param gamma gamma parameter
#' @param mu mu parameter
#' @export
pgigamma = Vectorize(.base_pgigamma)

.base_qgigamma <- function(p,alpha,beta,gamma=1,mu=0){
  force(p);force(alpha);force(beta);force(gamma);force(mu)
  pf = function(f){
    force(f)
    pgigamma(f,alpha=alpha,beta=beta,gamma=gamma,mu=mu)-p
  }

  upper = 2
  while(pf(upper)<0){
    upper=upper*2
  }
  
  uniroot(pf,c(0,upper))$root
}


#' quantile function for generalized inverse gamma
#'
#' @param p vector of probabilities
#' @param alpha alpha parameter
#' @param beta beta parameter
#' @param gamma gamma parameter
#' @param mu mu parameter
#' @export

qgigamma = Vectorize(.base_qgigamma)
