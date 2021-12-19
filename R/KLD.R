#' @title Kullback-Leibler divergence(KLD)
#' @description Calculate the Kullback-Leibler divergence between two probability distributions.
#' @param px discrete probability distributions
#' @param py discrete probability distributions
#' @param base the logarithmic base, defaults to \code{e}=exp(1)
#' @return forward Kullback-Leibler divergence and 
#' backward Kullback-Leibler divergence of discrete probability distributions p and q
#' @examples
#' \dontrun{
#' p <- c(105,24,10,2,120,56)
#' q <- c(1,4,8,15,200,78)
#' KLD(p, q)
#' KLD(p, q, base = exp(1))
#' KLD(p, q, base = 2)
#' }
#' @useDynLib StatComp21062
#' @export
KLD <- function(px, py, base = exp(1))
{
  ### Judge
  if(!is.vector(px)) px <- as.vector(px)
  if(!is.vector(py)) py <- as.vector(py)
  n1 <- length(px)
  n2 <- length(py)
  if(!identical(n1, n2))
    stop("px and py must have the same length.")
  if(any(!is.finite(px)) || any(!is.finite(py)))
    stop("px and py must be finite values.")
  if(any(px < 0) || any(py < 0))
    stop("px and py must be nonnegative values.")
  ### Calculate the Kullback-Leibler divergence(KLD)
  px <- px / sum(px)
  py <- py / sum(py)
  kld.forward <- px * (log(px, base = base)-log(py, base = base))
  kld.backward <- py * (log(py, base = base)-log(px, base = base))
  KLD.forward <- sum(kld.forward)
  KLD.backward <- sum(kld.backward)
  result <- list(KLD.forward = KLD.forward,
                 KLD.backward = KLD.backward)
  return(result)
}
