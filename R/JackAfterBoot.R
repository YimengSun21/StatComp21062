#' @title Jackknife-after-Bootstrap
#' @description Calculate the bias and standard deviation estimate for each leave-one-out Bootstrap sample.
#' @param data the data as a vector
#' @param func the function which applied to data
#' @param B the number of replicates in Bootstrap
#' @return the bias and standard deviation estimate for both Bootstrap and Jackknife-after-Bootstrap methods
#' @examples
#' \dontrun{
#' data("cars")
#' data <- cars$speed
#' JackAfterBoot(data, mean, 1000)
#' JackAfterBoot(data, median, 1000)
#' }
#' @importFrom stats sd
#' @useDynLib StatComp21062
#' @export
JackAfterBoot <- function(data, func = NULL, B){
  n <- length(data)
  theta <- func(data)
  theta_b <- numeric(B)
  mat <- matrix(0, nrow = B, ncol = n)
  ### Bootstrap
  for(i in 1:B){
    sam <- sample(1:n, size = n, replace = TRUE)
    mat[i, ] <- sam
    dat <- data[sam]
    theta_b[i] <- func(dat)
  }
  ### Jackknife After Bootstrap
  se_jack <- numeric(n)
  theta_j <- numeric(n)
  for (i in 1:n) {
    keep <- (1:B)[apply(mat, MARGIN = 1, FUN = function(k) {!any(k == i)})]
    theta_j[i] <- mean(theta_b[keep])
    se_jack[i] <- sd(theta_b[keep])
  }
  se_boot <- sd(theta_b)
  bias_boot <- mean(theta_b) - theta
  se_jackafterboot <- sqrt((n-1) * mean((se_jack - mean(se_jack))^2))
  bias_jackafterboot <- (n-1) * (mean(theta_j) - theta)
  result <- list(se_Bootstarp = se_boot,
                 se_JackAfterBoot = se_jackafterboot,
                 bias_Bootstarp = bias_boot,
                 bias_JackAfterBoot = bias_jackafterboot)
  return(result)
}
