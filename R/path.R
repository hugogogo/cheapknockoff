#' Multiple knockoff path
#'
#' This function generates a path of selected variables using multiple knockoff
#' given the test statistics (kappa, tau)
#'
#' @param kappa A \code{p} vector of test statistics, with kappa_i = 1 indicating the original variable winning
#' @param tau A \code{p} vector of test statistics, showing the manitude/importance of the variable
#'
#' @return An list of selected variable sets
#'
#' @examples
#' library(cheapknockoff)
#' set.seed(123)
#' n <- 100
#' p <- 30
#' x <- matrix(data = rnorm(n * p), nrow = n, ncol = p)
#' y <- x[, 1] - 2 * x[, 2] + rnorm(n)
#' omega <- c(2, 9, sample(seq(2, 9), size = 28, replace = TRUE))
#' # construct multiple knockoffs
#' X_k <- multiple_knockoff_Gaussian(X = x, mu = rep(0, p), Sigma = diag(1, p), omega = omega)
#' # compute knockoff statistics
#' stat <- cheapknockoff::stat_glmnet_coef(X = x, X_k = X_k, y = y, omega = omega)
#' # yield the path of selected variables
#' path <- cheapknockoff::generate_path(kappa = stat$kappa, tau = stat$tau)
#' @export

generate_path <- function(kappa, tau){
  p <- length(kappa)

  # input check
  stopifnot(length(kappa) == length(tau))

  # now tau[ord] is in non-increasing order
  ord <- order(tau, decreasing = TRUE)

  # al and kap are the permutations of omega and kappa, respectively, in the order ord
  kp <- kappa[ord]

  path <- list()
  for (k in seq(p)){
    path[[k]] <- ord[which(kp[1:k] == 1)]
  }
  return(path)
}
