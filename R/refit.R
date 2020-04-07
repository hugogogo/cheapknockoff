#' Refitting a path of selected variables from the mknockoff filter and compute prediction on the newdata
#'
#' @param path A path of selected variables set
#' @param x A \code{n}-by-\code{p} matrix of original variables
#' @param y A \code{n} vector of response
#' @param newdata A new data set of features that the predictions are computed
#' @param family The conditional distribution of y given X. See the family option for \code{glmnet}.
#'
#' @import glmnet
#'
#' @return a list consisting the refitted model and the corresponding prediction for each model on the path returned by \code{mk_path}.
#' @examples
#' library(mknockoff)
#' set.seed(123)
#' n <- 100
#' p <- 30
#' x <- matrix(data = rnorm(n * p), nrow = n, ncol = p)
#' y <- x[, 1] - 2 * x[, 2] + rnorm(n)
#' omega <- c(2, 9, sample(seq(2, 9), size = 28, replace = TRUE))
#' # construct multiple knockoffs
#' X_k <- knockoff_Gaussian(X = x, mu = rep(0, p), Sigma = diag(1, p), omega = omega)
#' # compute knockoff statistics
#' stat <- mknockoff::stat_glmnet_coef(X = x, X_k = X_k, y = y, omega = omega)
#' # yield the path of selected variables
#' path <- mknockoff::mk_path(kappa = stat$kappa, tau = stat$tau)
#' # refit and compute predictions on the original data
#' refit_output <- refit(path, x, y, x)
#'
#' @importFrom stats coef glm predict.glm rnorm
#' @export
refit <- function(path, x, y, newdata, family = "gaussian") {
  n_mod <- length(path)
  n <- nrow(newdata)
  pred <- matrix(NA, n, n_mod)

  result <- list()

  lam <- ifelse(n <= n_mod, 1e-3, 1e-6)

  for(i in seq(n_mod)){
    idx <- path[[i]]

    # added a tiny ridge
    if (length(idx) < 2){
      dat_tr <- data.frame(y = y, z = x[, idx])
      mod <- glm(formula = y ~ ., data = dat_tr, family = family)
      dat_pr <- data.frame(z = newdata[, idx])
      pred <- as.numeric(predict.glm(mod, newdata = dat_pr, type = "response"))
    }
    else{
      # glmnet does not work when there is only one column in the design matrix
      mod <- glmnet(x = x[, idx], y = y, family = family, alpha = 0, lambda = lam)
      pred <- as.numeric(predict.glmnet(object = mod, newx = newdata[, idx], type = "response"))
    }

    # store results
    result$mod[[i]] <- mod
    result$pred[[i]] <- pred
  }
  result$path <- path
  return(result)
}
