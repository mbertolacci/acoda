#' @importFrom coda mcpar
#' @importFrom coda spectrum0.ar
#' @importFrom coda thin
#' @importFrom stats end
#' @importFrom stats frequency
#' @importFrom stats quantile
#' @importFrom stats start
#' @importFrom stats time
#' @importFrom stats ts
#' @importFrom stats var
#' @importFrom stats window

#' @export
'[.mcmca' <- function(x, i, ..., drop = missing(i)) {
  y <- NextMethod('[')

  if (missing(i)) {
    return(mcmca(y, start = start(x), thin = thin(x)))
  }
  return(y)
}

#' Markov Chain Monte Carlo Objects for Arrays
#' @param data an array of MCMC output
#' @param start the iteration number of the first observation
#' @param end the iteration number of the last observation
#' @param thin the thinning interval between consecutive observations
#' @export
mcmca <- function(data = NA, start = 1, end = numeric(0), thin = 1) {
  if (is.vector(data)) {
    n_given <- length(data)
  } else {
    n_given <- dim(data)[1]
  }

  if (missing(end)) {
    end <- start + (n_given - 1) * thin
  } else if (missing(start)) {
    start <- end - (n_given - 1) * thin
  }

  n_actual <- floor((end - start) / thin + 1.0)
  if (n_given < n_actual) {
    stop('Start, end and thin incompatible with data')
  }
  end <- start + (n_actual - 1) * thin
  if (n_actual < n_given) {
    data <- index_array(data, 1, 1 : n_actual)
  }

  attr(data, 'mcpar') <- c(start, end, thin)
  class(data) <- 'mcmca'
  return(data)
}

#' @export
print.mcmca <- function(x, ...) {
  x.orig <- x
  cat(sprintf(
    paste0(
      'MCMC array output:\n',
      'Start = %d\n',
      'End = %s\n',
      'Thinning interval = %d\n'
      ),
    start(x), end(x), thin(x)
  ))
  attr(x, 'mcpar') <- NULL
  attr(x, 'class') <- NULL
  NextMethod('print', ...)
  invisible(x.orig)
}

#' @export
start.mcmca <- function(x, ...) {
  mcpar(x)[1]
}

#' @export
end.mcmca <- function(x, ...) {
  mcpar(x)[2]
}

#' @export
frequency.mcmca <- function(x, ...) {
  1 / thin.mcmca(x)
}

#' @export
thin.mcmca <- function(x, ...) {
  mcpar(x)[3]
}

#' @export
time.mcmca <- function(x, ...) {
  ts(
    seq(
      from = start(x), to = end(x), by = thin(x)
    ),
    start = start(x),
    end = end(x),
    deltat = thin(x)
  )
}

#' Dimensions of MCMC Array objects
#' @param x An mcmca object
#' @export
nvariables <- function(x) {
  if (!is.array(x)) {
    return(1)
  }
  n_dims <- length(dim(x))
  return(prod(dim(x)[2 : n_dims]))
}

#' @describeIn nvariables returns the number of iterations
#' @export
niterations <- function(x) {
  if (!is.array(x)) {
    return(length(x))
  }
  return(dim(x)[1])
}

#' @export
window.mcmca <- function(x, start = numeric(0), end = numeric(0),
                         thin = numeric(0), ...) {
  ts_eps <- getOption('ts.eps')
  x_start <- start(x)
  x_end <- end(x)
  x_thin <- thin(x)

  if (missing(thin)) {
    thin <- x_thin
  } else if (thin %% x_thin != 0) {
    thin <- x_thin
    warning('thin value not changed')
  }

  x_time <- as.vector(time(x))
  if (missing(start)) {
    start <- x_start
  } else if (length(start) != 1) {
    stop('bad value for start')
  } else if (start < x_start) {
    start <- x_start
    warning('start value not changed')
  }
  if (missing(end)) {
    end <- x_end
  } else if (length(end) != 1) {
    stop('bad value for end')
  } else if (end > x_end) {
    end <- x_end
    warning('end value not changed')
  }
  if (start > end) {
    stop('start cannot be after end')
  }

  if (all(abs(x_time - start) > abs(start) * ts_eps)) {
    start <- x_time[(x_time > start) & ((start + x_thin) > x_time)]
  }
  if (all(abs(end - x_time) > abs(end) * ts_eps)) {
    end <- x_time[(x_time < end) & ((end - x_thin) < x_time)]
  }
  use <- 1 : niterations(x)
  use <- use[
    use >= trunc((start - x_start) / x_thin + 1.5) &
    use <= trunc((end - x_start) / x_thin + 1.5) &
    (use - trunc((start - x_start) / x_thin + 1.5)) %% (thin %/% x_thin) == 0
  ]
  y <- index_array(x, 1, use)

  return(mcmca(y, start = start, end = end, thin = thin))
}

safespec0 <- function (x) {
  result <- try(spectrum0.ar(x)$spec)
  if (class(result) == 'try-error') {
    result <- NA
  }
  result
}

#' @export
summary.mcmca <- function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                          ...) {
  x <- object
  if (!is.array(x)) {
    # Convert vectors to matrix to allow use of dim()
    x <- mcmca(matrix(x, ncol = 1))
  }
  n_dims <- length(dim(x))
  x_stats <- array(0, dim = c(dim(x)[2 : n_dims], 4))
  dimnames(x_stats) <- vector('list', n_dims)
  if (!is.null(dimnames(x))) {
    dimnames(x_stats)[1 : (n_dims - 1)] <- dimnames(x)[2 : n_dims]
  }
  dimnames(x_stats)[[n_dims]] <- c('Mean', 'SD', 'Naive SE', 'Time-series SE')
  index_array(x_stats, n_dims, 1) <- apply(x, 2 : n_dims, mean)
  index_array(x_stats, n_dims, 2) <- sqrt(apply(x, 2 : n_dims, var))
  index_array(x_stats, n_dims, 3) <- (
    sqrt(apply(x, 2 : n_dims, var) / niterations(x))
  )
  index_array(x_stats, n_dims, 4) <- (
    sqrt(apply(x, 2 : n_dims, safespec0) / niterations(x))
  )

  x_quantiles <- apply(x, 2 : n_dims, quantile, quantiles)
  # Put the quantiles on the columns
  if (n_dims == 2) {
    x_quantiles <- t(x_quantiles)
  } else {
    x_quantiles <- aperm(x_quantiles, c(2, 1, 3 : n_dims))
  }

  output <- list(
    statistics = x_stats,
    quantiles = x_quantiles,
    start = start(x),
    end = end(x),
    thin = thin(x)
  )
  class(output) <- 'summary.mcmca'
  return(output)
}

#' @export
print.summary.mcmca <- function (x, digits = max(3, .Options$digits - 3), ...) {
  cat('\n', 'Iterations = ', x$start, ':', x$end, '\n', sep = '')
  cat('Thinning interval =', x$thin, '\n')
  cat('Sample size =', (x$end - x$start) / x$thin + 1, '\n')
  cat('\n1. Empirical mean and standard deviation for each variable,')
  cat('\n   plus standard error of the mean:\n\n')
  print(x$statistics, digits = digits, ...)
  cat('\n2. Quantiles for each variable:\n\n')
  print(x$quantiles, digits = digits, ...)
  cat('\n')
  invisible(x)
}
