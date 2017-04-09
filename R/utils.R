# As per http://stackoverflow.com/questions/14500707/select-along-one-of-n-dimensions-in-array
index_array <- function(x, dim, i, drop = FALSE) {
  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- i

  # Generate the call to [
  call <- as.call(c(
    list(as.name('['), quote(x)),
    indices,
    list(drop = drop)
  ))
  # Finally, evaluate it
  eval(call)
}

'index_array<-' <- function(x, dim, i, value) {
  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- i

  # Generate the call to [
  call <- as.call(c(
    list(as.name('[<-'), quote(x)),
    indices,
    list(quote(value))
  ))
  # Finally, evaluate it
  eval(call)
}
