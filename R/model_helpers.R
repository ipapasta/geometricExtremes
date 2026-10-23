# Internal model vocabulary. Numeric labels follow the revised manuscript.
directional_architecture <- function(model) {
  labels <- c(M1="decoupled", M2="shared-geometry", M3="residual-adjustment")
  if (length(model) != 1L || is.na(model) || !model %in% names(labels))
    stop("W.model must be M1 (decoupled), M2 (shared-geometry), or M3 (residual-adjustment).")
  unname(labels[[model]])
}

# Surface weights for the piecewise-linear finite-element basis.
fem_surface_weights <- function(mesh) {
  as.numeric(Matrix::colSums(fmesher::fm_fem(mesh, order=1)$c0))
}

# Exact integral of exp(piecewise-linear log intensity) on a degree-one circle
# or interval. Includes the cyclic seam and remains stable under constant shifts.
circle_log_integral <- function(log.L, mesh, limits=mesh$interval) {
  if (mesh$degree != 1L || length(log.L) != length(mesh$loc) ||
      any(!is.finite(log.L))) stop("Expected finite nodal log intensities on a degree-one mesh.")
  if (length(limits)!=2L || any(!is.finite(limits)) ||
      limits[1] > limits[2] || limits[1] < mesh$interval[1] ||
      limits[2] > mesh$interval[2]) stop("Invalid angular integration interval.")
  if (limits[1] == limits[2]) return(-Inf)
  x <- sort(unique(c(limits, mesh$loc[mesh$loc>limits[1] & mesh$loc<limits[2]])))
  y <- as.vector(fmesher::fm_basis(mesh, loc=x) %*% log.L)
  a <- head(y,-1); b <- tail(y,-1); delta <- abs(b-a)
  ratio <- rep(1,length(delta))
  nz <- delta > 0
  ratio[nz] <- -expm1(-delta[nz])/delta[nz]
  terms <- log(diff(x)) + pmax(a,b) + log(ratio)
  offset <- max(terms)
  offset + log(sum(exp(terms-offset)))
}
