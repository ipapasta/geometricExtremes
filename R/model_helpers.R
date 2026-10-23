# Internal model vocabulary. Numeric labels follow the revised manuscript.
directional_architecture <- function(model) {
  labels <- c(M1="decoupled", M2="shared-geometry", M3="residual-adjustment")
  if (length(model) != 1L || is.na(model) || !model %in% names(labels))
    stop("W.model must be M1 (decoupled), M2 (shared-geometry), or M3 (residual-adjustment).")
  unname(labels[[model]])
}

# Identify the two GP likelihood implementations. "EGP" is experimental and
# fixes the extended-GP body-shape parameter kappa at one.
is_gp_excess_family <- function(family) {
  identical(family, "GP") || identical(family, "EGP")
}

# Stable conversion from the fitted alpha-quantile to the GP scale:
# sigma = q_alpha * xi / {(1-alpha)^(-xi)-1}. The direct expression is 0/0
# at xi=0 and loses precision near zero.
gp_quantile_to_scale <- function(q_alpha, xi, alpha) {
  if (length(xi) != 1L || !is.finite(xi) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("xi must be finite and alpha must lie in (0,1).")
  q_alpha <- as.numeric(q_alpha)
  if (any(!is.finite(q_alpha)) || any(q_alpha <= 0))
    stop("GP quantiles must be finite and positive.")
  a <- -log1p(-alpha)
  factor <- if (abs(xi) < 1e-8) {
    # xi/expm1(a*xi) = 1/a - xi/2 + a*xi^2/12 + O(xi^4)
    1/a - xi/2 + a*xi^2/12
  } else {
    xi/expm1(a*xi)
  }
  q_alpha * factor
}

# Extract xi by semantic name. This remains correct for EGP fits even if other
# hyperparameters are present or INLA changes their ordering.
extract_gp_tail <- function(fit) {
  tab <- fit$summary.hyperpar
  if (is.null(tab) || !nrow(tab)) stop("The fitted model has no GP tail summary.")
  idx <- grep("^Tail parameter for (the )?(e)?gp observations(\\[[0-9]+\\])?$",
              rownames(tab), ignore.case=TRUE)
  if (length(idx) != 1L)
    stop("Could not uniquely identify the GP tail parameter in the INLA fit.")
  list(mean=unname(tab[idx,"mean"]), summary=tab[idx,,drop=FALSE],
       marginal=fit$marginals.hyperpar[[rownames(tab)[idx]]])
}

# INLA control for the experimental signed-GP likelihood. EGP with kappa=1 is
# exactly GP, while permitting xi below zero.
egp_kappa1_control <- function(alpha, interval, prior_rate) {
  tail <- INLA::inla.models()$likelihood$egp$hyper$theta1
  list(
    control.link=list(model="quantile", quantile=alpha),
    hyper=list(
      tail=list(initial=tail$to.theta(0, interval=interval), fixed=FALSE,
                prior="pc.egptail", param=c(prior_rate, interval)),
      shape=list(initial=0, fixed=TRUE)
    )
  )
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
