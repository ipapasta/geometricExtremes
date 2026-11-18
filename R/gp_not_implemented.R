# Placeholders for the radial GP/EGP implementation.
#
# The GP extension (commit 05d3309) calls these helpers, but their definitions
# were never committed. Until they are restored, set.options() rejects
# excess.dist.fam = "GP"/"EGP", and these stubs give a clear error if a GP fit
# reaches them by another route. Restoring the GP code means deleting this file.

gp_not_implemented <- function(family = "GP/EGP") {
  stop(sprintf("excess.dist.fam = '%s' is not implemented yet; use excess.dist.fam = 'E'.",
               paste(family, collapse = "/")), call. = FALSE)
}

gp_store_draws      <- function(...) gp_not_implemented()
gp_sample_posterior <- function(...) gp_not_implemented()
gp_return_sets      <- function(...) gp_not_implemented()
gp_uniform_ball     <- function(...) gp_not_implemented()
gp_box_probability  <- function(...) gp_not_implemented()
gp_radial_pp        <- function(...) gp_not_implemented()
gp_directional_pp   <- function(...) gp_not_implemented()
