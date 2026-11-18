# geometricExtremes (development version)

- `prob_estimation()` gains a `theta.post` argument (d = 2,
  `excess.dist.fam = "E"`), implementing eq. (11) of the directional extremal
  index manuscript. Given posterior samples of theta from `fit_theta()` (fitted
  on the same `fit_Qq()` object, with `N.theta >= N.GW`), the radial integral
  is weighted by theta at the sampled directions, pairing posterior sample j of
  theta with posterior sample j of the geometric model for each sample of Q_q.
  It returns posterior samples of the probability `P` that a time step is the
  radial maximum of a cluster falling in the set, and of the return period
  `T = 1/P` in time steps. It warns if `fitted.mod` was not fitted to cluster
  radial maxima. With `theta.post = NULL` (default) the output is unchanged.

- `fit_GL()` gains a `radial.idx` argument (d = 2, `excess.dist.fam = "E"`).
  Passing the output of `get_radial_maxima()` fits the radial likelihood to
  the cluster radial maxima only, while the directional likelihood still uses
  `W.data` (all exceedances, or all angles, which target the same directional
  distribution for a quantile set). It requires `W.model = "M1"`; `"M3"` is
  allowed with a warning (the field zeta is shared by both likelihoods) and
  `"M2"` is rejected. With `radial.idx = NULL` (default) the fit is unchanged.

- New functions `get_filtered_Z()`, `fit_theta()`, `plot_theta()` and
  `plot_mean_cluster_size()`, ported from TemporalGeometricClustering, for the
  directional extremal index theta(w) and the mean cluster size 1/theta(w).
  The model is unchanged. Posterior samples of theta are now grouped by
  posterior sample of Q_q, `fit_theta()` uses one INLA thread and can return
  the fitted inlabru objects, the code uses the fmesher API and
  `generate(x = )` (the original no longer ran on current inlabru), and the
  plotting functions no longer rely on global variables. Plotting is
  available for d = 2 only.

- New functions `get_exceedances()`, `get_clusters_runs()` and
  `get_radial_maxima()`, ported from TemporalGeometricClustering. For each
  posterior sample of Q_q (or its posterior mean), they return the exceedance
  times, the runs clusters of exceedances, and the radial maximum of each
  cluster. Compared with the original code, the radial maximum is taken over
  the exceedances of a cluster only (with a direction-dependent threshold, a
  non-exceedance inside a cluster can have a larger radius), and
  `use.mean.Qq = TRUE` is now honoured.

- `excess.dist.fam = "GP"` and `"EGP"` now stop with a "not implemented yet"
  error in `set.options()`. The GP extension referenced helper functions whose
  definitions were never committed, so GP/EGP fits failed in `fit_GL()`.
  Placeholders live in `R/gp_not_implemented.R`; deleting that file and the
  guard in `set.options()` restores the GP path once the helpers are recovered.

- `return_set_2d_isotropic()` and `return_set_3d_isotropic()` now warn, once
  per return period, when posterior draws give a negative isotropic radius in
  some direction. The 2D function previously printed draw indices to the
  console; the 3D function did not check at all.

- `set.options()` now uses the default values of the prior arguments
  (`QR.prior.*`, `zeta.prior.*`, `phi.prior.*`), `mesh.knots.2d` and
  `mesh.res.3d` when they are not supplied. Previously it stopped with
  "Specify value of ..." even though defaults were declared.

- `set.configs()` now defaults to `file.nm = "DataName"`. The previous default
  referred to an options list that is not available inside the function, so
  calling `set.configs()` without `file.nm` always failed.

- The error for a missing `N.GW` in `set.options()` now names `N.GW`
  (previously it said `N.Qq`).

- Replaced the deprecated `INLA::inla.mesh.projector()` and
  `INLA::inla.mesh.1d()` with `fmesher::fm_basis()` and `fmesher::fm_mesh_1d()`
  throughout. The projection matrices are identical; this only removes the
  deprecation warnings.

- Fixed direction sampling in 3D (`sample_QGW_posterior()` for d = 3). It
  called `inlabru::sample.lgcp(..., ignore.CRS = TRUE)`, which makes inlabru
  treat the sphere mesh as planar and return 2D points; with current inlabru it
  also failed outright on the returned `sf` object. Sampling now keeps the
  sphere CRS, uses unit-sphere areas (`R = 1`) and returns unit vectors. `sf`
  is added to Imports.

- Documented the `q` and `n.integration` arguments of
  `X_to_uniform_on_Ball()` and its 2D/3D variants, and gave the exported
  2D/3D helper functions proper help-page titles.

# geometricExtremes 0.2.0

## Model definitions and identifiability

- Model labels now follow the revised manuscript:
  - **M1 (decoupled):** $\mathcal L=\mathcal B$;
  - **M2 (shared-geometry):** $\mathcal L=\mathcal G$;
  - **M3 (residual-adjustment):** $\mathcal L=\mathcal B\mathcal G$.

  Earlier package versions used M1 for the shared-geometry model and M2 for the
  decoupled model; M3 is unchanged.

- In M1 and M3, the residual directional field is identified by the
  surface-weighted finite-element constraint
  $\int_{\mathbb S^{d-1}}\phi(\boldsymbol w)\,d\boldsymbol w=0.$
  while the auxiliary Poisson intercept `Int.W` remains unconstrained. The
  constraint uses finite-element integration weights on the circle for $d=2$
  and on the sphere for $d=3$.

- The 2D tail fit now reuses the mesh of the fitted quantile set
  $\mathcal Q_q$, so threshold and tail quantities share the same
  finite-element basis ordering.

## Signed radial GP model

- Added `excess.dist.fam = "EGP"`. This uses INLA's EGP likelihood with
  $\kappa=1$, giving a GP excess model with one global signed shape parameter
  $\xi\in(-0.5,0.5)$ by default. The admissible interval and PC-prior rate can
  be changed through `egp.xi.interval` and `egp.prior.rate`.

- Posterior draws of $\xi$, the radial scale set and the directional field are
  now sampled jointly. The fit stores draw-specific GP scales in `Sigma` and
  retains `G` as a compatibility alias. The posterior summary and marginal
  distribution of $\xi$ are stored in `xi.summary` and `xi.marginal`.

- Conversion from the fitted conditional GP quantile to the radial scale is
  evaluated stably near $\xi=0$, where the GP converges continuously to the
  exponential distribution.

- Added GP threshold stability and finite-support handling for negative
  $\xi$. 
- Added support-aware functionality for:
  - canonical and isotropic probability and return sets in dimensions two and
    three;
  - conditional simulation from the fitted GP tail;
  - radial and directional PP calculations;
  - the uniform-ball transformation and its existing diagnostics;
  - rectangular event probabilities through the existing 2D interface.

- Isotropic boundaries are required to remain outside the fitted quantile set
  $\mathcal Q_q$. The package reports an informative error when a requested
  level would require extrapolation below the fitted tail threshold.

- Saved GP fits from earlier package versions do not contain paired draws of
  $\xi$ and the latent fields. They must be refitted before using the new GP
  prediction, return-set or diagnostic methods.

- General event-probability integration in three dimensions is not yet
  implemented; this is unchanged from the previous package interface.

## Directional integration and sampling

- Directional probabilities on $\mathbb S^1$ now integrate the exponential
  of the fitted piecewise-linear log-intensity analytically on each mesh
  segment. This is exact for the fitted finite-element representation, while
  the representation itself remains mesh-dependent.

- Corrected periodic seam handling at $-\pi$ and $\pi$, rejection sampling
  on angular intervals containing no internal mesh nodes, and use of the
  circle's surface measure $2\pi$.

- Directional probabilities, isotropic sets, directional-density plots and the
  uniform-ball transformation now use a consistent directional normalisation.

## Plotting

- `plot_G()` now uses the same 2D geometric layout for radial exponential and
  radial GP scale sets: scaled observations, angular direction labels,
  uncertainty bands and the $[-1,1]^2$ standard-Laplace reference box are
  displayed consistently.

- For GP/EGP fits, the plotted object is the radial scale set $\Sigma$. Under
  the radial exponential model, this scale set is identified with the geometric
  limit set $\mathcal G$.

## Compatibility and robustness

- Updated posterior sampling to the current `inlabru::generate(x = ...)` API,
  modernised likelihood construction with `inlabru::bru_obs()`, and updated the
  3D mesh and mapper calls for current `fmesher` and `inlabru` versions.

- Added checks for non-finite observations, zero vectors, thresholds with no
  radial exceedances, inadmissible GP support and below-threshold isotropic
  boundaries.

- Added automated consistency tests for model labels, finite-element
  constraints, directional normalisation, periodic integration, circular
  sampling and the common E/GP plotting layout.
