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
  $\int_{\mathbb S^{d-1}}\phi(\bm w)\,d\bm w=0.$
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
  $\xi$. Calculations enforce
  $
  1+\xi\{r-r_{\mathcal Q_q}(\bm w)\}/r_{\Sigma}(\bm w)>0.
  $
  and do not silently clip observations or posterior draws at the finite radial
  endpoint.

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
