# geometricExtremes 0.2.0

## Model definitions and identifiability

- Model labels now follow the revised manuscript:
  - **M1 (decoupled):** \(\mathcal L=\mathcal B\);
  - **M2 (shared-geometry):** \(\mathcal L=\mathcal G\);
  - **M3 (residual-adjustment):** \(\mathcal L=\mathcal B\mathcal G\).

  Earlier package versions used M1 for the shared-geometry model and M2 for the
  decoupled model; M3 is unchanged.

- In M1 and M3, the residual directional field is now identified by imposing
  the surface-weighted finite-element constraint
  \[
  \int_{\mathbb S^{d-1}}\phi(\bm w)\,d\bm w=0,
  \]
  while leaving the auxiliary Poisson intercept `Int.W` unconstrained. The
  constraint uses finite-element integration weights on the circle for \(d=2\)
  and on the sphere for \(d=3\).

- The 2D tail fit now reuses the mesh of the fitted quantile set
  \(\mathcal Q_q\), ensuring that threshold and tail quantities have the same
  finite-element basis ordering.

## Radial excess models

- Added the experimental option `excess.dist.fam = "EGP"`. It uses INLA's EGP
  likelihood with \(\kappa=1\), which gives the GP family while permitting a
  global signed GP shape parameter \(\xi\). The fitted object stores the
  posterior summary and marginal distribution of \(\xi\).

- Conversion from the fitted conditional GP quantile to the radial scale is
  evaluated with a numerically stable formula near \(\xi=0\), where the GP
  converges continuously to the exponential distribution.

- Probability and return-set methods for `EGP` remain unavailable in this
  release. Supporting \(\xi<0\) requires the radial GP stability map and
  enforcement of the directional support
  \[
  1+\xi\{\lVert\bm x\rVert-r_{\mathcal Q_q}(\bm w)\}/r_{\mathcal G}(\bm w)>0.
  \]
  The current downstream implementations use the radial exponential
  transformation and therefore cannot simply be reused when \(\xi\ne0\).

## Directional integration and sampling

- Directional probabilities on \(\mathbb S^1\) now integrate the exponential
  of the fitted piecewise-linear log-intensity analytically on every mesh
  segment. This is exact for the finite-element representation, although it
  does not remove mesh approximation error.

- Corrected periodic seam handling at \(-\pi\) and \(\pi\), rejection sampling
  on angular intervals containing no internal mesh nodes, and the use of the
  circle's surface measure \(2\pi\).

- The same directional normalisation is now used consistently in directional
  probabilities, isotropic probability and return sets, directional-density
  plots, and the transformation to the unit ball.

## Plotting

- `plot_G()` now uses the same 2D geometric layout for radial exponential and
  radial GP scale sets: scaled observations, angular direction labels,
  uncertainty bands, and the \([-1,1]^2\) standard-Laplace reference box are
  displayed consistently.

- The plotted object for a GP/EGP fit should be interpreted as the radial scale
  set. Under the radial exponential model used in the manuscript, this scale
  set is identified with the geometric limit set \(\mathcal G\).

## Compatibility and robustness

- Updated posterior sampling to the current `inlabru::generate(x = ...)` API,
  modernised likelihood construction with `inlabru::bru_obs()`, and updated
  the 3D mesh and mapper calls for current `fmesher` and `inlabru` versions.

- Added input checks for non-finite observations, zero vectors, and threshold
  draws with no radial exceedances.

- Added automated consistency tests for:
  - the M1/M2/M3 definitions;
  - the surface-weighted finite-element constraint;
  - normalisation of the directional density;
  - periodic integration and circular sampling;
  - the common E and GP plotting layout.
