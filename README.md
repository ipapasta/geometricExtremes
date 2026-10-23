# Geometric multivariate extreme value analysis 

> A Bayesian approach to geometric extreme value analysis based on:
>
> Papastathopoulos, I., De Monte, L., Campbell, R. and Rue, H. (2023) 'Statistical inference for radially-stable generalized Pareto distributions and return level-sets in geometric extremes', [arXiv](https://arxiv.org/abs/2310.06130) preprint.

<p align="center"><img src="/figures/GaussCop_LapMargins.gif" width="60%" height="60%"/> </p>



## Installation
``` r
# install.packages("remotes")
remotes::install_github("ipapasta/geometricExtremes")
```

## Directional models and migration

From version 0.0.0.9001, model labels follow the revised manuscript:

| Label | Name | Directional body |
|---|---|---|
| M1 | decoupled | L = B |
| M2 | shared-geometry | L = G |
| M3 | residual-adjustment | L = BG |

**Breaking label change:** older package versions used M1 for shared geometry
and M2 for the decoupled model. Swap M1/M2 when reproducing old analyses;
M3 is unchanged. Do not reinterpret historical saved objects or filenames.
New tail fits record `model.convention = "manuscript-v1"` and `model.name`.

M1 and M3 constrain the finite-element surface integral of `phi` to zero,
while estimating `Int.W` freely. The implementation uses `phi = -d*log(r_B)`;
`phi.prior.sigma` describes the phi scale, not the log-radial scale. Existing
prior hyperparameters are retained. Constraints change the prior and may change
posterior uncertainty; they do not by themselves guarantee calibrated intervals.

Use `fit_GL(..., return_fitted_obj=TRUE, bru_options=list(...))` to retain the
underlying fits and investigate solver/prior sensitivity. Tail fitting and
posterior sampling default to one thread. `log.L` contains unnormalised
log-intensities, not normalised directional densities.

## Examples

Examples and explanations are found [here](EXAMPLES.md).

## Citing $`\texttt{geometricExtremes}`$

Please cite this package using this [BibTeX](Citation.bib) file.

## Reporting a bug

Please report any bugs or comments [here](https://github.com/ipapasta/geometricExtremes/issues).

## References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-Papastathopoulos2023" class="csl-entry">

Papastathopoulos, I., De Monte, L., Campbell, R. and Rue, L. (2023) Statistical inference for radially-stable generalized Pareto distributions and return level-sets in geometric extremes, arXiv preprint. <https://arxiv.org/abs/2310.06130>. 

</div>

<div id="ref-Rueetal2009" class="csl-entry">

Rue, H., Martino, S. and Chopin, N. (2009), Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 71: 319-392. <https://doi.org/10.1111/j.1467-9868.2008.00700.x>.

</div>

<div id="ref-BachletAl2019" class="csl-entry">

Bachl F, Lindgren F, Borchers D, Illian J (2019). inlabru: an R package for Bayesian spatial modelling from ecological survey data. Methods in Ecology and Evolution, 10, 760–766. <https://doi.org/10.1111/2041-210X.13168>.

</div>

<div id="ref-BolinLindgren2018" class="csl-entry">

Bolin D, Lindgren F (2018). Calculating Probabilistic Excursion Sets and Related Quantities Using excursions. Journal of Statistical Software, 86(5), 1–20. <https://doi.org/10.18637/jss.v086.i05>.

</div>

</div>
