# Examples of $`\texttt{geometricExtremes}`$ functionalities

This file contains example-code to run the main functions of the $`\texttt{geometricExtremes}`$ package. 
As of now, the functions are available for 2- and 3-dimensional modelling,
with examples respectively in Sections [Examples 2D](#2d) and [Examples 3D](#3d).

<a id="2d"></a>
# Examples 2D

In this section, we detail how to perform Bayesian geometric inference for exceedances of $`\mathcal{Q}_q`$. The procedure to fit exponential exceedances is described [here](#2d-exp), and the one to fit generalised Pareto exceedances is described [here](#2d-GP).

<a id="2d-exp"></a>

## Exponential exceedances

### Simulated data and transformation to Laplace marginal distributions

``` r
rm(list=ls())
library(mvtnorm)
library(rmutil)
library(geometricExtremes)

# Sample from a bivariate Gaussian copula in Laplace margins
set.seed(44)
N.X <- 5000
mu <- c(0,0)
Sigma <- rbind(c(1,0.8),
               c(0.8,1))

X.N <- rmvnorm(N.X, mu, Sigma) # \boldsymbol{X}_N
``` 

Here, the marginal distributions of the vector `X.N` are known and so is the true transformation to Laplace margins. When the marginal distributions are unknown, one can use the $`\texttt{geometricExtremes}`$ function toLaplaceMargins as follows:

``` r
# Upper/lower row of upper/lower quantiles above/below which to fit GPD tails
q <- rbind(c(0.90,0.90),
           c(0.10,0.10))

# Perform transformation of marginal distributions of X.N to Laplace
LapTransf <- toLaplaceMargins(X.N,q)
LapTransf$X.L # Observed data transformed to Laplace margins
```

### Model fitting

The model fitting procedure involves specifying options and saving configurations respectively through the functions set.options and set.configs.

``` r 
# Set fitting options, see ?set.options for description of variables
options <- set.options(X                = LapTransf$X.L,
                       excess.dist.fam  = "E",
                       W.data           = "ExcOnly",
                       W.model          = "M3",
		       use.mean.Qq      = FALSE
                       q                = 0.9, 
                       N.Qq             = 20,
                       N.GW             = 50,
                       QR.prior.range   = c(1, 0.999),
                       QR.prior.sigma   = c(1, 0.8),
                       zeta.prior.range = c(0.5, 0.999),
                       zeta.prior.sigma = c(5, 0.8),
                       phi.prior.range  = c(0.5, 0.999),
                       phi.prior.sigma  = c(5, 0.8),
                       mesh.knots.2d    = seq(-pi,pi,by=pi/400),
                       seed             = 0L)

# Set fitting configurations
config <- set.configs(save.path = "path/to/output/folder/", # Path of folder to save fitted objects if save == T
                      file.nm   = paste0("GaussCop_LapMargins_",
		                         options$excess.dist.fam,"_",
					 options$W.model,"_",
					 options$W.data),
                      save      = FALSE, # Save fitted objects to save.path if save == T
                      progress  = FALSE) # Save progression in .txt file in save.path if progress == T
```

To obtain posterior realisations from $`\mathcal{Q}_q`$ and $`\mathcal{G}`$, one simply runs:

``` r
# Fit the quantile set Q_q
fitted.Qq <- fit_Qq(LapTransf$X.L,options,config,return_fitted_obj=F)

# Fit the sets G and L
fitted.mod <- fit_GL(fitted.Qq,config)
```
Running the lines below will plot the estimated $`\mathcal{Q}_q`$, $`\mathcal{G}`$ and $`\mathcal{W}`$  sets with simultaneous predictive intervals.

``` r
par(mfrow=c(1,3),mar=c(2,2,0,0),mgp=c(2.6,0.8,0),pty="s")
plot_Qq(fitted.Qq,xlim=c(-10,10),ylim=c(-10,10),by=4)
plot_G(fitted.mod)
plot_W(fitted.mod)
```

<p align="center"><img src="/figures/Plot_Qq_G_W.png" width="80%" height="80%"/> </p>

### Probability estimation

``` r
# Define rectangular region of interest
x_B <- c(3,5); y_B <- c(3,5)
B <- expand.grid(x=x_B,y=y_B)
S_B <- range(atan2(y=B[,2],x=B[,1]))

# Sample from the posterior distribution of G and L to obtain samples on S_B
post.samp <- sample_QGW_posterior_2d(fitted.mod,N.w=1000,S_B=S_B,transf.G=F)

# Perform Bayesian probability estimation
ps <- prob_estimation(fitted.mod,post.samp,x_B,y_B)
median(ps[,2])
```

### Return level-sets estimation

``` r
par(mfrow=c(1,2),mar=c(2,2,0,0),mgp=c(2.6,0.8,0),pty="s")
return_sets <- return_set(fitted.mod,q.prime=c(0.95,0.99,0.999),include.Qq = FALSE)
plot_X_t(fitted.mod,return_sets,plt="boundary",xylim=c(-10,10))
plot_X_t(fitted.mod,return_sets,plt="set",xylim=c(-10,10))
```

<p align="center"><img src="/figures/Ret_sets.png" width="70%" height="70%"/> </p>

``` r
par(mfrow=c(1,2),mar=c(2,2,0,0),mgp=c(2.6,0.8,0),pty="s")
return_sets <- return_set(fitted.mod,q.prime=c(0.95,0.99,0.999),include.Qq = FALSE, LapTransf = LapTransf)
plot_X_t(fitted.mod,return_sets,plt="boundary",xylim=c(-4,4))
plot_X_t(fitted.mod,return_sets,plt="set",xylim=c(-4,4))
```


<p align="center"><img src="/figures/Ret_sets_orig_margins.png" width="70%" height="70%"/> </p>


### Measures of extremal dependence

#### Conditional Extremes parameter $`\alpha`$

Below is a function to obtain posterior samples from the Conditional Extremes parameter $`\alpha`$.

``` r
alphas <- alpha_posterior(fitted.mod)
mean(alphas[,1])
mean(alphas[,2])
```

#### Coefficient of residual tail dependence $`\eta`$

Below is a function to obtain posterior samples from the coefficient of residual tail dependence $`\eta`$.

``` r
etas <- eta_posterior(fitted.mod)
mean(etas)
```

### Diagnotics $`K_B`$ and $`K_C`$

Below is a function to assess the quality of fit of the model.

``` r
# Transform exceedances of Q_q^\star to the unit ball
U_on_ball <- X_to_uniform_on_Ball(fitted.mod)

# Thin U_on_ball to obtain same number of observation per threshold sample 
thin_U <- thin_U_on_ball(U_on_ball)

# Create the credible envelopes from a truly uniform sample
d    <- ncol(fitted.mod$X)
N    <- nrow(thin_U[[1]][[1]])
s <- seq(0, 2, len=30)
K.CI_B <- K_envelope(n=N, d=d, M=1000, subset="ball", s=s)
phi <- seq(0, pi, len=30)
K.CI_C <- K_envelope(n=N, d=d, M=1000, subset="ball", phi=phi)

# Estimate the K functions from the ball geometry
K.hat_B <- K_hat(thin_U,"ball")

# Estimate the K functions from the spherical cone geometry
K.hat_C <- K_hat(thin_U,"sph cone")
```

Below are plots of the estimated $`K_B`$ and $`K_C`$ (in light grey) with 
credible envelope (in red).
``` r
par(mfrow=c(1,2),pty="s",mar=c(4.5,4.5,1,4.5))
matplot(s, K.hat_B, type="l", lty=1, col="grey60", lwd=1, xlab="distance (Ball)", ylab = "K",frame.plot=F)
matplot(s, t(K.CI_B$envelope), type="l", lty=1, col="red", lwd=2, add=TRUE)

matplot(phi, K.hat_C, type="l", lty=1, col="grey60", lwd=1, xlab="distance (Spherical cone)", ylab = "K",frame.plot=F)
matplot(phi, t(K.CI_C$envelope), type="l", lty=1, col="red", lwd=2, add=TRUE)
```
<p align="center"><img src="/figures/K_diagnostics.png" width="70%" height="70%"/> </p>

<a id="2d-GP"></a>

## Generalised Pareto exceedances

### Simulated data in heavy-tailed margins

``` r
rm(list=ls())
library(QRM)
library(geometricExtremes)

# Sample from a bivariate student-t copula in student-t margins
set.seed(44)
N.X <- 5000
df <- 3
Sigma <- rbind(c(1,-0.7),
               c(-0.7,1))

X.t <- qt(rcopula.t(N.X, df=df, Sigma=Sigma),df=df)
```

The model fitting procedure involves specifying options and saving configurations respectively through the functions set.options and set.configs.
A key difference with the previous section on exponential exceedances, is that we now require the additional parameter $`\alpha`$. It arises through the inla parameterisation of the generalised Pareto fitting, see the [inla GPd documentation](https://inla.r-inla-download.org/r-inla.org/doc/likelihood/genPareto.pdf).

### Model fitting

``` r
# Set fitting options, see ?set.options for description of variables
options <- set.options(X                = X.t,
                       excess.dist.fam  = "GP",
                       W.data           = "ExcOnly",
                       W.model          = "M3",
                       q                = 0.9,
		       use.mean.Qq      = FALSE,
		       alpha            = 0.5, # GP specific parameter
                       N.Qq             = 20,
                       N.GW             = 50,
                       QR.prior.range   = c(1, 0.999),
                       QR.prior.sigma   = c(1, 0.8),
                       zeta.prior.range = c(0.5, 0.999),
                       zeta.prior.sigma = c(5, 0.8),
                       phi.prior.range  = c(0.5, 0.999),
                       phi.prior.sigma  = c(5, 0.8),
                       mesh.knots.2d    = seq(-pi,pi,by=pi/400),
                       seed             = 0L)

# Set fitting configurations
config <- set.configs(save.path = "path/to/output/folder/", # Path of folder to save fitted objects if save == T
                      file.nm   = paste0("tCop_tMargins_",
                                         options$excess.dist.fam,"_",
                                         options$W.model,"_",
                                         options$W.data),
                      save      = FALSE, # Save fitted objects to save.path if save == T
                      progress  = FALSE) # Save progression in .txt file in save.path if progress == T
```

To obtain posterior realisations from $`\mathcal{Q}_q`$ and $`\mathcal{G}`$, one again simply runs:

``` r
# Fit the quantile set Q_q
fitted.Qq <- fit_Qq(X.t,options,config,return_fitted_obj=F)

# Fit the sets G and L
fitted.mod <- fit_GL(fitted.Qq,config)
```
Running the lines below will plot the estimated $`\mathcal{Q}_q`$, $`\mathcal{G}`$, and $`\mathcal{W}`$ sets with simultaneous predictive intervals.

``` r
par(mfrow=c(1,3),mar=c(2,2,0,0),mgp=c(2.6,0.8,0),pty="s")
plot_Qq(fitted.Qq,xlim=c(-10,10),ylim=c(-10,10),by=4)
plot_G(fitted.mod)
plot_W(fitted.mod)
```

<p align="center"><img src="/figures/Plot_Qq_G_W_gp.png" width="80%" height="80%"/> </p>

### Probability estimation

``` r
# Define rectangular region of interest
x_B <- c(-15,-10); y_B <- c(10,15)
B <- expand.grid(x=x_B,y=y_B)
S_B <- range(atan2(y=B[,2],x=B[,1]))

# Sample from the posterior distribution of G and L to obtain samples on S_B
post.samp <- sample_QGW_posterior_2d(fitted.mod,N.w=1000,S_B=S_B,transf.G=F)

# Perform Bayesian probability estimation
ps <- prob_estimation(fitted.mod,post.samp,x_B,y_B)
median(ps[,2])
```


<a id="3d"></a>

# Examples 3D

## Exponential exceedances

### Simulated data and transformation to Laplace marginal distributions

``` r
rm(list=ls())
library(mvtnorm)
library(rmutil)
library(geometricExtremes)

set.seed(44)
mu <- c(0,0,0)
Sigma <- rbind(c(1,0.8,0.7), c(0.8,1,0.6),c(0.7,0.6,1))
X.N <- rmvnorm(5000, mu, Sigma)
```

``` r
# Upper/lower row of upper/lower quantiles above/below which to fit GPD tails
q <- rbind(c(0.90,0.90,0.90),
           c(0.10,0.10,0.10))

# Perform transformation of marginal distributions of X.N to Laplace
LapTransf <- toLaplaceMargins(X.N,q)
LapTransf$X.L # Observed data transformed to Laplace margins
```


### Model fitting

``` r
options <- set.options(X                = LapTransf$X.L,               
                       excess.dist.fam  = "E",       
                       W.data           = "ExcOnly",                
                       W.model          = "M3",    
                       use.mean.Qq      = FALSE,             
                       q                = 0.9,              
                       N.Qq             = 20,
                       N.GW             = 50,
                       QR.prior.range   = c(2, 0.999),
                       QR.prior.sigma   = c(1, 0.8),
                       zeta.prior.range = c(1, 0.999),
                       zeta.prior.sigma = c(2, 0.8),
                       phi.prior.range  = c(1, 0.999),
                       phi.prior.sigma  = c(1, 0.8),
                       mesh.res.3d      = 20,
                       seed             = 0L)


config <- set.configs(save.path = "path/to/folder/",
                      file.nm   = paste0("GaussCop3d_LapMargins_",options$excess.dist.fam,"_",options$W.model,"_",options$W.data),
                      save      = FALSE,
                      progress  = FALSE)
```

To obtain posterior realisations from $`\mathcal{Q}_q`$ and $`\mathcal{G}`$, one runs:

``` r
# Fit the quantile set Q_q
fitted.Qq <- fit_Qq(LapTransf$X.L,options,config,return_fitted_obj=F)

# Fit the sets G and L
fitted.mod <- fit_GL(fitted.Qq,config)
```

To plot the mean posterior sets $`\mathcal{Q}_q`$ and $`\mathcal{G}`$, one runs: 

``` r
plot_Qq(fitted.Qq,surface3d="mean")
plot_G(fitted.mod,surface3d="mean")
```

To obtain the lower and upper surfaces of the ($`1-\alpha`$)-simultaneous predictive intervals of $`\mathcal{Q}_q`$ and $`\mathcal{G}`$, it suffices to set, for instance:

``` r
plot_Qq(fitted.Qq,alpha=0.05,surface3d="lower")
plot_G(fitted.mod,alpha=0.05,surface3d="upper")
```


## Generalised Pareto exceedances

Available soon.

