# Examples of $`\texttt{geometricExtremes}`$ functionalities

This file contains example-code for running the main functions of the $`\texttt{geometricExtremes}`$ package. As of now, the functions are available for 2- and 3- dimensional modelling, with examples respectively in Sections [2D](#2d) and [3D](#3d).

<a id="2d"></a>
# Example 2D

## Simulated data and transformation to Laplace marginal distributions

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

Here, the marginal distributions of the vector $`\boldsymbol{X}_N`$ are known and so is the true transformation to Laplace margins. When the marginal distributions are unknown, one can use the $`\texttt{geometricExtremes}`$ function toLaplaceMargins as follows:

``` r
# Upper/lower row of upper/lower quantiles above/below which to fit GPD tails
q <- rbind(c(0.90,0.90),
           c(0.10,0.10))

# Perform transformation of marginal distributions of X.N to Laplace
X <- toLaplaceMargins(X.N,q)
X$X.L # Observed data transformed to Laplace margins
```

## Model fitting

The model fitting procedure involves specifying options and saving configurations respectively through the functions set.options and set.configs.

``` r 
library(geometricExtremes)

# Set fitting options, see ?set.options for description of variables
options <- set.options(X                = X$X.L,
                       excess.dist.fam  = "E",
                       W.data           = "ExcOnly",
                       W.model          = "M3",
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
                      save      = TRUE, # Save fitted objects to save.path if save == T
                      progress  = TRUE) # Save progression in .txt file in save.path if progress == T
```
To obtain posterior realisations from $`\mathcal{Q}_q`$ and $`\mathcal{G}`$, one simply runs:

``` r
# Fit the quantile set Q_q
fitted.Qq <- fit_Qq(X$X.L,options,config,return_fitted_obj=F)

# Fit the sets G and L
fitted.mod <- fit_GL(fitted.Qq,config)
```
Running the lines below will plot the estimated $`\mathcal{Q}_q`$ and $`\mathcal{G}`$ sets with simultaneous predictive intervals.

``` r
par(mfrow=c(1,2),mar=c(2,2,0,0),mgp=c(2.6,0.8,0),pty="s")
plot_Qq(fitted.Qq,cex.pts=0.4,cex.axis=1.4,xlim=c(-8,8),ylim=c(-8,8),by=4)
plot_G(fitted.mod,by=4)
```

<p align="center"><img src="/figures/Plot_Qq_G.png" width="70%" height="70%"/> </p>

## Probability estimation

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

## Return level-sets estimation

``` r
par(mfrow=c(1,2),mar=c(2,2,0,0),mgp=c(2.6,0.8,0),pty="s")
list_ret_sets <- return_set(fitted.mod,q.prime=c(0.95,0.99,0.999))
plot_X_t(fitted.mod,list_ret_sets,plt="boundary",xylim=c(-10,10))
plot_X_t(fitted.mod,list_ret_sets,plt="set",xylim=c(-10,10))
```

<p align="center"><img src="/figures/Ret_sets.png" width="70%" height="70%"/> </p>


## Measures of extremal dependence

### Conditional Extremes parameter $`\alpha`$

Below is a function to obtain posterior samples from the Conditional Extremes parameter $`\alpha`$.

``` r
alphas <- alpha_posterior(fitted.mod)
mean(alphas[,1])
mean(alphas[,2])
```

### Coefficient of residual tail dependence $`\eta`$

Below is a function to obtain posterior samples from the coefficient of residual tail dependence $`\eta`$.

``` r
etas <- eta_posterior(fitted.mod)
mean(etas)
```

<a id="3d"></a>

# Example 3D
``` r
rm(list=ls())
library(mvtnorm)
library(rmutil)

mu <- c(0,0,0)
Sigma <- rbind(c(1,0.8,0.7), c(0.8,1,0.6),c(0.7,0.6,1))
X <- rmutil::qlaplace(pnorm(rmvnorm(5000, mu, Sigma),0,1))

#####################
### Model fitting ###
#####################

options <- set.options(X                = X,               
                       excess.dist.fam  = "E",       
                       W.data           = "ExcOnly",                
                       W.model          = "M3",    
                       use.mean.Qq      = TRUE,             
                       q                = 0.9,              
                       N.Qq             = 1,
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

fitted.Qq <- fit_Qq(X,options,config,return_fitted_obj=F)
plot_Qq(fitted.Qq,cex.pts=0.4,cex.axis=1.4,xlim=c(-8,8),ylim=c(-8,8),by=4)

fitted.mod <- fit_GL(fitted.Qq,config)
plot_G(fitted.mod,surface3d="mean",surf.col = "red")
```

