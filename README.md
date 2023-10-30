# Geometric multivariate extreme value analysis 

> A Bayesian approach to geometric extreme value analysis based on:
>
> I. Papastathopoulos, L. De Monte, R. Campbell, and H. Rue. (2023) 'Statistical inference for radially-stable generalized Pareto distributions and return level-sets in geometric extremes', [arXiv](https://arxiv.org/abs/2310.06130) preprint.

<p align="center"><img src="/GaussCop_LapMargins.gif" width="60%" height="60%"/> </p>



## Installation
``` r
# install.packages("remotes")
remotes::install_github("ipapasta/NAMEofPACKAGE")
```
## Example 2D
``` r 
rm(list=ls())
library(mvtnorm)
library(rmutil)

# Sample from a bivariate Gaussian copula in Laplace margins
set.seed(44)
mu <- c(0,0)
Sigma <- rbind(c(1,0.8), c(0.8,1))
X <- rmutil::qlaplace(pnorm(mvtnorm::rmvnorm(N.X, mu, Sigma),0,1))

#####################
### Model fitting ###
#####################

# Set fitting options, see ?set.options for description of variables
options <- set.options(X                = X,
                       data.nm          = "GaussCop_LaplaceMargins",
                       excess.dist.fam  = "E",
                       W.data           = "ExcOnly",
                       W.model          = "M3",
                       q                = 0.9,
                       alpha            = 0.5,
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
                      file.nm   = paste0(options$data.nm,"_",
		                         options$excess.dist.fam,"_",
					 options$W.model,"_",
					 options$W.data),
                      save      = TRUE, # Save fitted objects to save.path if save == T
                      progress  = TRUE) # Save progression in .txt file in save.path if progress == T

par(mfrow=c(1,2),mar=c(4.1,4.1,0.5,1.1),mgp=c(2.6,0.8,0),pty="s")

# Fit the quantile set Q_q
fitted.Qq <- fit_Qq(X,options,config,return_fitted_obj=F)
plot_Qq(fitted.Qq,cex.pts=0.4,cex.axis=1.4,xlim=c(-8,8),ylim=c(-8,8),by=4)

# Fit the sets G and L
fitted.mod <- fit_GL(fitted.Qq,config)
plot_G(fitted.mod,by=4)

##############################
### Probability estimation ###
##############################

# Define rectangular region of interest
x_B <- c(3,5); y_B <- c(3,5)
B <- expand.grid(x=x_B,y=y_B)
S_B <- range(atan2(y=B[,2],x=B[,1]))

# Sample from the posterior distribution of G and L to obtain samples on S_B
post.samp <- sample_QGW_posterior_2d(fitted.mod,N.w=1000,S_B=S_B,transf.G=F)

# Perform Bayesian probability estimation
ps <- prob_estimation(fitted.mod,post.samp,x_B,y_B)
median(ps[,2])

####################################
### Return level-sets estimation ###
####################################

list_ret_sets <- return_set(fitted.mod,q.prime=c(0.95,0.99,0.999))
plot_boundary_X_t(fitted.mod,list_ret_sets,xylim=c(-10,10))

```

## Example 3D
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
                       data.nm          = "Gauss",               
                       excess.dist.fam  = "E",       
                       W.data           = "ExcOnly",                
                       W.model          = "M3",    
                       use.mean.Qq      = TRUE,             
                       q                = 0.9,
                       alpha            = 0.5,                   
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
                      file.nm   = paste0(options$data.nm,"_",options$excess.dist.fam,"_",options$W.model,"_",options$W.data),
                      save      = FALSE,
                      progress  = FALSE)

fitted.Qq <- fit_Qq(X,options,config,return_fitted_obj=F)
plot_Qq(fitted.Qq,cex.pts=0.4,cex.axis=1.4,xlim=c(-8,8),ylim=c(-8,8),by=4)

fitted.mod <- fit_GL(fitted.Qq,config)
plot_G(fitted.mod,surface3d="mean",surf.col = "red")
```

