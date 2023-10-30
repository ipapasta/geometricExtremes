#' Set options used by functions fit_Qq and fit_GL
#'
#' @param X n by d observations matrix.
#' @param excess.dist.fam Exponential ("E") or generalised Pareto ("GP") excesses.
#' @param W.model "M1", "M2", or "M3", as indicated in Papastathopoulos et al. (2023).
#' @param W.data Angles to include in likelihood: all angles ("AllW") or only those of exceedances of Q_q ("ExcOnly").
#' @param use.mean.Qq Boolean: use mean of Qq to fit G and L if TRUE, use samples from Qq otherwise.
#' @param q Quantile of the quantile set Q_q to fit.
#' @param alpha Quantile of the generalised Pareto to use in fitting parameterisation, alpha in (0,1).
#' @param N.Qq Number of quantile sets to sample from the posterior distribution of Q_q. Not used if use.mean.Qq is TRUE.
#' @param N.GW Number of G and W sets to sample from the posterior distribution of G and W.
#' @param QR.prior.range A length 2 vector, with (range0,Prange) specifying that P(rho < rho_0)=p_rho, where rho is the spatial range of the random field Q_q. If Prange is NA, then range0 is used as a fixed range value.
#' @param QR.prior.sigma A length 2 vector, with (sigma0,Psigma) specifying that P(sigma > sigma_0)=p_sigma, where sigma is the marginal standard deviation of the field Q_q. If Psigma is NA, then sigma0 is used as a fixed range value.
#' @param zeta.prior.range A length 2 vector, with (range0,Prange) specifying that P(rho < rho_0)=p_rho, where rho is the spatial range of the random field G. If Prange is NA, then range0 is used as a fixed range value.
#' @param zeta.prior.sigma A length 2 vector, with (sigma0,Psigma) specifying that P(σ > σ_0)=p_σ, where σ is the marginal standard deviation of the field B. If Psigma is NA, then sigma0 is used as a fixed range value.
#' @param phi.prior.range A length 2 vector, with (range0,Prange) specifying that P(rho < rho_0)=p_rho, where rho is the spatial range of the random field B. If Prange is NA, then range0 is used as a fixed range value.
#' @param phi.prior.sigma A length 2 vector, with (sigma0,Psigma) specifying that P(σ > σ_0)=p_σ, where σ is the marginal standard deviation of the field B. If Psigma is NA, then sigma0 is used as a fixed range value.
#' @param mesh.knots.2d Sequence of knots in a subset of of the interval -pi to pi.
#' @param mesh.res.3d Subdivision resolution for a semi-regular spherical triangulation with equidistant points along equidistant latitude bands.
#' @param seed Random number generator seed passed on to INLA::inla.posterior.sample.
#'
#' @return A list of options.
#' @export
#'
#' @examples
#' X <- rmutil::qlaplace(cbind(runif(1000),runif(1000),runif(1000)))
#' options <- set.options(X                = X,
#'                        excess.dist.fam  = "E",
#'                        W.data           = "ExcOnly",
#'                        W.model          = "M3",
#'                        use.mean.Qq      = "FALSE",
#'                        q                = 0.9,
#'                        alpha            = 0.5,
#'                        N.Qq             = 20,
#'                        N.GW             = 50,
#'                        QR.prior.range   = c(1, 0.999),
#'                        QR.prior.sigma   = c(1, 0.8),
#'                        zeta.prior.range = c(0.5, 0.999),
#'                        zeta.prior.sigma = c(5, 0.8),
#'                        phi.prior.range  = c(0.5, 0.999),
#'                        phi.prior.sigma  = c(5, 0.8),
#'                        mesh.knots.2d    = seq(-pi,pi,by=pi/400),
#'                        mesh.res.3d      = 20,
#'                        seed             = 0L)
set.options <- function(X,excess.dist.fam,W.model,W.data,use.mean.Qq,q,alpha,N.Qq,N.GW,
                        QR.prior.range = c(1, 0.999),QR.prior.sigma   = c(1, 0.8),
                        zeta.prior.range = c(0.5, 0.999),zeta.prior.sigma = c(5, 0.8),
                        phi.prior.range  = c(0.5, 0.999),phi.prior.sigma  = c(5, 0.8),
                        mesh.knots.2d    = seq(-pi,pi,by=pi/400),mesh.res.3d      = 20,
                        seed             = 0L){
  if(missing(X)){
    stop("Specify value of X.")
  }
  if(missing(excess.dist.fam)){
    stop("Specify value of excess.dist.fam: 'E' or 'GP'.")
  }else if(!(excess.dist.fam %in% c("E","GP"))){
    stop("Choose family for {R-r_Q(w) | R>r_Q(w), W=w} from: 'E', 'GP'")
  }
  if(excess.dist.fam=="GP"){
    stop("Functionality available very soon.")
  }

  if(missing(W.model)){
    stop("Specify value of W.data: 'M1', 'M2' or 'M3'.")
  }else if(!(W.model %in% c("M1","M2","M3"))){
    stop("Choose model for angles from: 'M1', 'M2', 'M3'")
  }

  if(missing(W.data)){
    stop("Specify value of W.data: 'ExcOnly' or 'AllW'.")
  }else if(!(W.data %in% c("AllW","ExcOnly"))){
    stop("Choose data for angles from: 'AllW', 'ExcOnly'")
  }
  if(missing(use.mean.Qq)){
    stop("Specify if use.mean.Qq is TRUE or FALSE.")
  }
  if(missing(q)){
    stop("Specify value of q in (0,1).")
  }
  if(missing(alpha)){
    stop("Specify value of alpha in (0,1).")
  }
  if(missing(N.Qq)){
    stop("Specify value of N.Qq.")
  }
  if(missing(N.GW)){
    stop("Specify value of N.Qq.")
  }
  if(missing(QR.prior.range)){
    stop("Specify value of QR.prior.range.")
  }
  if(missing(QR.prior.sigma)){
    stop("Specify value of QR.prior.sigma.")
  }
  if(missing(zeta.prior.range)){
    stop("Specify value of zeta.prior.range.")
  }
  if(missing(zeta.prior.sigma)){
    stop("Specify value of zeta.prior.sigma.")
  }
  if(missing(phi.prior.range)){
    stop("Specify value of phi.prior.range.")
  }
  if(missing(phi.prior.sigma)){
    stop("Specify value of phi.prior.sigma.")
  }
  if(ncol(X)==2){
    if(missing(mesh.knots.2d)){
      stop("Specify value of mesh.knots.2d.")
    }else if(max(mesh.knots.2d)>pi){
      stop("mesh.knots.2d must take values in [-pi,pi].")
    }else if(min(mesh.knots.2d)< -pi){
      stop("mesh.knots.2d must take values in [-pi,pi].")
    }
  }else if(ncol(X)==3){
    if(missing(mesh.res.3d)){
      stop("Specify value of mesh.res.3d.")
    }else if(mesh.res.3d<1){
      stop("mesh.res.3d must be greater than 1.")
    }
  }else{
    stop("X must be an n by d matrix with d=2 or d=3.")
  }

  return(list(excess.dist.fam  = excess.dist.fam,
              W.model          = W.model,
              W.data           = W.data,
              use.mean.Qq      = use.mean.Qq,
              q                = q,
              alpha            = alpha,
              N.Qq             = N.Qq,
              N.GW             = N.GW,
              QR.prior.range   = QR.prior.range,
              QR.prior.sigma   = QR.prior.sigma,
              zeta.prior.range = zeta.prior.range,
              zeta.prior.sigma = zeta.prior.sigma,
              phi.prior.range  = phi.prior.range,
              phi.prior.sigma  = phi.prior.sigma,
              mesh.knots.2d    = mesh.knots.2d,
              mesh.res.3d      = mesh.res.3d,
              seed             = seed,
              N.X              = nrow(X)))
}

#' Set configurations used by functions fit_Qq and fit_GL
#'
#' @param save.path Path to folder where .RData objects are to be saved.
#' @param file.nm Name of the .RData file.
#' @param save Boolean: TRUE to save .RData files, FALSE otherwise.
#' @param progress Boolean: TRUE to save progression in a .txt file, FALSE otherwise.
#'
#' @return A list of configurations.
#' @export
#'
#' @examples
#' config <- set.configs(save.path = "path/to/folder/",
#'                       file.nm   = "file_name", # default based on options if no argument given
#'                       save      = TRUE,
#'                       progress  = TRUE)
set.configs <- function(save.path = "",
                        file.nm   = paste0("DataName","_",options$excess.dist.fam,"_",options$W.model,"_",options$W.data),
                        save      = FALSE,
                        progress  = FALSE){
  if(save==TRUE & save.path==""){
    stop("Specify a path to save the <file.nm>.RData of fitted objects.")
  }
  if(!inherits(save,"logical")){
    stop("Argument save must be set to TRUE or FALSE.")
  }
  if(!inherits(progress,"logical")){
    stop("Argument progress must be set to TRUE or FALSE.")
  }
  config <- list(save.path = save.path,
                 file.nm   = file.nm,
                 save      = save,
                 progress  = progress)
}

#' Transform matrix of observed data to matrix with Laplace marginal distributions
#'
#' @param X n by d data matrix.
#' @param q 2 by d matrix. The top/bottom row denote the upper/lower quantiles above/below which to fit a GPd. NA for empirical.
#'
#'
#' @import evd
#' @import stats
#' @return A list including the transformed data and the marginal transformation parameters.
#' @export
#'
#' @examples n <- 500
#' X <- cbind(rnorm(n,0,1),rnorm(n,0,3),rnorm(n,0,5),rnorm(n,0,7))
#' q <- cbind(c(0.95,0.05),c(NA,0.05),c(0.95,NA),c(NA,NA))
#' LapTransf <- toLaplaceMargins(X,q)
toLaplaceMargins <- function(X,q){
  X <- as.matrix(X)
  d <- ncol(X)
  LapTransf <- list(q     = q,
                    mu    = matrix(NA,2,d),
                    sigma = matrix(NA,2,d),
                    xi    = matrix(NA,2,d),
                    X     = X,
                    X.L   = matrix(NA,ncol=d,nrow=nrow(X)))
  for(i in 1:d){
    q <- LapTransf$q[,i]
    LapTransf$mu[,i] <- quantile(X[,i],q)
    if(!is.na(q[1])){
      Y.up <- X[X[,i]>LapTransf$mu[1,i],i]-LapTransf$mu[1,i]

      gpd.fit <- fpot(Y.up, 0, model ="gpd")
      LapTransf$sigma[1,i] <- gpd.fit$estimate[1]
      LapTransf$xi[1,i]    <- gpd.fit$estimate[2]

    }
    if(!is.na(q[2])){
      Y.low <- -(X[X[,i]<LapTransf$mu[2,i],i]-LapTransf$mu[2,i])
      gpd.fit <- fpot(Y.low, 0, model ="gpd")
      LapTransf$sigma[2,i] <- gpd.fit$estimate[1]
      LapTransf$xi[2,i]    <- gpd.fit$estimate[2]

    }

    if(sum(is.na(q))==0){
      X.mid <- X[X[,i]>LapTransf$mu[2,i] & X[,i]<LapTransf$mu[1,i],i]
      F_X <- stats::ecdf(X.mid)
      n <- length(X.mid)

      pX <- function(x){
        if(x<LapTransf$mu[2,i]){
          return(q[2]*(1-evd::pgpd(-x,-LapTransf$mu[2,i],LapTransf$sigma[2,i],LapTransf$xi[2,i])))
        }else if(x>LapTransf$mu[2,i] & x <LapTransf$mu[1,i]){
          return(q[2] + (q[1]-q[2])*F_X(x)*n/(n+1))
        }else{
          return(q[1] + (1-q[1])*evd::pgpd(x,LapTransf$mu[1,i],LapTransf$sigma[1,i],LapTransf$xi[1,i]))
        }
      }

      U <- sapply(X[,i],function(xx) pX(xx))
      LapTransf$X.L[,i] <- rmutil::qlaplace(U)
    }else if(sum(is.na(q))==2){
      F_X <- stats::ecdf(X[,i])
      n <- length(X[,i])
      U <- F_X(X[,i])*n/(n+1)
      LapTransf$X.L[,i] <- rmutil::qlaplace(U)
    }else if(is.na(q)[1]){
      X.up <- X[X[,i]>LapTransf$mu[2,i],i]
      F_X <- stats::ecdf(X.up)
      n <- length(X.up)
      pX <- function(x){
        if(x<LapTransf$mu[2,i]){
          return(q[2]*(1-evd::pgpd(-x,-LapTransf$mu[2,i],LapTransf$sigma[2,i],LapTransf$xi[2,i])))
        }else if(x>LapTransf$mu[2,i]){
          return(q[2] + (1-q[2])*F_X(x)*n/(n+1))
        }
      }
      U <- sapply(X[,i],function(xx) pX(xx))
      LapTransf$X.L[,i] <- rmutil::qlaplace(U)
    }else{
      X.low <- X[X[,i]<LapTransf$mu[1,i],i]
      F_X <- stats::ecdf(X.low)
      n <- length(X.low)
      pX <- function(x){
        if(x<LapTransf$mu[1,i]){
          return(q[1]*F_X(x)*n/(n+1))
        }else{
          return(q[1] + (1-q[1])*evd::pgpd(x,LapTransf$mu[1,i],LapTransf$sigma[1,i],LapTransf$xi[1,i]))
        }
      }

      U <- sapply(X[,i],function(xx) pX(xx))
      LapTransf$X.L[,i] <- rmutil::qlaplace(U)
    }
  }
  LapTransf
}

#' Transform matrix of data with Laplace marginal distributions to matrix with original marginal distributions
#'
#' @param X.L n by d matrix of observations in Laplace margins, same dimensions as LapTransf$X.
#' @param LapTransf Object returned by the function toLaplaceMargins.
#'
#' @import rmutil
#' @import evd
#' @return An n by d matrix of observations in original marginal distributions.
#' @export
#'
#' @examples n <- 500
#' X <- rnorm(n,0,1)
#' q <- cbind(c(0.95,0.05))
#' LapTransf <- toLaplaceMargins(X,q)
#' X.L <- rmutil::qlaplace(runif(100))
#' X.Orig <- toOriginalMargins(X.L,LapTransf)
toOriginalMargins <- function(X.L,LapTransf){
  X.L <- as.matrix(X.L)
  d <- ncol(X.L)
  X.obs <- LapTransf$X
  if(d!=ncol(X.obs)){
    stop("X and X.L do not have same number of columns.")
  }
  X.orig <- matrix(NA,nrow(X.L),ncol(X.L))
  for(i in 1:d){
    q <- LapTransf$q[,i]
    if(sum(is.na(q))==0){
      X.mid.obs <- X.obs[X.obs[,i]>=LapTransf$mu[2,i] & X.obs[,i]<=LapTransf$mu[1,i],i]
      quantileX <- function(u){
        if(u<q[2]){
          u.low <- 1-u/q[2]
          return(-evd::qgpd(u.low,-LapTransf$mu[2,i],LapTransf$sigma[2,i],LapTransf$xi[2,i]))
        }else if(u>=q[2] & u<=q[1]){
          u.mid <- (u-q[2])/(q[1]-q[2])
          return(quantile(X.mid.obs,u.mid))
        }else{
          u.up <- (u-q[1])/(1-q[1])
          return(evd::qgpd(u.up,LapTransf$mu[1,i],LapTransf$sigma[1,i],LapTransf$xi[1,i]))
        }
      }
      n <- length(X.L[,i])
      X.orig[,i] <- sapply(rmutil::plaplace(X.L[,i]),function(uu) quantileX(uu)*(n+1)/n)
    }else if(sum(is.na(q))==2){
      U <- rmutil::plaplace(X.L[,i])
      n <- length(X.L[,i])
      X.orig[,i] <- sapply(U,function(uu) quantile(X.obs[,i],uu*(n+1)/n))

    }else if(is.na(q)[1]){
      X.up.obs <- X.obs[X.obs[,i]>=LapTransf$mu[2,i],i]
      n <- length(X.up.obs)
      quantileX <- function(u){
        if(u<q[2]){
          u.low <- 1-u/q[2]
          return(-evd::qgpd(u.low,-LapTransf$mu[2,i],LapTransf$sigma[2,i],LapTransf$xi[2,i]))
        }else if(u>=q[2]){
          u.mid <- (u-q[2])/(1-q[2])*(n+1)/n
          return(quantile(X.up.obs,u.mid))
        }
      }
      X.orig[,i] <- sapply(rmutil::plaplace(X.L[,i]),function(uu) quantileX(uu))
    }else{
      X.low.obs <- X.obs[X.obs[,i]<=LapTransf$mu[1,i],i]
      n <- length(X.low.obs)
      quantileX <- function(u){
        if(u<=q[1]){
          u.mid <- u/q[1]*(n+1)/n
          return(quantile(X.low.obs,u.mid))
        }else{
          u.up <- (u-q[1])/(1-q[1])
          return(evd::qgpd(u.up,LapTransf$mu[1,i],LapTransf$sigma[1,i],LapTransf$xi[1,i]))
        }
      }
      X.orig[,i] <- sapply(rmutil::plaplace(X.L[,i]),function(uu) quantileX(uu))
    }
  }
  X.orig
}

#' Obtain realisations from the posterior distribution of the quantile set Q_q
#'
#' @param X n by d data matrix.
#' @param options Object returned by the function set.options.
#' @param config Object returned by the function set.configs.
#' @param return_fitted_obj Boolean: TRUE to return fitted inlabru quantile regression object, FALSE otherwise.
#'
#' @import excursions
#' @return A list of posterior realisations of Q_q and other useful parameters.
#' @export
#'
#' @examples
fit_Qq <- function(X,options,config,return_fitted_obj=FALSE){
  if(ncol(X)==2){
    fitted.Qq <- fit_Qq_2d(X,options,config,return_fitted_obj=return_fitted_obj)
  }else if(ncol(X)==3){
    fitted.Qq <- fit_Qq_3d(X,options,config,return_fitted_obj=return_fitted_obj)
  }else{
    return("X must be an n by p matrix, with p in {2,3}")
  }
  fitted.Qq
}

#' Obtain realisations from the posterior distribution of the sets G and L
#'
#' @param fitted.Qq Object returned by the function fit_Qq.
#' @param config Object returned by the function set.configs.
#'
#' @return A list of posterior realisations of G and L and other useful parameters.
#' @export
#'
#' @examples
fit_GL <- function(fitted.Qq,config){
  if(ncol(fitted.Qq$X)==2){
    fitted.mod <- fit_GL_2d(fitted.Qq,config)
  }else if(ncol(fitted.Qq$X)==3){
    fitted.mod <- fit_GL_3d(fitted.Qq,config)
  }else{
    return("X must be an n by p matrix, with p in {2,3}")
  }
  fitted.mod
}

#' Sample angles from W and obtain g and G evaluated at these angles.
#'
#' @param fitted.mod  Object returned from the function fit_GL.
#' @param N.w         Number of angles to sample.
#' @param S_B         Subset of the 1-sphere S from -pi to pi. NA for S.
#' @param transf.G    Boolean, use transf.G if TRUE, G otherwise.
#'
#' @import INLA
#' @import inlabru
#' @return
#' @export
#'
#' @examples
sample_QGW_posterior <- function(fitted.mod,N.w,S_B=NA,transf.G=FALSE){
  if(ncol(fitted.mod$X)==2){
    post.samp <- sample_QGW_posterior_2d(fitted.mod,N.w,S_B=S_B,transf.G=transf.G)
  }else if(ncol(fitted.mod$X)>2){
    stop("Not yet implemented for d > 2.")
  }else{
    return("X must be an n by p matrix, with p = 2.")
  }
  post.samp
}

#' Probability estimation
#'
#' @param fitted.mod Object returned from the function fit_GL.
#' @param post.sample Object returned from the function sample_QGW_posterior.
#' @param x_B x-component of the rectangular set B of interest.
#' @param y_B y-component of the rectangular set B of interest.
#'
#' @import INLA
#' @import inlabru
#' @import evd
#' @return Vector of posterior samples for the probability that a new observation falls in B.
#' @export
#'
#' @examples
prob_estimation <- function(fitted.mod,post.sample,x_B,y_B){
  if(x_B[1]<0 & x_B[2] >0 & y_B[1]<0 & y_B[2] >0){
    stop("The origin 0 cannot be in the set B.")
  }
  # if(fitted.mod$options$excess.dist.fam=="GP"){
  #   stop("GP exceedances probabilities not yet implemented")
  # }
  if(ncol(fitted.mod$X)==2){
    p.post <- prob_estimation_2d(fitted.mod,post.sample,x_B,y_B)
  }else if(ncol(fitted.mod$X)>2){
    stop("Not yet implemented for d > 2.")
  }else{
    return("X must be an n by p matrix, with p = 2.")
  }
  p.post
}


#' Obtain realisations from the radial function of the return level-set X_t.
#'
#' @param fitted.mod Object returned from the function fit_GL.
#' @param t Return period t <= 1/(1-q) for the return level-set X_t.
#' @param q.prime Probability of the canonical return level-set level-set X_{1/(1-q)}.
#'
#' @return Matrix of realisations from the radial function of the return level-set X_t.
#' @export
#'
#' @examples
return_set <- function(fitted.mod,t=NA,q.prime=NA){
  if(inherits(t,"logical") & inherits(q.prime,"logical")){
    stop("Specify a value for t or q.prime.")
  }else if(!inherits(t,"logical") & !inherits(q.prime,"logical")){
    stop("Specify a value for t or for q.prime, not both.")
  }else if(!inherits(q.prime,"logical")){
    if(min(q.prime)<=fitted.mod$options$q){
      stop("q.prime must be greater than q.")
    }
    t <- sort(1/(1-q.prime))
  }else{
    if(min(1-1/t)<=fitted.mod$options$q){
      stop("t must be less than 1/(1-q).")
    }
  }

  if(ncol(fitted.mod$X)==2){
    ret_set <- return_set_2d(fitted.mod,t=t)
  }else if(ncol(fitted.mod$X)>2){
    if(length(t)>1){
      stop("Specifiy only one value of t or q.prime for 3d plots.")
    }
    ret_set <- return_set_3d(fitted.mod,t=t)
  }else{
    return("X must be an n by p matrix, with p = 2.")
  }
}

#' Obtain posterior samples from the conditional extremes alpha parameter
#'
#' @param fitted.mod Object returned from the function fit_GL.
#'
#' @return Matrix of posterior samples from alpha_{2|1} and alpha_{2|1}
#' @export
#'
#' @examples
alpha_posterior <- function(fitted.mod){
  if(ncol(fitted.mod$X)==2){
    alphas <- alpha_posterior_2d(fitted.mod)
  }else if(ncol(fitted.mod$X)==3){
    stop("Not yet implemented for d>2.")
  }
  alphas
}

#' Obtain posterior samples from the coefficient of residual tail dependence eta
#'
#' @param fitted.mod Object returned from the function fit_GL.
#'
#' @return Vector of posterior samples from eta.
#' @export
#'
#' @examples
eta_posterior <- function(fitted.mod){
  if(ncol(fitted.mod$X)==2){
    etas <- eta_posterior_2d(fitted.mod)
  }else if(ncol(fitted.mod$X)==3){
    stop("Not yet implemented for d>2.")
  }
  etas
}

#' Plot the posterior mean and (1-alpha)-predictive intervals for Q_q
#'
#' @param fitted.Qq Object returned by the function fit_Qq.
#' @param alpha Value in (0,1) for the (1-alpha)-simultaneous predictive interval of Q_q.
#' @param surface3d Surface to plot between the mean ("mean"), and simultaneous predictive interval bounds ("lower_alpha" or "upper_alpha").
#' @param cex.pts Size of points representing exceedances of Q_q.
#' @param cex.axis Size of axes' labels.
#' @param xlim Plotting range on x-axis.
#' @param ylim Plotting range on y-axis.
#' @param by Length of intervals between axis ticks and displayed values.
#'
#' @return None.
#' @export
#'
#' @examples
plot_Qq <- function(fitted.Qq,alpha=0.05,surface3d="mean",cex.pts=0.4,cex.axis=1.4,xlim=c(0,0),ylim=c(0,0),by=2){
  if(ncol(fitted.Qq$X)==2){
    plot_Qq_2d(fitted.Qq,alpha,cex.pts=0.4,cex.axis=1.4,xlim=c(0,0),ylim=c(0,0),by=2)
  }else if(ncol(fitted.Qq$X)==3){
    plot_Qq_3d(fitted.Qq,surface3d,alpha,xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3]))
  }
}

#' Plot the posterior mean and (1-alpha)-predictive intervals for G
#'
#' @param fitted.mod Object returned by the function fit_GL.
#' @param alpha Value in (0,1) for the (1-alpha)-simultaneous predictive interval of G.
#' @param surface3d Surface (3d) to plot between the mean ("mean"), and simultaneous predictive interval bounds ("lower_alpha" or "upper_alpha").
#' @param surf.col Color of surface (3d).
#' @param cex.pts Size of points representing exceedances of Q_q.
#' @param cex.axis Size of axes' labels.
#' @param xlim Plotting range on x-axis.
#' @param ylim Plotting range on y-axis.
#' @param by Length of intervals between axis ticks and displayed values.
#'
#' @return None.
#' @export
#'
#' @examples
plot_G <- function(fitted.mod,alpha=0.05,surface3d="mean",surf.col="grey",cex.pts=0.4,cex.axis=1.4,xlim=c(0,0),ylim=c(0,0),by=2){
  if(ncol(fitted.mod$X)==2){
    plot_G_2d(fitted.mod,alpha,mean_med=TRUE,cex.txt=1.4,cex.pts=0.6,transf.G = FALSE,cex.axis=1.4,by=2)
  }else if(ncol(fitted.mod$X)==3){
    plot_G_3d(fitted.mod,alpha,surface=surface3d,xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3]),surf.col=surf.col)
  }
}

#' Plot the posterior mean and (1-alpha)-predictive intervals for W
#'
#' @param fitted.mod Object returned by the function fit_GL.
#' @param alpha Value in (0,1) for the (1-alpha)-simultaneous predictive interval of W.
#' @param xylim Plotting range on the x and y axes.
#' @param main Title of figure.
#' @param mid.gap Radius of the middle circle of zero density.
#' @param txt.gap Gap between angle marks and axes' lines.
#' @param cex.txt Size of angle marks.
#' @param exconly Boolean: TRUE if circular histogram should be based only on exceedances, FALSE for all angles.
#'
#' @return None.
#' @export
#'
#' @examples
plot_W <- function(fitted.mod,alpha=0.05,xylim,main="",mid.gap=0.1, txt.gap=0.02,cex.txt=1.4,exconly=FALSE){
  if(ncol(fitted.mod$X)==2){
    plot_W_2d(fitted.mod,alpha,xylim,main="",mid.gap=0.1, txt.gap=0.02,cex.txt=1.4,exconly=exconly)
  }else if(ncol(fitted.mod$X)==3){
    return("Currently not implemented for 3d.")
  }
}

#' Plot the boundary of the return set X_t and its simultaneous predictive intervals
#'
#' @param fitted.mod Object returned by the function fit_GL.
#' @param list_ret_sets Object returned by the function return_set.
#' @param plt Plot return level-set ("set") or boundary ("boundary").
#' @param xylim Plotting limits of x and y axes.
#' @param xyzlim Plotting limits of x, y, and y axes.
#' @param xlab Label of x-axis.
#' @param ylab Label of y-axis.
#' @param zlab Label of z-axis.
#' @param by Length of intervals between axis ticks and displayed values.
#' @param cex.pts Size of observed data points.
#' @param cex.axis Size of text on axes.
#'
#' @return None
#' @export
#'
#' @examples
plot_X_t <- function(fitted.mod,list_ret_sets,plt = "boundary",xylim=c(0,0),xyzlim=c(0,0),xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3]),by=2,cex.pts=0.4,cex.axis=1.4){
  if(!(plt %in% c("boundary","set"))){
    stop("list_ret_sets must take value 'boundary' or 'set'.")
  }
  if(ncol(fitted.mod$X)==2){
    if(plt=="boundary"){
      plot_return_bdry_2d(fitted.mod,list_ret_sets,cex.pts=cex.pts,cex.axis=cex.axis,xlim=xylim,ylim=xylim,by=by)
    }else if(plt=="set"){
      plot_return_sets_2d(fitted.mod,list_ret_sets,cex.pts=cex.pts,cex.axis=cex.axis,xlim=xylim,ylim=xylim,by=by)
    }
  }else if(ncol(fitted.mod$X)==3){
    if(plt=="boundary"){
      plot_return_bdry_3d(fitted.mod,list_ret_sets,cex.pts=cex.pts,cex.axis=cex.axis,xyzlim=xyzlim,xlab=xlab,ylab=ylab,zlab=zlab)
    }else if(plt=="set"){
      stop("Not available for 3d.")
    }
  }
}
