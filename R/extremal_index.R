#' Exceedances of the quantile set Q_q
#'
#' @param fitted.Qq Object returned by the function fit_Qq (or fit_GL).
#' @param use.mean.Qq Boolean: use the posterior mean of Q_q if TRUE, the posterior samples of Q_q otherwise.
#'
#' @return A list with, for each posterior sample of Q_q (or for its posterior mean), the indices of the observations exceeding Q_q.
#' @export
get_exceedances <- function(fitted.Qq,use.mean.Qq=FALSE){
  X <- fitted.Qq$X
  R <- sqrt(rowSums(X^2))
  if(ncol(X)==2){
    W <- atan2(y=X[,2], x=X[,1])
  }else if(ncol(X)==3){
    W <- X/R
  }else{
    stop("X must be a n x d matrix with d=2 or d=3.")
  }
  A.Qq <- fmesher::fm_basis(fitted.Qq$mesh,loc=W)
  if(use.mean.Qq==TRUE){
    thres.Qq <- list(as.vector(A.Qq %*% fitted.Qq$Qq.mean$mean))
  }else{
    thres.Qq <- lapply(seq_len(ncol(fitted.Qq$Qq)),function(i) as.vector(A.Qq %*% fitted.Qq$Qq[,i]))
  }
  lapply(thres.Qq,function(thres) which(R>thres))
}

runs_clusters <- function(exc.ind,m){
  if(length(exc.ind)==0){
    return(list())
  }
  unname(split(exc.ind,cumsum(c(1,diff(exc.ind)>m))))
}

#' Clusters of exceedances of Q_q via the runs method
#'
#' @param fitted.Qq Object returned by the function fit_Qq (or fit_GL).
#' @param m Runs length: exceedances separated by more than m time steps belong to different clusters.
#' @param use.mean.Qq Boolean: use the posterior mean of Q_q if TRUE, the posterior samples of Q_q otherwise.
#'
#' @return A list with, for each posterior sample of Q_q (or for its posterior mean), a list of clusters, each given by the indices of its exceedances.
#' @export
get_clusters_runs <- function(fitted.Qq,m,use.mean.Qq=FALSE){
  exc.ind <- get_exceedances(fitted.Qq,use.mean.Qq=use.mean.Qq)
  lapply(exc.ind,runs_clusters,m=m)
}

#' Radial maxima of the clusters of exceedances of Q_q
#'
#' @param fitted.Qq Object returned by the function fit_Qq (or fit_GL).
#' @param m Runs length: exceedances separated by more than m time steps belong to different clusters.
#' @param use.mean.Qq Boolean: use the posterior mean of Q_q if TRUE, the posterior samples of Q_q otherwise.
#'
#' @return A list with, for each posterior sample of Q_q (or for its posterior mean), the indices of the exceedance with the largest radius in each cluster.
#' @export
get_radial_maxima <- function(fitted.Qq,m,use.mean.Qq=FALSE){
  R <- sqrt(rowSums(fitted.Qq$X^2))
  clusters <- get_clusters_runs(fitted.Qq,m,use.mean.Qq=use.mean.Qq)
  lapply(clusters,function(cl) vapply(cl,function(inds) inds[which.max(R[inds])],integer(1)))
}

#' Direction-tagged binary variables for the directional extremal index
#'
#' @param fitted.Qq Object returned by the function fit_Qq (or fit_GL).
#' @param m Runs length.
#' @param use.mean.Qq Boolean: use the posterior mean of Q_q if TRUE, the posterior samples of Q_q otherwise.
#'
#' @return A list with, for each posterior sample of Q_q (or for its posterior mean), a matrix whose first column is the binary variable z (1 if none of the m observations following an exceedance exceed Q_q) and whose other columns are the direction of the exceedance (angle for d=2, Cartesian coordinates on the sphere for d=3).
#' @export
get_filtered_Z <- function(fitted.Qq,m,use.mean.Qq=FALSE){
  X <- fitted.Qq$X
  if(ncol(X)==2){
    filtered_Z <- get_filtered_Z_2d(fitted.Qq,m,use.mean.Qq=use.mean.Qq)
  }else if(ncol(X)==3){
    filtered_Z <- get_filtered_Z_3d(fitted.Qq,m,use.mean.Qq=use.mean.Qq)
  }else{
    stop("X must be a n x d matrix with d=2 or d=3.")
  }
  filtered_Z
}

filtered_Z_indicators <- function(exc.ind,m){
  n.Z <- max(length(exc.ind)-m,0)
  Z <- rep(0,n.Z)
  for(j in seq_len(n.Z)){
    k <- exc.ind[j]
    Z[j] <- ifelse(sum(c((k+1):(k+m)) %in% exc.ind[(j+1):(j+m)])==0,1,0)
  }
  Z
}

get_filtered_Z_2d <- function(fitted.Qq,m,use.mean.Qq=FALSE){
  exc.ind <- get_exceedances(fitted.Qq,use.mean.Qq=use.mean.Qq)
  W <- atan2(y=fitted.Qq$X[,2], x=fitted.Qq$X[,1])
  lapply(exc.ind,function(inds){
    Z <- filtered_Z_indicators(inds,m)
    cbind(Z,W[inds[seq_along(Z)]])
  })
}

get_filtered_Z_3d <- function(fitted.Qq,m,use.mean.Qq=FALSE){
  exc.ind <- get_exceedances(fitted.Qq,use.mean.Qq=use.mean.Qq)
  R <- sqrt(rowSums(fitted.Qq$X^2))
  W <- fitted.Qq$X/R
  lapply(exc.ind,function(inds){
    Z <- filtered_Z_indicators(inds,m)
    cbind(Z,W[inds[seq_along(Z)],,drop=FALSE])
  })
}

#' Fit the directional extremal index theta
#'
#' @param filtered_Z Object returned by the function get_filtered_Z.
#' @param options List with elements mesh.knots.2d (d=2) or mesh.res.3d (d=3), zeta.prior.range, zeta.prior.sigma, N.theta and seed.
#' @param return_fitted_obj Boolean: TRUE to return the fitted inlabru objects, FALSE otherwise.
#'
#' @return A list with the mesh and, for each element of filtered_Z, the posterior samples of theta at the mesh nodes (a matrix with N.theta columns), or the posterior summary of theta if N.theta is 0.
#' @export
fit_theta <- function(filtered_Z,options,return_fitted_obj=FALSE){
  if(ncol(filtered_Z[[1]])==2){
    theta <- fit_theta_2d(filtered_Z,options,return_fitted_obj)
  }else if(ncol(filtered_Z[[1]])==4){
    theta <- fit_theta_3d(filtered_Z,options,return_fitted_obj)
  }else{
    stop("The argument 'filtered_Z' is not a valid object returned by the 'get_filtered_Z' function.")
  }
  theta
}

fit_theta_2d <- function(filtered_Z,options,return_fitted_obj=FALSE){
  mesh  <- fmesher::fm_mesh_1d(sort(options$mesh.knots.2d),boundary="cyclic",degree=1)

  circular.Matern_zeta  <- inla.spde2.pcmatern(mesh, prior.range = options$zeta.prior.range,
                                               prior.sigma = options$zeta.prior.sigma)

  cmp1 <- Z ~ Intercept(1) + f(W, model=circular.Matern_zeta)

  theta.hat <- list(mesh=mesh,theta=list())
  for(i in seq_along(filtered_Z)){
    message(paste0("Fitting theta(w) for Qq ",i,"/",length(filtered_Z),"."))
    Z <- filtered_Z[[i]][,1]
    W <- filtered_Z[[i]][,2]

    fit1 <- bru(cmp1, family = "binomial",Ntrials=1 ,data = data.frame(Z=Z,W=W),
                options = list(num.threads=1))
    if(return_fitted_obj) theta.hat$fit.obj[[i]] <- fit1

    if(options$N.theta==0){
      theta.hat$theta[[i]] <- predict(fit1,data.frame(W=mesh$loc),
                                      ~1-1/(1+exp(Intercept+f)))
    }else{
      theta.hat$theta[[i]] <- generate(n.samples = options$N.theta,
                                       x         = fit1,
                                       newdata   = data.frame(W=mesh$loc),
                                       formula   = ~1-1/(1+exp(Intercept+f)),
                                       seed      = options$seed)
    }
  }
  theta.hat
}

fit_theta_3d <- function(filtered_Z,options,return_fitted_obj=FALSE){
  mesh  <- fmesher::fm_rcdt_2d_inla(globe=options$mesh.res.3d, crs=fmesher::fm_CRS("sphere"))

  spherical.Matern  <- inla.spde2.pcmatern(mesh,
                                           prior.range = options$zeta.prior.range,
                                           prior.sigma = options$zeta.prior.sigma)

  cmp1 <- Z ~ Intercept(1) + f(cbind(x,y,z),
                               model=spherical.Matern,
                               mapper=bru_mapper(mesh))

  newdata <- data.frame(x=mesh$loc[,1],y=mesh$loc[,2],z=mesh$loc[,3])
  theta.hat <- list(mesh=mesh,theta=list())
  for(i in seq_along(filtered_Z)){
    message(paste0("Fitting theta(w) for Qq ",i,"/",length(filtered_Z),"."))
    Z <- filtered_Z[[i]][,1]
    W <- filtered_Z[[i]][,2:4,drop=FALSE]

    fit1 <- bru(cmp1, family = "binomial",Ntrials=1,
                data = data.frame(Z=Z,x=W[,1],y=W[,2],z=W[,3]),
                options = list(num.threads=1))
    if(return_fitted_obj) theta.hat$fit.obj[[i]] <- fit1

    if(options$N.theta==0){
      theta.hat$theta[[i]] <- predict(fit1,newdata,~1-1/(1+exp(Intercept+f)))
    }else{
      theta.hat$theta[[i]] <- generate(n.samples = options$N.theta,
                                       x         = fit1,
                                       newdata   = newdata,
                                       formula   = ~1-1/(1+exp(Intercept+f)),
                                       seed      = options$seed)
    }
  }
  theta.hat
}

#' Plot the posterior mean and (1-alpha)-credible intervals of theta
#'
#' @param filtered_Z Object returned by the function get_filtered_Z.
#' @param theta_hat Object returned by the function fit_theta, with N.theta > 0.
#' @param plot_filt_Z Boolean: plot the binary variables z of the first element of filtered_Z if TRUE.
#' @param sim_conf Boolean: simultaneous credible intervals if TRUE, marginal otherwise.
#' @param alpha Value in (0,1) for the (1-alpha)-credible intervals.
#'
#' @return None.
#' @export
plot_theta <- function(filtered_Z,theta_hat,plot_filt_Z=TRUE,sim_conf=TRUE,alpha=0.05){
  if(ncol(filtered_Z[[1]])==2){
    plot_theta_2d(filtered_Z,theta_hat,plot_filt_Z=plot_filt_Z,
                  sim_conf=sim_conf,alpha=alpha)
  }else if(ncol(filtered_Z[[1]])==4){
    stop("Currently not implemented for 3d.")
  }else{
    stop("The argument 'filtered_Z' is not a valid object returned by the 'get_filtered_Z' function.")
  }
}

#' Plot the posterior mean and (1-alpha)-credible intervals of the mean cluster size 1/theta
#'
#' @param filtered_Z Object returned by the function get_filtered_Z.
#' @param theta_hat Object returned by the function fit_theta, with N.theta > 0.
#' @param sim_conf Boolean: simultaneous credible intervals if TRUE, marginal otherwise.
#' @param alpha Value in (0,1) for the (1-alpha)-credible intervals.
#'
#' @return None.
#' @export
plot_mean_cluster_size <- function(filtered_Z,theta_hat,sim_conf=TRUE,alpha=0.05){
  if(ncol(filtered_Z[[1]])==2){
    plot_mean_cluster_size_2d(filtered_Z,theta_hat,
                              sim_conf=sim_conf,alpha=alpha)
  }else if(ncol(filtered_Z[[1]])==4){
    stop("Currently not implemented for 3d.")
  }else{
    stop("The argument 'filtered_Z' is not a valid object returned by the 'get_filtered_Z' function.")
  }
}

theta_samples <- function(theta_hat){
  if(!all(vapply(theta_hat$theta,is.matrix,logical(1))))
    stop("theta_hat must contain posterior samples (N.theta > 0).")
  do.call(cbind,theta_hat$theta)
}

plot_theta_2d <- function(filtered_Z,theta_hat,plot_filt_Z=TRUE,sim_conf=TRUE,alpha=0.05){
  W <- filtered_Z[[1]][,2]
  Z <- filtered_Z[[1]][,1]
  R <- Z + 0.25
  Zs <- cbind(R*cos(W),R*sin(W))

  theta.samp <- theta_samples(theta_hat)

  par(mfcol=c(1,1),mar=c(4,4,3,1),mgp=c(3,1,0),pty="s")
  theta.mean <- apply(theta.samp,1,mean)
  if(sim_conf){
    excurs <- simconf.mc(samples = theta.samp,alpha = alpha)
    theta.low <- excurs$a
    theta.high <- excurs$b
  }else{
    theta.low <- apply(theta.samp,1,function(x) quantile(x,prob=alpha/2))
    theta.high <- apply(theta.samp,1,function(x) quantile(x,prob=1-alpha/2))
  }

  w.seq <- seq(0,2*pi,by=0.01)
  Z.0 <- cbind(0.25*cos(w.seq),0.25*sin(w.seq))
  Z.1 <- cbind(1.25*cos(w.seq),1.25*sin(w.seq))
  plot(Z.1[,1],Z.1[,2],col="red",
       main=expression(theta),xlab="",ylab="",
       xaxt="n",yaxt="n",bty="n",type="l",lty="dashed")
  lines(Z.0[,1],Z.0[,2],col="red",lty="dashed")

  if(plot_filt_Z==TRUE){
    points(Zs,pch=16)
  }
  W.xy <- cbind(cos(theta_hat$mesh$loc),sin(theta_hat$mesh$loc))
  lines((theta.mean+0.25)*W.xy[,1],(theta.mean+0.25)*W.xy[,2])
  lines((theta.low+0.25)*W.xy[,1],(theta.low+0.25)*W.xy[,2],lty="dotted")
  lines((theta.high+0.25)*W.xy[,1],(theta.high+0.25)*W.xy[,2],lty="dotted")
}

plot_mean_cluster_size_2d <- function(filtered_Z,theta_hat,sim_conf=TRUE,alpha=0.05){
  par(mfcol=c(1,1),mar=c(4,4,3,1),mgp=c(3,1,0),pty="s")
  mean_clust_size <- 1/theta_samples(theta_hat)
  mean_clust_size.mean <- apply(mean_clust_size,1,mean)
  if(sim_conf){
    excurs <- simconf.mc(samples = mean_clust_size,alpha = alpha)
    mean_clust_size.low <- excurs$a
    mean_clust_size.high <- excurs$b
  }else{
    mean_clust_size.low <- apply(mean_clust_size,1,function(x) quantile(x,prob=alpha/2))
    mean_clust_size.high <- apply(mean_clust_size,1,function(x) quantile(x,prob=1-alpha/2))
  }

  mesh.loc <- theta_hat$mesh$loc
  mean_clust_size.high.xy <- cbind(mean_clust_size.high*cos(mesh.loc),
                                   mean_clust_size.high*sin(mesh.loc))
  mean_clust_size.low.xy <- cbind(mean_clust_size.low*cos(mesh.loc),
                                  mean_clust_size.low*sin(mesh.loc))
  mean_clust_size.mean.xy <- cbind(mean_clust_size.mean*cos(mesh.loc),
                                   mean_clust_size.mean*sin(mesh.loc))
  plot(mean_clust_size.high.xy[,1],mean_clust_size.high.xy[,2],
       xlab=expression(X[1]),ylab=expression(X[2]),
       xaxt="n",yaxt="n",bty="n",type="l",lty="dotted")
  axis(1)
  axis(2)
  lines(mean_clust_size.low.xy[,1],mean_clust_size.low.xy[,2],lty="dotted")
  lines(mean_clust_size.mean.xy[,1],mean_clust_size.mean.xy[,2])
}
