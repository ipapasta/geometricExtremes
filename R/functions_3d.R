#' Title
#'
#' @param X
#' @param options
#' @param config
#' @param return_fitted_obj
#'
#' @import INLA
#' @import inlabru
#' @import pracma
#' @return
#' @noRd
#'
#' @examples
fit_Qq_3d <- function(X,options,config,return_fitted_obj=F){

  if(config$progress){
    prog <- paste0("Begin quantile regression")
    write.table(prog, file = paste0(config$save.path,"Progression.txt"),append=TRUE)
  }

  # Convert from cartesian to spherical coordinates
  polar.coord    <- cart2sph(as.matrix(X))
  theta <- polar.coord[,1]
  phi   <- polar.coord[,2]
  R     <- polar.coord[,3]
  cart.coord.on.sphere <- sph2cart(cbind(theta, phi, rep(1, nrow(X))))
  dfRphitheta <- data.frame(R=R, theta=theta, phi=phi,
                            x=cart.coord.on.sphere[,1],
                            y=cart.coord.on.sphere[,2],
                            z=cart.coord.on.sphere[,3],
                            X1=X[,1], X2=X[,2], X3=X[,3])

  ## Definition of mesh and prior
  mesh.globe     <- inla.mesh.create(globe=options$mesh.res.3d, crs=fm_CRS("sphere"))
  spherical.spde  <- inla.spde2.pcmatern(mesh.globe,
                                         prior.range = options$QR.prior.range,
                                         prior.sigma = options$QR.prior.sigma)

  cmp.spherical   <- R ~ Intercept(1) + spde3(cbind(x, y, z), model=spherical.spde, mapper=bru_mapper(mesh.globe, indexed=TRUE))

  fit.Qq <- bru(cmp.spherical, data = as.data.frame(dfRphitheta),
                family="gamma",
                options=list(num.threads=8, verbose = F), #safe=TRUE,
                control.family = list(control.link = list(model="quantile",
                                                          quantile=options$q)))

  if(config$progress){
    prog <- paste0("QR done")
    write.table(prog, file = paste0(config$save.path,"Progression.txt"),append=TRUE)
  }

  mesh.globe.df  <- data.frame(x=mesh.globe$loc[,1], y=mesh.globe$loc[,2], z=mesh.globe$loc[,3]) #%>% arrange(x, y, z)

  Qq.mean  <- predict(fit.Qq,newdata=mesh.globe.df, ~ exp(Intercept+spde3),seed = options$seed)
  Qq.post.samp  <- generate(n.samples = options$N.Qq,
                            object    = fit.Qq,
                            newdata   = mesh.globe.df,#data.frame(x=X[,1], y=X[,2]),
                            formula   = ~ exp(Intercept+spde3),
                            seed      = options$seed)

  prop.exc <- rep(NA,options$N.Qq)

  data.unit.sphere <- sph2cart(cbind(as.matrix(dfRphitheta[,2:3]),1))
  A.Qq <- inla.mesh.projector(mesh=mesh.globe,loc=data.unit.sphere)$proj$A

  for(i in 1:options$N.Qq){
    Qq.thetaphi <- as.vector(A.Qq %*% Qq.post.samp[,i])
    is.excess  <- ifelse(dfRphitheta$R - Qq.thetaphi > 0, 1, 0)
    prop.exc[i]       <- sum(is.excess)/options$N.X
  }

  if(return_fitted_obj==F){
    fitted.Qq <- list(X         = X,
                      mesh      = mesh.globe,
                      options    = options,
                      Qq      = Qq.post.samp,
                      Qq.mean  = Qq.mean,
                      prop.exc  = prop.exc)
  }else{
    fitted.Qq <- list(X          = X,
                      mesh       = mesh.globe,
                      options     = options,
                      Qq       = Qq.post.samp,
                      Qq.mean   = Qq.mean,
                      prop.exc   = prop.exc,
                      QR.fit.obj = fit.Qq)
  }
  if(config$save){
    filenm <- paste0(config$file.nm,"_Qq.RData")
    save(fitted.Qq,file = paste0(config$save.path,filenm))
  }
  return(fitted.Qq)
}

#' Title
#'
#' @param fitted.Qq
#' @param config
#'
#' @import INLA
#' @import inlabru
#' @import pracma
#' @import sp
#' @return
#' @noRd
#'
#' @examples
fit_GL_3d <- function(fitted.Qq,config){
  X <- fitted.Qq$X
  options <- fitted.Qq$options

  polar.coord    <- cart2sph(as.matrix(X))
  theta <- polar.coord[,1]
  phi   <- polar.coord[,2]
  R     <- polar.coord[,3]
  cart.coord.on.sphere <- sph2cart(cbind(theta, phi, rep(1, nrow(X))))
  dfRphitheta <- data.frame(R=R, theta=theta, phi=phi,
                            x=cart.coord.on.sphere[,1],
                            y=cart.coord.on.sphere[,2],
                            z=cart.coord.on.sphere[,3],
                            X1=X[,1], X2=X[,2], X3=X[,3])

  ## Definition of mesh and prior
  mesh.globe     <- fitted.Qq$mesh
  mesh.globe.df  <- data.frame(x=mesh.globe$loc[,1], y=mesh.globe$loc[,2], z=mesh.globe$loc[,3])

  fitted.mod <- list(X         = X,
                     mesh      = mesh.globe,
                     options    = options,
                     Qq      = fitted.Qq$Qq,
                     Qq.mean  = fitted.Qq$Qq.mean,
                     prop.exc  = fitted.Qq$prop.exc,
                     log.L = list(),
                     g         = list(),
                     G         = list())

  spherical.spde.zeta  <- inla.spde2.pcmatern(mesh.globe,
                                              prior.range = options$zeta.prior.range,
                                              prior.sigma = options$zeta.prior.sigma)
  spherical.spde.phi   <- inla.spde2.pcmatern(mesh.globe,
                                              prior.range = options$phi.prior.range,
                                              prior.sigma = options$phi.prior.sigma)

  if(options$W.model=="M1"){
    cmp <- ~ 0 + Int.W(1) + Int.G(1) + zeta(coordinates, model = spherical.spde.zeta)
  }else if(options$W.model=="M2" | options$W.model=="M3"){
    cmp <- ~ 0 + Int.W(1) + Int.G(1) + zeta(coordinates, model = spherical.spde.zeta) + phi(coordinates, model = spherical.spde.phi)
  }

  data.unit.sphere <- sph2cart(cbind(as.matrix(dfRphitheta[,2:3]),1))
  A.Qq <- inla.mesh.projector(mesh=mesh.globe,loc=data.unit.sphere)$proj$A

  if(options$use.mean.Qq==TRUE){
    N.Qq <- 1
  }else{
    N.Qq <- options$N.Qq
  }
  for(i in 1:N.Qq){
    message(paste0("Estimating G and W for Qq ",i,"/",options$N.Qq))
    if(config$progress){
      prog <- paste0("Estimating G and W for Qq ",i,"/",options$N.Qq)
      write.table(prog, file = paste0(config$save.path,"Progression.txt"),append=TRUE)
    }

    if(options$use.mean.Qq==TRUE){
      Qq.thetaphi <- as.vector(A.Qq %*% fitted.Qq$Qq.mean$mean)
    }else{
      Qq.thetaphi <- as.vector(A.Qq %*% fitted.Qq$Qq[,i])
    }

    dfRphitheta$exc  <- dfRphitheta$R - Qq.thetaphi
    dfRphitheta.exc    <- dfRphitheta[dfRphitheta$exc > 0,] #dfRphitheta %>% dplyr::filter(exc > 0)

    ## To fit a point pattern on the sphere with inlabru, we need to
    ## create a SpatialPointsDataFrame and store the Cartesian coordinates
    ## (x,y,z) as "coordinates".

    if(options$W.data=="AllW"){
      df.data.W <- SpatialPointsDataFrame(coords = dfRphitheta[,c("x", "y", "z")],
                                          data=data.frame(exc=dfRphitheta$exc))
      k <- nrow(df.data.W)
    }else if(options$W.data=="ExcOnly"){
      df.data.W <- SpatialPointsDataFrame(coords = dfRphitheta.exc[,c("x", "y", "z")],
                                          data=data.frame(exc=dfRphitheta.exc$exc))
      k <- nrow(df.data.W)
    }
    df.data.exc <- SpatialPointsDataFrame(coords = dfRphitheta.exc[,c("x", "y", "z")],
                                          data=data.frame(exc=dfRphitheta.exc$exc))

    alpha <- options$alpha
    if(options$excess.dist.fam=="E"){
      if(options$W.model=="M1"){
        form.W  <- coordinates  ~ Int.W + log(k) - 3*Int.G - 3* zeta
      }else if(options$W.model=="M2"){
        form.W  <- coordinates  ~ Int.W + log(k) - phi
      }else if(options$W.model=="M3"){
        form.W  <- coordinates  ~ Int.W + log(k) - 3*Int.G - 3* zeta - phi
      }
      form.G   <- exc ~ Int.G + zeta
      lik_G <- like("exponential", formula=form.G, data=df.data.exc)
    }else if(options$excess.dist.fam=="GP"){
      if(options$W.model=="M1"){
        form.W  <- coordinates ~ 0 + Int.W + 3*Int.G + 3*zeta + log(k)
      }else if(options$W.model=="M2"){
        form.W  <- coordinates ~ 0 + Int.W - phi + log(k)
      }else if(options$W.model=="M3"){
        form.W  <- coordinates ~ 0 + Int.W + 3*Int.G + 3*zeta - phi + log(k)
      }
      form.G  <- exc ~ 0 + Int.G + zeta
      lik_G <- like("gp", formula=form.G, data=df.data.exc,
                    control.family = list(control.link = list(quantile = alpha)))
    }

    lik_W  <- like("cp", formula = form.W,
                   data = df.data.W,
                   domain = list(coordinates=mesh.globe),
                   ips = ipoints(mesh.globe))

    fit_GW <- bru(cmp, lik_W, lik_G,
                  options = list(bru_max_iter = 1,
                                 num.threads  = 8,
                                 verbose      = FALSE))

    if(config$progress){
      prog <- "Fitting done"
      write.table(prog, file = paste0(config$save.path,"Progression.txt"),append=TRUE)
    }


    if(options$excess.dist.fam=="E"){
      # Generate intensities and g
      if(options$W.model=="M1"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW,
                               newdata   = mesh.globe.df, #SpatialPointsDataFrame(coords=mesh.globe.df,data=data.frame(rep(NA,nrow(mesh.globe.df)))),
                               formula   = ~list(Int.W + log(k) - 3*Int.G - 3*zeta,
                                                 exp(Int.G+zeta),
                                                 exp(-Int.G-zeta)),
                               seed      = options$seed)
      }else if(options$W.model=="M2"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW,
                               newdata   = mesh.globe.df, #SpatialPointsDataFrame(coords=mesh.globe.df,data=data.frame(rep(NA,nrow(mesh.globe.df)))),
                               formula   = ~list(Int.W + log(k) - phi,
                                                 exp(Int.G+zeta),
                                                 exp(-Int.G-zeta)),
                               seed      = options$seed)
      }else if(options$W.model=="M3"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW,
                               newdata   = mesh.globe.df, #SpatialPointsDataFrame(coords=mesh.globe.df,data=data.frame(rep(NA,nrow(mesh.globe.df)))),
                               formula   = ~list(Int.W + log(k) - 3*Int.G - 3*zeta -phi,
                                                 exp(Int.G+zeta),
                                                 exp(-Int.G-zeta)),
                               seed      = options$seed)
      }
    }else if(options$excess.dist.fam=="GP"){
      xi.est <- fit_GW$summary.hyperpar$mean[1]
      if(options$W.model=="M1"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW ,
                               newdata   = mesh.globe.df,
                               formula   = ~list(Int.W + 3*Int.G + 3*zeta + log(k),
                                                 ((1-alpha)^(-xi.est)-1)/(xi.est*exp(Int.G+zeta)),
                                                 xi.est*exp(Int.G+zeta)/((1-alpha)^(-xi.est)-1)),
                               seed      = options$seed)
      }else if(options$W.model=="M2"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW ,
                               newdata   = mesh.globe.df,
                               formula   = ~list(Int.W - phi + log(k),
                                                 ((1-alpha)^(-xi.est)-1)/(xi.est*exp(Int.G+zeta)),
                                                 xi.est*exp(Int.G+zeta)/((1-alpha)^(-xi.est)-1)),
                               seed      = options$seed)
      }else if(options$W.model=="M3"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW ,
                               newdata   = mesh.globe.df,
                               formula   = ~list(Int.W + 3*Int.G + 3*zeta - phi + log(k),
                                                 ((1-alpha)^(-xi.est)-1)/(xi.est*exp(Int.G+zeta)),
                                                 xi.est*exp(Int.G+zeta)/((1-alpha)^(-xi.est)-1),
                                                 phi,
                                                 Int.W,
                                                 log(k)),
                               seed      = options$seed)
        fitted.mod$phi[[i]]         <- lapply(post_joint,function(xx) xx[[4]])
        fitted.mod$Int.W[[i]]         <- lapply(post_joint,function(xx) xx[[5]])
        fitted.mod$log.k[[i]]         <- lapply(post_joint,function(xx) xx[[6]])
      }
    }


    if(config$progress){
      prog <- "Sampling done"
      write.table(prog, file = paste0(config$save.path,"Progression.txt"),append=TRUE)
    }

    # w <- lapply(post_joint,function(xx) sample.lgcp(mesh.globe, xx[[1]], ignore.CRS=TRUE)@coords)
    # fitted.mod$w[[i]] <- w
    #
    # fitted.mod$g[[i]] <- list()
    # for(j in 1:options$N.GW){
    #   A <- inla.mesh.projector(mesh=mesh.globe,loc=w[[j]])$proj$A
    #   fitted.mod$g[[i]][[j]] <- as.vector(A %*% post_joint[[j]][[2]])
    # }
    #
    fitted.mod$log.L[[i]] <- lapply(post_joint,function(xx) xx[[1]])
    fitted.mod$g[[i]]         <- lapply(post_joint,function(xx) xx[[2]])
    fitted.mod$G[[i]]         <- lapply(post_joint,function(xx) xx[[3]])

    if(config$save){
      filenm <- paste0(config$file.nm,".RData")
      save(fitted.mod,file = paste0(config$save.path,filenm))
    }
  }
  fitted.mod
}

#' Title
#'
#' @param fitted.mod
#' @param t
#' @param q.prime
#'
#' @import excursions
#' @return
#' @noRd
#'
#' @examples
return_set_3d <- function(fitted.mod,alpha=0.05,t){
  ret_set_list <- list(t=t)
  for(k in 1:length(t)){
    K <- log(t[k]*(1-fitted.mod$options$q))
    n.mesh <- nrow(fitted.mod$mesh$loc)
    n.samp <- fitted.mod$options$N.Qq*fitted.mod$options$N.GW
    post.ret.set <- matrix(NA,nrow=n.samp,ncol=n.mesh)
    cnt <- 1
    for(i in 1:fitted.mod$options$N.Qq){
      for(j in 1:fitted.mod$options$N.GW){
        ret.set <- fitted.mod$Qq[,i] + K*fitted.mod$G[[i]][[j]]
        post.ret.set[cnt,] <- ret.set
        cnt <- cnt+1
      }
    }

    excurs <- simconf.mc(samples = t(post.ret.set),alpha = 1-alpha)

    ret_set_list[[k+1]] <- list(samp  = post.ret.set,
                                mean  = apply(post.ret.set,2,mean),
                                lower = excurs$a,
                                upper = excurs$b)
  }
  return(ret_set_list)
}

#' Title
#'
#' @param fitted.Qq
#' @param xlab
#' @param ylab
#' @param zlab
#' @param alpha
#' @param surface
#'
#' @import rgl
#' @import pracma
#' @return
#' @noRd
#'
#' @examples
plot_Qq_3d <- function(fitted.Qq,surface="mean",alpha=0.05,xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3])){

  X <- fitted.Qq$X
  mesh.globe <- fitted.Qq$mesh
  options <- fitted.Qq$options

  polar.coord    <- cart2sph(as.matrix(X))
  theta <- polar.coord[,1]
  phi   <- polar.coord[,2]
  R     <- polar.coord[,3]
  cart.coord.on.sphere <- sph2cart(cbind(theta, phi, rep(1, nrow(X))))
  dfRphitheta <- data.frame(R=R, theta=theta, phi=phi,
                            x=cart.coord.on.sphere[,1],
                            y=cart.coord.on.sphere[,2],
                            z=cart.coord.on.sphere[,3],
                            X1=X[,1], X2=X[,2], X3=X[,3])

  data.unit.sphere <- sph2cart(cbind(as.matrix(dfRphitheta[,2:3]),1))
  A.Qq <- inla.mesh.projector(mesh=mesh.globe,loc=data.unit.sphere)$proj$A

  exc.inds <- c()
  for(i in 1:options$N.Qq){
    Qq.thetaphi <- as.vector(A.Qq %*% fitted.Qq$Qq[,i])
    exc.inds  <- c(exc.inds,which(dfRphitheta$R - Qq.thetaphi > 0))
  }
  exc.inds <- unique(exc.inds)

  X.exc <- X[exc.inds,]

  Qqs <- fitted.Qq$Qq

  if(surface=="mean"){
    partial.Qq <- apply(Qqs,1,mean)
  }else if(surface=="lower"){
    excurs <- simconf.mc(samples = Qqs,alpha = 1-alpha)
    partial.Qq <- excurs$a #apply(all.Gs,2,function(xx) quantile(xx,0.025))
  }else if(surface=="upper"){
    excurs <- simconf.mc(samples = Qqs,alpha = 1-alpha)
    partial.Qq <- excurs$b #apply(all.Gs,2,function(xx) quantile(xx,0.975))
  }

  N <- 100
  locs.polar     <- as.matrix(data.frame(expand.grid(theta=seq(-pi, pi, len=N), phi = seq(-pi/2, pi/2, len=N)), r=rep(1, 10*10)))
  locs.cartesian <- sph2cart(locs.polar)

  A <- inla.mesh.projector(mesh=fitted.Qq$mesh,loc=locs.cartesian)$proj$A
  Qq.on.grid <- as.vector(A %*% partial.Qq)
  Qq.m <- matrix(Qq.on.grid, nrow=N, ncol=N)

  ## plot estimated unit level set
  plot3d(X.exc,xlab=xlab,ylab=ylab,zlab=zlab)
  surface3d(x=locs.cartesian[,1]*Qq.m,
            y=locs.cartesian[,2]*Qq.m,
            z=locs.cartesian[,3]*Qq.m, alpha=.3, col="gray")
  axes3d(edges="bbox")
}

#' Title
#'
#' @param fitted.mod
#' @param surface
#' @param xlab
#' @param ylab
#' @param zlab
#' @param alpha
#' @param surf.col
#'
#' @import rgl
#' @import pracma
#' @return
#' @noRd
#'
#' @examples
plot_G_3d <- function(fitted.mod,alpha=0.05,surface,xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3]),surf.col="grey"){

  X <- fitted.mod$X
  mesh.globe <- fitted.mod$mesh
  options <- fitted.mod$options

  polar.coord    <- cart2sph(as.matrix(X))
  theta <- polar.coord[,1]
  phi   <- polar.coord[,2]
  R     <- polar.coord[,3]
  cart.coord.on.sphere <- sph2cart(cbind(theta, phi, rep(1, nrow(X))))
  dfRphitheta <- data.frame(R=R, theta=theta, phi=phi,
                            x=cart.coord.on.sphere[,1],
                            y=cart.coord.on.sphere[,2],
                            z=cart.coord.on.sphere[,3],
                            X1=X[,1], X2=X[,2], X3=X[,3])

  data.unit.sphere <- sph2cart(cbind(as.matrix(dfRphitheta[,2:3]),1))
  A.Qq <- inla.mesh.projector(mesh=mesh.globe,loc=data.unit.sphere)$proj$A

  exc.inds <- c()
  for(i in 1:length(fitted.mod$G)){
    Qq.thetaphi <- as.vector(A.Qq %*% fitted.mod$Qq[,i])
    exc.inds  <- c(exc.inds,which(dfRphitheta$R - Qq.thetaphi > 0))
  }
  exc.inds <- unique(exc.inds)
  X.exc <- X[exc.inds,]/log(nrow(X)/2)


  all.Gs <- matrix(NA,length(fitted.mod$G)*options$N.GW,nrow(mesh.globe$loc))
  cnt <- 1
  for(i in 1:length(fitted.mod$G)){
    for(j in 1:length(fitted.mod$G[[1]])){
      all.Gs[cnt,] <- fitted.mod$G[[i]][[j]]
      cnt <- cnt + 1
    }
  }

  if(surface=="mean"){
    partial.G <- apply(all.Gs,2,mean)
  }else if(surface=="lower"){
    excurs <- simconf.mc(samples = t(all.Gs),alpha = 1-alpha)
    partial.G <- excurs$a #apply(all.Gs,2,function(xx) quantile(xx,0.025))
  }else if(surface=="upper"){
    excurs <- simconf.mc(samples = t(all.Gs),alpha = 1-alpha)
    partial.G <- excurs$b #apply(all.Gs,2,function(xx) quantile(xx,0.975))
  }

  N <- 100
  locs.polar     <- as.matrix(data.frame(expand.grid(theta=seq(-pi, pi, len=N), phi = seq(-pi/2, pi/2, len=N)), r=rep(1, 10*10)))
  locs.cartesian <- sph2cart(locs.polar)

  A <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=locs.cartesian)$proj$A
  partial.G.on.grid <- as.vector(A %*% partial.G)
  partial.G.m <- matrix(partial.G.on.grid, nrow=N, ncol=N)

  ## plot estimated unit level set
  plot3d(X/log(nrow(X)/2),xlab=xlab,ylab=ylab,zlab=zlab,
         xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),zlim=c(-1.2,1.2),cex.lab=3)
  surface3d(x=locs.cartesian[,1]*partial.G.m,
            y=locs.cartesian[,2]*partial.G.m,
            z=locs.cartesian[,3]*partial.G.m, alpha=.3, col=surf.col)
  axes3d(edges="bbox")
}

#' Title
#'
#' @param fitted.mod
#' @param list_ret_sets
#' @param surface
#' @param cex.pts
#' @param cex.axis
#' @param xyzlim
#' @param xlab
#' @param ylab
#' @param zlab
#'
#' @import rgl
#' @import pracma
#' @import excursions
#' @return
#' @noRd
#'
#' @examples
plot_return_bdry_3d <- function(fitted.mod,list_ret_sets,surface="mean",cex.pts=0.4,cex.axis=1.4,xyzlim=c(0,0),xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3])){

  t <- rev(sort(list_ret_sets$t))

  if(sum(xyzlim==c(0,0))==2){
    plot3d(fitted.mod$X,xlab=xlab,ylab=ylab,zlab=zlab)
  }else{
    plot3d(fitted.mod$X,xlab=xlab,ylab=ylab,zlab=zlab,
           xlim=xyzlim,ylim=xyzlim,zlim=xyzlim)
  }

  post.ret.set <- list_ret_sets[[2]]
  if(surface=="mean"){
    partial.G <- post.ret.set$mean
  }else if(surface=="lower"){
    partial.G <- post.ret.set$lower
  }else if(surface=="upper"){
    partial.G <- post.ret.set$upper
  }

  N <- 100
  locs.polar     <- as.matrix(data.frame(expand.grid(theta=seq(-pi, pi, len=N), phi = seq(-pi/2, pi/2, len=N)), r=rep(1, 10*10)))
  locs.cartesian <- sph2cart(locs.polar)

  A <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=locs.cartesian)$proj$A
  partial.G.on.grid <- as.vector(A %*% partial.G)
  partial.G.m <- matrix(partial.G.on.grid, nrow=N, ncol=N)

  surface3d(x=locs.cartesian[,1]*partial.G.m,
            y=locs.cartesian[,2]*partial.G.m,
            z=locs.cartesian[,3]*partial.G.m, alpha=.3, col="gray")
  axes3d(edges="bbox")
}

#' Title
#'
#' @param w
#' @param x
#' @param y
#' @param z
#'
#' @import pracma
#' @return
#' @noRd
#'
#' @examples
get_rinfsup_3d <- function(w,x,y,z){
  theta <- w[1]; phi <- w[2]

  s2c1 <- sph2cart(c(theta,phi,1))
  R.x <- x/s2c1[1]
  R.y <- y/s2c1[2]
  R.z <- z/s2c1[3]

  r.x <- cbind(x,
               R.x*s2c1[2],
               R.x*s2c1[3],
               R.x)
  r.y <- cbind(R.y*s2c1[1],
               y,
               R.y*s2c1[3],
               R.y)
  r.z <- cbind(R.z*s2c1[1],
               R.z*s2c1[3],
               z,
               R.z)
  r.infsup <- apply(rbind(r.x,r.y,r.z),1,
                    function(xx){
                      ifelse(xx[1]>=x[1]&xx[1]<=x[2]&xx[2]>=y[1]&xx[2]<=y[2]&xx[3]>=z[1]&xx[3]<=z[2],
                             xx[4],0)
                    })
  if(sum(r.infsup)==0){
    return(c(0,0))
  }else{
    return(sort(r.infsup[r.infsup>0]))
  }
}
