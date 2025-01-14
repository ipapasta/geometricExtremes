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
    prog <- paste0("Begin quantile regression: ", config$file.nm)
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
                   domain = list(coordinates=mesh.globe))#,ips = fm_int(domain=mesh.globe))

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
#' @param alpha
#' @param conf
#' @param t
#'
#' @import excursions
#' @return
#' @noRd
#'
#' @examples
return_set_3d <- function(fitted.mod,alpha=0.05,conf="sim",t){
  ret_set_list <- list(pars=list(t=t,
                                 marginals = "Laplace"),
                       X=fitted.mod$X)

  for(k in 1:length(t)){
    K <- log(t[k]*(1-fitted.mod$options$q))
    n.mesh <- nrow(fitted.mod$mesh$loc)
    n.samp <- fitted.mod$options$N.Qq*fitted.mod$options$N.GW
    post.ret.set <- matrix(NA,nrow=n.samp,ncol=n.mesh)
    cnt <- 1
    for(i in 1:fitted.mod$options$N.Qq){
      for(j in 1:fitted.mod$options$N.GW){
        message("Q_q: ",toString(i),"/",fitted.mod$options$N.Qq," and G,W: ", toString(j),"/",fitted.mod$options$N.GW)
        ret.set <- fitted.mod$Qq[,i] + K*fitted.mod$G[[i]][[j]]
        post.ret.set[cnt,] <- ret.set
        cnt <- cnt+1
      }
    }

    excurs <- simconf.mc(samples = t(post.ret.set),alpha = alpha)
    if(conf=="sim"){
      low <- excurs$a
      up <- excurs$b
    }else if(conf=="marg"){
      low <- excurs$a.marg
      up <- excurs$b.marg
    }

    ret_set_list[[k+2]] <- list(samp  = post.ret.set,
                                mean  = fitted.mod$mesh$loc*apply(post.ret.set,2,mean),
                                lower = fitted.mod$mesh$loc*low,
                                upper = fitted.mod$mesh$loc*up)
  }
  return(ret_set_list)
}


#' Title
#'
#' @param fitted.mod
#' @param alpha
#' @param conf
#' @param t
#'
#' @return
#' @export
#'
#' @examples
return_set_3d_isotropic <- function(fitted.mod,alpha=0.05,conf="sim",t){
  ret_set_list <- list(pars=list(t=t,
                                 marginals = "Laplace"),
                       X=fitted.mod$X)

  get_integrating_constant_3d <- function(mesh,log.f,n.MC.samp){
    Us <- t(apply(matrix(rnorm(3*n.MC.samp),ncol=3),1,function(x) x/sqrt(sum(x^2))))

    A <- inla.mesh.projector(mesh=mesh,loc=Us)$proj$A
    lambda_at_Us <- exp(as.vector(A %*%log.f))
    int_const <- 4*pi*mean(lambda_at_Us)
    int_const
  }

  fW.uniform  <- 1/(4*pi)
  # ret_set_list <- return_set_3d(fitted.mod,alpha,conf,t)
  #
  # plot_return_bdry_3d(fitted.mod,ret_set_list,"mean")

  for(k in 1:length(t)){
    K <- log(t[k]*(1-fitted.mod$options$q))
    n.mesh <- nrow(fitted.mod$mesh$loc)
    n.samp <- fitted.mod$options$N.Qq*fitted.mod$options$N.GW
    post.ret.set <- matrix(NA,nrow=n.samp,ncol=n.mesh)
    cnt <- 1
    for(i in 1:fitted.mod$options$N.Qq){
      for(j in 1:fitted.mod$options$N.GW){
        message("Q_q: ",toString(i),"/",fitted.mod$options$N.Qq," and G,W: ", toString(j),"/",length(fitted.mod$log.L[[1]]))
        C <- get_integrating_constant_3d(fitted.mod$mesh,fitted.mod$log.L[[i]][[j]],100000)
        f_W <- exp(fitted.mod$log.L[[i]][[j]])/C
        ret.set <- fitted.mod$Qq[,i] + K*fitted.mod$G[[i]][[j]]
        ret.set.omni <- ret.set + log(f_W/fW.uniform) * fitted.mod$G[[i]][[j]]
        post.ret.set[cnt,] <- ret.set.omni
        cnt <- cnt+1
      }
    }

    excurs <- simconf.mc(samples = t(post.ret.set),alpha = alpha)
    if(conf=="sim"){
      low <- excurs$a
      up <- excurs$b
    }else if(conf=="marg"){
      low <- excurs$a.marg
      up <- excurs$b.marg
    }

    ret_set_list[[k+2]] <- list(samp  = post.ret.set,
                                mean  = fitted.mod$mesh$loc*apply(post.ret.set,2,mean),
                                lower = fitted.mod$mesh$loc*low,
                                upper = fitted.mod$mesh$loc*up)
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
#' @param conf
#'
#' @import rgl
#' @import pracma
#' @return
#' @noRd
#'
#' @examples
plot_Qq_3d <- function(fitted.Qq,surface="mean",alpha=0.05,conf="marg",xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3])){

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

  # exc.inds <- c()
  # for(i in 1:options$N.Qq){
  #   Qq.thetaphi <- as.vector(A.Qq %*% fitted.Qq$Qq[,i])
  #   exc.inds  <- c(exc.inds,which(dfRphitheta$R - Qq.thetaphi > 0))
  # }
  # exc.inds <- unique(exc.inds)
  #
  # X.exc <- X[exc.inds,]

  Qqs <- fitted.Qq$Qq

  if(surface=="mean"){
    partial.Qq <- apply(Qqs,1,mean)
  }else if(surface=="lower"){
    excurs <- simconf.mc(samples = Qqs,alpha = alpha)
    if(conf=="sim"){
      partial.Qq <- excurs$a
    }else if(conf=="marg"){
      partial.Qq <- excurs$a.marg
    }
  }else if(surface=="upper"){
    excurs <- simconf.mc(samples = Qqs,alpha = alpha)
    if(conf=="sim"){
      partial.Qq <- excurs$b
    }else if(conf=="marg"){
      partial.Qq <- excurs$b.marg
    }
  }

  partial.Qq.at.W <- as.vector(A.Qq %*% partial.Qq)
  exc.inds <- R > partial.Qq.at.W
  X.exc <- X[exc.inds,]
  X.non.exc <- X[!exc.inds,]

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
  points3d(X.non.exc,col="grey40")
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
#' @param conf
#'
#' @import rgl
#' @import pracma
#' @return
#' @noRd
#'
#' @examples
plot_G_3d <- function(fitted.mod,alpha=0.05,conf="marg",surface,xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3]),surf.col="grey"){

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
    excurs <- simconf.mc(samples = t(all.Gs),alpha = alpha)
    if(conf=="sim"){
      partial.G <- excurs$a
    }else if(conf=="marg"){
      partial.G <- excurs$a.marg
    }
  }else if(surface=="upper"){
    excurs <- simconf.mc(samples = t(all.Gs),alpha = alpha)
    if(conf=="sim"){
      partial.G <- excurs$b
    }else if(conf=="marg"){
      partial.G <- excurs$b.marg
    }
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
plot_return_bdry_3d <- function(fitted.mod,list_ret_sets,surface="mean",cex.pts=0.4,cex.axis=1.4,xyzlim=c(0,0),xlab=expression(X[1]),ylab=expression(X[2]),zlab=expression(X[3]),col="grey"){

  t <- rev(sort(list_ret_sets$pars$t))

  xlim <- range(fitted.mod$X[,1])
  ylim <- range(fitted.mod$X[,2])
  zlim <- range(fitted.mod$X[,3])
  if(sum(xyzlim==c(0,0))==2){

    plot3d(NA,xlab=xlab,ylab=ylab,zlab=zlab,
           xlim=xlim,ylim=ylim,zlim=zlim)
  }else{
    plot3d(NA,xlab=xlab,ylab=ylab,zlab=zlab,
           xlim=xyzlim,ylim=xyzlim,zlim=xyzlim)
  }

  post.ret.set <- list_ret_sets[[3]]

  if(list_ret_sets$pars$marginals=="Laplace"){
    if(surface=="mean"){
      partial.G <- apply(post.ret.set$mean,1,function(xx) sqrt(sum(xx^2)))
    }else if(surface=="lower"){
      partial.G <- apply(post.ret.set$lower,1,function(xx) sqrt(sum(xx^2)))
    }else if(surface=="upper"){
      partial.G <- apply(post.ret.set$upper,1,function(xx) sqrt(sum(xx^2)))
    }

    N <- 100
    locs.polar     <- as.matrix(data.frame(expand.grid(theta=seq(-pi, pi, len=N), phi = seq(-pi/2, pi/2, len=N)), r=rep(1, 10*10)))
    locs.cartesian <- sph2cart(locs.polar)

    A <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=locs.cartesian)$proj$A
    partial.G.on.grid <- as.vector(A %*% partial.G)
    partial.G.m <- matrix(partial.G.on.grid, nrow=N, ncol=N)

    surface3d(x=locs.cartesian[,1]*partial.G.m,
              y=locs.cartesian[,2]*partial.G.m,
              z=locs.cartesian[,3]*partial.G.m, alpha=.3, col=col)

    R <- apply(fitted.mod$X,1,function(xx) sqrt(sum(xx^2)))
    W <- fitted.mod$X/R
    A.obs <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=W)$proj$A
    r_Qq_at_w_obs  <- as.vector(A.obs %*%partial.G)

    X.exc <- fitted.mod$X[R>r_Qq_at_w_obs,]
    points3d(X.exc[,1],X.exc[,2],X.exc[,3],col="black")
    X.nonexc <- fitted.mod$X[!(R>r_Qq_at_w_obs),]
    points3d(X.nonexc[,1],X.nonexc[,2],X.nonexc[,3],col="grey60")


  }else{
    if(surface=="mean"){
      lines3d(post.ret.set$mean,col="red")
    }else if(surface=="lower"){
      lines3d(post.ret.set$lower,col="red")
    }else if(surface=="upper"){
      lines3d(post.ret.set$upper,col="red")
    }
  }

  # if(){
  #   lines3d(ret_sets[[3]]$mean,col="red")
  # }

  # partial.G.m <- matrix(partial.G, nrow=N, ncol=N)
  # surface3d(x=fitted.mod$mesh$loc[,1]*partial.G.m,
  #           y=fitted.mod$mesh$loc[,2]*partial.G.m,
  #           z=fitted.mod$mesh$loc[,3]*partial.G.m, alpha=.3, col="gray")
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

#' Title
#'
#' @param fitted.mod
#'
#' @return
#' @export
#'
#' @examples
X_to_uniform_on_Ball_3d <- function(fitted.mod){

  get_integrating_constant_3d <- function(mesh,log.f,n.MC.samp){
    Us <- t(apply(matrix(rnorm(3*n.MC.samp),ncol=3),1,function(x) x/sqrt(sum(x^2))))

    A <- inla.mesh.projector(mesh=mesh,loc=Us)$proj$A
    lambda_at_Us <- exp(as.vector(A %*%log.f))
    int_const <- 4*pi*mean(lambda_at_Us)
    int_const
  }

  X <- fitted.mod$X

  U_on_Ball <- list()

  ## Convert from cartesian to polar coordinates
  d <- dim(X)[2]
  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- X/R

  A.obs <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=W)$proj$A
  if(fitted.mod$options$use.mean.Qq==TRUE){
    N.Qq <- 1
  }else{
    N.Qq <- fitted.mod$options$N.Qq
  }

  fW.uniform  <- 1/(4*pi)

  for(i in 1:N.Qq){

    U_on_Ball[[i]] <- list()

    if(fitted.mod$options$use.mean.Qq==TRUE){
      r_Qq <- fitted.mod$Qq.mean$mean
    }else{
      r_Qq <- fitted.mod$Qq[,i]
    }
    r_Qq_at_w_obs  <- as.vector(A.obs %*%r_Qq)

    for(j in 1:length(fitted.mod$log.L[[1]])){
      message("Q_q: ",toString(i),"/",N.Qq," and G,W: ", toString(j),"/",length(fitted.mod$log.L[[1]]))
      C <- get_integrating_constant_3d(fitted.mod$mesh,fitted.mod$log.L[[i]][[j]],100000)
      f_W_at_wobs <- exp(as.vector(A.obs %*%fitted.mod$log.L[[i]][[j]]))/C

      r_G_at_wobs <- as.vector(A.obs %*% fitted.mod$G[[i]][[j]])

      location    <- r_Qq_at_w_obs + log(f_W_at_wobs/fW.uniform) * r_G_at_wobs
      Rtilde  <-  ((R-location)/r_G_at_wobs)
      ind.pos <- which(Rtilde > 0)

      Rtilde.pos.unif <- pexp(Rtilde[ind.pos])
      XX <- cbind(Rtilde.pos.unif^(1/d)*W[ind.pos,1],
                  Rtilde.pos.unif^(1/d)*W[ind.pos,2],
                  Rtilde.pos.unif^(1/d)*W[ind.pos,3])

      U_on_Ball[[i]][[j]] <- XX
    }
  }
  return(U_on_Ball)
}
