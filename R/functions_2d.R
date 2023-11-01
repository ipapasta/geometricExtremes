#' Obtain posterior samples from the radial function of quantile set, r_{Q_q}, through Bayesian gamma quantile regression.
#'
#' @param X                 matrix (n by d) of observations.
#' @param options            hyper-parameters for the model fits.
#' @param config            save configurations.
#' @param return_fitted_obj boolean to return the fitted Q_q.
#'
#' @import INLA
#' @import inlabru
#' @return Posterior samples from the radial function of quantile set, r_{Q_q}.
#' @noRd
#'
#' @examples
fit_Qq_2d <- function(X,options,config,return_fitted_obj=F){

  if(config$progress){
    prog <- paste0("Begin quantile regression")
    write.table(prog, file = paste0(config$save.path,"Progression.txt"))
  }

  ## Convert from cartesian to polar coordinates
  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- atan2(y=X[,2], x=X[,1])
  dfRW <- data.frame(R=R, W=W)

  ## Definition of mesh and prior
  if(sum(range(options$mesh.knots.2d)==c(-pi,pi))==2){
    mesh             <- inla.mesh.1d(sort(options$mesh.knots.2d),
                                     boundary="cyclic",degree=1)
  }else{
    mesh             <- inla.mesh.1d(sort(options$mesh.knots.2d),
                                     boundary="neumann",degree=1)
  }

  circular.Matern  <- inla.spde2.pcmatern(mesh, prior.range = options$QR.prior.range,
                                          prior.sigma = options$QR.prior.sigma)

  cmp.Qq  <- R ~  Intercept(1) + spde1(W, model=circular.Matern, mapper=bru_mapper(mesh, indexed=TRUE))
  fit.Qq <- bru(cmp.Qq,
                data = as.data.frame(dfRW),
                family="gamma",
                options=list(num.threads=8, verbose = FALSE, safe=TRUE),
                control.family = list(control.link = list(model="quantile", quantile=options$q)))

  Qq.post.samp <- generate(n.samples = options$N.Qq,
                           object    = fit.Qq,
                           newdata   = data.frame(W=mesh$loc),#data.frame(x=X[,1], y=X[,2]),
                           formula   = ~ exp(Intercept+spde1),
                           seed      = options$seed)

  prop.exc <- rep(NA,options$N.Qq)
  A.Qq <- inla.mesh.projector(mesh=mesh,loc=dfRW$W)$proj$A
  for(i in 1:options$N.Qq){
    dfRW$thres.gamma  <- as.vector(A.Qq %*%Qq.post.samp[,i])
    is.excess         <- ifelse(dfRW$R - dfRW$thres.gamma > 0 , 1, 0)

    prop.exc[i]       <- sum(is.excess)/options$N.X
  }

  Qq.mean      <- predict(fit.Qq,
                          data.frame(W=mesh$loc),
                          ~ exp(Intercept+spde1),
                          seed      = options$seed)

  fitted.Qq <- list(X         = X,
                    mesh      = mesh,
                    options    = options,
                    Qq      = Qq.post.samp,
                    Qq.mean  = Qq.mean,
                    prop.exc  = prop.exc)
  if(return_fitted_obj==T){
    fitted.Qq$QR.fit.obj <- fit.Qq
  }

  if(config$save){
    filenm <- paste0(config$file.nm,"_Qq.RData")
    save(fitted.Qq,file = paste0(config$save.path,filenm))
  }

  if(config$progress){
    prog <- paste0("Quantile regression done")
    write.table(prog, file = paste0(config$save.path,"Progression.txt"),append=TRUE)
  }
  return(fitted.Qq)
}

#' Obtain posterior samples from the radial function of the limit set G and the angle set W
#'
#' @param fitted.Qq returned object from the function fit.Q_q.2d.
#' @param config      save configurations.
#'
#' @import INLA
#' @import inlabru
#' @import utils
#' @return
#' @noRd
#'
#' @examples
fit_GL_2d <- function(fitted.Qq,config){

  X <- fitted.Qq$X
  options <- fitted.Qq$options

  ## Convert from cartesian to polar coordinates
  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- atan2(y=X[,2], x=X[,1])
  dfRW <- data.frame(R=R, W=W)

  ## Definition of mesh and prior
  if(sum(range(options$mesh.knots.2d)==c(-pi,pi))==2){
    mesh             <- inla.mesh.1d(sort(options$mesh.knots.2d),
                                     boundary="cyclic",degree=1)
  }else{
    mesh             <- inla.mesh.1d(sort(options$mesh.knots.2d),degree=1)
  }

  Qq.post.samp <- fitted.Qq$Qq
  if(options$excess.dist.fam=="E"){
    fitted.mod <- list(X         = X,
                       mesh     = mesh,
                       options  = options,
                       Qq       = fitted.Qq$Qq,
                       Qq.mean  = fitted.Qq$Qq.mean,
                       prop.exc = fitted.Qq$prop.exc,
                       g        = list(),
                       G        = list(),
                       log.L    = list(),
                       phi      = list(),
                       Int.W    = list(),
                       log.k    = list())
  }else{
    fitted.mod <- list(X        = X,
                       mesh     = mesh,
                       options  = options,
                       Qq       = fitted.Qq$Qq,
                       Qq.mean  = fitted.Qq$Qq.mean,
                       prop.exc = fitted.Qq$prop.exc,
                       G        = list(),
                       xi       = list(),
                       log.L    = list(),
                       phi      = list(),
                       Int.W    = list(),
                       log.k    = list())
  }

  ## Model definition
  circular.Matern_zeta  <- inla.spde2.pcmatern(mesh, prior.range = options$zeta.prior.range,
                                               prior.sigma = options$zeta.prior.sigma)

  # B.sigma <- cbind(0, 1, 0, mesh$loc)
  # B.range <- cbind(0, 0, 1, mesh$loc)
  # circular.Matern_zeta  <- rspde.matern(mesh = mesh,B.sigma = B.sigma,B.range = B.range,
  #                                       parameterization = "matern")

  if(options$W.model=="M1"){
    cmp.GW <- ~ 0 + Int.G(1) + Int.W(1) + zeta(W, model = circular.Matern_zeta)
  }else if(options$W.model=="M2" | options$W.model=="M3"){
    weights.domain     <- ipoints(domain=mesh, samplers=c(-pi,pi))
    locs               <- weights.domain$x
    A.constr           <- inla.spde.make.A(mesh=mesh, loc=locs, weights=weights.domain$weight*499/(2*pi),#weights.domain$weight,
                                           block=rep(1, nrow(weights.domain)))
    A.constr <- matrix(c(0.5,rep(1,mesh$n-1)),nrow=1)
    circular.Matern_phi  <- inla.spde2.pcmatern(mesh, prior.range = options$phi.prior.range,
                                                prior.sigma = options$phi.prior.sigma,
                                                extraconstr=list(A=as.matrix(A.constr,nrow=1), e=matrix(0, ncol = 1)))
    cmp.GW <- ~ 0 + Int.G(1) + Int.W(1) + zeta(W, model = circular.Matern_zeta) + phi(W, model = circular.Matern_phi)
  }

  A.Qq <- inla.mesh.projector(mesh=mesh,loc=dfRW$W)$proj$A
  if(options$use.mean.Qq==TRUE){
    N.Qq <- 1
  }else{
    N.Qq <- options$N.Qq
  }
  for(i in 1:N.Qq){
    message(paste0("Estimating G and W for Qq ",i,"/",N.Qq))
    if(config$progress){
      prog <- paste0("Estimating G and W for Qq ",i,"/",options$N.Qq)
      write.table(prog, file = paste0(config$save.path,"Progression.txt"),append=TRUE)
    }

    if(options$use.mean.Qq==TRUE){
      dfRW$thres.gamma  <- as.vector(A.Qq %*% fitted.mod$Qq.mean$mean)
    }else{
      dfRW$thres.gamma  <- as.vector(A.Qq %*%Qq.post.samp[,i])
    }

    dfRW$is.excess        <- ifelse(dfRW$R - dfRW$thres.gamma > 0 , 1, 0)
    dfRW$excess.gamma     <- pmax(dfRW$R - dfRW$thres.gamma, 0)

    dfRW.exc.only         <- dfRW[dfRW$excess.gamma > 0,] #dfRW %>% dplyr::filter(excess.gamma > 0)
    dfRW$excess.gamma     <- ifelse(dfRW$excess.gamma>0,dfRW$excess.gamma, NA)

    if(options$W.data=="AllW"){
      log.k <- log(nrow(dfRW))
    }else if(options$W.data=="ExcOnly"){
      log.k <- log(nrow(dfRW.exc.only))
    }

    alpha <- options$alpha
    if(options$excess.dist.fam=="E"){
      if(options$W.model=="M1"){
        form.W  <- W ~ 0 + Int.W - 2*Int.G - 2*zeta + log.k
      }else if(options$W.model=="M2"){
        form.W  <- W ~ 0 + Int.W - phi + log.k
      }else if(options$W.model=="M3"){
        form.W  <- W ~ 0 + Int.W - 2*Int.G - 2*zeta - phi + log.k
      }
      form.G <- excess.gamma ~ 0 + Int.G + zeta
      lik_G <- like("exponential", formula=form.G, data=dfRW.exc.only)
    }else if(options$excess.dist.fam=="GP"){
      if(options$W.model=="M1"){
        form.W  <- W ~ 0 + Int.W + 2*Int.G + 2*zeta + log.k
      }else if(options$W.model=="M2"){
        form.W  <- W ~ 0 + Int.W - phi + log.k
      }else if(options$W.model=="M3"){
        form.W  <- W ~ 0 + Int.W + 2*Int.G + 2*zeta - phi + log.k
      }
      form.G  <- excess.gamma ~ 0 + Int.G + zeta
      lik_G <- like("gp", formula=form.G, data=dfRW.exc.only,
                    control.family = list(control.link = list(quantile = alpha)))
    }


    if(options$W.data=="AllW"){
      lik_W  <- like("cp", formula=form.W, data=dfRW, domain=list(W=mesh))
    }else if(options$W.data=="ExcOnly"){
      lik_W  <- like("cp", formula=form.W, data=dfRW.exc.only, domain=list(W=mesh))
    }


    fit_GW <- bru(cmp.GW, lik_W, lik_G,
                  options = list(bru_max_iter = 1,
                                 num.threads=8,
                                 verbose = FALSE))

    # Generate intensities and g
    if(options$excess.dist.fam=="E"){
      if(options$W.model=="M1"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW ,
                               newdata   = data.frame(W=mesh$loc),
                               formula   = ~list(Int.W - 2*Int.G - 2*zeta + log.k,
                                                 exp(Int.G+zeta),
                                                 exp(-Int.G-zeta)),
                               seed      = options$seed)
      }else if(options$W.model=="M2"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW ,
                               newdata   = data.frame(W=mesh$loc),
                               formula   = ~list(Int.W - phi + log.k,
                                                 exp(Int.G+zeta),
                                                 exp(-Int.G-zeta),
                                                 phi,
                                                 Int.W,
                                                 log.k),
                               seed      = options$seed)
      }else if(options$W.model=="M3"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW,
                               newdata   = data.frame(W=mesh$loc),
                               formula   = ~list(Int.W - 2*Int.G - 2*zeta - phi + log.k,
                                                 exp(Int.G+zeta),
                                                 exp(-Int.G-zeta),
                                                 phi,
                                                 Int.W,
                                                 log.k),
                               seed      = options$seed)
        fitted.mod$phi[[i]]         <- lapply(post_joint,function(xx) xx[[4]])
        fitted.mod$Int.W[[i]]         <- lapply(post_joint,function(xx) xx[[5]])
        fitted.mod$log.k[[i]]         <- lapply(post_joint,function(xx) xx[[6]])
      }

      fitted.mod$log.L[[i]] <- lapply(post_joint,function(xx) xx[[1]])
      fitted.mod$g[[i]]         <- lapply(post_joint,function(xx) xx[[2]])
      fitted.mod$G[[i]]   <- lapply(post_joint,function(xx) xx[[3]])
    }else if(options$excess.dist.fam=="GP"){
      xi.est <- fit_GW$summary.hyperpar$mean[1]
      if(options$W.model=="M1"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW ,
                               newdata   = data.frame(W=mesh$loc),
                               formula   = ~list(Int.W + 2*Int.G + 2*zeta + log.k,
                                                 xi.est*exp(Int.G+zeta)/((1-alpha)^(-xi.est)-1)),
                               seed      = options$seed)
      }else if(options$W.model=="M2"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW ,
                               newdata   = data.frame(W=mesh$loc),
                               formula   = ~list(Int.W - phi + log.k,
                                                 xi.est*exp(Int.G+zeta)/((1-alpha)^(-xi.est)-1)),
                               seed      = options$seed)
      }else if(options$W.model=="M3"){
        post_joint <- generate(n.samples = options$N.GW,
                               object    = fit_GW ,
                               newdata   = data.frame(W=mesh$loc),
                               formula   = ~list(Int.W + 2*Int.G + 2*zeta - phi + log.k,
                                                 xi.est*exp(Int.G+zeta)/((1-alpha)^(-xi.est)-1),
                                                 phi,
                                                 Int.W,
                                                 log.k),
                               seed      = options$seed)
        fitted.mod$phi[[i]]         <- lapply(post_joint,function(xx) xx[[3]])
        fitted.mod$Int.W[[i]]         <- lapply(post_joint,function(xx) xx[[4]])
        fitted.mod$log.k[[i]]         <- lapply(post_joint,function(xx) xx[[5]])
      }

      fitted.mod$log.L[[i]] <- lapply(post_joint,function(xx) xx[[1]])
      fitted.mod$G[[i]]     <- lapply(post_joint,function(xx) xx[[2]])
      fitted.mod$xi[[i]]    <- xi.est
    }

    if(config$save){
      filenm <- paste0(config$file.nm,"_n",nrow(X),"_q",options$q,".RData")
      save(fitted.mod,file = paste0(config$save.path,filenm))
    }
  }
  if(config$save){
    filenm <- paste0(config$file.nm,"_n",nrow(X),"_q",options$q,".RData")
    save(fitted.mod,file = paste0(config$save.path,filenm))
  }
  fitted.mod
}

#' Sample angles from the LGCP on an interval of the 1-sphere
#'
#' @param n           number of observations w.
#' @param log.L   natural logarithm of the unnormalised density of angles W.
#' @param mesh        INLA object for the cyclic mesh on the 1-sphere.
#' @param limits.W_B  vector of lower and upper limits of angles (polar) for sampling region.
#' @param fitted.mod
#'
#' @import INLA
#' @return Set of n observed angles from log.L.
#' @noRd
#'
#' @examples
sample_LGCP_acc_rej_B <- function(n,fitted.mod,log.L,mesh,limits.W_B){
  mesh <- fitted.mod$mesh

  mesh.len <- length(mesh$loc)

  mesh.mod <- c(mesh$loc,abs(mesh$loc[1]))

  L.mod <- exp(c(log.L,abs(log.L[1]))+log(0.01))

  A             <- inla.mesh.projector(mesh=mesh,loc=limits.W_B)$proj$A
  L.at.limits <- exp(as.vector(A %*%log.L)+log(0.01))

  min.mesh.B <- min(which(mesh.mod>limits.W_B[1]))
  max.mesh.B <- max(which(mesh.mod<limits.W_B[2]))

  log.L.B <- c(L.at.limits[1],
               L.mod[min.mesh.B:max.mesh.B],
               L.at.limits[2])
  mesh.B      <- c(limits.W_B[1],
                   mesh$loc[min.mesh.B:max.mesh.B],
                   limits.W_B[2])

  inds <- cbind(c(1:(length(mesh.B)-1)),c(2:length(mesh.B)))
  step.fn <- apply(inds,1,function(x){max(log.L.B[x[1]],log.L.B[x[2]])})

  ps <- step.fn*diff(mesh.B)/sum(step.fn*diff(mesh.B))

  n.by.mesh.int <- rmultinom(1,n,ps)

  new.ws <- unlist(apply(inds,1,
                         function(x){
                           runif(n.by.mesh.int[x[1]],mesh.B[x[1]],mesh.B[x[2]])
                         }))
  prop.dist <- unlist(apply(inds,1,function(x) rep(step.fn[x[1]],each=n.by.mesh.int[x[1]])))
  acc.ref <- runif(n,0,prop.dist)

  A <- inla.mesh.projector(mesh=mesh,loc=new.ws)$proj$A
  log.L.new.ws <- exp(as.vector(A %*%log.L)+log(0.01))

  new.ws[acc.ref<log.L.new.ws]
}

#' Calculate the probability that W falls in an interval of the 1-sphere
#'
#' @param log.L   natural logarithm of the unnormalised density of angles W.
#' @param mesh        INLA object for the cyclic mesh on the 1-sphere.
#' @param limits.W_B  vector of lower and upper limits of angles (polar) for sampling region.
#'
#' @import INLA
#' @import inlabru
#' @return probability that W falls in the range limits.W_B.
#' @export
#'
#' @examples
get_P_W <- function(log.L,mesh,limits.W_B){

  mesh.len    <- length(mesh$loc)
  mesh.mod    <- c(mesh$loc,abs(mesh$loc[1]))
  L.mod <- exp(c(log.L,abs(log.L[1]))+log(0.01))

  A             <- inla.mesh.projector(mesh=mesh,loc=limits.W_B)$proj$A
  L.at.limits <- exp(as.vector(A %*%log.L)+log(0.01))

  min.mesh.B <- min(which(mesh.mod>limits.W_B[1]))
  max.mesh.B <- max(which(mesh.mod<limits.W_B[2]))

  log.L.B <- c(L.at.limits[1],
               L.mod[min.mesh.B:max.mesh.B],
               L.at.limits[2])
  mesh.B      <- c(limits.W_B[1],
                   mesh$loc[min.mesh.B:max.mesh.B],
                   limits.W_B[2])

  h <- diff(mesh.B)
  integrand <- sapply(c(1:length(h)),
                      function(x){
                        (log.L.B[x[1]]+log.L.B[x[1]+1])*h[x[1]]/2
                      })

  h <- diff(mesh.mod)
  int.const <- sapply(c(1:mesh.len),
                      function(x){
                        (L.mod[x[1]]+L.mod[x[1]+1])*h[x[1]]/2
                      })
  sum(integrand)/sum(int.const)
}

#' Sample n angles from the LGCP on an interval of the 1-sphere
#'
#' @param n        number of observations w.
#' @param log.L   natural logarithm of the unnormalised density of angles W.
#' @param mesh        INLA object for the cyclic mesh on the 1-sphere.
#' @param limits.W_B  vector of lower and upper limits of angles (polar) for sampling region.
#' @param fitted.mod
#'
#' @import INLA
#' @import inlabru
#' @return
#' @noRd
#'
#' @examples
sample_LGCP_fix <- function(n,log.L,fitted.mod,mesh,limits.W_B=NA){
  if(!inherits(limits.W_B,"logical")){
    new.w <- sample_LGCP_acc_rej_B(n,fitted.mod,log.L,mesh,limits.W_B)
    while(length(new.w)<n){
      new.w <- c(new.w,
                 sample_LGCP_acc_rej_B(n-length(new.w),fitted.mod,log.L,mesh,limits.W_B))
    }
    P.W_B <- get_P_W(log.L,
                     mesh,
                     limits.W_B)
    list(w=new.w,P.W_B=P.W_B)
  }else{
    exp.num.samp <- numer_int(exp(log.L),mesh$loc)
    mod.exp.num.samp <- log(n/exp.num.samp)
    new.w <- sample.lgcp(mesh,log.L+mod.exp.num.samp)$x
    while(length(new.w)<n){
      new.w <- c(new.w,
                 sample.lgcp(mesh,log.L+mod.exp.num.samp)$x)
    }
    new.w[1:n]
  }
}

#' Title
#'
#' @param fitted.mod  returned object from the function fit_GL_2d.
#' @param N.w         number of angles to sample.
#' @param S_B         Subset of the 1-sphere S from -pi to pi. NA for S.
#' @param transf.G    boolean, use transf.G if TRUE, G otherwise.
#'
#' @import INLA
#' @import inlabru
#' @return
#' @noRd
#'
#' @examples
sample_QGW_posterior_2d <- function(fitted.mod,N.w,S_B=NA,transf.G=F){

  mesh <- fitted.mod$mesh
  N.Qq <- length(fitted.mod$G)
  N.GW <- length(fitted.mod$G[[1]])

  if(fitted.mod$options$excess.dist.fam=="E"){
    post.samp <- list(w     = list(),
                      P.W_B = list(),
                      g.at.w   = list(),
                      Qq.at.w = list())
  }else{
    post.samp <- list(w     = list(),
                      P.W_B = list(),
                      G.at.w   = list(),
                      Qq.at.w = list())
  }

  W_B    <- S_B
  for(i in 1:N.Qq){
    message(paste0("Sample W for threshold ",i,"/",N.Qq))
    samp.w <- lapply(fitted.mod$log.L[[i]],function(xx) sample_LGCP_fix(N.w,xx,fitted.mod,mesh,W_B))

    if(!inherits(S_B,"logical")){
      post.samp$w[[i]]     <- w <- lapply(samp.w,function(xx) xx$w)
      post.samp$P.W_B[[i]] <- lapply(samp.w,function(xx) xx$P.W_B)
    }else{
      post.samp$w[[i]] <- w <- samp.w
    }

    if(fitted.mod$options$excess.dist.fam=="E"){
      post.samp$g.at.w[[i]] <- list()
    }else{
      post.samp$G.at.w[[i]] <- list()
    }

    post.samp$Qq.at.w[[i]] <- list()
    for(j in 1:N.GW){
      A <- inla.mesh.projector(mesh=mesh,loc=w[[j]])$proj$A
      if(fitted.mod$options$excess.dist.fam=="E"){
        if(transf.G==F){
          post.samp$g.at.w[[i]][[j]] <- as.vector(A %*% fitted.mod$g[[i]][[j]])
        }else{
          post.samp$g.at.w[[i]][[j]] <- 1/as.vector(A %*% fitted.mod$G_T[[i]][[j]])
        }
      }else{
        post.samp$G.at.w[[i]][[j]] <- as.vector(A %*% fitted.mod$G[[i]][[j]])
      }

      post.samp$Qq.at.w[[i]][[j]]  <- as.vector(A %*% fitted.mod$Qq[,i])
    }
  }
  post.samp
}

#' Title
#'
#' @param w
#' @param x
#' @param y
#'
#' @return
#' @noRd
#'
#' @examples
get_rinfsup <- function(w,x,y){
  if(x[2]==Inf && y[2]==Inf){
    r.x <- x[1]/cos(w)
    r.y <- y[1]/sin(w)
    return(c(max(r.x,r.y),Inf))
  }
  r.x <- x/cos(w)
  yy <- r.x*sin(w)
  ind <- which(yy>=y[1] & yy <= y[2])
  if(length(ind)>0){
    r.y <- y/sin(w)
    xx <- r.y*cos(w)
    r_infsup <- c(r.x[ind],r.y[which(xx>=x[1] & xx <= x[2])])
    return(sort(r_infsup))
  }else{
    return(c(0,0))
  }
}

#' Title
#'
#' @param fitted.mod
#' @param post.sample
#' @param x
#' @param y
#'
#' @import INLA
#' @import inlabru
#' @import evd
#' @return
#' @noRd
#'
#' @examples
prob_estimation_2d <- function(fitted.mod,post.sample,x,y){

  excess.dist.fam <- fitted.mod$options$excess.dist.fam
  x <- sort(x); y <- sort(y)

  N.Qq <- length(fitted.mod$G)
  N.GW <- length(fitted.mod$G[[1]])
  X    <- fitted.mod$X

  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- atan2(y=X[,2], x=X[,1])

  r_infsup <- t(sapply(W,get_rinfsup,x=x,y=y))
  ind.W.B <- apply(r_infsup,1,function(xx){min(xx)>0})

  post.pred <- rep(0,N.Qq*N.GW)
  post.pred.correct <- rep(0,N.Qq*N.GW)
  cnt <- 0
  if(length(post.sample$P.W_B)==0){
    for(i in 1:N.Qq){
      message(paste0("Probabilities for threshold Q_q ",i,"/",N.Qq))

      A <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=W)$proj$A
      Qq  <- as.vector(A %*% fitted.mod$Qq[,i])

      ind.R.Bi <- R> r_infsup[,1] & R < Qq

      p.Bi <- sum(ind.R.Bi & ind.W.B)/dim(X)[1]

      xi <- fitted.mod$xi
      for(j in 1:N.GW){
        cnt <- cnt+1
        w     <- post.sample$w[[i]][[j]]
        Qq.at.w <- post.sample$Qq.at.w[[i]][[j]]

        n <- length(w)
        r <- t(sapply(w,get_rinfsup,x=x,y=y))

        if(excess.dist.fam=="E"){
          g.at.w   <- post.sample$g.at.w[[i]][[j]]
          pars <- cbind(w,Qq.at.w,r,g.at.w)#; pars <- pars[order(w),];
          pars <- pars[pars[,3]>0,]
          p_inf <- pexp(pars[,3]-pars[,2],rate=pars[,5])
          p_sup <- pexp(pars[,4]-pars[,2],rate=pars[,5])
        }else{
          G.at.w   <- post.sample$G.at.w[[i]][[j]]
          pars <- cbind(w,Qq.at.w,r,G.at.w)#; pars <- pars[order(w),];
          pars <- pars[pars[,3]>0,]
          p_inf <- pgpd(pars[,3]-pars[,2],0,scale=pars[,5],shape=xi[[i]])
          p_sup <- pgpd(pars[,4]-pars[,2],0,scale=pars[,5],shape=xi[[i]])
        }
        p <- sum(p_sup-p_inf)/n
        post.pred.correct[cnt] <- p*(1-fitted.mod$options$q) + p.Bi
        post.pred[cnt] <- p*(1-fitted.mod$options$q)
      }
    }
  }else{
    for(i in 1:N.Qq){
      message(paste0("Probabilities for threshold Q_q ",i,"/",N.Qq))

      A <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=W)$proj$A
      Qq  <- as.vector(A %*% fitted.mod$Qq[,i])

      ind.R.Bi <- R> r_infsup[,1] & R < Qq

      p.Bi <- sum(ind.R.Bi & ind.W.B)/dim(X)[1]

      xi <- fitted.mod$xi
      for(j in 1:N.GW){
        cnt <- cnt+1
        w     <- post.sample$w[[i]][[j]]
        Qq.at.w <- post.sample$Qq.at.w[[i]][[j]]
        P.W_B <- post.sample$P.W_B[[i]][[j]]

        r <- t(sapply(w,get_rinfsup,x=x,y=y))

        if(excess.dist.fam=="E"){
          g.at.w   <- post.sample$g.at.w[[i]][[j]]
          pars <- cbind(w,Qq.at.w,r,g.at.w)#; pars <- pars[order(w),];
          pars <- pars[pars[,3]>0,]
          p_inf <- pexp(pars[,3]-pars[,2],rate=pars[,5])
          p_sup <- pexp(pars[,4]-pars[,2],rate=pars[,5])
        }else{
          G.at.w   <- post.sample$G.at.w[[i]][[j]]
          pars <- cbind(w,Qq.at.w,r,G.at.w)#; pars <- pars[order(w),];
          pars <- pars[pars[,3]>0,]
          p_inf <- pgpd(pars[,3]-pars[,2],0,scale=pars[,5],shape=xi[[i]])
          p_sup <- pgpd(pars[,4]-pars[,2],0,scale=pars[,5],shape=xi[[i]])
        }
        p <- mean(p_sup-p_inf)#sum(p_inf-p_sup)/n
        post.pred.correct[cnt] <- p*P.W_B*(1-fitted.mod$options$q) + p.Bi
        post.pred[cnt] <- p*P.W_B*(1-fitted.mod$options$q)
      }
    }
  }
  cbind(post.pred,post.pred.correct)
}

#' Title
#'
#' @param fitted.mod
#' @param u
#' @param conditioning.marg
#' @param N.w
#' @param transf.G
#'
#' @import rmutil
#' @return
#' @noRd
#'
#' @examples
chi_posterior <- function(fitted.mod,u,conditioning.marg=1,N.w=5000,transf.G=F){

  N.Qq        <- fitted.mod$options$N.Qq
  N.GW <- fitted.mod$options$N.GW
  X            <- fitted.mod$X
  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- atan2(y=X[,2], x=X[,1])

  S_B <- c(0,pi/2)
  post.samp.marg.x <- sample_QGW_posterior(fitted.mod,N.w=N.w,S_B=S_B,transf.G=transf.G)

  # Get P(W in [0,pi/2]) for each i,j
  P.W_B.sq <- list()
  for(i in 1:fitted.mod$options$N.Qq){
    P.W_B.sq[[i]] <- list()
    for(j in 1:fitted.mod$options$N.GW){
      P.W_B.sq[[i]][[j]] <- get_P_W(fitted.mod$log.L[[i]][[j]],
                                    fitted.mod$mesh,limits.W_B=S_B)
    }
  }

  chis <- list() # to be returned

  for(k in 1:length(u)){
    message(paste0("Chi(u) estimation for u=",u[k],"."))
    x.sq <- c(u[k],Inf); y.sq <- c(u[k],Inf)
    r_infsup.sq <- t(sapply(W,get_rinfsup,x=x.sq,y=y.sq))
    ind.W.B.sq <- apply(r_infsup.sq,1,function(xx){min(xx)>0})

    if(conditioning.marg==1){
      x.marg <- c(u[k],Inf); y.marg <- c(-Inf,Inf)
    }else if(conditioning.marg==2){
      x.marg <- c(-Inf,Inf); y.marg <- c(u[k],Inf)
    }
    r_infsup.marg <- t(sapply(W,get_rinfsup,x=x.marg,y=y.marg))
    ind.W.B.marg <- apply(r_infsup.marg,1,function(xx){min(xx)>0})

    post.pred <- rep(0,N.Qq*N.GW)
    post.pred.correct <- rep(0,N.Qq*N.GW)
    post.pred.true.marg <- rep(0,N.Qq*N.GW)
    cnt <- 0

    for(i in 1:N.Qq){
      message(paste0("Chi(u) estimation for threshold Q_q ",i,"/",N.Qq))
      A <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=W)$proj$A
      Qq  <- as.vector(A %*% fitted.mod$Qq[,i])

      ind.R.Bi.sq <- R> r_infsup.sq[,1] & R < Qq
      p.Bi.sq <- sum(ind.R.Bi.sq & ind.W.B.sq)/dim(X)[1]

      ind.R.Bi.marg <- R> r_infsup.marg[,1] & R < Qq
      p.Bi.marg <- sum(ind.R.Bi.marg & ind.W.B.marg)/dim(X)[1]

      samp.p.Bi.sq <- rbinom(N.GW,fitted.mod$options$N.X,p.Bi.sq)/fitted.mod$options$N.X
      samp.p.Bi.marg <- rbinom(N.GW,fitted.mod$options$N.X,p.Bi.marg)/fitted.mod$options$N.X


      for(j in 1:N.GW){
        cnt <- cnt+1
        w     <- post.samp.marg.x$w[[i]][[j]]
        g.at.w   <- post.samp.marg.x$g[[i]][[j]]
        Qq.at.w <- post.samp.marg.x$Qq.at.w[[i]][[j]]
        P.W_B <- post.samp.marg.x$P.W_B[[i]][[j]]
        P.W_B.sq.ij <- P.W_B.sq[[i]][[j]]

        ind.sq <- w>0 & w<pi/2
        w.sq <- w[ind.sq]
        r.sq <- t(sapply(w.sq,get_rinfsup,x=x.sq,y=y.sq))
        pars <- cbind(w.sq,Qq.at.w[ind.sq],r.sq,g.at.w[ind.sq])#; pars <- pars[order(w.sq),];
        pars <- pars[pars[,3]>0,]
        p_inf.sq <- 1-pexp(pars[,3]-pars[,2],rate=pars[,5])
        p.sq <- mean(p_inf.sq)

        r <- t(sapply(w,get_rinfsup,x=x.marg,y=y.marg))
        pars <- cbind(w,Qq.at.w,r,g.at.w); pars <- pars[order(w),]; pars <- pars[pars[,3]>0,]
        p_inf <- 1-pexp(pars[,3]-pars[,2],rate=pars[,5])
        p.marg <- mean(p_inf)

        numer <- p.sq*P.W_B.sq.ij*(1-fitted.mod$options$q)
        denom <- p.marg*P.W_B*(1-fitted.mod$options$q)

        post.pred.correct[cnt] <- (numer + samp.p.Bi.sq[j])/(denom + samp.p.Bi.marg[j])
        post.pred[cnt] <- numer/denom
        post.pred.true.marg[cnt] <- (numer + p.Bi.sq)/(1-rmutil::plaplace(u[k]))
      }
    }
    chis[[k]] <- cbind(chis=post.pred,
                       chis.corrected = post.pred.correct,
                       chis.true.marg = post.pred.true.marg)
  }
  chis
}

#' Title
#'
#' @param fitted.mod
#' @param t
#' @param q.prime
#'
#' @return
#' @noRd
#'
#' @examples
return_set_2d <- function(fitted.mod,alpha=0.05,t){
  ret_set_list <- list(t=t)
  for(k in 1:length(t)){
    K <- log(t[k]*(1-fitted.mod$options$q))
    n.mesh <- length(fitted.mod$mesh$loc)
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
#' @param pts
#' @param w
#'
#' @return
#' @noRd
#'
#' @examples
numer_int <- function(pts,w){
  delta <- diff(w)
  summ <- 0
  for(i in 1:(length(pts)-1)){
    summ <- summ + delta[i]*(pts[i]+pts[i+1])/2
  }
  summ
}

#' Title
#'
#' @param fitted.mod
#'
#' @return
#' @noRd
#'
#' @examples
posterior_alpha_a <- function(fitted.mod){
  alphas <- matrix(NA,length(fitted.mod$G)*length(fitted.mod$G[[1]]),2)
  cnt <- 0
  for(i in 1:length(fitted.mod$g)){
    for(j in 1:length(fitted.mod$g[[1]])){
      cnt <- cnt+1
      xy <- pol2cart(cbind(fitted.mod$mesh$loc,fitted.mod$G[[i]][[j]]))
      ind_max_1 <- which.max(xy[,1]); ind_max_2 <- which.max(xy[,2])
      alphas[cnt,] <- c(xy[ind_max_1,2],xy[ind_max_2,1])
    }
  }
  alphas
}

#' Title
#'
#' @param fitted.mod
#'
#' @return
#' @noRd
#'
#' @examples
alpha_posterior_2d <- function(fitted.mod){
  alphas <- matrix(NA,length(fitted.mod$G)*length(fitted.mod$G[[1]]),2)
  cnt <- 0
  for(i in 1:length(fitted.mod$g)){
    for(j in 1:length(fitted.mod$g[[1]])){
      cnt <- cnt+1
      xy <- pol2cart(cbind(fitted.mod$mesh$loc,fitted.mod$G[[i]][[j]]))
      ind_max_1 <- which.max(xy[,1]); ind_max_2 <- which.max(xy[,2])

      wr_max1 <- cart2pol(xy[ind_max_1,])
      r_1 <- 1/cos(wr_max1[1])
      alpha_1 <- pol2cart(c(wr_max1[1],r_1))[2]

      wr_max2 <- cart2pol(xy[ind_max_2,])
      r_2 <- 1/sin(wr_max2[1])
      alpha_2 <- pol2cart(c(wr_max2[1],r_2))[1]
      alphas[cnt,] <- c(alpha_1,alpha_2)
    }
  }
  alphas
}

#' Title
#'
#' @param fitted.mod
#'
#' @return
#' @noRd
#'
#' @examples
dG_transformation <- function(fitted.mod){
  G <- fitted.mod$G
  mesh <- list()
  cnt <- 0
  for(i in 1:length(fitted.mod$g)){
    mesh[[i]] <- list()
    for(j in 1:length(fitted.mod$g[[1]])){
      cnt <- cnt+1
      xy <- pol2cart(cbind(fitted.mod$mesh$loc,fitted.mod$G[[i]][[j]]))
      ind_max_1 <- which.max(xy[,1]); ind_max_2 <- which.max(xy[,2])

      if(ind_max_1!=ind_max_2){
        pt1 <- xy[ind_max_1,]; pt2 <- xy[ind_max_2,]
        wr_max1 <- cart2pol(pt1)
        r_1 <- 1/cos(wr_max1[1])

        wr_max2 <- cart2pol(pt2)
        r_2 <- 1/sin(wr_max2[1])

        S <- cbind(pt1,pt2)
        M <- matrix(c(r_1/wr_max1[2],0,0,r_2/wr_max2[2]),2,2)

        Tr <- S%*%M%*%solve(S)

        G.ij <- cart2pol(t(Tr %*% t(xy)))
        G.ij <- G.ij[order(G.ij[,1]),]

        G.ij <- rbind(c(-pi,mean(G.ij[c(1,nrow(G.ij)),2])),
                      G.ij,
                      c(pi,mean(G.ij[c(1,nrow(G.ij)),2])))

        mesh <- inla.mesh.1d(G.ij[,1],
                             boundary="cyclic",degree=1)
        A.G.ij <- inla.mesh.projector(mesh=mesh,loc=fitted.mod$mesh$loc)$proj$A
        G.ij.mesh <- as.vector(A.G.ij %*% G.ij[-nrow(G.ij),2])
        G[[i]][[j]] <- G.ij.mesh
      }else{
        pt1 <- xy[ind_max_1,]; pt2 <- xy[ind_max_2,]
        wr_max1 <- cart2pol(pt1)
        r_1 <- 1/cos(wr_max1[1])

        wr_max2 <- cart2pol(pt2)
        r_2 <- 1/sin(wr_max2[1])

        # Find vector orthogonal to unique eigenvector since pt1=pt2
        find_vec <- function(x){pt1%*%c(1,x)}

        pt2 <- c(1,uniroot(find_vec,c(-20,20))$root)
        pt2 <- pt2/sqrt(sum(pt2^2))

        S <- cbind(pt1,pt2)
        M <- matrix(c(min(r_1/wr_max1[2],r_2/wr_max2[2]),0,0,1),2,2)

        Tr <- S%*%M%*%solve(S)

        G.ij <- cart2pol(t(Tr %*% t(xy)))
        G.ij <- G.ij[order(G.ij[,1]),]

        # Project observed points of $ back onto the mesh
        G.ij <- rbind(c(-pi,mean(G.ij[c(1,nrow(G.ij)),2])),
                      G.ij,
                      c(pi,mean(G.ij[c(1,nrow(G.ij)),2])))

        mesh <- inla.mesh.1d(G.ij[,1],
                             boundary="cyclic",degree=1)
        A.G.ij <- inla.mesh.projector(mesh=mesh,loc=fitted.mod$mesh$loc)$proj$A
        G.ij.mesh <- as.vector(A.G.ij %*% G.ij[-nrow(G.ij),2])
        G[[i]][[j]] <- G.ij.mesh
      }
    }
  }
  G
}

#' Title
#'
#' @param fitted.mod
#'
#' @return
#' @noRd
#'
#' @examples
eta_posterior_2d <- function(fitted.mod){
  etas <- rep(NA,length(fitted.mod$G)*length(fitted.mod$G[[1]]))
  A <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=pi/4)$proj$A
  cnt <- 1
  for(i in 1:length(fitted.mod$G)){
    for(j in 1:length(fitted.mod$G[[1]])){
      r <- as.vector(A %*%fitted.mod$G[[i]][[j]])
      etas[cnt] <- r*cos(pi/4)
      cnt <- cnt + 1
    }
  }
  etas
}

#' Title
#'
#' @param fitted.Qq
#' @param cex.pts
#' @param cex.axis
#' @param xlim
#' @param ylim
#' @param by
#' @param alpha
#'
#' @import excursions
#' @import grDevices
#' @import graphics
#' @return
#' @noRd
#'
#' @examples
plot_Qq_2d <- function(fitted.Qq,alpha,cex.pts=0.4,cex.axis=1.4,xlim=c(0,0),ylim=c(0,0),by=2){

  mean.Qq <- apply(fitted.Qq$Qq,1,mean,na.rm=T)
  meann <- pol2cart(cbind(fitted.Qq$mesh$loc,mean.Qq))

  if(sum(xlim==c(0,0) & ylim==c(0,0))==2){
    plot(fitted.Qq$X,col="grey35",
         pch=16,cex=cex.pts,cex.lab=1.6,cex.axis=cex.axis,bty="n",xlab="",ylab="")
  }else{
    plot(fitted.Qq$X,col="grey35",xaxt="n",yaxt="n",xlim=xlim,ylim=ylim,
         pch=16,cex=cex.pts,cex.lab=1.6,cex.axis=cex.axis,bty="n",xlab="",ylab="")
    axis(1, at = seq(xlim[1], xlim[2], by = by),cex.axis=cex.axis)
    axis(2, at = seq(ylim[1], ylim[2], by = by),cex.axis=cex.axis)
  }

  if(ncol(fitted.Qq$Qq)==1){
    col <- col2rgb("grey30")/255
    polygon(meann,col=rgb(1,1,1),border=rgb(col[1], col[2], col[3], alpha = 0.3))
  }else{
    excurs <- simconf.mc(samples = fitted.Qq$Qq,alpha = 1-alpha)
    low <- pol2cart(cbind(fitted.Qq$mesh$loc,excurs$a))
    upp <- pol2cart(cbind(fitted.Qq$mesh$loc,excurs$b))

    col <- col2rgb("grey30")/255
    polygon(upp,col=rgb(col[1], col[2], col[3], alpha = 0.3),
            border=rgb(col[1], col[2], col[3], alpha = 0.3))
    polygon(low,col=rgb(1,1,1),border=rgb(col[1], col[2], col[3], alpha = 0.3))
    lines(meann)
  }

  X <- fitted.Qq$X
  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- atan2(y=X[,2], x=X[,1])
  ind.exc <- c()
  for(i in 1:ncol(fitted.Qq$Qq)){
    A.Qq <- inla.mesh.projector(mesh=fitted.Qq$mesh,loc=W)$proj$A
    R.W <- as.vector(A.Qq %*%fitted.Qq$Qq[,i])
    ind.exc <- c(ind.exc,which(R>R.W))
  }
  ind.exc <- unique(ind.exc)
  points(fitted.Qq$X[ind.exc,],col="grey35",pch=16,cex=cex.pts)
}

#' Title
#'
#' @param fitted.mod
#' @param mean_med
#' @param cex.txt
#' @param cex.pts
#' @param transf.G
#' @param alpha
#' @param cex.axis
#' @param by
#'
#' @import excursions
#' @import grDevices
#' @import graphics
#' @return
#' @noRd
#'
#' @examples
plot_G_2d <- function(fitted.mod,alpha,mean_med="mean",cex.txt=1.4,cex.pts=0.6,transf.G = F,cex.axis=1.4,by=2){
  if(transf.G==T){
    G <- fitted.mod$G_T
  }else{
    G <- fitted.mod$G
  }
  n.samp <- length(G)*length(G[[1]])
  all.est.G <- matrix(NA,nrow=length(fitted.mod$mesh$loc),ncol=n.samp)
  count <- 0
  for(i in 1:length(G)){
    est.G <- G[[i]]
    for(j in 1:length(G[[1]])){
      count <- count +1
      # lines(pol2cart(cbind(fitted.mod$mesh$loc,est.G[[j]])),col="grey50")
      all.est.G[,count] <- est.G[[j]]
    }
  }
  emp.quants <- apply(all.est.G,1,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
  mean.G <- apply(all.est.G,1,mean,na.rm=T)

  excurs <- simconf.mc(samples = all.est.G,alpha = 1-alpha)
  low <- pol2cart(cbind(fitted.mod$mesh$loc,excurs$a))
  upp <- pol2cart(cbind(fitted.mod$mesh$loc,excurs$b))

  if(fitted.mod$options$excess.dist.fam=="E"){
    plot(fitted.mod$X/log(nrow(fitted.mod$X)/2),col="grey35",
         xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),pch=16,cex=cex.pts,
         cex.lab=1.6,cex.axis=1.2,xaxt='n',yaxt='n',bty="n",xlab="",ylab="")

    col <- col2rgb("grey30")/255
    polygon(upp,col=rgb(col[1], col[2], col[3], alpha = 0.3),
            border=rgb(col[1], col[2], col[3], alpha = 0.3))
    polygon(low,col=rgb(1,1,1),border=rgb(col[1], col[2], col[3], alpha = 0.3))
    lines(c(-1.2,1.2),c(0,0),col="black")
    lines(c(0,0),c(-1.2,1.2),col="black")
    text(c(-1.30,1.30,0,0),c(0,0,-1.30,1.30),cex=cex.txt,
         c(expression(pi),0,expression(3*pi/2),expression(pi/2)),col="black")
    points(fitted.mod$X/log(nrow(fitted.mod$X)/2),col="grey35",pch=16,cex=cex.pts)
    if(mean_med=="mean"){
      lines(pol2cart(cbind(fitted.mod$mesh$loc,mean.G)))
    }else{
      lines(pol2cart(cbind(fitted.mod$mesh$loc,emp.quants[2,])))
    }

    segments(-1,1,1,1,lty="dotted",col="grey50")
    segments(-1,1,-1,-1,lty="dotted",col="grey50")
    segments(-1,-1,1,-1,lty="dotted",col="grey50")
    segments(1,-1,1,1,lty="dotted",col="grey50")
  }else if(fitted.mod$options$excess.dist.fam=="GP"){

    if(mean_med=="mean"){
      mean_med <- pol2cart(cbind(fitted.mod$mesh$loc,mean.G))
    }else{
      mean_med <- pol2cart(cbind(fitted.mod$mesh$loc,emp.quants[2,]))
    }

    lims <- c(floor(min(mean_med)),ceiling(max(mean_med)))*1.1

    plot(mean_med,col="grey35",xaxt="n",yaxt="n",xlim=lims,ylim=lims,
         pch=16,cex=cex.pts,cex.lab=1.6,cex.axis=cex.axis,bty="n",xlab="",ylab="")
    axis(1, at = seq(lims[1], lims[2], by = by),cex.axis=cex.axis)
    axis(2, at = seq(lims[1], lims[2], by = by),cex.axis=cex.axis)

    col <- col2rgb("grey30")/255
    polygon(upp,col=rgb(col[1], col[2], col[3], alpha = 0.3),
            border=rgb(col[1], col[2], col[3], alpha = 0.3))
    polygon(low,col=rgb(1,1,1),border=rgb(col[1], col[2], col[3], alpha = 0.3))
  }
}

#' Title
#'
#' @param fitted.mod
#' @param xylim
#' @param main
#' @param mid.gap
#' @param txt.gap
#' @param cex.txt
#' @param exconly
#' @param alpha
#'
#' @import pracma
#' @import excursions
#' @import grDevices
#' @import graphics
#' @return
#' @noRd
#'
#' @examples
plot_W_2d <- function(fitted.mod,alpha,xylim,main="",mid.gap=0.1, txt.gap=0.02,cex.txt=1.4,exconly=F){

  X <- fitted.mod$X
  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- atan2(y=X[,2], x=X[,1])

  if(exconly==T){
    A.Qq <- inla.mesh.projector(mesh=fitted.mod$mesh,loc=W)$proj$A
    R.W <- as.vector(A.Qq %*% fitted.mod$Qq.mean$mean)
    ind.exc <- R-R.W>0
    X <- X[ind.exc,]
    W <- atan2(y=X[,2], x=X[,1])
  }

  h <- hist(W,breaks=seq(-pi,pi,length.out=60),plot = F)

  mesh <- fitted.mod$mesh$loc
  m <- matrix(NA,length(fitted.mod$G)*length(fitted.mod$G[[1]]),length(mesh))
  cnt <- 1
  for(i in 1:length(fitted.mod$G)){
    for(j in 1:length(fitted.mod$G[[1]])){
      f_W <- exp(fitted.mod$log.L[[i]][[j]])
      intt <- numer_int(f_W,mesh)
      m[cnt,] <- f_W/intt
      cnt <- cnt + 1
    }
  }

  # low <- pol2cart(cbind(mesh,apply(m,2,quantile,probs=c(0.025))+mid.gap))
  # up <- pol2cart(cbind(mesh,apply(m,2,quantile,probs=c(0.975))+mid.gap))

  excurs <- simconf.mc(samples = t(m),alpha = 1-alpha)
  low <- pol2cart(cbind(mesh,excurs$a+mid.gap))
  up <- pol2cart(cbind(mesh,excurs$b+mid.gap))

  if(xylim==0){
    xylim <- max(c(excurs$b,h$density))
  }

  lines.gap <- 0.15
  plot(NA,NA,xlim=c(-xylim-0.05-mid.gap-lines.gap,xylim+0.05+mid.gap+lines.gap),
       ylim=c(-xylim-0.05-mid.gap-lines.gap,xylim+0.05+mid.gap+lines.gap),
       xaxt='n',yaxt='n',bty="n",xlab="",ylab="")
  title(main,line=-1)

  col <- col2rgb("grey30")/255
  polygon(up,col=rgb(col[1], col[2], col[3], alpha = 0.3),
          border=rgb(col[1], col[2], col[3], alpha = 0.3))
  polygon(low,col=rgb(1,1,1),border=rgb(col[1], col[2], col[3], alpha = 0.3))

  ws <- c(seq(-pi,pi,by=0.001),-pi)
  for(i in 1:ceiling(xylim*10)){
    if(i==1){
      col <- "black"
    }else{
      col<-"lightgrey"
    }
    lines(i*cos(ws)/10,i*sin(ws)/10,col=col,lwd=0.7)
  }

  lines(c(-xylim-lines.gap,xylim+lines.gap),c(0,0),col="black")
  lines(c(0,0),c(-xylim-lines.gap,xylim+lines.gap),col="black")
  text(c(-xylim-txt.gap-lines.gap,xylim+txt.gap+lines.gap,0,0),
       c(0,0,-xylim-0.005-txt.gap-lines.gap,xylim+0.005+txt.gap+lines.gap),
       cex=cex.txt,
       c(expression(pi),0,expression(3*pi/2),expression(pi/2)),col="black")

  for(i in 1:length(h$density)){
    pt1 <- pol2cart(c(h$breaks[i],h$density[i]+mid.gap))
    pt2 <- pol2cart(c(h$breaks[i+1],h$density[i]+mid.gap))
    strt1 <- pol2cart(c(h$breaks[i],mid.gap))
    strt2 <- pol2cart(c(h$breaks[i+1],mid.gap))
    segments(strt1[1],strt1[2],pt1[1],pt1[2],col="grey40",lwd=1)
    segments(strt2[1],strt2[2],pt2[1],pt2[2],col="grey40",lwd=1)
    segments(pt2[1],pt2[2],pt1[1],pt1[2],col="grey40",lwd=1)
  }
  post_mean <- pol2cart(cbind(mesh,apply(m,2,mean)+mid.gap))
  lines(post_mean,lwd=1.2)
}

#' Title
#'
#' @param fitted.mod
#' @param t
#' @param list_ret_sets
#' @param cex.axis
#' @param xlim
#' @param ylim
#' @param by
#'
#' @import grDevices
#' @import graphics
#' @return
#' @noRd
#'
#' @examples
plot_return_bdry_2d <- function(fitted.mod,list_ret_sets,cex.pts=0.4,cex.axis=1.4,xlim=c(0,0),ylim=c(0,0),by=2){

  t <- rev(sort(list_ret_sets$t))

  if(sum(xlim==c(0,0) & ylim==c(0,0))==2){
    plot(fitted.mod$X,col="grey35",
         pch=16,cex=cex.pts,cex.lab=1.6,cex.axis=cex.axis,bty="n",xlab="",ylab="")
  }else{
    plot(fitted.mod$X,col="grey35",xaxt="n",yaxt="n",xlim=xlim,ylim=ylim,
         pch=16,cex=cex.pts,cex.lab=1.6,cex.axis=cex.axis,bty="n",xlab="",ylab="")
    axis(1, at = seq(xlim[1], xlim[2], by = by),cex.axis=cex.axis)
    axis(2, at = seq(ylim[1], ylim[2], by = by),cex.axis=cex.axis)
  }

  alpha.col <- 0.4
  cols <- rev(seq(col2rgb("grey20")[1],col2rgb("grey80")[1],length.out=length(t)+1)/255)
  for(i in length(t):1){
    post.ret.set <- list_ret_sets[[i+1]]#return_set_2d(fitted.mod,t=t[i])
    low <- pol2cart(cbind(fitted.mod$mesh$loc,post.ret.set$lower))
    up <- pol2cart(cbind(fitted.mod$mesh$loc,post.ret.set$upper))
    polygon(up,col=rgb(cols[i], cols[i], cols[i], alpha = alpha.col),
            border=rgb(cols[i], cols[i], cols[i], alpha = alpha.col))
    polygon(low,col=rgb(1,1,1),border=rgb(cols[i], cols[i], cols[i], alpha = alpha.col))

  }
  points(fitted.mod$X,col="grey35",pch=16,cex=cex.pts)

  if(ncol(fitted.mod$Qq)==1){
    mean_Qq <- apply(fitted.mod$Qq,1,mean,na.rm=T)
    mean_Qq<- pol2cart(cbind(fitted.mod$mesh$loc,mean_Qq))
    col <- cols[length(cols)]
    polygon(mean_Qq,col=rgb(1,1,1),border=rgb(col, col, col, alpha = alpha.col))
  }else{
    excurs <- simconf.mc(samples = fitted.mod$Qq,alpha = 0.95)#,u=0,type = "=")
    low <- pol2cart(cbind(fitted.mod$mesh$loc,excurs$a))
    up <- pol2cart(cbind(fitted.mod$mesh$loc,excurs$b))

    col <- cols[length(cols)]
    polygon(up,col=rgb(col, col, col, alpha = alpha.col),border=rgb(col, col, col, alpha = alpha.col))
    polygon(low,col=rgb(1,1,1),border=rgb(col, col, col, alpha = alpha.col))
  }
}

#' Title
#'
#' @param fitted.mod
#' @param list_ret_sets
#' @param cex.pts
#' @param cex.axis
#' @param xlim
#' @param ylim
#' @param by
#'
#' @import grDevices
#' @import graphics
#' @return
#' @noRd
#'
#' @examples
plot_return_sets_2d <- function(fitted.mod,list_ret_sets,cex.pts=0.4,cex.axis=1.4,xlim=c(0,0),ylim=c(0,0),by=2){

  t <- rev(sort(list_ret_sets$t))

  if(sum(xlim==c(0,0) & ylim==c(0,0))==2){
    plot(fitted.mod$X,col="grey35",
         pch=16,cex=cex.pts,cex.lab=1.6,cex.axis=cex.axis,bty="n",xlab="",ylab="")
  }else{
    plot(fitted.mod$X,col="grey35",xaxt="n",yaxt="n",xlim=xlim,ylim=ylim,
         pch=16,cex=cex.pts,cex.lab=1.6,cex.axis=cex.axis,bty="n",xlab="",ylab="")
    axis(1, at = seq(xlim[1], xlim[2], by = by),cex.axis=cex.axis)
    axis(2, at = seq(ylim[1], ylim[2], by = by),cex.axis=cex.axis)
  }

  alpha.col <- 0.4
  cols <- rev(seq(col2rgb("grey20")[1],col2rgb("grey80")[1],length.out=length(t)+1)/255)

  lim <- 40
  big.rect <- rbind(c(-lim,-lim),c(lim,-lim),c(lim,lim),c(-lim,lim))
  col <- cols[length(cols)]
  polygon(big.rect,col=rgb(cols[1], cols[1], cols[1], alpha = alpha.col),
          border=rgb(col, col, col, alpha = alpha.col))

  for(i in length(t):1){
    post.ret.set <- list_ret_sets[[i+1]]#return_set_2d(fitted.mod,t=t[i])
    mean_post <- post.ret.set$mean
    mean_post <- pol2cart(cbind(fitted.mod$mesh$loc,mean_post))

    polygon(mean_post,col="white",border="white")
    col <- cols[length(cols)-i+1]
    polygon(mean_post,col=rgb(col, col, col, alpha = alpha.col),
            border=rgb(col, col, col, alpha = alpha.col))
  }
  points(fitted.mod$X,col="grey35",pch=16,cex=cex.pts)

  mean_Qq <- apply(fitted.mod$Qq,1,mean,na.rm=T)
  mean_Qq<- pol2cart(cbind(fitted.mod$mesh$loc,mean_Qq))

  col <- cols[length(cols)]
  # polygon(mean_Qq,col=rgb(col, col, col, alpha = alpha.col),border=rgb(col, col, col, alpha = alpha.col))
  polygon(mean_Qq,col="white",border=col)
  # lines(mean_Qq)
}

#' Title
#'
#' @param fitted.mod
#' @param n.samp
#'
#' @import stats
#' @return
#' @noRd
#'
#' @examples
pp_R_w_2d <- function(fitted.mod,n.samp){

  options <- fitted.mod$options
  mesh <- fitted.mod$mesh
  X <- fitted.mod$X
  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- atan2(y=X[,2], x=X[,1])
  dfRW <- data.frame(R=R, W=W)

  min.Qq <- apply(fitted.mod$Qq,1,min)
  A.Qq <- inla.mesh.projector(mesh=mesh,loc=dfRW$W)$proj$A
  dfRW$thres.gamma      <- as.vector(A.Qq %*%min.Qq)
  dfRW$is.excess        <- ifelse(dfRW$R - dfRW$thres.gamma > 0 , 1, 0)
  dfRW.exc <- dfRW[dfRW$is.excess==1,]
  A.Qq.exc <- inla.mesh.projector(mesh=mesh,loc=dfRW.exc$W)$proj$A

  Qq.at.w <- matrix(NA,nrow(dfRW.exc),options$N.Qq)
  g.at.w   <- list()
  for(i in 1:options$N.Qq){
    Qq.at.w[,i] <- as.vector(A.Qq.exc %*%fitted.mod$Qq[,i])
    g.at.w[[i]] <- matrix(NA,nrow(dfRW.exc),options$N.GW)
    for(j in 1:options$N.GW){
      g.at.w[[i]][,j] <- as.vector(A.Qq.exc %*%fitted.mod$g[[i]][[j]])
    }
  }

  ps <- rep(NA,nrow(dfRW.exc))
  samp.p.p <- matrix(NA,options$N.GW*n.samp,options$N.Qq)
  for(k in 1:nrow(dfRW.exc)){
    for(i in 1:options$N.Qq){
      samp.p.p[,i] <- Qq.at.w[k,i] + c(sapply(g.at.w[[i]][k,],function(xx) rexp(n.samp,xx)))
    }
    F.R_w <- stats::ecdf(c(samp.p.p))
    ps[k] <- F.R_w(dfRW.exc$R[k])
  }
  ps
}

#' Title
#'
#' @param fitted.mod
#' @param n.samp
#' @param thresh.num
#'
#' @import stats
#' @return
#' @noRd
#'
#' @examples
pp_W_2d <- function(fitted.mod,n.samp,thresh.num=NA){
  mesh <- fitted.mod$mesh
  X <- fitted.mod$X
  R <- apply(X, 1, function(x) sqrt(sum(x^2)))
  W <- atan2(y=X[,2], x=X[,1])
  dfRW <- data.frame(R=R, W=W)

  if(is.na(thresh.num)){

    min.Qq <- apply(fitted.mod$Qq,1,min)
    A.Qq <- inla.mesh.projector(mesh=mesh,loc=dfRW$W)$proj$A
    dfRW$thres.gamma      <- as.vector(A.Qq %*%min.Qq)
    dfRW$is.excess        <- ifelse(dfRW$R - dfRW$thres.gamma > 0 , 1, 0)
    dfRW.exc <- dfRW[dfRW$is.excess==1,]

    ws <- unlist(lapply(fitted.mod$log.L,
                        function(yy) lapply(yy,
                                            function(xx) sample_LGCP_fix(n.samp,xx,fitted.mod,mesh))))
    F.W <- stats::ecdf(ws)

  }else{
    A.Qq <- inla.mesh.projector(mesh=mesh,loc=dfRW$W)$proj$A
    dfRW$thres.gamma      <- as.vector(A.Qq %*%fitted.mod$Qq[,thresh.num])
    dfRW$is.excess        <- ifelse(dfRW$R - dfRW$thres.gamma > 0 , 1, 0)
    dfRW.exc <- dfRW[dfRW$is.excess==1,]

    ws <- unlist(lapply(fitted.mod$log.L[[thresh.num]],function(xx) sample_LGCP_fix(n.samp,xx,fitted.mod,mesh)))
    F.W <- stats::ecdf(ws)
  }
  F.W(dfRW.exc$W)
}
