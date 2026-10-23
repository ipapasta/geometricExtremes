test_that('model names follow the manuscript', {
  expect_identical(directional_architecture('M1'),'decoupled')
  expect_identical(directional_architecture('M2'),'shared-geometry')
  expect_identical(directional_architecture('M3'),'residual-adjustment')
  expect_error(directional_architecture('M4'),'W.model')
})

test_that('experimental EGP helpers support negative xi stably', {
  expect_true(is_gp_excess_family('GP'))
  expect_true(is_gp_excess_family('EGP'))
  expect_false(is_gp_excess_family('E'))

  alpha <- .5
  q <- c(.2,1,3)
  at_zero <- gp_quantile_to_scale(q,0,alpha)
  expect_equal(at_zero,q/-log1p(-alpha),tolerance=1e-13)
  for(xi in c(-1e-10,1e-10,-.2,.2)) {
    direct <- q*xi/expm1(-log1p(-alpha)*xi)
    expect_equal(gp_quantile_to_scale(q,xi,alpha),direct,tolerance=1e-10)
    expect_true(all(gp_quantile_to_scale(q,xi,alpha)>0))
  }

  mock <- list(
    summary.hyperpar=rbind(
      'Tail parameter for the gp observations[2]'=c(mean=.1,sd=.01),
      'Range for zeta'=c(mean=1,sd=.1)),
    marginals.hyperpar=list(
      'Tail parameter for the gp observations[2]'=cbind(c(0,.2),c(1,1))))
  expect_equal(extract_gp_tail(mock)$mean,.1)
  rownames(mock$summary.hyperpar)[1] <- 'Tail parameter for egp observations[2]'
  names(mock$marginals.hyperpar)[1] <- rownames(mock$summary.hyperpar)[1]
  expect_equal(extract_gp_tail(mock)$mean,.1)

  ctl <- egp_kappa1_control(.5,c(-.5,.5),5)
  expect_identical(ctl$hyper$shape$fixed,TRUE)
  expect_equal(ctl$hyper$shape$initial,0)
  expect_equal(ctl$hyper$tail$param,c(5,-.5,.5))
  expect_error(gp_quantile_to_scale(1,0,1),'alpha')
})

test_that('experimental EGP options are explicit and validated', {
  X <- matrix(rnorm(40),ncol=2)
  args <- list(X=X,excess.dist.fam='EGP',W.model='M1',W.data='ExcOnly',
    use.mean.Qq=TRUE,q=.8,alpha=.5,N.Qq=1,N.GW=1,
    QR.prior.range=c(1,.5),QR.prior.sigma=c(1,.5),
    zeta.prior.range=c(1,.5),zeta.prior.sigma=c(1,.5),
    phi.prior.range=c(1,.5),phi.prior.sigma=c(1,.5),
    mesh.knots.2d=seq(-pi,pi,length.out=9),mesh.res.3d=2)
  expect_warning(op <- do.call(set.options,args),'Experimental EGP')
  expect_identical(op$excess.dist.fam,'EGP')
  expect_equal(op$egp.xi.interval,c(-.5,.5))
  args$egp.xi.interval <- c(0,.5)
  expect_error(suppressWarnings(do.call(set.options,args)),'straddle zero')
})

test_that('native constraints use surface weights including the cyclic seam', {
  meshes <- list(
    fmesher::fm_mesh_1d(seq(-pi,pi,length.out=33),boundary='cyclic'),
    fmesher::fm_mesh_1d(c(-pi,-2.8,-1.3,-.2,.1,1.7,pi),boundary='cyclic'),
    fmesher::fm_rcdt_2d_inla(globe=2,crs=fmesher::fm_CRS('sphere')))
  for(mesh in meshes) {
    b <- as.numeric(Matrix::colSums(fmesher::fm_fem(mesh,order=1)$c0))
    spde <- INLA::inla.spde2.pcmatern(mesh,prior.range=c(1,.5),
                                     prior.sigma=c(1,.5),constr=TRUE)
    A <- as.matrix(spde$f$extraconstr$A)
    expect_equal(nrow(A),1L)
    expect_equal(as.numeric(A)/sum(A),b/sum(b),tolerance=1e-10)
    expect_true(all(b>0))
    expect_true(all(spde$f$extraconstr$e==0))
  }
  b <- as.numeric(Matrix::colSums(fmesher::fm_fem(meshes[[1]],order=1)$c0))
  expect_equal(sum(b),2*pi,tolerance=1e-10)
  expect_equal(b,rep(2*pi/32,32),tolerance=1e-10)
})

test_that('circle normalisation is exact, periodic, and shift-invariant', {
  m <- fmesher::fm_mesh_1d(c(-pi,-2.8,-1.3,-.2,.1,1.7,pi),boundary='cyclic')
  ell <- sin(m$loc) + .3*cos(2*m$loc)
  expect_equal(circle_log_integral(rep(0,m$n),m),log(2*pi),tolerance=1e-12)
  expect_equal(get_P_W(rep(0,m$n),m,c(-pi,0)),.5,tolerance=1e-12)
  expect_equal(get_P_W(ell,m,c(-.1,.05)),get_P_W(ell+1000,m,c(-.1,.05)),tolerance=1e-10)
  expect_equal(get_P_W(ell,m,c(-pi,pi)),1)
  expect_equal(get_P_W(ell,m,c(0,0)),0)
  Z <- circle_log_integral(ell,m)
  # Independent adaptive quadrature, split at every knot.
  x <- c(m$loc,pi)
  value <- sum(vapply(seq_len(length(x)-1),function(j)
    integrate(function(w) exp(as.vector(fmesher::fm_basis(m,loc=w)%*%ell)-Z),
      x[j],x[j+1],rel.tol=1e-10)$value,numeric(1)))
  expect_equal(value,1,tolerance=1e-9)
  expect_error(get_P_W(ell,m,c(1,-1)),'Invalid')
})

test_that('E and GP scale-set plots use the same geometric layout', {
  make_fit <- function(family) {
    mesh <- fmesher::fm_mesh_1d(seq(-pi,pi,length.out=17),boundary='cyclic')
    set.seed(11)
    G <- replicate(20,exp(.05*sin(mesh$loc)+rnorm(mesh$n,0,.01)),simplify=FALSE)
    list(X=matrix(rnorm(400),ncol=2),mesh=mesh,G=list(G),
         options=list(excess.dist.fam=family))
  }
  tf <- tempfile(fileext='.pdf'); grDevices::pdf(tf)
  on.exit({grDevices::dev.off();unlink(tf)},add=TRUE)
  e <- plot_G(make_fit('E'),conf='marg')
  gp <- plot_G(make_fit('EGP'),conf='marg')
  expect_equal(e$data.scale,gp$data.scale)
  expect_equal(e$xlim,gp$xlim)
  expect_equal(e$ylim,gp$ylim)
  expect_true(e$show.reference.box && gp$show.reference.box)
})

test_that('angular rejection sampling handles an interval without internal knots', {
  m <- fmesher::fm_mesh_1d(seq(-pi,pi,length.out=17),boundary='cyclic')
  set.seed(19)
  z <- sample_LGCP_fix(100,rep(-3,m$n),list(mesh=m),m,c(.01,.02))
  expect_length(z$w,100)
  expect_true(all(z$w>=.01 & z$w<=.02))
  expect_equal(z$P.W_B,.01/(2*pi),tolerance=1e-10)
})
