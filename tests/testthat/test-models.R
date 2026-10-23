test_that('model names follow the manuscript', {
  expect_identical(directional_architecture('M1'),'decoupled')
  expect_identical(directional_architecture('M2'),'shared-geometry')
  expect_identical(directional_architecture('M3'),'residual-adjustment')
  expect_error(directional_architecture('M4'),'W.model')
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

test_that('angular rejection sampling handles an interval without internal knots', {
  m <- fmesher::fm_mesh_1d(seq(-pi,pi,length.out=17),boundary='cyclic')
  set.seed(19)
  z <- sample_LGCP_fix(100,rep(-3,m$n),list(mesh=m),m,c(.01,.02))
  expect_length(z$w,100)
  expect_true(all(z$w>=.01 & z$w<=.02))
  expect_equal(z$P.W_B,.01/(2*pi),tolerance=1e-10)
})
