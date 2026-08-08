context("Tests using binary dataset")

test_that("Binary", {
  testthat::skip_on_cran()
  data(binary)
  N <- length(binary$Y)
  k <- NROW(binary$X)
  nvars <- as.integer(N*k + k)
  start <- rnorm(nvars) ## random starting values
  priors <- list(inv.Sigma = rWishart(1,k+5,10 * diag(k))[,,1],
                 inv.Omega = 10 * diag(k))

  m <- list(list(hs=binary.hess, method="Sparse", precond=0),
            list(hs=NULL, method="BFGS", precond=0),
            list(hs=NULL, method="SR1", precond=0),
            list(hs=binary.hess, method="Sparse", precond=1)
            )

  for (meth in m) {

    if (!(Sys.info()[['sysname']] == 'sunos' & meth$method %in% c('BFGS', 'SR1'))) {


    opt0 <- trust.optim(start,
                        fn=binary.f,
                        gr=binary.grad,
                        hs=meth$hs,
                        method=meth$method,
                        control=list(
                          start.trust.radius=5,
                          stop.trust.radius = 1e-9,
                          prec=1e-5,
                          cg.tol=1e-5,
                          report.precision=1L,
                          preconditioner=meth$precond,
                          report.freq=20L,
                          report.header.freq=10,
                          maxit=5000L,
                          function.scale.factor= -1,
                          report.level=0
                        ),
                        data=binary, priors=priors
                        )

    norm_gr <- sqrt(sum(opt0$gradient ^ 2))

    expect_equal(norm_gr, 0,  tolerance=.001)
    expect_match(opt0$status, "Success")
    expect_match(opt0$method, meth$method)

    }

  }

  ## opt <- trust.optim(start, fn=binary.f,
  ##                    gr = binary.grad,
  ##                    hs = binary.hess,
  ##                    method = "Sparse",
  ##                    control = list(
  ##                      start.trust.radius=5,
  ##                      stop.trust.radius = 1e-7,
  ##                      prec=1e-7,
  ##                      report.precision=1L,
  ##                      maxit=500L,
  ##                      preconditioner=1L,
  ##                      function.scale.factor=-1
  ##                  ),
  ##                  data=binary, priors=priors
  ##                  )

  })
