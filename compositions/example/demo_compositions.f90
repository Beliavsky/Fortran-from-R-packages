program demo_compositions
  use compositions
  implicit none
  real(dp), allocatable :: x(:,:),cov(:,:),center(:)
  type(fit_dirichlet_result) :: fit
  call rng_seed(2026)
  x=rdirichlet(2000,[2.0_dp,3.0_dp,5.0_dp])
  fit=fit_dirichlet(x,max_iter=200,tol=1e-8_dp)
  call compositional_covariance(x,cov,center)
  print '(a,3f10.4)', 'Dirichlet alpha: ',fit%alpha
  print '(a,3f10.4)', 'Aitchison center: ',center
  cov=variation_matrix(x)
  print '(a,f10.4)', 'variation(1,2): ',cov(1,2)
end program
