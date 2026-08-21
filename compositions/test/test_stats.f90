program test_stats
  use compositions
  implicit none
  real(dp), allocatable :: x(:,:),cov(:,:),center(:),rcov(:,:),rcenter(:),y(:,:),pred(:,:)
  real(dp) :: alpha(3),design(4000,2),z(4000,2),btrue(2,2)
  type(fit_dirichlet_result) :: fr
  type(pca_result) :: pr
  type(compositional_lm_result) :: lm
  type(normal_location_result) :: nt
  integer :: i
  call rng_seed(2468)
  alpha=[2.0_dp,3.0_dp,5.0_dp]; x=rdirichlet(5000,alpha)
  fr=fit_dirichlet(x,max_iter=200,tol=1e-8_dp)
  if(maxval(abs(fr%alpha-alpha))>0.35_dp) then; print *,fr%alpha; error stop 'fit_dirichlet'; end if
  call compositional_covariance(x,cov,center)
  if(abs(sum(center)-1.0_dp)>1e-12_dp) error stop 'center closure'
  pr=compositional_pca(x)
  if(size(pr%eigenvalues)/=2) error stop 'pca rank'
  nt=acomp_normal_location_one_sample(x)
  if(nt%p_value<0.0_dp.or.nt%p_value>1.0_dp) error stop 'location pvalue'

  ! Exact ILR regression recovery.
  do i=1,size(design,1)
    design(i,1)=1.0_dp; design(i,2)=real(i-2000,dp)/2000.0_dp
  end do
  btrue=reshape([0.2_dp,-0.1_dp,0.4_dp,0.25_dp],[2,2])
  z=matmul(design,btrue); y=ilr_inv_rows(z)
  lm=compositional_lm_fit(design,y)
  if(.not.lm%ok) error stop 'lm fit status'
  if(maxval(abs(lm%coefficients-btrue))>1e-10_dp) error stop 'lm coefficients'
  pred=compositional_lm_predict(lm,design(1:5,:))
  if(maxval(abs(pred-y(1:5,:)))>1e-10_dp) error stop 'lm prediction'

  ! Exercise supplied robustbase dependency through robust covariance.
  call robust_compositional_covariance(x(1:500,:),rcov,rcenter,alpha=0.5_dp)
  if(any(rcov/=rcov)) error stop 'robust covariance NaN'
  if(abs(sum(rcenter)-1.0_dp)>1e-10_dp) error stop 'robust center closure'
  print *, 'test_stats: PASS'
end program
