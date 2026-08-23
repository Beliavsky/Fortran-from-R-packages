program basic
  use mixspe, only : dp, rpe, spe_model, em_fit
  implicit none
  real(dp), allocatable :: a(:,:), b(:,:), x(:,:)
  real(dp) :: s(2,2)
  type(spe_model) :: fit
  s=0.0_dp; s(1,1)=1.0_dp; s(2,2)=1.0_dp
  a=rpe(100,2.0_dp,[0.0_dp,0.0_dp],s)
  b=rpe(100,3.0_dp,[3.0_dp,0.0_dp],s)
  allocate(x(200,2)); x(1:100,:)=a; x(101:200,:)=b
  call em_fit(x,2,'EIIV',fit,max_iter=100,tol=1.0e-5_dp)
  print '(a,f12.4)', 'logLik = ',fit%loglik
  print '(a,f12.4)', 'BIC    = ',fit%bic
  print '(a,2f10.4)', 'pi     = ',fit%pi
end program
