program test_rsadd
  use relsurv
  implicit none
  real(dp)::x(4,1),off(4),h(4),ds(4),dur(4)
  integer::event(4)
  type(rsadd_result)::fit
  x=1.0_dp; off=0.0_dp; h=0.0_dp; ds=0.0_dp; dur=1.0_dp; event=[1,0,1,0]
  call rsadd_ml_rows(x,off,event,h,ds,dur,fit,maxiter=100,tol=1.0e-12_dp)
  if(.not.fit%converged) then; print *, 'FAIL convergence'; error stop 1; end if
  if(abs(fit%coef(1)-log(0.5_dp))>1.0e-8_dp) then
    print *, 'FAIL coefficient',fit%coef(1),log(0.5_dp); error stop 1
  end if
  print *, 'test_rsadd: PASS'
end program
