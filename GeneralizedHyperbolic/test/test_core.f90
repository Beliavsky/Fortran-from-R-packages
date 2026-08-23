program test_core
  use generalized_hyperbolic
  implicit none
  real(dp) :: p, q, m, v
  real(dp) :: x(200)
  real(dp) :: a(3), b(3), g1(5), g2(5)
  integer :: fails
  type(dist_fit) :: fit
  real(dp) :: gx(30)
  fails=0

  p=pskewlap(0.0_dp,0.0_dp,1.0_dp,1.0_dp)
  if(abs(p-0.5_dp)>1e-12_dp) fails=fails+1
  q=qskewlap(0.2_dp,0.0_dp,1.0_dp,2.0_dp)
  if(abs(pskewlap(q,0.0_dp,1.0_dp,2.0_dp)-0.2_dp)>1e-10_dp) fails=fails+1

  m=gig_mean(1.0_dp,1.0_dp,0.5_dp)
  v=gig_var(1.0_dp,1.0_dp,0.5_dp)
  if(.not.(m>0 .and. v>0)) fails=fails+1
  q=qgig(0.5_dp,1.0_dp,1.0_dp,0.5_dp)
  if(abs(pgig(q,1.0_dp,1.0_dp,0.5_dp)-0.5_dp)>2e-3_dp) fails=fails+1

  a=[1.0_dp,4.0_dp,0.7_dp]
  call gig_change_pars(1,4,a,b)
  call gig_change_pars(4,1,b,a)
  if(maxval(abs(a-[1.0_dp,4.0_dp,0.7_dp]))>1e-10_dp) fails=fails+1

  m=ghyp_mean(0.0_dp,1.0_dp,2.0_dp,0.3_dp,1.0_dp)
  v=ghyp_var(0.0_dp,1.0_dp,2.0_dp,0.3_dp,1.0_dp)
  if(.not.(v>0)) fails=fails+1
  q=qghyp(0.5_dp,0.0_dp,1.0_dp,2.0_dp,0.3_dp,1.0_dp)
  if(abs(pghyp(q,0.0_dp,1.0_dp,2.0_dp,0.3_dp,1.0_dp)-0.5_dp)>1e-2_dp) fails=fails+1
  if(abs(dhyperb(0.2_dp,0.0_dp,1.0_dp,2.0_dp,0.3_dp)-dghyp(0.2_dp,0.0_dp,1.0_dp,2.0_dp,0.3_dp,1.0_dp))>1e-12_dp) fails=fails+1
  if(abs(dnig(0.2_dp,0.0_dp,1.0_dp,2.0_dp,0.3_dp)-dghyp(0.2_dp,0.0_dp,1.0_dp,2.0_dp,0.3_dp,-0.5_dp))>1e-12_dp) fails=fails+1

  g1=[0.2_dp,1.3_dp,2.2_dp,0.4_dp,0.7_dp]
  call ghyp_change_pars(1,5,g1,g2)
  call ghyp_change_pars(5,1,g2,g1)
  if(maxval(abs(g1-[0.2_dp,1.3_dp,2.2_dp,0.4_dp,0.7_dp]))>1e-9_dp) fails=fails+1

  call rskewlap(x,0.0_dp,1.0_dp,1.0_dp)
  if(any(.not.(x<huge(1.0_dp)))) fails=fails+1

  call rskewlap(gx,0.2_dp,1.0_dp,1.0_dp)
  ! Smoke-test likelihood optimizer on positive GIG-like values.
  gx=abs(gx)+0.2_dp
  call gig_fit(gx,fit,maxit=15)
  if(.not.(fit%loglik > -huge(1.0_dp))) fails=fails+1

  if(fails==0)then
    print *, 'test_core: PASS'
  else
    print *, 'test_core: FAIL', fails
    error stop 1
  endif
end program
