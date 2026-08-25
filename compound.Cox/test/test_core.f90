program test_core
  use compound_cox
  implicit none
  real(dp) :: t(8), pi(8), x1(8,2)
  integer :: d(8)
  type(cg_result) :: cg
  type(univariate_result) :: ur
  type(cg_test_result) :: ct
  real(dp), allocatable :: xp(:,:), xt(:,:)

  t = [1.0_dp,3.0_dp,5.0_dp,4.0_dp,7.0_dp,8.0_dp,10.0_dp,13.0_dp]
  d = [1,0,0,1,1,0,1,0]
  call cg_clayton(t,d,18.0_dp,cg)
  call assert_close(cg%surv(3),0.62624616_dp,2.0e-7_dp,'Clayton reference')
  call assert_close(cg%surv(7),0.12500003_dp,2.0e-7_dp,'Clayton tail')
  call cg_gumbel(t,d,2.0_dp,cg)
  call assert_close(cg%surv(3),0.64716821_dp,2.0e-7_dp,'Gumbel reference')
  call cg_frank(t,d,9.0_dp,cg)
  call assert_close(cg%surv(3),0.65242443_dp,2.0e-7_dp,'Frank reference')
  call assert_close(cg%tau,0.6367259208_dp,2.0e-7_dp,'Frank tau')

  x1(:,1) = [0.2_dp,-0.1_dp,0.4_dp,0.1_dp,0.7_dp,-0.2_dp,0.9_dp,0.0_dp]
  x1(:,2) = [1.0_dp,0.4_dp,0.2_dp,0.3_dp,-0.1_dp,-0.5_dp,-0.6_dp,-0.8_dp]
  call uni_score(t,d,x1,ur)
  call assert_close(ur%beta(1),0.7631160572_dp,2.0e-9_dp,'score beta 1')
  call assert_close(ur%beta(2),3.6649699213_dp,2.0e-9_dp,'score beta 2')
  call assert_close(ur%z(1),0.6379617783_dp,2.0e-9_dp,'score z 1')
  call assert_close(ur%z(2),2.7493636469_dp,2.0e-9_dp,'score z 2')
  if (any(.not.(ur%p >= 0.0_dp .and. ur%p <= 1.0_dp))) error stop 'uni_score p range'
  call uni_wald(t,d,x1,ur)
  if (any(.not.(ur%p >= 0.0_dp .and. ur%p <= 1.0_dp))) error stop 'uni_wald p range'

  pi = [8.0_dp,7.0_dp,6.0_dp,5.0_dp,4.0_dp,3.0_dp,2.0_dp,1.0_dp]
  call cg_test(t,d,pi,18.0_dp,ct,copula='clayton',nperm=20)
  if (ct%rmst_good < 0.0_dp .or. ct%rmst_poor < 0.0_dp) error stop 'cg_test rmst'
  if (ct%p_value < 0.0_dp .or. ct%p_value > 1.0_dp) error stop 'cg_test p'

  call x_pathway(4000,6,2,2,xp,rho1=0.2_dp,rho2=0.8_dp)
  call assert_near_zero(maxval(abs(sum(xp,dim=1)/real(size(xp,1),dp))),0.08_dp,'x_pathway means')
  call x_tag(4000,6,2,xt,s=1)
  call assert_near_zero(maxval(abs(sum(xt,dim=1)/real(size(xt,1),dp))),0.08_dp,'x_tag means')

  print '(a)', 'test_core: PASS'
contains
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol)then
    print *,trim(msg),x,y
    error stop 1
    end if
  end subroutine
  subroutine assert_near_zero(x,tol,msg)
    real(dp),intent(in)::x,tol
    character(len=*),intent(in)::msg
    if(abs(x)>tol)then
    print *,trim(msg),x
    error stop 1
    end if
  end subroutine
end program test_core
