program test_em_parity
  use bzinb
  use test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  integer, parameter :: n=80
  integer :: x(n),y(n),fail,i
  integer :: xx(1),yy(1),ff(1)
  real(dp) :: p(9),pb(9),h,fd,pp(9),pm(9),ll0
  type(em_expectation_result) :: er
  type(bnb_fit_result) :: fb
  type(bzinb_fit_result) :: fz
  fail=0

  p=[1.5_dp,0.8_dp,1.1_dp,0.7_dp,1.0_dp,0.55_dp,0.15_dp,0.20_dp,0.10_dp]
  call bzinb_expectation(2,1,1,p,er,.true.,.false.)
  call assert_close('expectation loglik',er%expt(1), &
    bzinb_logpmf(2,1,p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9)),3.0e-7_dp,fail)
  do i=1,5
    h=1.0e-6_dp*(1.0_dp+abs(p(i)))
    pp=p;pm=p;pp(i)=pp(i)+h;pm(i)=pm(i)-h
    fd=(bzinb_logpmf(2,1,pp(1),pp(2),pp(3),pp(4),pp(5),pp(6),pp(7),pp(8),pp(9))- &
        bzinb_logpmf(2,1,pm(1),pm(2),pm(3),pm(4),pm(5),pm(6),pm(7),pm(8),pm(9)))/(2.0_dp*h)
    call assert_close('source score a/b',er%score(i),fd,2.0e-6_dp,fail)
  end do
  ! Upstream em.cpp/expt.cpp uses a constrained-mixture score convention in
  ! which a positive-positive pair has score 1 for the p1 coordinate.
  call assert_close('source p1 score convention',er%score(6),1.0_dp,2.0e-7_dp,fail)
  call assert_close('source p2 score convention',er%score(7),0.0_dp,2.0e-7_dp,fail)
  call assert_close('source p3 score convention',er%score(8),0.0_dp,2.0e-7_dp,fail)

  call set_bzinb_seed(7291)
  call rbnb_sample(n,1.3_dp,0.9_dp,1.0_dp,0.8_dp,1.1_dp,x,y)
  pb=[1.2_dp,1.0_dp,1.0_dp,0.9_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
  ll0=loglik_bnb(x,y,pb(1:5))
  fb=fit_bnb_em(x,y,maxiter=700,initial=pb(1:5),tol=2.0e-6_dp)
  call assert_true('bnb em improves',fb%loglik>=ll0-1.0e-7_dp,fail)
  call assert_true('bnb source information',fb%covariance_ok.and.all(fb%se>=0.0_dp),fail)
  call assert_close('bnb returned historical max',fb%loglik,maxval(fb%trajectory(2:)),2.0e-8_dp,fail)

  call rbzinb_sample(n,1.4_dp,0.8_dp,1.0_dp,0.75_dp,1.05_dp,p(6:9),x,y)
  ll0=loglik_bzinb(x,y,p)
  fz=fit_bzinb_em(x,y,maxiter=900,initial=p,tol=5.0e-6_dp)
  call assert_true('bzinb em improves',fz%loglik>=ll0-1.0e-7_dp,fail)
  call assert_close('bzinb p sum',sum(fz%param(6:9)),1.0_dp,2.0e-12_dp,fail)
  call assert_close('bzinb returned historical max',fz%loglik,maxval(fz%trajectory(2:)),2.0e-8_dp,fail)
  call assert_true('bzinb source information',fz%covariance_ok.and.all(fz%se>=0.0_dp),fail)
  fz=fit_bzinb_direct(x,y,maxiter=140,initial=p,tol=1.0e-5_dp)
  call assert_true('direct fallback finite',all(fz%param(1:5)>0.0_dp).and.ieee_is_finite(fz%loglik),fail)
  call assert_close('direct fallback p sum',sum(fz%param(6:9)),1.0_dp,2.0e-12_dp,fail)
  fz=fit_bzinb_em(x,y,maxiter=900,initial=p,tol=5.0e-6_dp)
  call assert_close('p4 variance mapping',fz%covariance(9,9),sum(fz%covariance(6:8,6:8)),2.0e-8_dp,fail)

  x=0;y=0
  fb=fit_bnb(x,y,maxiter=10)
  call assert_close('bnb all-zero shortcut',fb%loglik,0.0_dp,0.0_dp,fail)
  fz=fit_bzinb(x,y,maxiter=10)
  call assert_close('bzinb all-zero p1',fz%param(6),1.0_dp,0.0_dp,fail)
  call assert_close('bzinb all-zero p4',fz%param(9),0.0_dp,0.0_dp,fail)

  xx=[0];yy=[0];ff=[1]
  call bzinb_expectation_vec(xx,yy,ff,p,er,.true.,.false.)
  call assert_true('zero pair finite information',all(ieee_is_finite(er%information)),fail)

  if(fail==0)then
    print *,'test_em_parity: PASS'
  else
    error stop 1
  end if
end program test_em_parity
