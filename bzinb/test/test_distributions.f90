program test_distributions
  use bzinb
  use test_support
  implicit none
  integer :: fail,x,y
  real(dp)::s,p(4),a
  fail=0
  call assert_close('bp00',bp_pmf(0,0,1.0_dp,2.0_dp,3.0_dp),exp(-6.0_dp),1e-13_dp,fail)
  s=0
  do x=0,25;do y=0,25;s=s+bp_pmf(x,y,0.7_dp,1.2_dp,0.8_dp);end do;end do
  call assert_close('bp mass',s,1.0_dp,2e-10_dp,fail)
  p=[0.5_dp,0.2_dp,0.2_dp,0.1_dp]
  s=0
  do x=0,60;do y=0,60;s=s+bzinb_pmf(x,y,1.5_dp,0.8_dp,1.1_dp,0.7_dp,1.0_dp,p(1),p(2),p(3),p(4));end do;end do
  call assert_close('bzinb mass',s,1.0_dp,2e-8_dp,fail)
  call assert_close('bzipA identity',bzip_a_pmf(1,2,0.5_dp,1.0_dp,1.5_dp,0.0_dp), &
                    bp_pmf(1,2,0.5_dp,1.0_dp,1.5_dp),1e-14_dp,fail)
  a=1.25_dp
  call assert_close('idigamma',digamma_fn(inverse_digamma(digamma_fn(a))),digamma_fn(a),1e-10_dp,fail)
  if(fail==0)then;print *,'test_distributions: PASS';else;error stop 1;end if
end program
