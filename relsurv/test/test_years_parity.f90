program test_years_parity
  use relsurv, only : dp, years_result, yl2013_result, years_difference, years_yl2013
  implicit none
  real(dp)::t(3),s(3),nr(3),ne(3),fp(3),bh(2,3),bp(2,3),ph(3)
  type(years_result)::yd
  type(yl2013_result)::y13
  t=[1.0_dp,2.0_dp,3.0_dp];s=[0.8_dp,0.6_dp,0.6_dp];nr=[10.0_dp,8.0_dp,6.0_dp]
  ne=[2.0_dp,2.0_dp,0.0_dp];fp=[0.1_dp,0.2_dp,0.25_dp]
  bh(1,:)=[0.18_dp,0.38_dp,0.40_dp];bh(2,:)=[0.22_dp,0.42_dp,0.44_dp]
  bp(1,:)=[0.09_dp,0.19_dp,0.24_dp];bp(2,:)=[0.11_dp,0.21_dp,0.26_dp]
  call years_difference(t,s,nr,ne,fp,yd,scale=1.0_dp,boot_failure=bh,boot_pop_failure=bp)
  call assert_close(yd%observed_area(3),0.6_dp,1.0e-12_dp,'observed area')
  call assert_close(yd%population_area(3),0.3_dp,1.0e-12_dp,'population area')
  call assert_close(yd%estimate(3),0.3_dp,1.0e-12_dp,'years difference')
  if(yd%standard_error(3)<=0.0_dp)error stop 'bootstrap years se'
  ph=[0.01_dp,0.01_dp,0.01_dp]
  call years_yl2013(t,s,nr,ne,ph,y13,scale=1.0_dp)
  if(y13%estimate(3)<=0.0_dp)error stop 'yl2013 estimate'
  print *, 'test_years_parity: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol;character(len=*),intent(in)::msg
    if(abs(a-b)>tol)then;print *,'FAIL ',msg,a,b;error stop 1;end if
  end subroutine assert_close
end program test_years_parity
