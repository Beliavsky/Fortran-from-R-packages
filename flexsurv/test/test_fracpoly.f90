program test_fracpoly
  use flexsurv
  implicit none
  real(dp)::x(2),p(2),b(2,2),db(2,2)
  integer::fails
  fails=0;x=[2.0_dp,3.0_dp];p=[0.0_dp,0.0_dp]
  call bfp(x,p,b);call dbfp(x,p,db)
  if(maxval(abs(b(:,1)-log(x)))>1e-13_dp)fails=fails+1
  if(maxval(abs(b(:,2)-log(x)**2))>1e-13_dp)fails=fails+1
  if(maxval(abs(db(:,1)-1.0_dp/x))>1e-13_dp)fails=fails+1
  if(maxval(abs(db(:,2)-2.0_dp*log(x)/x))>1e-13_dp)fails=fails+1
  if(fails>0)error stop 1
  print *,'test_fracpoly: PASS'
end program test_fracpoly
