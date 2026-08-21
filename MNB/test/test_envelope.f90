program test_envelope
  use mnb, only : dp,mnb_envelope_result,envelope_mnb,set_mnb_seed
  implicit none
  integer,parameter::n=5,mi=2,nn=n*mi
  real(dp)::x(nn,2),y(nn),start(3)
  type(mnb_envelope_result)::e
  integer::i
  x(:,1)=1.0_dp;do i=1,nn;x(i,2)=real(mod((i-1)/mi,2),dp);end do
  y=[1.0_dp,0.0_dp,2.0_dp,1.0_dp,0.0_dp,1.0_dp,3.0_dp,2.0_dp,1.0_dp,1.0_dp]
  start=[2.0_dp,0.0_dp,0.1_dp];call set_mnb_seed(99)
  e=envelope_mnb(start,y,x,n,mi,3,3)
  if(size(e%lower)/=nn)error stop 1
  if(any(e%lower>e%upper))error stop 2
  print *, 'test_envelope: PASS'
end program
