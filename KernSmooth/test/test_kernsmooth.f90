program test_kernsmooth
  use kernsmooth_mod
  implicit none
  real(dp) :: x(7), y(7), c(21), g(101), d(101), fit(51), gg(51)
  real(dp) :: bw, dx, err, hll
  real(dp) :: xy(5,2), g1(31), g2(31), f2(31,31)
  integer :: i
  x = [-2.0_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,2.0_dp]
  y = 1.0_dp + 2.0_dp*x
  call linbin(x,-2.0_dp,2.0_dp,c)
  call check(abs(sum(c)-7.0_dp)<1e-12_dp,'linbin mass')
  call bkde(x,0.5_dp,g,d,'normal',[-3.0_dp,3.0_dp])
  dx=g(2)-g(1)
  call check(abs(sum(d)*dx-1.0_dp)<2e-2_dp,'bkde normalization')
  call check(minval(d)>=0.0_dp,'bkde nonnegative')
  call locpoly(x,y,0.9_dp,gg,fit,drv=0,degree=1,range_x=[-1.5_dp,1.5_dp])
  err=maxval(abs(fit-(1.0_dp+2.0_dp*gg)))
  call check(err<1e-8_dp,'locpoly linear recovery')
  bw=dpih(x,level=0)
  call check(bw>0.0_dp .and. bw<10.0_dp,'dpih positive')
  bw=dpik(x,level=0)
  call check(bw>0.0_dp .and. bw<10.0_dp,'dpik positive')
  hll=dpill(x,y)
  call check(hll>0.0_dp .and. hll<10.0_dp,'dpill positive')
  xy(:,1)=[-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp]
  xy(:,2)=[-1.0_dp,0.5_dp,0.0_dp,-0.5_dp,1.0_dp]
  call bkde2d(xy,[0.5_dp,0.5_dp],g1,g2,f2)
  call check(minval(f2)>=0.0_dp,'bkde2d nonnegative')
  call check(abs(sum(f2)*(g1(2)-g1(1))*(g2(2)-g2(1))-1.0_dp)<3e-2_dp,'bkde2d normalization')
  print *, 'test_kernsmooth: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine
end program
