program smooth_example
  use locfit
  implicit none
  integer,parameter::n=21,m=9
  real(dp)::x(n,1),y(n),xe(m,1)
  type(locfit_options)::opt
  type(locfit_result)::fit
  integer::i
  do i=1,n
    x(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
    y(i)=sin(x(i,1))+0.08_dp*cos(7.0_dp*x(i,1))
  end do
  do i=1,m
    xe(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/real(m-1,dp)
  end do
  opt%nn=0.55_dp
  opt%degree=2
  opt%kernel=wtcub
  call locfit_fit(x,y,xe,fit,opt)
  write(*,'(a)') '       x             fit             se       status'
  do i=1,m
    write(*,'(3f16.8,2x,i0)') xe(i,1),fit%fit(i),fit%se(i),fit%status(i)
  end do
end program smooth_example
