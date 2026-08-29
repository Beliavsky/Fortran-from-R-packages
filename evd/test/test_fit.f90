program test_fit
   use evd
   implicit none
   integer, parameter :: n=120
   real(dp) :: x(n), p, start(3), dat(n,2), bstart(7)
   integer :: i, fails
   type(evd_fit_t) :: fg, fb
   fails=0
   do i=1,n
      p=(real(i,dp)-0.5_dp)/real(n,dp)
      x(i)=qgev(p,0.5_dp,1.4_dp,0.10_dp)
   end do
   start=[0.4_dp,1.2_dp,0.05_dp]
   fg=fit_gev(x,start,maxit=250)
   call check(fg%nll<1.0e90_dp .and. fg%estimate(2)>0.0_dp,'GEV fit finite',fails)
   call check(abs(fg%estimate(1)-0.5_dp)<0.25_dp,'GEV fitted location',fails)
   call check(abs(fg%estimate(2)-1.4_dp)<0.30_dp,'GEV fitted scale',fails)

   do i=1,n
      p=(real(i,dp)-0.5_dp)/real(n,dp)
      dat(i,1)=qgev(p,0.0_dp,1.0_dp,0.0_dp)
      dat(i,2)=qgev(mod(real(37*i,dp),real(n,dp))/real(n,dp)*0.98_dp+0.01_dp,0.0_dp,1.0_dp,0.0_dp)
   end do
   bstart=[0.9_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp]
   fb=fit_bvevd(dat,'log',bstart,maxit=120)
   call check(fb%nll<1.0e90_dp .and. fb%estimate(1)>0.0_dp .and. fb%estimate(1)<1.0_dp,'bivariate fit smoke',fails)
   if(fails>0) then
      write(*,'(a,i0)') 'test_fit: FAIL ',fails
      error stop 1
   end if
   write(*,'(a)') 'test_fit: PASS'
contains
   subroutine check(ok,msg,fails)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      integer,intent(inout)::fails
      if(.not.ok) then
      write(*,'(a)') 'FAIL: '//msg
      fails=fails+1
      end if
   end subroutine
end program
