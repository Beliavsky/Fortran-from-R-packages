program test_pkpd
   use rmutil
   implicit none
   real(dp), allocatable :: y(:), y2(:)
   real(dp) :: p1(5), times(4), ind(4)
   times=[0.1_dp,0.4_dp,0.8_dp,1.5_dp]
   p1=[log(10.0_dp),log(2.0_dp),log(0.4_dp),log(0.8_dp),log(0.3_dp)]
   y=mu1_0o1c(p1(1:2),times,2.0_dp,0.5_dp)
   call check(maxval(y)<1.0_dp .and. minval(y)>0.0_dp,"mu1_0o1c positive")
   y=mu1_1o1c(p1(1:3),times,2.0_dp)
   call check(minval(y)>0.0_dp,"mu1_1o1c positive")
   y=mu1_1o2c(p1(1:4),times,2.0_dp)
   call check(all(abs(y)<100.0_dp),"mu1_1o2c finite")
   ind=[1.0_dp,0.0_dp,1.0_dp,0.0_dp]
   y=mu2_0o1c(p1,times,ind,1.5_dp,0.5_dp)
   call check(all(y>=0.0_dp),"mu2_0o1c nonnegative")
   y2=mu2_0o1cfp([p1,0.2_dp],times,ind,1.5_dp,0.5_dp)
   call check(all(y2>=0.0_dp),"mu2_0o1cfp nonnegative")
   print *, "test_pkpd: PASS"
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok; character(*),intent(in)::msg
      if(.not.ok)then; print *,"FAIL: ",trim(msg); error stop 1; end if
   end subroutine check
end program test_pkpd
