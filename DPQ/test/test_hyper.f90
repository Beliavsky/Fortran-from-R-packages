program test_hyper
   use r_compat, only: dp, phyper, dhyper
   use dpq_hyper
   implicit none
   real(dp)::pe
   pe=phyper(4.0_dp,12,18,8)
   call check_close(phyper_r(4.0_dp,12.0_dp,18.0_dp,8.0_dp),pe,2.0e-14_dp,'phyper_r')
   call check_close(phyper_r2(4.0_dp,12.0_dp,18.0_dp,8.0_dp),pe,2.0e-14_dp,'phyper_r2')
   call check_close(pdhyper(4.0_dp,12.0_dp,18.0_dp,8.0_dp),pe/dhyper(4.0_dp,12,18,8),2.0e-13_dp,'pdhyper')
   if (.not.(phyper_appr_as152(4.0_dp,12.0_dp,18.0_dp,8.0_dp) >= 0.0_dp .and. &
      phyper_appr_as152(4.0_dp,12.0_dp,18.0_dp,8.0_dp) <= 1.0_dp)) error stop 'AS152 range'
   if (.not.(phyper_ibeta(4.0_dp,12.0_dp,18.0_dp,8.0_dp) >= 0.0_dp .and. &
      phyper_ibeta(4.0_dp,12.0_dp,18.0_dp,8.0_dp) <= 1.0_dp)) error stop 'Ibeta range'
   print '(a)', 'test_hyper: PASS'
contains
   subroutine check_close(got,want,tol,label)
      real(dp),intent(in)::got,want,tol
      character(len=*),intent(in)::label
      if (.not.(abs(got-want) <= tol*max(1.0_dp,abs(want)))) then
         write(*,'(a,2es24.15)') trim(label)//' failed: ',got,want
         error stop 1
      end if
   end subroutine
end program
