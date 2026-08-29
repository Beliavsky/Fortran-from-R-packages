program test_toms1006
   use r_compat, only: dp, r_lgamma, pgamma
   use dpq_toms1006
   implicit none
   type(dltgamma_result)::z
   real(dp)::want
   call check_close(lgamma_p11(3.2_dp),r_lgamma(3.2_dp),3.0e-14_dp,'lgammaP11')
   z=dltgamma_inc(0.0_dp,2.0_dp,1.0_dp,3.0_dp)
   want=exp(r_lgamma(3.0_dp))*pgamma(2.0_dp,3.0_dp,1.0_dp)
   call check_close(z%value(),want,3.0e-14_dp,'dltgammaInc positive mu')
   z=dltgamma_inc(0.0_dp,1.0_dp,-1.0_dp,1.0_dp)
   call check_close(z%value(),exp(1.0_dp)-1.0_dp,3.0e-14_dp,'dltgammaInc negative mu')
   print '(a)', 'test_toms1006: PASS'
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
