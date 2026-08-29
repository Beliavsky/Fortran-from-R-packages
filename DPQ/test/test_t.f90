program test_t
   use r_compat, only: dp, pt, dt
   use dpq_t
   implicit none
   real(dp)::q
   call check_close(pnt_r(1.2_dp,7.5_dp,0.8_dp),0.6346261683484926_dp,3.0e-11_dp,'nct cdf')
   call check_close(dnt_jkbf(1.2_dp,7.5_dp,0.8_dp),0.33818993986527696_dp,3.0e-11_dp,'nct pdf')
   q=qnt_r(0.73_dp,7.5_dp,0.8_dp)
   call check_close(q,1.506689456775855_dp,5.0e-11_dp,'nct quantile')
   call check_close(pnt_r(1.2_dp,7.5_dp,0.0_dp),pt(1.2_dp,7.5_dp),2.0e-14_dp,'central CDF identity')
   call check_close(dnt_jkbf(1.2_dp,7.5_dp,0.0_dp),dt(1.2_dp,7.5_dp),2.0e-14_dp,'central PDF identity')
   if (.not.(pnt_jw39(1.2_dp,7.5_dp,0.8_dp) > 0.0_dp .and. pnt_jw39(1.2_dp,7.5_dp,0.8_dp) < 1.0_dp)) &
      error stop 'JW39 range'
   if (.not.(dt_wv(1.2_dp,7.5_dp,0.8_dp) > 0.0_dp)) error stop 'dtWV range'
   print '(a)', 'test_t: PASS'
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
