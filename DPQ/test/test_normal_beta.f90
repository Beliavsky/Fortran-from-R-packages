program test_normal_beta
   use r_compat, only: dp, pbeta
   use dpq_normal_beta
   implicit none
   real(dp)::p,x

   p=1.0e-6_dp
   x=qnorm_r(p)
   call check_close(0.5_dp*erfc(-x/sqrt(2.0_dp)),p,2.0e-12_dp,'qnorm_r inversion')
   p=0.73_dp
   x=qbeta_r(p,2.5_dp,4.0_dp)
   call check_close(pbeta(x,2.5_dp,4.0_dp),p,2.0e-12_dp,'qbeta_r inversion')
   call check_close(pbeta_rv1(0.37_dp,2.5_dp,4.0_dp),pbeta(0.37_dp,2.5_dp,4.0_dp),2.0e-14_dp,'pbetaRv1')
   call check_close(pnbeta_as310(0.37_dp,2.5_dp,4.0_dp,0.0_dp),pbeta(0.37_dp,2.5_dp,4.0_dp),2.0e-13_dp,'pnbeta ncp=0')
   if (.not.(pnorm_u_s53(5.0_dp) > 0.0_dp .and. pnorm_u_s53(5.0_dp) < 1.0_dp)) error stop 'pnormU range'
   if (.not.(qbeta_appr(0.5_dp,10.0_dp,12.0_dp) > 0.0_dp .and. qbeta_appr(0.5_dp,10.0_dp,12.0_dp) < 1.0_dp)) &
      error stop 'qbetaAppr range'
   print '(a)', 'test_normal_beta: PASS'
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
