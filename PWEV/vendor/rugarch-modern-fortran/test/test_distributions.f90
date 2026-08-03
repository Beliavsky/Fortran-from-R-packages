program test_distributions
   use rugarch
   implicit none
   real(dp), parameter :: probs(5)=[0.01_dp,0.10_dp,0.50_dp,0.90_dp,0.99_dp]
   real(dp) :: q, p
   integer :: i

   do i=1,size(probs)
      q=qstd(probs(i),0.0_dp,1.0_dp,7.0_dp)
      p=pstd(q,0.0_dp,1.0_dp,7.0_dp)
      call assert_close(p,probs(i),2.0e-6_dp,'std round trip')
      q=qsged(probs(i),0.0_dp,1.0_dp,1.6_dp,1.3_dp)
      p=psged(q,0.0_dp,1.0_dp,1.6_dp,1.3_dp)
      call assert_close(p,probs(i),3.0e-6_dp,'sged round trip')
      q=qjsu(probs(i),0.0_dp,1.0_dp,0.5_dp,1.8_dp)
      p=pjsu(q,0.0_dp,1.0_dp,0.5_dp,1.8_dp)
      call assert_close(p,probs(i),2.0e-7_dp,'jsu round trip')
   end do
   print '(a)', 'distribution tests passed'
contains
   subroutine assert_close(a,b,tol,label)
      real(dp), intent(in) :: a,b,tol
      character(len=*), intent(in) :: label
      if (abs(a-b)>tol) then
         print *,trim(label),a,b
         error stop 'assertion failed'
      end if
   end subroutine assert_close
end program test_distributions
