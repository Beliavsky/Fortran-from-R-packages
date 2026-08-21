program test_vectors
   use bridgedist, only : dp, dbridge, pbridge, qbridge, dbridge_recycle, pbridge_recycle, qbridge_recycle
   implicit none
   real(dp), allocatable :: got(:)
   real(dp) :: x(5), phi(2), p(3)
   integer :: i, fails

   fails = 0
   x = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   phi = [0.25_dp, 0.75_dp]
   p = [0.1_dp, 0.5_dp, 0.9_dp]

   got = dbridge_recycle(x, phi)
   do i = 1, size(got)
      call check(got(i), dbridge(x(i), phi(1 + modulo(i - 1, 2))), 1.0e-14_dp, 'density recycle')
   end do

   got = pbridge_recycle(x, phi)
   do i = 1, size(got)
      call check(got(i), pbridge(x(i), phi(1 + modulo(i - 1, 2))), 1.0e-14_dp, 'cdf recycle')
   end do

   got = qbridge_recycle(p, phi)
   do i = 1, size(got)
      call check(got(i), qbridge(p(1 + modulo(i - 1, 3)), phi(1 + modulo(i - 1, 2))), 1.0e-14_dp, 'quantile recycle')
   end do

   if (fails /= 0) error stop 1
   print '(a)', 'test_vectors: PASS'

contains

   subroutine check(a, b, rtol, label)
      real(dp), intent(in) :: a, b, rtol
      character(len=*), intent(in) :: label
      if (abs(a - b) > rtol * max(1.0_dp, abs(b))) then
         print '(a,2(1x,es24.16))', trim(label)//' FAIL', a, b
         fails = fails + 1
      end if
   end subroutine check

end program test_vectors
