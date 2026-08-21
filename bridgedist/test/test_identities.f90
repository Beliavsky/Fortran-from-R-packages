program test_identities
   use bridgedist, only : dp, dbridge, pbridge, qbridge
   implicit none
   real(dp), parameter :: phis(4) = [0.1_dp, 0.3_dp, 0.65_dp, 0.95_dp]
   real(dp), parameter :: probs(7) = [1.0e-8_dp, 0.01_dp, 0.2_dp, 0.5_dp, 0.8_dp, 0.99_dp, 1.0_dp - 1.0e-8_dp]
   real(dp) :: x, h, deriv
   integer :: i, j, fails

   fails = 0
   do i = 1, size(phis)
      do j = 1, size(probs)
         x = qbridge(probs(j), phis(i))
         call check(pbridge(x, phis(i)), probs(j), 2.0e-10_dp, 'cdf-quantile')
      end do
      do j = -4, 4
         x = real(j, dp)
         call check(qbridge(pbridge(x, phis(i)), phis(i)), x, 2.0e-10_dp, 'quantile-cdf')
         call check(pbridge(-x, phis(i)), 1.0_dp - pbridge(x, phis(i)), 2.0e-13_dp, 'symmetry')
         h = 1.0e-5_dp
         deriv = (pbridge(x + h, phis(i)) - pbridge(x - h, phis(i))) / (2.0_dp * h)
         call check(deriv, dbridge(x, phis(i)), 2.0e-8_dp, 'cdf derivative')
      end do
   end do

   if (fails /= 0) error stop 1
   print '(a)', 'test_identities: PASS'

contains

   subroutine check(got, expected, rtol, label)
      real(dp), intent(in) :: got, expected, rtol
      character(len=*), intent(in) :: label
      if (abs(got - expected) > rtol * max(1.0_dp, abs(expected))) then
         print '(a,2(1x,es24.16))', trim(label)//' FAIL', got, expected
         fails = fails + 1
      end if
   end subroutine check

end program test_identities
