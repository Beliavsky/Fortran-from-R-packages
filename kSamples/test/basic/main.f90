program test_basic
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use ksamples, only : dp, conv, djt, pjt, qjt
   implicit none
   real(dp), allocatable :: x(:), p(:), d(:), cdf(:), quantile(:)
   real(dp) :: x1(2), p1(2), x2(2), p2(2), support(5), probs(4)
   integer :: nn(2)

   x1 = [0.0_dp, 1.0_dp]
   p1 = [0.5_dp, 0.5_dp]
   x2 = [0.0_dp, 1.0_dp]
   p2 = [0.5_dp, 0.5_dp]
   call conv(x1, p1, x2, p2, x, p)
   call assert_true(size(x) == 3, 'conv support size')
   call assert_close(x(1), 0.0_dp, 1.0e-14_dp, 'conv x(1)')
   call assert_close(x(2), 1.0_dp, 1.0e-14_dp, 'conv x(2)')
   call assert_close(x(3), 2.0_dp, 1.0e-14_dp, 'conv x(3)')
   call assert_close(p(1), 0.25_dp, 1.0e-14_dp, 'conv p(1)')
   call assert_close(p(2), 0.50_dp, 1.0e-14_dp, 'conv p(2)')
   call assert_close(p(3), 0.25_dp, 1.0e-14_dp, 'conv p(3)')

   nn = [2, 2]
   support = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   call djt(support, nn, d)
   call assert_close(sum(d), 1.0_dp, 2.0e-14_dp, 'djt normalization')
   call assert_close(d(1), 1.0_dp/6.0_dp, 2.0e-14_dp, 'djt lower endpoint')
   call assert_close(d(3), 2.0_dp/6.0_dp, 2.0e-14_dp, 'djt center mass')
   call assert_close(d(5), d(1), 2.0e-14_dp, 'djt symmetry')

   call pjt(support, nn, cdf)
   call assert_close(cdf(3), 4.0_dp/6.0_dp, 2.0e-14_dp, 'pjt center cdf')
   call assert_close(cdf(5), 1.0_dp, 2.0e-14_dp, 'pjt upper endpoint')

   probs = [0.0_dp, 1.0_dp/6.0_dp, 0.5_dp, 1.1_dp]
   call qjt(probs, nn, quantile)
   call assert_true(.not. ieee_is_finite(quantile(1)) .and. quantile(1) < 0.0_dp, 'qjt negative infinity')
   call assert_close(quantile(2), 0.0_dp, 1.0e-14_dp, 'qjt first positive quantile')
   call assert_close(quantile(3), 2.0_dp, 1.0e-14_dp, 'qjt median')
   call assert_true(ieee_is_nan(quantile(4)), 'qjt invalid probability NaN')

   print '(a)', 'basic kSamples tests passed'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed value being checked.
      real(dp), intent(in) :: expected !! Reference value for the check.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference.
      character(len=*), intent(in) :: label !! Human-readable test label shown on failure.
      if (abs(actual - expected) > tolerance) then
         write (*, '(a,2(1x,es24.16))') trim(label)//' failed:', actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Logical condition that must hold.
      character(len=*), intent(in) :: label !! Human-readable test label shown on failure.
      if (.not. condition) then
         write (*, '(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_basic
