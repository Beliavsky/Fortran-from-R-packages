! SPDX-License-Identifier: GPL-2.0-or-later
program test_shapes
   use icsnp, only : dp, icsnp_ok, tyler_shape, duembgen_shape, duembgen_shape_wt, &
      symm_huber, symm_huber_wt, HP1_shape, HR_Mest, spatial_sign, &
      location_scatter_result, spatial_sign_result
   implicit none
   integer, parameter :: n = 80, p = 3
   real(dp) :: x(n, p), weights(n), zero(p)
   real(dp), allocatable :: shape(:,:)
   type(location_scatter_result) :: robust
   type(spatial_sign_result) :: signs
   integer :: i, status, iterations

   do i = 1, n
      x(i, 1) = sin(0.37_dp * real(i, dp)) + 0.2_dp * cos(0.11_dp * real(i, dp))
      x(i, 2) = 0.5_dp * x(i, 1) + cos(0.23_dp * real(i, dp))
      x(i, 3) = -0.3_dp * x(i, 1) + 0.4_dp * x(i, 2) + sin(0.17_dp * real(i, dp))
      weights(i) = 0.5_dp + real(mod(i, 7), dp) / 7.0_dp
   end do
   zero = 0.0_dp

   call tyler_shape(x, shape, status, iterations)
   call check(status == icsnp_ok .and. symmetric_positive_diag(shape), 'tyler_shape')
   call duembgen_shape(x, shape, status, iterations)
   call check(status == icsnp_ok .and. symmetric_positive_diag(shape), 'duembgen_shape')
   call duembgen_shape_wt(x, weights, shape, status, iterations)
   call check(status == icsnp_ok .and. symmetric_positive_diag(shape), 'duembgen_shape_wt')
   call symm_huber(x, shape, status, iterations)
   call check(status == icsnp_ok .and. symmetric_positive_diag(shape), 'symm_huber')
   call symm_huber_wt(x, weights, shape, status, iterations)
   call check(status == icsnp_ok .and. symmetric_positive_diag(shape), 'symm_huber_wt')

   call spatial_sign(x, signs, center=zero, estimate_center=.false., estimate_shape=.true.)
   call check(signs%status == icsnp_ok .and. size(signs%signs, 1) == n, 'spatial_sign')
   call check(maxval(sqrt(sum(signs%signs * signs%signs, dim=2))) <= 1.0_dp + 1.0e-10_dp, &
      'spatial_sign norms')

   call HR_Mest(x, robust, maxiter=100)
   call check(robust%status == icsnp_ok .and. size(robust%center) == p, 'HR_Mest')
   call HP1_shape(x, shape, status, location=robust%center)
   call check(status == icsnp_ok .and. symmetric_positive_diag(shape), 'HP1_shape')
   print '(a)', 'test_shapes: PASS'
contains
   logical function symmetric_positive_diag(a)
      real(dp), intent(in) :: a(:,:)
      symmetric_positive_diag = size(a, 1) == p .and. size(a, 2) == p .and. &
         maxval(abs(a - transpose(a))) < 1.0e-8_dp .and. all(diagonal(a) > 0.0_dp)
   end function symmetric_positive_diag

   pure function diagonal(a) result(d)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: d(min(size(a, 1), size(a, 2)))
      integer :: j
      do j = 1, size(d)
         d(j) = a(j, j)
      end do
   end function diagonal

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // trim(label)
         error stop 1
      end if
   end subroutine check
end program test_shapes
