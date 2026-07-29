! SPDX-License-Identifier: MIT
program test_distribution
   use, intrinsic :: iso_fortran_env, only : int64
   use vasicekfit, only : dp, vasicek_density, vasicek_cdf, vasicek_quantile, random_vasicek
   implicit none
   real(dp), parameter :: x(5) = [0.01_dp, 0.03_dp, 0.05_dp, 0.10_dp, 0.20_dp]
   real(dp), parameter :: density_ref(5) = [ &
      26.380177556686796_dp, 16.789064887976860_dp, 6.946711528052643_dp, &
       0.7472374032107312_dp, 0.012219066333581687_dp ]
   real(dp), parameter :: cdf_ref(5) = [ &
      0.1511644505055207_dp, 0.6198971909679910_dp, 0.8444772584241776_dp, &
      0.9822643577526721_dp, 0.9996900194957168_dp ]
   real(dp), parameter :: prob(4) = [0.01_dp, 0.5_dp, 0.95_dp, 0.99_dp]
   real(dp), parameter :: quantile_ref(4) = [ &
      0.002907989337226529_dp, 0.023709946567140396_dp, &
      0.0757510379137456_dp, 0.11370042836916827_dp ]
   real(dp), allocatable :: sample(:)
   real(dp) :: value, integral, h, kappa(2), factors(2)
   integer :: i, n
   logical :: ok

   do i = 1, size(x)
      call assert_close(vasicek_density(x(i), 0.03_dp, 0.10_dp), density_ref(i), 2.0e-11_dp)
      call assert_close(vasicek_cdf(x(i), 0.03_dp, 0.10_dp), cdf_ref(i), 2.0e-12_dp)
      call assert_close(vasicek_quantile(cdf_ref(i), 0.03_dp, 0.10_dp), x(i), 2.0e-11_dp)
   end do
   do i = 1, size(prob)
      call assert_close(vasicek_quantile(prob(i), 0.03_dp, 0.10_dp), quantile_ref(i), 2.0e-11_dp)
   end do

   call assert_close(vasicek_density(0.0_dp, 0.03_dp, 0.10_dp), 0.0_dp, 0.0_dp)
   call assert_close(vasicek_density(1.0_dp, 0.03_dp, 0.10_dp), 0.0_dp, 0.0_dp)
   call assert_close(vasicek_cdf(0.05_dp, 0.03_dp, 0.10_dp, lower_tail=.false.) + &
      vasicek_cdf(0.05_dp, 0.03_dp, 0.10_dp), 1.0_dp, 2.0e-15_dp)
   call assert_close(vasicek_cdf(0.05_dp, 0.03_dp, 0.10_dp, log_probability=.true.), &
      log(vasicek_cdf(0.05_dp, 0.03_dp, 0.10_dp)), 2.0e-15_dp)

   n = 20000
   h = (1.0_dp - 2.0e-7_dp) / real(n,dp)
   integral = 0.0_dp
   do i = 0, n
      value = vasicek_density(1.0e-7_dp + real(i,dp) * h, 0.03_dp, 0.10_dp)
      if (i == 0 .or. i == n) then
         integral = integral + value
      else if (mod(i,2) == 0) then
         integral = integral + 2.0_dp * value
      else
         integral = integral + 4.0_dp * value
      end if
   end do
   integral = integral * h / 3.0_dp
   call assert_close(integral, 1.0_dp, 2.0e-5_dp)

   kappa = [0.1_dp, -0.05_dp]
   factors = [1.5_dp, 0.5_dp]
   value = vasicek_density(0.05_dp, 0.03_dp, 0.10_dp, kappa, factors, ok=ok)
   call assert_true(ok .and. value > 0.0_dp)

   allocate(sample(100000))
   call random_vasicek(sample, 0.04_dp, 0.08_dp, seed=42_int64, ok=ok)
   call assert_true(ok)
   call assert_close(sum(sample) / real(size(sample),dp), 0.04_dp, 0.002_dp)
   call assert_true(all(sample > 0.0_dp) .and. all(sample < 1.0_dp))

   value = vasicek_density(0.05_dp, 0.0_dp, 0.1_dp, ok=ok)
   call assert_true(.not. ok)

   print '(a)', 'test_distribution: PASS'

contains

   subroutine assert_close(actual, expected, tolerance)
      real(dp), intent(in) :: actual, expected, tolerance
      if (abs(actual - expected) > tolerance + tolerance * abs(expected)) then
         print '(a,3es25.16)', 'mismatch: ', actual, expected, abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true

end program test_distribution
