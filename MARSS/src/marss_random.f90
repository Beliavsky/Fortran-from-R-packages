! SPDX-License-Identifier: GPL-2.0-only
module marss_random
   use marss_kinds, only : dp, i8
   use r_linalg, only : symmetric_eigen
   implicit none
   private
   public :: marss_rng
   public :: rng_seed
   public :: rng_uniform
   public :: rng_normal
   public :: sample_mvn

   type :: marss_rng
      integer(i8) :: state = 12345_i8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type marss_rng

contains

   pure subroutine rng_seed(rng, seed)
      type(marss_rng), intent(out) :: rng !! Generator state to initialize.
      integer(i8), intent(in) :: seed !! Positive integer seed, zero is replaced by one.

      rng%state = modulo(abs(seed), 2147483646_i8) + 1_i8
      rng%has_spare = .false.
      rng%spare = 0.0_dp
   end subroutine rng_seed

   real(dp) function rng_uniform(rng) result(u)
      type(marss_rng), intent(inout) :: rng !! Generator state advanced by one Park-Miller draw.
      integer(i8), parameter :: modulus = 2147483647_i8
      integer(i8), parameter :: multiplier = 16807_i8

      rng%state = modulo(multiplier * rng%state, modulus)
      if (rng%state == 0_i8) rng%state = 1_i8
      u = real(rng%state, dp) / real(modulus, dp)
   end function rng_uniform

   real(dp) function rng_normal(rng) result(z)
      type(marss_rng), intent(inout) :: rng !! Generator state advanced as needed for a standard-normal draw.
      real(dp) :: r
      real(dp) :: theta
      real(dp) :: u1
      real(dp) :: u2
      real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp

      if (rng%has_spare) then
         z = rng%spare
         rng%has_spare = .false.
         return
      end if
      u1 = max(rng_uniform(rng), tiny(1.0_dp))
      u2 = rng_uniform(rng)
      r = sqrt(-2.0_dp * log(u1))
      theta = two_pi * u2
      z = r * cos(theta)
      rng%spare = r * sin(theta)
      rng%has_spare = .true.
   end function rng_normal

   subroutine sample_mvn(mean, covariance, rng, draw, info)
      real(dp), intent(in) :: mean(:) !! Mean vector of the multivariate normal distribution.
      real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-semidefinite covariance matrix.
      type(marss_rng), intent(inout) :: rng !! Generator state advanced by the required normal draws.
      real(dp), intent(out) :: draw(:) !! Generated multivariate-normal vector with size equal to mean.
      integer, intent(out) :: info !! Zero on success, nonzero if the eigendecomposition fails or shapes mismatch.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp), allocatable :: z(:)
      integer :: i
      integer :: n

      n = size(mean)
      if (size(covariance, 1) /= n .or. size(covariance, 2) /= n .or. size(draw) /= n) then
         info = -1
         return
      end if
      call symmetric_eigen(covariance, values, vectors, info)
      if (info /= 0) return
      if (any(values < -1.0e-10_dp)) then
         info = -2
         return
      end if
      allocate(z(n))
      do i = 1, n
         z(i) = rng_normal(rng) * sqrt(max(values(i), 0.0_dp))
      end do
      draw = mean + matmul(vectors, z)
      info = 0
   end subroutine sample_mvn

end module marss_random
