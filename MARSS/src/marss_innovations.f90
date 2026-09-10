! SPDX-License-Identifier: GPL-2.0-only
module marss_innovations
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use marss_kinds, only : dp, i8
   use marss_types, only : marss_model, marss_kf_result, marss_innov_boot_result
   use marss_kalman, only : marss_kfss
   use marss_parameters, only : marss_b_at, marss_u_at, marss_z_at, marss_a_at
   use marss_random, only : marss_rng, rng_seed, rng_uniform
   use r_linalg, only : symmetric_eigen
   implicit none
   private
   public :: marss_innovations_boot

contains

   subroutine marss_innovations_boot(model, nboot, min_index, seed, result, info)
      type(marss_model), intent(in) :: model !! Complete-data time-invariant fitted model supplying B,U,Z,A,x0 and observations.
      integer, intent(in) :: nboot !! Number of innovations bootstrap series to generate, must be positive.
      integer, intent(in) :: min_index !! Number of leading innovations not resampled when the series has more than five steps.
      integer(i8), intent(in) :: seed !! Deterministic seed for bootstrap innovation-index resampling.
      type(marss_innov_boot_result), intent(out) :: result !! Bootstrap states and observations, with replicate in dimension three.
      integer, intent(out) :: info !! Zero on success, positive values identify validation or eigendecomposition failures.
      type(marss_kf_result) :: kf
      type(marss_rng) :: rng
      real(dp), allocatable :: sigma_sqrt(:, :, :)
      real(dp), allocatable :: std_innov(:, :)
      real(dp), allocatable :: bks(:, :, :)
      real(dp), allocatable :: states_ext(:, :)
      real(dp), allocatable :: e(:, :)
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      real(dp), allocatable :: root(:, :)
      real(dp), allocatable :: invroot(:, :)
      real(dp) :: at(size(model%a))
      real(dp) :: bt(size(model%b, 1), size(model%b, 2))
      real(dp) :: ut(size(model%u))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      integer, allocatable :: sample_index(:)
      integer :: b
      integer :: eig_info
      integer :: i
      integer :: j
      integer :: m
      integer :: min_use
      integer :: n
      integer :: t
      integer :: tt

      info = 0
      if (nboot < 1) then
         info = 1
         return
      end if
      if (min_index < 0) then
         info = 2
         return
      end if
      if (.not. allocated(model%y)) then
         info = 3
         return
      end if
      do t = 1, size(model%y, 2)
         do i = 1, size(model%y, 1)
            if (ieee_is_nan(model%y(i, t))) then
               info = 4
               return
            end if
         end do
      end do
      call marss_kfss(model, kf)
      if (.not. kf%ok) then
         info = 100 + kf%info
         return
      end if

      n = size(model%y, 1)
      m = size(model%b, 1)
      tt = size(model%y, 2)
      min_use = min_index
      if (tt <= 5) min_use = 1
      min_use = min(min_use, max(tt - 1, 0))
      allocate(sigma_sqrt(n, n, tt), std_innov(n, tt), bks(m, n, tt))
      do t = 1, tt
         call symmetric_eigen(kf%sigma(:, :, t), values, vectors, eig_info)
         if (eig_info /= 0 .or. any(values <= 0.0_dp)) then
            info = 200 + t
            return
         end if
         allocate(root(n, n), invroot(n, n))
         root = 0.0_dp
         invroot = 0.0_dp
         do i = 1, n
            do j = 1, n
               root(i, j) = sum(vectors(i, :) * sqrt(values) * vectors(j, :))
               invroot(i, j) = sum(vectors(i, :) / sqrt(values) * vectors(j, :))
            end do
         end do
         sigma_sqrt(:, :, t) = root
         std_innov(:, t) = matmul(invroot, kf%innov(:, t))
         bt = marss_b_at(model, t)
         bks(:, :, t) = matmul(matmul(bt, kf%gain(:, :, t)), root)
         deallocate(values, vectors, root, invroot)
      end do

      allocate(result%data(n, tt, nboot), result%states(m, tt, nboot))
      allocate(states_ext(m, tt + 1), e(n, tt), sample_index(tt))
      call rng_seed(rng, seed)
      do b = 1, nboot
         do t = 1, tt
            sample_index(t) = t
         end do
         if (tt > min_use) then
            do t = min_use + 1, tt
               sample_index(t) = min_use + 1 + int(rng_uniform(rng) * real(tt - min_use, dp))
               sample_index(t) = min(sample_index(t), tt)
            end do
         end if
         do t = 1, tt
            e(:, t) = std_innov(:, sample_index(t))
         end do
         states_ext(:, 1) = model%x0
         do t = 1, tt
            bt = marss_b_at(model, t)
            ut = marss_u_at(model, t)
            zt = marss_z_at(model, t)
            at = marss_a_at(model, t)
            result%data(:, t, b) = matmul(zt, states_ext(:, t)) + at + &
               matmul(sigma_sqrt(:, :, t), e(:, t))
            states_ext(:, t + 1) = matmul(bt, states_ext(:, t)) + ut + &
               matmul(bks(:, :, t), e(:, t))
         end do
         result%states(:, :, b) = states_ext(:, 2:tt + 1)
      end do
   end subroutine marss_innovations_boot

end module marss_innovations
