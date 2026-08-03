! SPDX-License-Identifier: GPL-2.0-or-later
module fints_apca
   use fints_kinds, only : dp
   use fints_status, only : fints_ok, fints_invalid_input, fints_singular
   use fints_types, only : apca_result
   use fints_linalg, only : symmetric_eigen, least_squares
   use fints_summary_mod, only : sample_mean
   implicit none
   private
   public :: apca

contains

   subroutine apca(x, number_factors, result)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: number_factors
      type(apca_result), intent(out) :: result
      real(dp), allocatable :: centered(:,:), omega(:,:), eigenvalues(:), eigenvectors(:,:)
      real(dp), allocatable :: design(:,:), beta(:), residuals(:), scales(:), scaled(:,:)
      real(dp) :: sse, total_ss, residual_variance
      integer :: observations, series, i, stat, dof

      result = apca_result()
      observations = size(x, 1)
      series = size(x, 2)
      dof = observations - number_factors - 1
      if (observations < 3 .or. series < 1 .or. number_factors < 1 .or. &
         number_factors > min(observations, series) .or. dof <= 0) then
         result%status = fints_invalid_input
         return
      end if

      allocate(centered(observations, series))
      do i = 1, series
         centered(:, i) = x(:, i) - sample_mean(x(:, i))
      end do
      omega = matmul(centered, transpose(centered)) / real(series, dp)
      call symmetric_eigen(omega, eigenvalues, eigenvectors, stat)
      if (stat /= fints_ok) then
         result%status = stat
         return
      end if

      allocate(design(observations, number_factors + 1), scales(series))
      design(:, 1) = 1.0_dp
      design(:, 2:) = eigenvectors(:, 1:number_factors)
      do i = 1, series
         call least_squares(design, x(:, i), beta, residuals, sse, stat)
         if (stat /= fints_ok) then
            result%status = stat
            return
         end if
         residual_variance = sse / real(dof, dp)
         if (residual_variance <= tiny(1.0_dp)) then
            result%status = fints_singular
            return
         end if
         scales(i) = 1.0_dp / sqrt(residual_variance)
      end do

      allocate(scaled(observations, series))
      do i = 1, series
         scaled(:, i) = x(:, i) * scales(i)
         scaled(:, i) = scaled(:, i) - sample_mean(scaled(:, i))
      end do
      omega = matmul(scaled, transpose(scaled)) / real(series, dp)
      call symmetric_eigen(omega, eigenvalues, eigenvectors, stat)
      if (stat /= fints_ok) then
         result%status = stat
         return
      end if

      allocate(result%eigenvalues(size(eigenvalues)))
      allocate(result%factors(observations, number_factors))
      allocate(result%loadings(series, number_factors))
      allocate(result%r_squared(series))
      result%eigenvalues = eigenvalues
      result%factors = eigenvectors(:, 1:number_factors)
      design(:, 2:) = result%factors
      do i = 1, series
         call least_squares(design, x(:, i), beta, residuals, sse, stat)
         if (stat /= fints_ok) then
            result%status = stat
            return
         end if
         result%loadings(i, :) = beta(2:)
         total_ss = sum((x(:, i) - sample_mean(x(:, i))) ** 2)
         if (total_ss > 0.0_dp) then
            result%r_squared(i) = max(0.0_dp, min(1.0_dp, 1.0_dp - sse / total_ss))
         else
            result%r_squared(i) = 0.0_dp
         end if
      end do
      result%status = fints_ok
   end subroutine apca

end module fints_apca
