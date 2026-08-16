! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_information
   use vgam_kinds, only : dp
   use vgam_optim, only : numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   abstract interface
      function scalar_objective(x) result(f)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function scalar_objective
   end interface

   public :: observed_information, score_outer_information
   public :: constrained_information, lift_constrained_covariance

contains

   subroutine observed_information(objective, theta, information, covariance, status)
      procedure(scalar_objective) :: objective
      real(dp), intent(in) :: theta(:)
      real(dp), allocatable, intent(out) :: information(:, :), covariance(:, :)
      integer, intent(out) :: status
      allocate(information(size(theta), size(theta)))
      call numerical_hessian(objective, theta, information)
      call invert_matrix(information, covariance, status)
   end subroutine observed_information

   subroutine score_outer_information(scores, information, weights, center)
      real(dp), intent(in) :: scores(:, :)
      real(dp), allocatable, intent(out) :: information(:, :)
      real(dp), intent(in), optional :: weights(:)
      logical, intent(in), optional :: center
      real(dp), allocatable :: w(:), s(:, :), mean_score(:)
      logical :: ctr
      integer :: n, p, i, j, k
      n = size(scores, 1); p = size(scores, 2)
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            allocate(information(0, 0)); return
         end if
         w = weights
      end if
      ctr = .false.; if (present(center)) ctr = center
      s = scores
      if (ctr .and. sum(w) > 0.0_dp) then
         mean_score = matmul(transpose(s), w)/sum(w)
         do i = 1, n
            s(i, :) = s(i, :) - mean_score
         end do
      end if
      allocate(information(p, p)); information = 0.0_dp
      do i = 1, n
         do j = 1, p
            do k = j, p
               information(j, k) = information(j, k) + w(i)*s(i, j)*s(i, k)
               information(k, j) = information(j, k)
            end do
         end do
      end do
   end subroutine score_outer_information

   subroutine constrained_information(full_information, constraint, free_information)
      real(dp), intent(in) :: full_information(:, :), constraint(:, :)
      real(dp), allocatable, intent(out) :: free_information(:, :)
      if (size(full_information, 1) /= size(full_information, 2) .or. &
          size(full_information, 1) /= size(constraint, 1)) then
         allocate(free_information(0, 0)); return
      end if
      free_information = matmul(transpose(constraint), matmul(full_information, constraint))
   end subroutine constrained_information

   subroutine lift_constrained_covariance(free_covariance, constraint, full_covariance)
      real(dp), intent(in) :: free_covariance(:, :), constraint(:, :)
      real(dp), allocatable, intent(out) :: full_covariance(:, :)
      if (size(free_covariance, 1) /= size(free_covariance, 2) .or. &
          size(free_covariance, 1) /= size(constraint, 2)) then
         allocate(full_covariance(0, 0)); return
      end if
      full_covariance = matmul(constraint, matmul(free_covariance, transpose(constraint)))
   end subroutine lift_constrained_covariance

end module vgam_information
