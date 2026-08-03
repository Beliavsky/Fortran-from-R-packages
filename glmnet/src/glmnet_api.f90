! SPDX-License-Identifier: GPL-2.0-only
module glmnet_api
   use glmnet_kinds, only : dp
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument
   use glmnet_types, only : glmnet_control_type, glmnet_path_result, glmnet_sparse_csc
   use glmnet_gaussian, only : fit_gaussian_path, fit_mgaussian_path
   use glmnet_glm, only : fit_binomial_path, fit_poisson_path
   use glmnet_multinomial, only : fit_multinomial_path, fit_multinomial_matrix_path
   use glmnet_cox, only : fit_cox_path
   use glmnet_utils, only : sparse_to_dense
   implicit none
   private
   public :: fit_glmnet, fit_glmnet_sparse, big_glm
   public :: fit_multinomial_path, fit_multinomial_matrix_path, fit_cox_path

   interface fit_glmnet
      module procedure fit_glmnet_vector
      module procedure fit_glmnet_matrix
   end interface fit_glmnet
contains
   subroutine fit_glmnet_vector(x, y, family, result, control, weights, offset, &
      lambda, penalty_factor, lower, upper, excluded)
      real(dp), intent(in) :: x(:,:), y(:)
      character(len=*), intent(in) :: family
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights(:), offset(:), lambda(:)
      real(dp), intent(in), optional :: penalty_factor(:), lower(:), upper(:)
      logical, intent(in), optional :: excluded(:)
      select case (trim(adjustl(family)))
      case ('gaussian')
         call fit_gaussian_path(x, y, result, control, weights, offset, lambda, &
            penalty_factor, lower, upper, excluded)
      case ('binomial')
         call fit_binomial_path(x, y, result, control, weights, offset, lambda, &
            penalty_factor, lower, upper, excluded)
      case ('poisson')
         call fit_poisson_path(x, y, result, control, weights, offset, lambda, &
            penalty_factor, lower, upper, excluded)
      case default
         result%status = glmnet_invalid_argument
         result%family = adjustl(family)
         result%nobs = size(x, 1)
         result%nvars = size(x, 2)
      end select
   end subroutine fit_glmnet_vector

   subroutine fit_glmnet_matrix(x, y, family, result, control, weights, offset, &
      lambda, penalty_factor, lower, upper, excluded)
      real(dp), intent(in) :: x(:,:), y(:,:)
      character(len=*), intent(in) :: family
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights(:), offset(:,:), lambda(:)
      real(dp), intent(in), optional :: penalty_factor(:), lower(:), upper(:)
      logical, intent(in), optional :: excluded(:)
      select case (trim(adjustl(family)))
      case ('mgaussian')
         call fit_mgaussian_path(x, y, result, control, weights, offset, lambda, &
            penalty_factor, lower, upper, excluded)
      case ('multinomial')
         call fit_multinomial_matrix_path(x, y, result, control, weights, offset, &
            lambda, penalty_factor, lower, upper, excluded)
      case default
         result%status = glmnet_invalid_argument
         result%family = adjustl(family)
         result%nobs = size(x, 1)
         result%nvars = size(x, 2)
      end select
   end subroutine fit_glmnet_matrix

   subroutine fit_glmnet_sparse(x, y, family, result, control, weights, offset, &
      lambda, penalty_factor, lower, upper, excluded)
      type(glmnet_sparse_csc), intent(in) :: x
      real(dp), intent(in) :: y(:)
      character(len=*), intent(in) :: family
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights(:), offset(:), lambda(:)
      real(dp), intent(in), optional :: penalty_factor(:), lower(:), upper(:)
      logical, intent(in), optional :: excluded(:)
      real(dp), allocatable :: dense(:,:)
      integer :: status
      call sparse_to_dense(x, dense, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call fit_glmnet_vector(dense, y, family, result, control, weights, offset, &
         lambda, penalty_factor, lower, upper, excluded)
   end subroutine fit_glmnet_sparse

   subroutine big_glm(x, y, family, lambda, result, control, weights, offset)
      real(dp), intent(in) :: x(:,:), y(:), lambda
      character(len=*), intent(in) :: family
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights(:), offset(:)
      real(dp) :: lambda_vector(1)
      lambda_vector(1) = max(lambda, 0.0_dp)
      call fit_glmnet_vector(x, y, family, result, control, weights, offset, &
         lambda=lambda_vector)
   end subroutine big_glm
end module glmnet_api
