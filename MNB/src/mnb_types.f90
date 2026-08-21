! SPDX-License-Identifier: GPL-2.0-or-later
module mnb_types
  use mnb_kinds, only : dp
  implicit none
  private
  public :: mnb_fit_result, mnb_residual_result, mnb_global_result
  public :: mnb_local_result, mnb_envelope_result

  type :: mnb_fit_result
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: hessian(:,:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: se(:), z(:), p_value(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: convergence = 1
  end type

  type :: mnb_residual_result
    real(dp), allocatable :: weighted(:)
    real(dp), allocatable :: standardized_weighted(:)
    real(dp), allocatable :: pearson(:)
    real(dp), allocatable :: standardized_pearson(:)
    real(dp), allocatable :: deviance(:)
    real(dp), allocatable :: leverage(:)
  end type

  type :: mnb_global_result
    real(dp), allocatable :: cook_distance(:)
    real(dp), allocatable :: likelihood_displacement(:)
  end type

  type :: mnb_local_result
    real(dp), allocatable :: direction(:)
    real(dp), allocatable :: total_curvature(:)
    real(dp) :: selected_eigenvalue = 0.0_dp
  end type

  type :: mnb_envelope_result
    real(dp), allocatable :: lower(:), mean(:), upper(:)
    real(dp), allocatable :: residual(:)
  end type
end module mnb_types
