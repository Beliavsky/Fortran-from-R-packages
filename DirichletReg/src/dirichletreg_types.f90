! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_types
  use dirichletreg_kinds, only : dp
  implicit none
  private
  public :: design_block, dirichletreg_model

  type :: design_block
    real(dp), allocatable :: x(:,:)
  end type design_block

  type :: dirichletreg_model
    logical :: alternative = .false.
    integer :: dims = 0
    integer :: base = 1
    integer :: npar = 0
    integer :: convergence = 1
    integer :: bfgs_iterations = 0
    integer :: newton_iterations = 0
    real(dp) :: loglik = 0.0_dp
    real(dp) :: aic = 0.0_dp
    real(dp) :: bic = 0.0_dp
    real(dp) :: nobs = 0.0_dp
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: hessian(:,:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: se(:)
    real(dp), allocatable :: alpha(:,:)
    real(dp), allocatable :: mu(:,:)
    real(dp), allocatable :: phi(:)
    integer, allocatable :: n_vars(:)
  end type dirichletreg_model

end module dirichletreg_types
