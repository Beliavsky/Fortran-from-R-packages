! SPDX-License-Identifier: GPL-3.0-only
module ao_types
  use ao_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: AO_PARTITION_SEQUENTIAL = 1
  integer, parameter, public :: AO_PARTITION_RANDOM     = 2
  integer, parameter, public :: AO_PARTITION_NONE       = 3
  integer, parameter, public :: AO_PARTITION_CUSTOM     = 4

  integer, parameter, public :: AO_BASE_BFGS        = 1
  integer, parameter, public :: AO_BASE_NELDER_MEAD = 2
  integer, parameter, public :: AO_BASE_NEWTON      = 3

  abstract interface
    function ao_objective_fn(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function ao_objective_fn

    subroutine ao_gradient_fn(x, gradient)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
    end subroutine ao_gradient_fn

    subroutine ao_hessian_fn(x, hessian)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hessian(:,:)
    end subroutine ao_hessian_fn

    function ao_norm_fn(x, y) result(value)
      import dp
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: value
    end function ao_norm_fn
  end interface

  type, public :: ao_block
    integer, allocatable :: index(:)
  end type ao_block

  type, public :: ao_real_vector
    real(dp), allocatable :: value(:)
  end type ao_real_vector

  type, public :: ao_options
    integer :: partition = AO_PARTITION_SEQUENTIAL
    real(dp) :: new_block_probability = 0.3_dp
    integer :: minimum_block_number = 1
    logical :: minimize = .true.
    integer :: iteration_limit = 1000000000
    real(dp) :: seconds_limit = huge(1.0_dp)
    real(dp) :: tolerance_value = 1.0e-6_dp
    real(dp) :: tolerance_parameter = 1.0e-6_dp
    integer :: tolerance_history = 1
    integer :: base_optimizer = AO_BASE_BFGS
    integer :: base_max_iterations = 10
    integer :: lbfgs_memory = 5
    real(dp) :: base_gradient_tolerance = 1.0e-8_dp
    real(dp) :: finite_difference_step = 1.0e-6_dp
    logical :: add_details = .true.
    logical :: verbose = .false.
    type(ao_block), allocatable :: custom_partition(:)
  end type ao_options

  type, public :: ao_history
    integer :: n = 0
    integer :: capacity = 0
    integer, allocatable :: iteration(:)
    real(dp), allocatable :: value(:)
    real(dp), allocatable :: parameter(:,:)
    logical, allocatable :: active(:,:)
    real(dp), allocatable :: seconds(:)
  end type ao_history

  type, public :: ao_result
    real(dp), allocatable :: estimate(:)
    real(dp) :: value = huge(1.0_dp)
    real(dp) :: seconds = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
    character(len=160) :: stopping_reason = 'not terminated yet'
    type(ao_history) :: details
  end type ao_result

  type, public :: ao_multi_result
    type(ao_result), allocatable :: process(:)
    integer :: best_process = 0
    real(dp), allocatable :: estimate(:)
    real(dp) :: value = huge(1.0_dp)
    real(dp) :: seconds = 0.0_dp
  end type ao_multi_result

  public :: ao_objective_fn, ao_gradient_fn, ao_hessian_fn, ao_norm_fn

end module ao_types
