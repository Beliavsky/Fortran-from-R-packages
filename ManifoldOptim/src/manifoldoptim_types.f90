! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
module manifoldoptim_types
  use manifoldoptim_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: MANI_EUCLIDEAN = 1
  integer, parameter, public :: MANI_SPHERE    = 2
  integer, parameter, public :: MANI_STIEFEL   = 3
  integer, parameter, public :: MANI_GRASSMANN = 4
  integer, parameter, public :: MANI_SPD       = 5
  integer, parameter, public :: MANI_LOWRANK   = 6
  integer, parameter, public :: MANI_ORTHGROUP = 7

  integer, parameter, public :: STATUS_SUCCESS = 0
  integer, parameter, public :: STATUS_MAXITER = 1
  integer, parameter, public :: STATUS_LINESEARCH = 2
  integer, parameter, public :: STATUS_BAD_INPUT = 3

  integer, parameter, public :: LINESEARCH_ARMIJO = 0
  integer, parameter, public :: LINESEARCH_WOLFE = 1
  integer, parameter, public :: LINESEARCH_STRONG_WOLFE = 2
  integer, parameter, public :: LINESEARCH_EXACT = 3
  integer, parameter, public :: LINESEARCH_INPUTFUN = 4

  type, public :: manifold_component
    integer :: kind = MANI_EUCLIDEAN
    integer :: n = 1
    integer :: m = 1
    integer :: p = 1
    integer :: numofmani = 1
    integer :: param_set = 1
  contains
    procedure :: block_length => component_block_length
    procedure :: total_length => component_total_length
  end type manifold_component

  type, public :: manifold_domain
    type(manifold_component), allocatable :: component(:)
  contains
    procedure :: length => domain_length
  end type manifold_domain

  abstract interface
    function line_search_callback(x, eta, initial_step, initial_slope) result(step)
      import dp
      real(dp), intent(in) :: x(:), eta(:), initial_step, initial_slope
      real(dp) :: step
    end function line_search_callback
  end interface

  type, public :: solver_options
    real(dp) :: tolerance = 1.0e-4_dp
    integer :: max_iteration = 1000
    integer :: debug = 0
    logical :: isconvex = .false.
    integer :: memory = 4
    integer :: max_linesearch = 30
    real(dp) :: initial_step = 1.0_dp
    real(dp) :: armijo = 1.0e-4_dp
    integer :: line_search = LINESEARCH_ARMIJO
    procedure(line_search_callback), pointer, nopass :: line_search_proc => null()
    real(dp) :: wolfe = 0.999_dp
    real(dp) :: max_step = 1000.0_dp
    real(dp) :: min_step = 1.0e-14_dp
    real(dp) :: eps_numerical_grad = 1.0e-6_dp
    real(dp) :: eps_numerical_hess = 1.0e-4_dp
    real(dp) :: trust_radius = 1.0_dp
    real(dp) :: max_trust_radius = 1000.0_dp
    real(dp) :: sr1_skip = 1.0e-8_dp
    real(dp) :: qn_nu = 1.0e-4_dp
    real(dp) :: qn_mu = 1.0_dp
    real(dp) :: broyden_phi = 1.0_dp
    character(len=16) :: cg_beta = 'HS'
  end type solver_options

  type, public :: solver_result
    real(dp), allocatable :: xopt(:)
    real(dp) :: fval = huge(1.0_dp)
    real(dp) :: normgf = huge(1.0_dp)
    real(dp) :: normgfgf0 = huge(1.0_dp)
    integer :: iter = 0
    integer :: num_obj_eval = 0
    integer :: num_grad_eval = 0
    integer :: nR = 0
    integer :: nV = 0
    integer :: nVp = 0
    integer :: nH = 0
    integer :: status = STATUS_BAD_INPUT
    character(len=160) :: message = ''
    real(dp), allocatable :: fun_series(:)
    real(dp), allocatable :: grad_series(:)
  end type solver_result

  abstract interface
    subroutine objective_callback(x, f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f
    end subroutine objective_callback

    subroutine gradient_callback(x, g)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
    end subroutine gradient_callback

    subroutine hessvec_callback(x, eta, hess_eta)
      import dp
      real(dp), intent(in) :: x(:), eta(:)
      real(dp), intent(out) :: hess_eta(:)
    end subroutine hessvec_callback
  end interface

  public :: objective_callback, gradient_callback, hessvec_callback, line_search_callback
  public :: make_component

contains

  pure function make_component(kind, n, m, p, numofmani, param_set) result(c)
    integer, intent(in) :: kind, n
    integer, intent(in), optional :: m, p, numofmani, param_set
    type(manifold_component) :: c
    c%kind = kind
    c%n = n
    if (present(m)) c%m = m
    if (present(p)) c%p = p
    if (present(numofmani)) c%numofmani = numofmani
    if (present(param_set)) c%param_set = param_set
    if (kind == MANI_SPHERE) then
      c%m = 1
      c%p = 1
    else if (kind == MANI_SPD .or. kind == MANI_ORTHGROUP) then
      c%m = n
      c%p = n
    end if
  end function make_component

  pure integer function component_block_length(self) result(nel)
    class(manifold_component), intent(in) :: self
    select case (self%kind)
    case (MANI_EUCLIDEAN)
      nel = self%n * self%m
    case (MANI_SPHERE)
      nel = self%n
    case (MANI_STIEFEL, MANI_GRASSMANN)
      nel = self%n * self%p
    case (MANI_SPD, MANI_ORTHGROUP)
      nel = self%n * self%n
    case (MANI_LOWRANK)
      nel = self%n * self%p + self%p * self%p + self%m * self%p
    case default
      nel = 0
    end select
  end function component_block_length

  pure integer function component_total_length(self) result(nel)
    class(manifold_component), intent(in) :: self
    nel = self%numofmani * self%block_length()
  end function component_total_length

  pure integer function domain_length(self) result(nel)
    class(manifold_domain), intent(in) :: self
    integer :: i
    nel = 0
    if (.not. allocated(self%component)) return
    do i = 1, size(self%component)
      nel = nel + self%component(i)%total_length()
    end do
  end function domain_length

end module manifoldoptim_types
