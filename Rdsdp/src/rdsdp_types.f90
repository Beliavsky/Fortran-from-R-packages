! Modern Fortran data model for Rdsdp/DSDP5.
! DSDP copyright/license: see licenses/DSDP-LICENSE.
module rdsdp_types
   use rdsdp_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: dsdp_sdp_block = 1
   integer, parameter, public :: dsdp_lp_block  = 2

   integer, parameter, public :: dsdp_data_dense   = 1
   integer, parameter, public :: dsdp_data_sparse  = 2
   integer, parameter, public :: dsdp_data_lowrank = 3

   integer, parameter, public :: dsdp_pdfeasible = 1
   integer, parameter, public :: dsdp_unknown    = 0
   integer, parameter, public :: dsdp_unbounded  = 3
   integer, parameter, public :: dsdp_infeasible = 4

   integer, parameter, public :: dsdp_converged = 0
   integer, parameter, public :: dsdp_max_iterations = 1
   integer, parameter, public :: dsdp_linesearch_failure = 2
   integer, parameter, public :: dsdp_singular_schur = 3
   integer, parameter, public :: dsdp_numeric_failure = 4

   type, public :: dsdp_sparse_sym
      integer :: n = 0
      integer :: nnz = 0
      ! Full symmetric COO representation.  Off-diagonal entries appear twice.
      integer, allocatable :: row(:)
      integer, allocatable :: col(:)
      real(dp), allocatable :: val(:)
   end type dsdp_sparse_sym

   type, public :: dsdp_lowrank_sym
      integer :: n = 0
      integer :: rank = 0
      ! A = sum_j coeff(j) * vec(:,j) * transpose(vec(:,j)).
      real(dp), allocatable :: coeff(:)
      real(dp), allocatable :: vec(:,:)
   end type dsdp_lowrank_sym

   type, public :: dsdp_block
      integer :: category = dsdp_sdp_block
      integer :: n = 0
      ! Dense storage retained for compatibility/reference mode.
      real(dp), allocatable :: c(:,:)
      real(dp), allocatable :: a(:,:,:)
      ! v0.2.0 abstract SDP data storage, modelled after DSDPDataMat.
      integer :: c_storage = dsdp_data_dense
      integer, allocatable :: a_storage(:)
      type(dsdp_sparse_sym) :: c_sparse
      type(dsdp_sparse_sym), allocatable :: a_sparse(:)
      type(dsdp_lowrank_sym) :: c_lowrank
      type(dsdp_lowrank_sym), allocatable :: a_lowrank(:)
      ! LP/nonnegative block: cdiag(n), adiag(n,m)
      real(dp), allocatable :: cdiag(:)
      real(dp), allocatable :: adiag(:,:)
   end type dsdp_block

   type, public :: dsdp_problem
      integer :: m = 0
      type(dsdp_block), allocatable :: block(:)
      real(dp), allocatable :: b(:)
   contains
      procedure :: nblocks => dsdp_problem_nblocks
      procedure :: barrier_dimension => dsdp_problem_barrier_dimension
   end type dsdp_problem

   type, public :: dsdp_primal_block
      integer :: category = dsdp_sdp_block
      integer :: n = 0
      real(dp), allocatable :: x(:,:)
      real(dp), allocatable :: xdiag(:)
   end type dsdp_primal_block

   type, public :: dsdp_control
      real(dp) :: gaptol = 1.0e-7_dp
      real(dp) :: pinfeastol = 1.0e-7_dp
      real(dp) :: rtol = 1.0e-8_dp
      integer :: maxiter = 80
      integer :: max_correctors = 20
      integer :: max_linesearch = 40
      real(dp) :: mu_factor = 0.25_dp
      real(dp) :: newton_tol = 5.0e-8_dp
      real(dp) :: armijo = 1.0e-4_dp
      real(dp) :: penalty = 1.0e6_dp
      real(dp) :: schur_regularization = 1.0e-12_dp
      integer :: print = 0
      ! v0.2.0 controls.
      logical :: use_sparse_data = .true.
      real(dp) :: sparse_density_threshold = 0.20_dp
      logical :: use_cg = .false.
      integer :: cg_threshold = 200
      integer :: cg_maxiter = 80
      real(dp) :: cg_tol = 1.0e-10_dp
      logical :: cg_fallback_direct = .true.
      ! v0.3.0 Schur-system controls.
      logical :: cg_matrix_free = .true.
      logical :: use_sparse_schur_factor = .true.
      integer :: sparse_schur_threshold = 80
      real(dp) :: sparse_schur_density_limit = 0.20_dp
      real(dp) :: sparse_schur_drop_tol = 0.0_dp
   end type dsdp_control

   type, public :: dsdp_solution
      type(dsdp_primal_block), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      integer :: stype = dsdp_unknown
      integer :: status = dsdp_numeric_failure
      integer :: iterations = 0
      integer :: newton_steps = 0
      integer :: line_search_backtracks = 0
      real(dp) :: dobj = 0.0_dp
      real(dp) :: pobj = 0.0_dp
      real(dp) :: r = huge(1.0_dp)
      real(dp) :: mu = huge(1.0_dp)
      real(dp) :: pstep = 0.0_dp
      real(dp) :: dstep = 0.0_dp
      real(dp) :: pnorm = huge(1.0_dp)
      real(dp) :: pinfeas = huge(1.0_dp)
      real(dp) :: relgap = huge(1.0_dp)
      real(dp) :: min_dual_eig = -huge(1.0_dp)
      ! v0.2.0 diagnostics.
      integer :: sdp_data_nnz = 0
      integer :: schur_assemblies = 0
      integer :: sparse_pair_evals = 0
      integer :: lowrank_pair_evals = 0
      integer :: dense_pair_evals = 0
      integer :: direct_schur_solves = 0
      integer :: cg_solves = 0
      integer :: cg_iterations = 0
      real(dp) :: schur_assembly_time = 0.0_dp
      real(dp) :: schur_solve_time = 0.0_dp
      ! v0.3.0 diagnostics.
      integer :: matrix_free_cg_solves = 0
      integer :: matrix_free_matvecs = 0
      integer :: sparse_factor_solves = 0
      integer :: sparse_factor_fallbacks = 0
      integer :: sparse_symbolic_analyses = 0
      integer :: sparse_numeric_factorizations = 0
      integer :: schur_matrix_nnz = 0
      integer :: schur_factor_nnz = 0
   end type dsdp_solution

contains

   pure integer function dsdp_problem_nblocks(this) result(n)
      class(dsdp_problem), intent(in) :: this
      if (allocated(this%block)) then
         n = size(this%block)
      else
         n = 0
      end if
   end function dsdp_problem_nblocks

   pure integer function dsdp_problem_barrier_dimension(this) result(nu)
      class(dsdp_problem), intent(in) :: this
      integer :: k
      nu = 1 ! residual variable r > 0
      if (.not. allocated(this%block)) return
      do k = 1, size(this%block)
         nu = nu + this%block(k)%n
      end do
   end function dsdp_problem_barrier_dimension

end module rdsdp_types
