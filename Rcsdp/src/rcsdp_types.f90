! Modern Fortran data model for Rcsdp/CSDP 6.1.1.
! Original software is distributed under the Common Public License 1.0.
module rcsdp_types
   use rcsdp_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: csdp_matrix = 1
   integer, parameter, public :: csdp_diag   = 2

   integer, parameter, public :: csdp_success          = 0
   integer, parameter, public :: csdp_primal_infeas    = 1
   integer, parameter, public :: csdp_dual_infeas      = 2
   integer, parameter, public :: csdp_partial_success  = 3
   integer, parameter, public :: csdp_max_iterations   = 4
   integer, parameter, public :: csdp_primal_edge      = 5
   integer, parameter, public :: csdp_dual_edge        = 6
   integer, parameter, public :: csdp_no_progress      = 7
   integer, parameter, public :: csdp_singular         = 8
   integer, parameter, public :: csdp_numeric_failure  = 9

   type, public :: csdp_block
      integer :: category = csdp_matrix
      integer :: n = 0
      real(dp), allocatable :: mat(:,:)
      real(dp), allocatable :: diag(:)
   end type csdp_block

   type, public :: csdp_block_matrix
      type(csdp_block), allocatable :: block(:)
   contains
      procedure :: order => block_matrix_order
      procedure :: nblocks => block_matrix_nblocks
   end type csdp_block_matrix

   type, public :: csdp_sparse_block
      integer :: blocknum = 0
      integer :: n = 0
      integer, allocatable :: i(:)
      integer, allocatable :: j(:)
      real(dp), allocatable :: v(:)
   contains
      procedure :: nnz => sparse_block_nnz
   end type csdp_sparse_block

   type, public :: csdp_constraint
      type(csdp_sparse_block), allocatable :: block(:)
   end type csdp_constraint

   type, public :: csdp_problem
      type(csdp_block_matrix) :: c
      type(csdp_constraint), allocatable :: a(:)
      real(dp), allocatable :: b(:)
      real(dp) :: constant_offset = 0.0_dp
   contains
      procedure :: order => problem_order
      procedure :: nconstraints => problem_nconstraints
   end type csdp_problem

   type, public :: csdp_control
      real(dp) :: axtol = 1.0e-8_dp
      real(dp) :: atytol = 1.0e-8_dp
      real(dp) :: objtol = 1.0e-8_dp
      real(dp) :: pinftol = 1.0e8_dp
      real(dp) :: dinftol = 1.0e8_dp
      integer :: maxiter = 100
      real(dp) :: minstepfrac = 0.90_dp
      real(dp) :: maxstepfrac = 0.97_dp
      real(dp) :: minstepp = 1.0e-8_dp
      real(dp) :: minstepd = 1.0e-8_dp
      logical :: usexzgap = .true.
      logical :: tweakgap = .false.
      logical :: affine = .false.
      integer :: printlevel = 1
      logical :: perturbobj = .true.
      logical :: fastmode = .false.
      logical :: use_sparse_schur = .true.
      logical :: use_fill_products = .true.
      real(dp) :: fill_density_limit = 0.01_dp
      logical :: use_schur_scaling = .true.
      logical :: use_lanczos_linesearch = .true.
      integer :: lanczos_threshold = 180
      integer :: lanczos_iterations = 30
   end type csdp_control

   type, public :: csdp_solution
      type(csdp_block_matrix) :: x
      type(csdp_block_matrix) :: z
      real(dp), allocatable :: y(:)
      real(dp) :: pobj = 0.0_dp
      real(dp) :: dobj = 0.0_dp
      real(dp) :: relgap = huge(1.0_dp)
      real(dp) :: pinfeas = huge(1.0_dp)
      real(dp) :: dinfeas = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = csdp_numeric_failure
      integer :: constraint_nnz = 0
      integer :: sparse_constraint_blocks = 0
      integer :: dense_constraint_blocks = 0
      integer :: schur_assemblies = 0
      integer :: schur_sparse_pairs = 0
      integer :: schur_dense_products = 0
      integer :: schur_refinements = 0
      real(dp) :: schur_seconds = 0.0_dp
      integer :: fill_nnz = 0
      integer :: fill_full_entries = 0
      integer :: fill_sparse_products = 0
      integer :: fill_dense_products = 0
      integer :: lanczos_linesearches = 0
      real(dp) :: schur_diagadd = 0.0_dp
   end type csdp_solution

contains

   pure integer function block_matrix_order(this) result(n)
      class(csdp_block_matrix), intent(in) :: this
      integer :: k
      n = 0
      if (.not. allocated(this%block)) return
      do k = 1, size(this%block)
         n = n + this%block(k)%n
      end do
   end function block_matrix_order

   pure integer function block_matrix_nblocks(this) result(n)
      class(csdp_block_matrix), intent(in) :: this
      if (allocated(this%block)) then
         n = size(this%block)
      else
         n = 0
      end if
   end function block_matrix_nblocks

   pure integer function sparse_block_nnz(this) result(n)
      class(csdp_sparse_block), intent(in) :: this
      if (allocated(this%v)) then
         n = size(this%v)
      else
         n = 0
      end if
   end function sparse_block_nnz

   pure integer function problem_order(this) result(n)
      class(csdp_problem), intent(in) :: this
      n = this%c%order()
   end function problem_order

   pure integer function problem_nconstraints(this) result(n)
      class(csdp_problem), intent(in) :: this
      if (allocated(this%b)) then
         n = size(this%b)
      else
         n = 0
      end if
   end function problem_nconstraints

end module rcsdp_types
