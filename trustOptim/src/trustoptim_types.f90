! SPDX-License-Identifier: MPL-2.0
module trustoptim_types
   use trustoptim_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: trust_success = 0
   integer, parameter, public :: trust_failure = -1
   integer, parameter, public :: trust_continue = -2
   integer, parameter, public :: trust_enoprog = 1
   integer, parameter, public :: trust_etolf = 3
   integer, parameter, public :: trust_etolx = 4
   integer, parameter, public :: trust_etolg = 5
   integer, parameter, public :: trust_emaxiter = 6
   integer, parameter, public :: trust_ebadlen = 7
   integer, parameter, public :: trust_enotsqr = 8
   integer, parameter, public :: trust_esing = 9
   integer, parameter, public :: trust_enomove = 10
   integer, parameter, public :: trust_failedcg = 11
   integer, parameter, public :: trust_moved = 12
   integer, parameter, public :: trust_contract = 13
   integer, parameter, public :: trust_expand = 14
   integer, parameter, public :: trust_unknown = 15
   integer, parameter, public :: trust_enegmove = 16

   integer, parameter, public :: trust_method_sr1 = 1
   integer, parameter, public :: trust_method_bfgs = 2
   integer, parameter, public :: trust_method_sparse = 3

   type, public :: sparse_symmetric_matrix
      integer :: n = 0
      integer :: nnz = 0
      integer, allocatable :: row(:)
      integer, allocatable :: col(:)
      real(dp), allocatable :: val(:)
   contains
      procedure :: clear => sparse_clear
      procedure :: matvec => sparse_matvec
      procedure :: to_dense => sparse_to_dense
      procedure :: set_from_dense => sparse_set_from_dense
   end type sparse_symmetric_matrix

   type, public :: trustoptim_control
      real(dp) :: start_trust_radius = 5.0_dp
      real(dp) :: stop_trust_radius = sqrt(epsilon(1.0_dp))
      real(dp) :: cg_tol = sqrt(epsilon(1.0_dp))
      real(dp) :: prec = sqrt(epsilon(1.0_dp))
      integer :: report_freq = 1
      integer :: report_level = 0
      integer :: report_precision = 6
      integer :: report_header_freq = 25
      integer :: maxit = 100
      real(dp) :: contract_factor = 0.5_dp
      real(dp) :: expand_factor = 3.0_dp
      real(dp) :: contract_threshold = 0.25_dp
      real(dp) :: expand_threshold_ap = 0.8_dp
      real(dp) :: expand_threshold_radius = 0.8_dp
      real(dp) :: function_scale_factor = 1.0_dp
      integer :: precond_refresh_freq = 1
      integer :: preconditioner = 0
      integer :: trust_iter = 2000
   end type trustoptim_control

   type, public :: trustoptim_result
      real(dp), allocatable :: solution(:)
      real(dp), allocatable :: gradient(:)
      type(sparse_symmetric_matrix) :: hessian
      real(dp) :: fval = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = trust_unknown
      real(dp) :: trust_radius = 0.0_dp
      integer :: nnz = 0
      integer :: method = 0
      integer :: hessian_update_method = 0
      integer :: last_cg_iterations = 0
      character(len=64) :: last_cg_reason = ''
   contains
      procedure :: status_message => trust_status_message
      procedure :: method_name => trust_method_name
   end type trustoptim_result

contains

   subroutine sparse_clear(this)
      class(sparse_symmetric_matrix), intent(inout) :: this
      this%n = 0
      this%nnz = 0
      if (allocated(this%row)) deallocate(this%row)
      if (allocated(this%col)) deallocate(this%col)
      if (allocated(this%val)) deallocate(this%val)
   end subroutine sparse_clear

   subroutine sparse_matvec(this, x, y)
      class(sparse_symmetric_matrix), intent(in) :: this
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(:)
      integer :: k, i, j

      y = 0.0_dp
      do k = 1, this%nnz
         i = this%row(k)
         j = this%col(k)
         y(i) = y(i) + this%val(k) * x(j)
         if (i /= j) y(j) = y(j) + this%val(k) * x(i)
      end do
   end subroutine sparse_matvec

   subroutine sparse_to_dense(this, a)
      class(sparse_symmetric_matrix), intent(in) :: this
      real(dp), intent(out) :: a(:,:)
      integer :: k, i, j

      a = 0.0_dp
      do k = 1, this%nnz
         i = this%row(k)
         j = this%col(k)
         a(i,j) = a(i,j) + this%val(k)
         if (i /= j) a(j,i) = a(j,i) + this%val(k)
      end do
   end subroutine sparse_to_dense

   subroutine sparse_set_from_dense(this, a, drop_tol)
      class(sparse_symmetric_matrix), intent(inout) :: this
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: drop_tol
      real(dp) :: tol
      integer :: n, i, j, k, count

      call this%clear()
      n = size(a,1)
      tol = 0.0_dp
      if (present(drop_tol)) tol = max(0.0_dp, drop_tol)
      count = 0
      do j = 1, n
         do i = j, n
            if (abs(a(i,j)) > tol) count = count + 1
         end do
      end do
      this%n = n
      this%nnz = count
      allocate(this%row(count), this%col(count), this%val(count))
      k = 0
      do j = 1, n
         do i = j, n
            if (abs(a(i,j)) > tol) then
               k = k + 1
               this%row(k) = i
               this%col(k) = j
               this%val(k) = a(i,j)
            end if
         end do
      end do
   end subroutine sparse_set_from_dense

   function trust_status_message(this) result(msg)
      class(trustoptim_result), intent(in) :: this
      character(len=64) :: msg

      select case (this%status)
      case (trust_success)
         msg = 'Success'
      case (trust_failure)
         msg = 'Unspecified failure'
      case (trust_continue)
         msg = 'Continuing'
      case (trust_enoprog)
         msg = 'Not making any progress'
      case (trust_etolf)
         msg = 'Cannot reach tolerance in F'
      case (trust_etolx)
         msg = 'Cannot reach tolerance in X'
      case (trust_etolg)
         msg = 'Radius below stop.trust.radius'
      case (trust_emaxiter)
         msg = 'Exceeded max iterations'
      case (trust_ebadlen)
         msg = 'Matrix/vector lengths not conformant'
      case (trust_enotsqr)
         msg = 'Matrix is not square'
      case (trust_esing)
         msg = 'Matrix is apparently singular'
      case (trust_enomove)
         msg = 'No proposed movement'
      case (trust_failedcg)
         msg = 'CG failed'
      case (trust_moved)
         msg = 'Continuing'
      case (trust_contract)
         msg = 'Continuing - TR contract'
      case (trust_expand)
         msg = 'Continuing - TR expand'
      case (trust_enegmove)
         msg = 'Negative predicted move'
      case default
         msg = 'Unspecified status'
      end select
   end function trust_status_message

   function trust_method_name(this) result(name)
      class(trustoptim_result), intent(in) :: this
      character(len=16) :: name
      select case (this%method)
      case (trust_method_sr1)
         name = 'SR1'
      case (trust_method_bfgs)
         name = 'BFGS'
      case (trust_method_sparse)
         name = 'Sparse'
      case default
         name = 'Unknown'
      end select
   end function trust_method_name

end module trustoptim_types
