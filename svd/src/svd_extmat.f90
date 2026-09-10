! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the computational extmat API from R package svd 0.5.8.
module svd_extmat
   use svd_kinds, only : dp
   use rspectra, only : linear_operator
   implicit none
   private

   integer, parameter, public :: svd_success = 0
   integer, parameter, public :: svd_invalid_shape = -1001
   integer, parameter, public :: svd_invalid_callback = -1002

   abstract interface
      subroutine extmat_matvec_callback(context, x, y)
         import :: dp
         class(*), intent(inout) :: context !! User-owned mutable state used by the matrix-vector callback.
         real(dp), intent(in) :: x(:) !! Input vector; its size is the operator column count for forward products.
         real(dp), intent(out) :: y(:) !! Output vector; its size is the operator row count for forward products.
      end subroutine extmat_matvec_callback
   end interface

   type, public :: extmat_callbacks
      !! Procedure-pointer bundle corresponding to the R `mul` and `tmul` closures.
      procedure(extmat_matvec_callback), pointer, nopass :: mul => null()
      procedure(extmat_matvec_callback), pointer, nopass :: tmul => null()
   end type extmat_callbacks

   type, extends(linear_operator), public :: extmat_operator
      !! Matrix-free real operator compatible with RSpectra's linear-operator API.
      type(extmat_callbacks) :: callbacks
      class(*), pointer :: context => null()
   contains
      procedure :: prod => extmat_prod
      procedure :: tprod => extmat_tprod
   end type extmat_operator

   public :: ematmul
   public :: extmat_ncol
   public :: extmat_nrow
   public :: extmat_vector
   public :: is_extmat
   public :: make_extmat
   public :: materialize_extmat

contains

   function make_extmat(callbacks, context, nrow, ncol, info) result(op)
      !! Constructs a matrix-free operator from forward/transpose callbacks and persistent user state.
      type(extmat_callbacks), intent(in) :: callbacks !! Forward and transpose callback procedure pointers.
      class(*), target, intent(inout) :: context !! Caller-owned state that must outlive every use of the returned operator.
      integer, intent(in) :: nrow !! Number of matrix rows; must be positive.
      integer, intent(in) :: ncol !! Number of matrix columns; must be positive.
      integer, intent(out), optional :: info !! Zero on success; negative when dimensions or callbacks are invalid.
      type(extmat_operator) :: op
      integer :: status

      status = svd_success
      if (nrow <= 0 .or. ncol <= 0) status = svd_invalid_shape
      if (.not. associated(callbacks%mul) .or. .not. associated(callbacks%tmul)) then
         status = svd_invalid_callback
      end if

      if (status == svd_success) then
         op%nrow = nrow
         op%ncol = ncol
         op%callbacks = callbacks
         op%context => context
      end if
      if (present(info)) info = status
   end function make_extmat

   pure logical function is_extmat(x) result(answer)
      !! Reports whether a polymorphic RSpectra operator is a valid translated `extmat` object.
      class(linear_operator), intent(in) :: x !! Operator to inspect; other RSpectra operator types return false.

      answer = .false.
      select type (x)
      type is (extmat_operator)
         answer = associated(x%callbacks%mul) .and. associated(x%callbacks%tmul) .and. associated(x%context)
      class default
         answer = .false.
      end select
   end function is_extmat

   pure integer function extmat_nrow(x) result(nrow)
      !! Returns the declared row count of a translated external matrix.
      type(extmat_operator), intent(in) :: x !! Matrix-free operator whose row count is requested.

      nrow = x%nrow
   end function extmat_nrow

   pure integer function extmat_ncol(x) result(ncol)
      !! Returns the declared column count of a translated external matrix.
      type(extmat_operator), intent(in) :: x !! Matrix-free operator whose column count is requested.

      ncol = x%ncol
   end function extmat_ncol

   subroutine ematmul(emat, v, y, transposed, info)
      !! Applies an external matrix or its transpose to one real vector.
      type(extmat_operator), intent(inout) :: emat !! Matrix-free operator to apply.
      real(dp), intent(in) :: v(:) !! Input vector; length is `ncol` normally or `nrow` when transposed.
      real(dp), allocatable, intent(out) :: y(:) !! Allocated result vector with the complementary operator dimension.
      logical, intent(in), optional :: transposed !! When true, apply the transpose; default is false.
      integer, intent(out), optional :: info !! Zero on success; negative for an invalid object or vector length.
      logical :: use_transpose
      integer :: status

      status = svd_success
      use_transpose = .false.
      if (present(transposed)) use_transpose = transposed

      if (.not. is_extmat(emat)) then
         allocate(y(0))
         status = svd_invalid_callback
      else if (use_transpose) then
         if (size(v) /= emat%nrow) then
            allocate(y(0))
            status = svd_invalid_shape
         else
            allocate(y(emat%ncol))
            call emat%tprod(v, y)
         end if
      else
         if (size(v) /= emat%ncol) then
            allocate(y(0))
            status = svd_invalid_shape
         else
            allocate(y(emat%nrow))
            call emat%prod(v, y)
         end if
      end if

      if (present(info)) info = status
   end subroutine ematmul

   subroutine materialize_extmat(emat, a, info)
      !! Materializes a matrix-free operator by applying it to canonical basis vectors.
      type(extmat_operator), intent(inout) :: emat !! Matrix-free operator to materialize.
      real(dp), allocatable, intent(out) :: a(:, :) !! Allocated dense matrix with shape `(nrow, ncol)`.
      integer, intent(out), optional :: info !! Zero on success; negative when the operator is invalid.
      real(dp), allocatable :: basis(:)
      integer :: j, status

      status = svd_success
      if (.not. is_extmat(emat)) then
         allocate(a(0, 0))
         status = svd_invalid_callback
         if (present(info)) info = status
         return
      end if

      allocate(a(emat%nrow, emat%ncol), basis(emat%ncol))
      do j = 1, emat%ncol
         basis = 0.0_dp
         basis(j) = 1.0_dp
         call emat%prod(basis, a(:, j))
      end do
      if (present(info)) info = status
   end subroutine materialize_extmat

   subroutine extmat_vector(emat, values, info)
      !! Materializes an external matrix and returns its column-major numeric vector.
      type(extmat_operator), intent(inout) :: emat !! Matrix-free operator to flatten in Fortran/R column-major order.
      real(dp), allocatable, intent(out) :: values(:) !! Allocated vector containing all matrix entries column by column.
      integer, intent(out), optional :: info !! Zero on success; negative when the external matrix is invalid.
      real(dp), allocatable :: a(:, :)
      integer :: status

      call materialize_extmat(emat, a, status)
      if (status /= svd_success) then
         allocate(values(0))
      else
         allocate(values(size(a)))
         values = reshape(a, [size(a)])
      end if
      if (present(info)) info = status
   end subroutine extmat_vector

   subroutine extmat_prod(self, x, y)
      !! Implements the forward product required by `linear_operator`.
      class(extmat_operator), intent(inout) :: self !! External matrix whose forward callback is invoked.
      real(dp), intent(in) :: x(:) !! Input vector with size equal to `self%ncol`.
      real(dp), intent(out) :: y(:) !! Output vector with size equal to `self%nrow`.

      if (.not. associated(self%context) .or. .not. associated(self%callbacks%mul)) then
         error stop "svd_extmat: invalid forward callback"
      end if
      if (size(x) /= self%ncol .or. size(y) /= self%nrow) then
         error stop "svd_extmat: non-conformable forward product"
      end if
      call self%callbacks%mul(self%context, x, y)
   end subroutine extmat_prod

   subroutine extmat_tprod(self, x, y)
      !! Implements the transpose product required by `linear_operator`.
      class(extmat_operator), intent(inout) :: self !! External matrix whose transpose callback is invoked.
      real(dp), intent(in) :: x(:) !! Input vector with size equal to `self%nrow`.
      real(dp), intent(out) :: y(:) !! Output vector with size equal to `self%ncol`.

      if (.not. associated(self%context) .or. .not. associated(self%callbacks%tmul)) then
         error stop "svd_extmat: invalid transpose callback"
      end if
      if (size(x) /= self%nrow .or. size(y) /= self%ncol) then
         error stop "svd_extmat: non-conformable transpose product"
      end if
      call self%callbacks%tmul(self%context, x, y)
   end subroutine extmat_tprod

end module svd_extmat
