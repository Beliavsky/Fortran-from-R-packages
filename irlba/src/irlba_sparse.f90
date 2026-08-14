! Computational translation of the R package irlba 2.3.7.
! Upstream package and native core: GPL-3 / GPL-3-or-later.
! See ../LICENSE and ../UPSTREAM.md for provenance and copyright details.
module irlba_sparse
  use irlba_kinds, only : dp
  implicit none
  private

  type, public :: csc_matrix
    integer :: nrow = 0
    integer :: ncol = 0
    integer, allocatable :: col_ptr(:)
    integer, allocatable :: row_ind(:)
    real(dp), allocatable :: value(:)
  contains
    procedure :: matvec => csc_matvec
    procedure :: tmatvec => csc_tmatvec
    procedure :: matmat => csc_matmat
    procedure :: tmatmat => csc_tmatmat
    procedure :: to_dense => csc_to_dense
  end type csc_matrix

  public :: csc_from_dense

contains

  function csc_from_dense(a, zero_tol) result(s)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: zero_tol
    type(csc_matrix) :: s
    real(dp) :: ztol
    integer :: i, j, k, nnz

    ztol = 0.0_dp
    if (present(zero_tol)) ztol = max(0.0_dp, zero_tol)
    s%nrow = size(a, 1)
    s%ncol = size(a, 2)
    nnz = count(abs(a) > ztol)
    allocate(s%col_ptr(s%ncol + 1), s%row_ind(nnz), s%value(nnz))
    k = 1
    s%col_ptr(1) = 1
    do j = 1, s%ncol
      do i = 1, s%nrow
        if (abs(a(i, j)) > ztol) then
          s%row_ind(k) = i
          s%value(k) = a(i, j)
          k = k + 1
        end if
      end do
      s%col_ptr(j + 1) = k
    end do
  end function csc_from_dense

  subroutine csc_matvec(self, x, y)
    class(csc_matrix), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: j, k

    if (size(x) /= self%ncol .or. size(y) /= self%nrow) error stop "csc_matvec: dimension mismatch"
    y = 0.0_dp
    do j = 1, self%ncol
      do k = self%col_ptr(j), self%col_ptr(j + 1) - 1
        y(self%row_ind(k)) = y(self%row_ind(k)) + self%value(k) * x(j)
      end do
    end do
  end subroutine csc_matvec

  subroutine csc_tmatvec(self, x, y)
    class(csc_matrix), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: j, k

    if (size(x) /= self%nrow .or. size(y) /= self%ncol) error stop "csc_tmatvec: dimension mismatch"
    y = 0.0_dp
    do j = 1, self%ncol
      do k = self%col_ptr(j), self%col_ptr(j + 1) - 1
        y(j) = y(j) + self%value(k) * x(self%row_ind(k))
      end do
    end do
  end subroutine csc_tmatvec


  subroutine csc_matmat(self, x, y)
    class(csc_matrix), intent(in) :: self
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: y(:, :)
    integer :: j, k

    if (size(x, 1) /= self%ncol .or. size(y, 1) /= self%nrow .or. &
        size(y, 2) /= size(x, 2)) error stop "csc_matmat: dimension mismatch"
    y = 0.0_dp
    do j = 1, self%ncol
      do k = self%col_ptr(j), self%col_ptr(j + 1) - 1
        y(self%row_ind(k), :) = y(self%row_ind(k), :) + self%value(k) * x(j, :)
      end do
    end do
  end subroutine csc_matmat

  subroutine csc_tmatmat(self, x, y)
    class(csc_matrix), intent(in) :: self
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(out) :: y(:, :)
    integer :: j, k

    if (size(x, 1) /= self%nrow .or. size(y, 1) /= self%ncol .or. &
        size(y, 2) /= size(x, 2)) error stop "csc_tmatmat: dimension mismatch"
    y = 0.0_dp
    do j = 1, self%ncol
      do k = self%col_ptr(j), self%col_ptr(j + 1) - 1
        y(j, :) = y(j, :) + self%value(k) * x(self%row_ind(k), :)
      end do
    end do
  end subroutine csc_tmatmat

  function csc_to_dense(self) result(a)
    class(csc_matrix), intent(in) :: self
    real(dp), allocatable :: a(:, :)
    integer :: j, k

    allocate(a(self%nrow, self%ncol))
    a = 0.0_dp
    do j = 1, self%ncol
      do k = self%col_ptr(j), self%col_ptr(j + 1) - 1
        a(self%row_ind(k), j) = a(self%row_ind(k), j) + self%value(k)
      end do
    end do
  end function csc_to_dense

end module irlba_sparse
