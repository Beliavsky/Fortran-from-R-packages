! Computational translation of the R package irlba 2.3.7.
! Upstream package and native core: GPL-3 / GPL-3-or-later.
! See ../LICENSE and ../UPSTREAM.md for provenance and copyright details.
module irlba_operator
  use irlba_kinds, only : dp
  use irlba_sparse, only : csc_matrix
  implicit none
  private

  type, abstract, public :: linear_operator
    integer :: nrow = 0
    integer :: ncol = 0
  contains
    procedure(matvec_iface), deferred :: matvec
    procedure(matvec_iface), deferred :: tmatvec
  end type linear_operator

  abstract interface
    subroutine matvec_iface(self, x, y)
      import :: linear_operator, dp
      class(linear_operator), intent(in) :: self
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(:)
    end subroutine matvec_iface
  end interface

  type, extends(linear_operator), public :: dense_operator
    real(dp), pointer :: a(:, :) => null()
  contains
    procedure :: matvec => dense_matvec
    procedure :: tmatvec => dense_tmatvec
  end type dense_operator

  type, extends(linear_operator), public :: csc_operator
    type(csc_matrix), pointer :: a => null()
  contains
    procedure :: matvec => sparse_matvec
    procedure :: tmatvec => sparse_tmatvec
  end type csc_operator

contains

  subroutine dense_matvec(self, x, y)
    class(dense_operator), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    y = matmul(self%a, x)
  end subroutine dense_matvec

  subroutine dense_tmatvec(self, x, y)
    class(dense_operator), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    y = matmul(transpose(self%a), x)
  end subroutine dense_tmatvec

  subroutine sparse_matvec(self, x, y)
    class(csc_operator), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    call self%a%matvec(x, y)
  end subroutine sparse_matvec

  subroutine sparse_tmatvec(self, x, y)
    class(csc_operator), intent(in) :: self
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    call self%a%tmatvec(x, y)
  end subroutine sparse_tmatvec

end module irlba_operator
