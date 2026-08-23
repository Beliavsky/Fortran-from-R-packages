! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module rspectra_operators
    use rspectra_kinds, only: dp
    implicit none
    private

    type, abstract, public :: linear_operator
        integer :: nrow = 0
        integer :: ncol = 0
    contains
        procedure(prod_iface), deferred :: prod
        procedure(tprod_iface), deferred :: tprod
    end type linear_operator

    abstract interface
        subroutine prod_iface(self, x, y)
            import :: linear_operator, dp
            class(linear_operator), intent(inout) :: self
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
        end subroutine prod_iface
        subroutine tprod_iface(self, x, y)
            import :: linear_operator, dp
            class(linear_operator), intent(inout) :: self
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
        end subroutine tprod_iface
    end interface

    type, extends(linear_operator), public :: dense_operator
        real(dp), allocatable :: a(:,:)
    contains
        procedure :: prod => dense_prod
        procedure :: tprod => dense_tprod
    end type dense_operator

    type, extends(linear_operator), public :: csr_operator
        integer, allocatable :: row_ptr(:)
        integer, allocatable :: col_ind(:)
        real(dp), allocatable :: val(:)
    contains
        procedure :: prod => csr_prod
        procedure :: tprod => csr_tprod
    end type csr_operator

    public :: make_dense_operator, make_csr_operator, csr_to_dense, is_symmetric_dense

contains

    function make_dense_operator(a) result(op)
        real(dp), intent(in) :: a(:,:)
        type(dense_operator) :: op
        op%nrow = size(a, 1)
        op%ncol = size(a, 2)
        op%a = a
    end function make_dense_operator

    function make_csr_operator(nrow, ncol, row_ptr, col_ind, val) result(op)
        integer, intent(in) :: nrow, ncol
        integer, intent(in) :: row_ptr(:), col_ind(:)
        real(dp), intent(in) :: val(:)
        type(csr_operator) :: op
        op%nrow = nrow
        op%ncol = ncol
        op%row_ptr = row_ptr
        op%col_ind = col_ind
        op%val = val
    end function make_csr_operator

    subroutine dense_prod(self, x, y)
        class(dense_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        y = matmul(self%a, x)
    end subroutine dense_prod

    subroutine dense_tprod(self, x, y)
        class(dense_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        y = matmul(transpose(self%a), x)
    end subroutine dense_tprod

    subroutine csr_prod(self, x, y)
        class(csr_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: i, j
        y = 0.0_dp
        do i = 1, self%nrow
            do j = self%row_ptr(i), self%row_ptr(i + 1) - 1
                y(i) = y(i) + self%val(j) * x(self%col_ind(j))
            end do
        end do
    end subroutine csr_prod

    subroutine csr_tprod(self, x, y)
        class(csr_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: i, j
        y = 0.0_dp
        do i = 1, self%nrow
            do j = self%row_ptr(i), self%row_ptr(i + 1) - 1
                y(self%col_ind(j)) = y(self%col_ind(j)) + self%val(j) * x(i)
            end do
        end do
    end subroutine csr_tprod

    function csr_to_dense(op) result(a)
        type(csr_operator), intent(in) :: op
        real(dp), allocatable :: a(:,:)
        integer :: i, j
        allocate(a(op%nrow, op%ncol))
        a = 0.0_dp
        do i = 1, op%nrow
            do j = op%row_ptr(i), op%row_ptr(i + 1) - 1
                a(i, op%col_ind(j)) = a(i, op%col_ind(j)) + op%val(j)
            end do
        end do
    end function csr_to_dense

    logical function is_symmetric_dense(a, tol) result(ans)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in), optional :: tol
        real(dp) :: t
        if (size(a, 1) /= size(a, 2)) then
            ans = .false.
            return
        end if
        t = sqrt(epsilon(1.0_dp))
        if (present(tol)) t = tol
        ans = maxval(abs(a - transpose(a))) <= t * max(1.0_dp, maxval(abs(a)))
    end function is_symmetric_dense

end module rspectra_operators
