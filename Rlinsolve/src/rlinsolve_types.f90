module rlinsolve_types
    use rlinsolve_kinds, only : dp
    implicit none
    private

    type, public :: sparse_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        integer, allocatable :: row(:)
        integer, allocatable :: col(:)
        real(dp), allocatable :: val(:)
    contains
        procedure :: nnz => sparse_nnz
    end type sparse_matrix

    type, public :: solver_result
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: errors(:)
        integer :: iter = 0
        integer :: flag = 0
        logical :: converged = .false.
    end type solver_result

    type, public :: solver_matrix_result
        real(dp), allocatable :: x(:, :)
        integer, allocatable :: iter(:)
        integer, allocatable :: flag(:)
        logical, allocatable :: converged(:)
        type(solver_result), allocatable :: column(:)
    end type solver_matrix_result

contains

    pure integer function sparse_nnz(self) result(n)
        class(sparse_matrix), intent(in) :: self !! Sparse matrix whose stored-entry count is requested.

        if (allocated(self%val)) then
            n = size(self%val)
        else
            n = 0
        end if
    end function sparse_nnz

end module rlinsolve_types
