! SPDX-License-Identifier: GPL-3.0-only
module ecos_matrixextra_adapter
    use matrix_sparse, only : matrix_csc => csc_matrix
    use matrixextra_types, only : coo_matrix
    use matrixextra_conversions, only : matrix_csc_from_coo => csc_from_coo
    use ecos, only : dp, ecos_csc_matrix, ecos_problem, ecos_dims, setup_problem_csc
    implicit none
    private
    public :: ecos_csc_from_matrix, ecos_csc_from_coo, setup_problem_matrixextra

contains

    subroutine ecos_csc_from_matrix(a,b)
        type(matrix_csc), intent(in) :: a
        type(ecos_csc_matrix), intent(out) :: b
        b%nrow=a%nrow; b%ncol=a%ncol
        allocate(b%colptr(size(a%col_ptr)),b%rowind(size(a%row_ind)),b%values(size(a%values)))
        b%colptr=a%col_ptr; b%rowind=a%row_ind; b%values=a%values
    end subroutine ecos_csc_from_matrix

    subroutine ecos_csc_from_coo(a,b,info)
        type(coo_matrix), intent(in) :: a
        type(ecos_csc_matrix), intent(out) :: b
        integer, intent(out), optional :: info
        type(matrix_csc) :: m
        integer :: ierr
        call matrix_csc_from_coo(a,m,ierr)
        if(ierr==0) call ecos_csc_from_matrix(m,b)
        if(present(info)) info=ierr
    end subroutine ecos_csc_from_coo

    subroutine setup_problem_matrixextra(prob,c,g,h,dims,a,b,bool_vars,int_vars,ierr)
        type(ecos_problem), intent(out) :: prob
        real(dp), intent(in) :: c(:),h(:)
        type(matrix_csc), intent(in) :: g
        type(ecos_dims), intent(in) :: dims
        type(matrix_csc), intent(in), optional :: a
        real(dp), intent(in), optional :: b(:)
        integer, intent(in), optional :: bool_vars(:),int_vars(:)
        integer, intent(out), optional :: ierr
        type(ecos_csc_matrix) :: ge,ae
        integer :: ie
        call ecos_csc_from_matrix(g,ge)
        if(present(a)) then
            call ecos_csc_from_matrix(a,ae)
            call setup_problem_csc(prob,c,ge,h,dims,ae,b,bool_vars,int_vars,ie)
        else
            call setup_problem_csc(prob,c,ge,h,dims,bool_vars=bool_vars,int_vars=int_vars,ierr=ie)
        end if
        if(present(ierr)) ierr=ie
    end subroutine setup_problem_matrixextra

end module ecos_matrixextra_adapter
