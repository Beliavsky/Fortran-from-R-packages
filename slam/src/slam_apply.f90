! SPDX-License-Identifier: GPL-2.0-only
!
! Statically typed analogues of slam's apply/crossapply routines.
module slam_apply
    use slam_kinds, only : dp
    use slam_stm, only : simple_triplet_matrix
    implicit none
    private

    abstract interface
        function vector_function(x) result(value)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function vector_function
        function two_vector_function(x,y) result(value)
            import dp
            real(dp), intent(in) :: x(:),y(:)
            real(dp) :: value
        end function two_vector_function
    end interface

    public :: colapply_simple_triplet_matrix
    public :: rowapply_simple_triplet_matrix
    public :: crossapply_simple_triplet_matrix
    public :: tcrossapply_simple_triplet_matrix
    public :: vector_function, two_vector_function

contains

    function colapply_simple_triplet_matrix(x, fun) result(out)
        type(simple_triplet_matrix), intent(in) :: x
        procedure(vector_function) :: fun
        real(dp), allocatable :: out(:)
        real(dp), allocatable :: work(:)
        integer :: c,k

        allocate(out(x%ncol),work(x%nrow))
        do c=1,x%ncol
            work=0.0_dp
            do k=1,x%nnz()
                if(x%j(k)==c) work(x%i(k))=x%v(k)
            end do
            out(c)=fun(work)
        end do
    end function colapply_simple_triplet_matrix

    function rowapply_simple_triplet_matrix(x,fun) result(out)
        type(simple_triplet_matrix), intent(in) :: x
        procedure(vector_function) :: fun
        real(dp), allocatable :: out(:)
        type(simple_triplet_matrix) :: tx
        tx=x%transpose()
        out=colapply_simple_triplet_matrix(tx,fun)
    end function rowapply_simple_triplet_matrix

    function crossapply_simple_triplet_matrix(x,y,fun) result(out)
        type(simple_triplet_matrix), intent(in) :: x,y
        procedure(two_vector_function) :: fun
        real(dp), allocatable :: out(:,:)
        real(dp), allocatable :: xd(:,:),yd(:,:)
        integer :: a,b

        if(x%nrow/=y%nrow) error stop "crossapply_simple_triplet_matrix: row counts do not conform"
        xd=x%to_dense(); yd=y%to_dense()
        allocate(out(x%ncol,y%ncol))
        do b=1,y%ncol
            do a=1,x%ncol
                out(a,b)=fun(xd(:,a),yd(:,b))
            end do
        end do
    end function crossapply_simple_triplet_matrix

    function tcrossapply_simple_triplet_matrix(x,y,fun) result(out)
        type(simple_triplet_matrix), intent(in) :: x,y
        procedure(two_vector_function) :: fun
        real(dp), allocatable :: out(:,:)
        type(simple_triplet_matrix) :: tx,ty
        tx=x%transpose(); ty=y%transpose()
        out=crossapply_simple_triplet_matrix(tx,ty,fun)
    end function tcrossapply_simple_triplet_matrix

end module slam_apply
