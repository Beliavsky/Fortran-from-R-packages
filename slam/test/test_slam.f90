program test_slam
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use slam
    implicit none

    call test_stm_core()
    call test_stm_arithmetic()
    call test_stm_indexing()
    call test_stm_nan()
    call test_ssa_core()
    call test_rollup()
    call test_apply()
    call test_io()
    print '(a)', 'All slam-fortran tests passed.'

contains

    subroutine test_stm_core()
        real(dp) :: a(3,2), b(3,2)
        real(dp), allocatable :: z(:,:),s(:)
        type(simple_triplet_matrix) :: x,y,t

        a=reshape([1.0_dp,0.0_dp,0.0_dp, 2.0_dp,1.0_dp,0.0_dp],[3,2])
        b=reshape([1.0_dp,2.0_dp,3.0_dp, 4.0_dp,5.0_dp,6.0_dp],[3,2])
        x=dense_to_stm(a)
        y=dense_to_stm(b)
        call assert_true(x%valid(),'stm valid')
        call assert_true(x%nnz()==3,'stm nnz')
        call assert_mat_close(x%to_dense(),a,'dense roundtrip')
        s=x%row_sums()
        call assert_vec_close(s,[3.0_dp,1.0_dp,0.0_dp],'row sums')
        s=x%col_sums()
        call assert_vec_close(s,[1.0_dp,3.0_dp],'col sums')
        s=x%row_means()
        call assert_vec_close(s,[1.5_dp,0.5_dp,0.0_dp],'row means')
        s=x%col_norms(2.0_dp)
        call assert_vec_close(s,[1.0_dp,sqrt(5.0_dp)],'col norms')
        t=x%transpose()
        call assert_mat_close(t%to_dense(),transpose(a),'transpose')
        t=x%reshape(2,3)
        call assert_mat_close(t%to_dense(),reshape(a,[2,3]),'stm reshape')
        z=stm_tcrossprod(x,y)
        call assert_mat_close(z,matmul(a,transpose(b)),'tcrossprod')
        z=stm_crossprod(x,y)
        call assert_mat_close(z,matmul(transpose(a),b),'crossprod')
        z=stm_matprod(x,y%transpose())
        call assert_mat_close(z,matmul(a,transpose(b)),'matprod')
        z=stm_matprod_dense(x,transpose(b))
        call assert_mat_close(z,matmul(a,transpose(b)),'stm dense matprod')
        z=dense_matprod_stm(transpose(b),x)
        call assert_mat_close(z,matmul(transpose(b),a),'dense stm matprod')
    end subroutine test_stm_core

    subroutine test_stm_arithmetic()
        real(dp) :: a(2,3),b(2,3)
        type(simple_triplet_matrix) :: x,y,z
        a=reshape([1.0_dp,0.0_dp,2.0_dp,0.0_dp,0.0_dp,3.0_dp],[2,3])
        b=reshape([0.0_dp,4.0_dp,2.0_dp,0.0_dp,5.0_dp,0.0_dp],[2,3])
        x=dense_to_stm(a); y=dense_to_stm(b)
        z=x+y
        call assert_mat_close(z%to_dense(),a+b,'stm add')
        z=x-y
        call assert_mat_close(z%to_dense(),a-b,'stm subtract')
        z=x*y
        call assert_mat_close(z%to_dense(),a*b,'stm hadamard')
        z=2.0_dp*x
        call assert_mat_close(z%to_dense(),2.0_dp*a,'stm scalar multiply')
        z=x**2.0_dp
        call assert_mat_close(z%to_dense(),a**2,'stm power')
        z=stm_row_scale(x,[2.0_dp,3.0_dp])
        b=a; b(1,:)=2*b(1,:); b(2,:)=3*b(2,:)
        call assert_mat_close(z%to_dense(),b,'row scale')
    end subroutine test_stm_arithmetic

    subroutine test_stm_indexing()
        real(dp) :: a(3,3),expected(2,2)
        type(simple_triplet_matrix) :: x,y
        a=reshape([1.0_dp,0.0_dp,2.0_dp, 0.0_dp,3.0_dp,0.0_dp, 4.0_dp,0.0_dp,5.0_dp],[3,3])
        x=dense_to_stm(a)
        y=x%extract([3,1],[3,2])
        expected=a([3,1],[3,2])
        call assert_mat_close(y%to_dense(),expected,'extract')
        call stm_set(x,[1,2,3],[1,2,3],[0.0_dp,7.0_dp])
        a(1,1)=0.0_dp; a(2,2)=7.0_dp; a(3,3)=0.0_dp
        call assert_mat_close(x%to_dense(),a,'subassign and recycling')
    end subroutine test_stm_indexing

    subroutine test_stm_nan()
        real(dp) :: a(2,2),nan
        real(dp), allocatable :: s(:),z(:,:)
        type(simple_triplet_matrix) :: x
        nan=ieee_value(0.0_dp,ieee_quiet_nan)
        a=reshape([1.0_dp,0.0_dp,nan,2.0_dp],[2,2])
        x=dense_to_stm(a)
        s=x%row_sums()
        call assert_true(ieee_is_nan(s(1)),'nan row sum propagates')
        s=x%row_sums(.true.)
        call assert_vec_close(s,[1.0_dp,2.0_dp],'nan row sum removed')
        z=stm_tcrossprod(x)
        call assert_mat_same_nan(z,matmul(a,transpose(a)),'nan tcrossprod dense fallback')
    end subroutine test_stm_nan

    subroutine test_ssa_core()
        integer :: dim3(3),coords(3,3)
        real(dp) :: flat(12)
        real(dp), allocatable :: back(:)
        type(simple_sparse_array) :: x,y,z
        type(simple_triplet_matrix) :: m

        dim3=[3,2,2]
        flat=0.0_dp; flat([1,5,12])=[2.0_dp,3.0_dp,4.0_dp]
        x=dense_flat_to_ssa(flat,dim3)
        call assert_true(x%valid(),'ssa valid')
        call assert_true(x%nnz()==3,'ssa nnz')
        back=x%to_dense_flat()
        call assert_vec_close(back,flat,'ssa dense roundtrip')
        y=x%permute([3,2,1])
        call assert_true(all(y%dim==[2,2,3]),'ssa permute dims')
        y=x%reshape([2,3,2])
        call assert_vec_close(y%to_dense_flat(),flat,'ssa reshape')
        y=extend_simple_sparse_array(x,-2)
        call assert_true(all(y%dim==[3,1,2,2]),'ssa extend')

        coords=reshape([1,2,3, 1,2,1, 1,1,2],[3,3])
        ! reshape above is column-major: rows are (1,1,1),(2,2,1),(3,1,2)
        z=make_ssa(coords,[1.0_dp,2.0_dp,3.0_dp],dim3)
        call ssa_set(z,coords(1:1,:),[5.0_dp])
        call assert_true(z%nnz()==3,'ssa set existing')
        y=abind_simple_sparse_array(x,x,-1)
        call assert_true(all(y%dim==[2,3,2,2]),'ssa abind new margin')

        m=ssa_to_stm(make_ssa(reshape([1,2,1,2],[2,2]),[2.0_dp,3.0_dp],[2,2]))
        call assert_mat_close(m%to_dense(),reshape([2.0_dp,0.0_dp,0.0_dp,3.0_dp],[2,2]),'ssa to stm')
        z=stm_to_ssa(m)
        call assert_true(z%valid(),'stm to ssa')
    end subroutine test_ssa_core

    subroutine test_rollup()
        real(dp) :: a(2,3),nan,expected(2,2)
        type(simple_triplet_matrix) :: x,y,m
        type(simple_sparse_array) :: s,t
        nan=ieee_value(0.0_dp,ieee_quiet_nan)
        a=reshape([1.0_dp,0.0_dp,0.0_dp,2.0_dp,1.0_dp,nan],[2,3])
        x=dense_to_stm(a)
        y=rollup_stm_sum(x,2,[1,2,1],na_rm=.true.,reduce_zeros=.true.)
        expected=reshape([2.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2])
        call assert_mat_close(y%to_dense(),expected,'stm rollup')
        s=stm_to_ssa(x)
        t=rollup_ssa_sum(s,2,[1,2,1],na_rm=.true.,reduce_zeros=.true.)
        m=ssa_to_stm(t)
        call assert_mat_close(m%to_dense(),expected,'ssa rollup')
    end subroutine test_rollup

    subroutine test_apply()
        real(dp) :: a(3,2)
        real(dp), allocatable :: s(:),z(:,:)
        type(simple_triplet_matrix) :: x
        a=reshape([1.0_dp,0.0_dp,2.0_dp, 3.0_dp,4.0_dp,0.0_dp],[3,2])
        x=dense_to_stm(a)
        s=colapply_simple_triplet_matrix(x,sum_fun)
        call assert_vec_close(s,[3.0_dp,7.0_dp],'col apply')
        s=rowapply_simple_triplet_matrix(x,sum_fun)
        call assert_vec_close(s,[4.0_dp,4.0_dp,2.0_dp],'row apply')
        z=crossapply_simple_triplet_matrix(x,x,dot_fun)
        call assert_mat_close(z,matmul(transpose(a),a),'cross apply')
    end subroutine test_apply

    subroutine test_io()
        real(dp) :: a(3,2)
        type(simple_triplet_matrix) :: x,y
        character(len=*), parameter :: cluto='test_slam_tmp.cluto'
        character(len=*), parameter :: mc='test_slam_tmp_mc'
        a=reshape([1.0_dp,0.0_dp,2.0_dp, 0.0_dp,3.0_dp,4.0_dp],[3,2])
        x=dense_to_stm(a)
        call write_stm_cluto(x,cluto)
        y=read_stm_cluto(cluto)
        call assert_mat_close(y%to_dense(),a,'CLUTO roundtrip')
        call write_stm_mc(x,mc)
        y=read_stm_mc(mc,'tfn')
        call assert_mat_close(y%to_dense(),transpose(a),'MC writer compatibility')
        call delete_file(cluto)
        call delete_file(mc//'_dim')
        call delete_file(mc//'_row_ccs')
        call delete_file(mc//'_col_ccs')
        call delete_file(mc//'_tfn_nz')
    end subroutine test_io

    function sum_fun(x) result(v)
        real(dp),intent(in)::x(:)
        real(dp)::v
        v=sum(x)
    end function sum_fun

    function dot_fun(x,y) result(v)
        real(dp),intent(in)::x(:),y(:)
        real(dp)::v
        v=sum(x*y)
    end function dot_fun

    subroutine assert_true(cond,msg)
        logical,intent(in)::cond
        character(len=*),intent(in)::msg
        if(.not.cond) then
            write(*,'(a)') 'FAIL: '//msg
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_vec_close(a,b,msg)
        real(dp),intent(in)::a(:),b(:)
        character(len=*),intent(in)::msg
        call assert_true(size(a)==size(b),msg//' size')
        if(size(a)>0) call assert_true(all(abs(a-b)<=1.0e-11_dp*(1.0_dp+abs(b))),msg)
    end subroutine assert_vec_close

    subroutine assert_mat_close(a,b,msg)
        real(dp),intent(in)::a(:,:),b(:,:)
        character(len=*),intent(in)::msg
        call assert_true(all(shape(a)==shape(b)),msg//' shape')
        if(size(a)>0) call assert_true(all(abs(a-b)<=1.0e-11_dp*(1.0_dp+abs(b))),msg)
    end subroutine assert_mat_close

    subroutine assert_mat_same_nan(a,b,msg)
        real(dp),intent(in)::a(:,:),b(:,:)
        character(len=*),intent(in)::msg
        integer::r,c
        call assert_true(all(shape(a)==shape(b)),msg//' shape')
        do c=1,size(a,2); do r=1,size(a,1)
            if(ieee_is_nan(b(r,c))) then
                call assert_true(ieee_is_nan(a(r,c)),msg//' nan')
            else
                call assert_true(abs(a(r,c)-b(r,c))<=1.0e-11_dp*(1.0_dp+abs(b(r,c))),msg)
            end if
        end do; end do
    end subroutine assert_mat_same_nan

    subroutine delete_file(name)
        character(len=*),intent(in)::name
        integer::u,ios
        open(newunit=u,file=name,status='old',iostat=ios)
        if(ios==0) close(u,status='delete')
    end subroutine delete_file

end program test_slam
