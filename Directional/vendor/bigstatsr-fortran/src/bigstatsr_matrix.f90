! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr_matrix
    use bigstatsr_kinds, only: dp
    use bigstatsr_fbm, only: fbm_real, fbm_code256
    implicit none
    private

    type, public :: colstats_result
        real(dp), allocatable :: sum(:)
        real(dp), allocatable :: var(:)
    end type colstats_result

    type, public :: scaling_result
        real(dp), allocatable :: center(:)
        real(dp), allocatable :: scale(:)
    end type scaling_result

    public :: big_colstats, big_scale, big_prod_vec, big_cprod_vec
    public :: big_prod_mat, big_cprod_mat, big_crossprod_self, big_tcrossprod_self
    public :: big_cor, big_counts_rows, big_counts_cols

contains

    function big_colstats(x, rows, cols) result(res)
        type(fbm_real), intent(in) :: x
        integer, intent(in), optional :: rows(:), cols(:)
        type(colstats_result) :: res
        integer, allocatable :: rr(:), cc(:)
        real(dp), allocatable :: col(:)
        real(dp) :: sx, sxx
        integer :: j, n
        call resolve_indices(x,rows,cols,rr,cc)
        n=size(rr)
        allocate(res%sum(size(cc)),res%var(size(cc)),col(n))
        do j=1,size(cc)
            call x%read_col(cc(j),col,rr)
            sx=sum(col)
            sxx=dot_product(col,col)
            res%sum(j)=sx
            if (n>1) then
                res%var(j)=(sxx-sx*sx/real(n,dp))/real(n-1,dp)
                if (res%var(j)<0.0_dp .and. abs(res%var(j))<100.0_dp*epsilon(1.0_dp)*max(1.0_dp,sxx)) &
                    res%var(j)=0.0_dp
            else
                res%var(j)=0.0_dp
            end if
        end do
    end function big_colstats

    function big_scale(x, center, scale, rows, cols) result(ms)
        type(fbm_real), intent(in) :: x
        logical, intent(in), optional :: center,scale
        integer, intent(in), optional :: rows(:),cols(:)
        type(scaling_result) :: ms
        type(colstats_result) :: st
        integer, allocatable :: rr(:),cc(:)
        logical :: do_center,do_scale
        call resolve_indices(x,rows,cols,rr,cc)
        do_center=.true.
        do_scale=.true.
        if (present(center)) do_center=center
        if (present(scale)) do_scale=scale
        if (.not.do_center) do_scale=.false.
        allocate(ms%center(size(cc)),ms%scale(size(cc)))
        if (do_center) then
            st=big_colstats(x,rr,cc)
            ms%center=st%sum/real(size(rr),dp)
            if (do_scale) then
                ms%scale=sqrt(max(st%var,0.0_dp))
            else
                ms%scale=1.0_dp
            end if
        else
            ms%center=0.0_dp
            ms%scale=1.0_dp
        end if
    end function big_scale

    function big_prod_vec(x,a,rows,cols,center,scale) result(y)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: a(:)
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: center(:),scale(:)
        real(dp), allocatable :: y(:)
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: col(:)
        real(dp) :: ctr,scl
        integer :: j
        call resolve_indices(x,rows,cols,rr,cc)
        if (size(a)/=size(cc)) error stop 'big_prod_vec: coefficient length mismatch'
        if (present(center)) then
            if (size(center)/=size(cc)) error stop 'big_prod_vec: center length mismatch'
        end if
        if (present(scale)) then
            if (size(scale)/=size(cc)) error stop 'big_prod_vec: scale length mismatch'
        end if
        allocate(y(size(rr)),col(size(rr)))
        y=0.0_dp
        do j=1,size(cc)
            ctr=0.0_dp
            scl=1.0_dp
            if (present(center)) ctr=center(j)
            if (present(scale)) scl=scale(j)
            if (scl<=0.0_dp) error stop 'big_prod_vec: nonpositive scale'
            if (abs(a(j))>0.0_dp) then
                call x%read_col(cc(j),col,rr)
                y=y+a(j)*(col-ctr)/scl
            end if
        end do
    end function big_prod_vec

    function big_cprod_vec(x,a,rows,cols,center,scale) result(y)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: a(:)
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: center(:),scale(:)
        real(dp), allocatable :: y(:)
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: col(:)
        real(dp) :: ctr,scl
        integer :: j
        call resolve_indices(x,rows,cols,rr,cc)
        if (size(a)/=size(rr)) error stop 'big_cprod_vec: vector length mismatch'
        allocate(y(size(cc)),col(size(rr)))
        do j=1,size(cc)
            ctr=0.0_dp
            scl=1.0_dp
            if (present(center)) ctr=center(j)
            if (present(scale)) scl=scale(j)
            if (scl<=0.0_dp) error stop 'big_cprod_vec: nonpositive scale'
            call x%read_col(cc(j),col,rr)
            y(j)=dot_product((col-ctr)/scl,a)
        end do
    end function big_cprod_vec

    function big_prod_mat(x,a,rows,cols,center,scale) result(y)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: center(:),scale(:)
        real(dp), allocatable :: y(:,:)
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: col(:)
        real(dp) :: ctr,scl
        integer :: j
        call resolve_indices(x,rows,cols,rr,cc)
        if (size(a,1)/=size(cc)) error stop 'big_prod_mat: row mismatch'
        allocate(y(size(rr),size(a,2)),col(size(rr)))
        y=0.0_dp
        do j=1,size(cc)
            ctr=0.0_dp
            scl=1.0_dp
            if (present(center)) ctr=center(j)
            if (present(scale)) scl=scale(j)
            call x%read_col(cc(j),col,rr)
            y=y+spread((col-ctr)/scl,2,size(a,2))*spread(a(j,:),1,size(rr))
        end do
    end function big_prod_mat

    function big_cprod_mat(x,a,rows,cols,center,scale) result(y)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: center(:),scale(:)
        real(dp), allocatable :: y(:,:)
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: col(:)
        real(dp) :: ctr,scl
        integer :: j
        call resolve_indices(x,rows,cols,rr,cc)
        if (size(a,1)/=size(rr)) error stop 'big_cprod_mat: row mismatch'
        allocate(y(size(cc),size(a,2)),col(size(rr)))
        do j=1,size(cc)
            ctr=0.0_dp
            scl=1.0_dp
            if (present(center)) ctr=center(j)
            if (present(scale)) scl=scale(j)
            call x%read_col(cc(j),col,rr)
            y(j,:)=matmul((col-ctr)/scl,a)
        end do
    end function big_cprod_mat

    function big_crossprod_self(x,rows,cols,center,scale) result(k)
        type(fbm_real), intent(in) :: x
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: center(:),scale(:)
        real(dp), allocatable :: k(:,:)
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: xi(:),xj(:)
        real(dp) :: ci,cj,si,sj
        integer :: i,j
        call resolve_indices(x,rows,cols,rr,cc)
        allocate(k(size(cc),size(cc)),xi(size(rr)),xj(size(rr)))
        k=0.0_dp
        do j=1,size(cc)
            call x%read_col(cc(j),xj,rr)
            cj=0.0_dp
            sj=1.0_dp
            if (present(center)) cj=center(j)
            if (present(scale)) sj=scale(j)
            xj=(xj-cj)/sj
            do i=1,j
                call x%read_col(cc(i),xi,rr)
                ci=0.0_dp
                si=1.0_dp
                if (present(center)) ci=center(i)
                if (present(scale)) si=scale(i)
                k(i,j)=dot_product((xi-ci)/si,xj)
                k(j,i)=k(i,j)
            end do
        end do
    end function big_crossprod_self

    function big_tcrossprod_self(x,rows,cols,center,scale) result(k)
        type(fbm_real), intent(in) :: x
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: center(:),scale(:)
        real(dp), allocatable :: k(:,:)
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: col(:)
        real(dp) :: ctr,scl
        integer :: j
        call resolve_indices(x,rows,cols,rr,cc)
        allocate(k(size(rr),size(rr)),col(size(rr)))
        k=0.0_dp
        do j=1,size(cc)
            ctr=0.0_dp
            scl=1.0_dp
            if (present(center)) ctr=center(j)
            if (present(scale)) scl=scale(j)
            call x%read_col(cc(j),col,rr)
            col=(col-ctr)/scl
            k=k+spread(col,2,size(rr))*spread(col,1,size(rr))
        end do
    end function big_tcrossprod_self

    function big_cor(x,rows,cols) result(k)
        type(fbm_real), intent(in) :: x
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), allocatable :: k(:,:)
        type(scaling_result) :: ms
        integer, allocatable :: rr(:),cc(:)
        call resolve_indices(x,rows,cols,rr,cc)
        ms=big_scale(x,.true.,.true.,rr,cc)
        ms%scale=ms%scale*sqrt(real(size(rr)-1,dp))
        k=big_crossprod_self(x,rr,cc,ms%center,ms%scale)
    end function big_cor

    function big_counts_rows(x,code_map,rows,cols) result(counts)
        type(fbm_code256), intent(in) :: x
        integer, intent(in) :: code_map(0:255)
        integer, intent(in), optional :: rows(:),cols(:)
        integer, allocatable :: counts(:,:)
        integer, allocatable :: rr(:),cc(:),col(:)
        integer :: kmax,i,j,k
        call resolve_code_indices(x,rows,cols,rr,cc)
        kmax=maxval(code_map)
        allocate(counts(kmax,size(rr)),col(x%nrow))
        counts=0
        do j=1,size(cc)
            call x%read_col(cc(j),col)
            do i=1,size(rr)
                k=code_map(col(rr(i)))
                if (k>0) counts(k,i)=counts(k,i)+1
            end do
        end do
    end function big_counts_rows

    function big_counts_cols(x,code_map,rows,cols) result(counts)
        type(fbm_code256), intent(in) :: x
        integer, intent(in) :: code_map(0:255)
        integer, intent(in), optional :: rows(:),cols(:)
        integer, allocatable :: counts(:,:)
        integer, allocatable :: rr(:),cc(:),col(:)
        integer :: kmax,i,j,k
        call resolve_code_indices(x,rows,cols,rr,cc)
        kmax=maxval(code_map)
        allocate(counts(kmax,size(cc)),col(x%nrow))
        counts=0
        do j=1,size(cc)
            call x%read_col(cc(j),col)
            do i=1,size(rr)
                k=code_map(col(rr(i)))
                if (k>0) counts(k,j)=counts(k,j)+1
            end do
        end do
    end function big_counts_cols

    subroutine resolve_indices(x,rows,cols,rr,cc)
        type(fbm_real), intent(in) :: x
        integer, intent(in), optional :: rows(:),cols(:)
        integer, allocatable, intent(out) :: rr(:),cc(:)
        integer :: i
        if (present(rows)) then
            rr=rows
        else
            allocate(rr(x%nrow))
            rr=[(i,i=1,x%nrow)]
        end if
        if (present(cols)) then
            cc=cols
        else
            allocate(cc(x%ncol))
            cc=[(i,i=1,x%ncol)]
        end if
        if (any(rr<1) .or. any(rr>x%nrow)) error stop 'row index out of bounds'
        if (any(cc<1) .or. any(cc>x%ncol)) error stop 'column index out of bounds'
    end subroutine resolve_indices

    subroutine resolve_code_indices(x,rows,cols,rr,cc)
        type(fbm_code256), intent(in) :: x
        integer, intent(in), optional :: rows(:),cols(:)
        integer, allocatable, intent(out) :: rr(:),cc(:)
        integer :: i
        if (present(rows)) then
            rr=rows
        else
            allocate(rr(x%nrow)); rr=[(i,i=1,x%nrow)]
        end if
        if (present(cols)) then
            cc=cols
        else
            allocate(cc(x%ncol)); cc=[(i,i=1,x%ncol)]
        end if
    end subroutine resolve_code_indices

end module bigstatsr_matrix
