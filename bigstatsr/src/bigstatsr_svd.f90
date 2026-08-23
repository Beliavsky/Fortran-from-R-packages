! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr_svd
    use bigstatsr_kinds, only: dp
    use bigstatsr_fbm, only: fbm_real
    use bigstatsr_matrix, only: scaling_result, big_scale, big_crossprod_self, &
        big_tcrossprod_self, big_prod_mat, big_cprod_mat
    use rspectra, only: linear_operator, svds, svds_opts, svds_result
    use rspectra_external, only: dsyev
    implicit none
    private

    type, public :: big_svd_result
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: u(:,:)
        real(dp), allocatable :: v(:,:)
        real(dp), allocatable :: center(:)
        real(dp), allocatable :: scale(:)
        integer :: niter = 0
        integer :: nops = 0
        integer :: info = 0
    end type big_svd_result

    type, extends(linear_operator) :: fbm_operator
        type(fbm_real), pointer :: x => null()
        integer, allocatable :: rows(:)
        integer, allocatable :: cols(:)
    contains
        procedure :: prod => fbm_op_prod
        procedure :: tprod => fbm_op_tprod
    end type fbm_operator

    public :: big_svd, big_random_svd, svd_predict

contains

    function big_random_svd(x,k,center,scale,rows,cols,tol,maxiter,ncv) result(out)
        type(fbm_real), target, intent(in) :: x
        integer, intent(in) :: k
        logical, intent(in), optional :: center,scale
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter,ncv
        type(big_svd_result) :: out
        type(fbm_operator) :: op
        type(scaling_result) :: ms
        type(svds_opts) :: opts
        type(svds_result) :: r
        integer, allocatable :: rr(:),cc(:)
        logical :: dc,ds
        integer :: i
        call resolve_indices(x,rows,cols,rr,cc)
        dc=.false.; ds=.false.
        if (present(center)) dc=center
        if (present(scale)) ds=scale
        ms=big_scale(x,dc,ds,rr,cc)
        if (any(ms%scale<=sqrt(epsilon(1.0_dp)))) then
            out%info=-10
            return
        end if
        op%x=>x
        op%rows=rr
        op%cols=cc
        op%nrow=size(rr)
        op%ncol=size(cc)
        opts=svds_opts()
        if (present(tol)) opts%tol=tol
        if (present(maxiter)) opts%maxitr=maxiter
        if (present(ncv)) opts%ncv=ncv
        r=svds(op,k,nu=k,nv=k,opts=opts,center=ms%center,scale=ms%scale)
        out%info=r%info
        out%niter=r%niter
        out%nops=r%nops
        out%center=ms%center
        out%scale=ms%scale
        if (allocated(r%d)) out%d=r%d
        if (allocated(r%u)) out%u=r%u
        if (allocated(r%v)) out%v=r%v
        if (allocated(out%d)) then
            do i=1,size(out%d)-1
                if (out%d(i)<out%d(i+1)) exit
            end do
        end if
    end function big_random_svd

    function big_svd(x,k,center,scale,rows,cols) result(out)
        type(fbm_real), intent(in) :: x
        integer, intent(in) :: k
        logical, intent(in), optional :: center,scale
        integer, intent(in), optional :: rows(:),cols(:)
        type(big_svd_result) :: out
        type(scaling_result) :: ms
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: gram(:,:),eig(:),tmp(:,:)
        logical :: dc,ds
        integer :: n,m,kk,info,i,src
        call resolve_indices(x,rows,cols,rr,cc)
        n=size(rr); m=size(cc); kk=min(k,min(n,m))
        if (kk<1) then
            out%info=-1
            return
        end if
        dc=.false.; ds=.false.
        if (present(center)) dc=center
        if (present(scale)) ds=scale
        ms=big_scale(x,dc,ds,rr,cc)
        if (any(ms%scale<=sqrt(epsilon(1.0_dp)))) then
            out%info=-10
            return
        end if
        out%center=ms%center
        out%scale=ms%scale
        if (m<=n) then
            gram=big_crossprod_self(x,rr,cc,ms%center,ms%scale)
            allocate(eig(m))
            call sym_eigen(gram,eig,info)
            if (info/=0) then
                out%info=info
                return
            end if
            allocate(out%d(kk),out%v(m,kk))
            do i=1,kk
                src=m-i+1
                out%d(i)=sqrt(max(0.0_dp,eig(src)))
                out%v(:,i)=gram(:,src)
            end do
            tmp=big_prod_mat(x,out%v,rr,cc,ms%center,ms%scale)
            allocate(out%u(n,kk))
            do i=1,kk
                if (out%d(i)>sqrt(epsilon(1.0_dp))) then
                    out%u(:,i)=tmp(:,i)/out%d(i)
                else
                    out%u(:,i)=0.0_dp
                end if
            end do
        else
            gram=big_tcrossprod_self(x,rr,cc,ms%center,ms%scale)
            allocate(eig(n))
            call sym_eigen(gram,eig,info)
            if (info/=0) then
                out%info=info
                return
            end if
            allocate(out%d(kk),out%u(n,kk))
            do i=1,kk
                src=n-i+1
                out%d(i)=sqrt(max(0.0_dp,eig(src)))
                out%u(:,i)=gram(:,src)
            end do
            tmp=big_cprod_mat(x,out%u,rr,cc,ms%center,ms%scale)
            allocate(out%v(m,kk))
            do i=1,kk
                if (out%d(i)>sqrt(epsilon(1.0_dp))) then
                    out%v(:,i)=tmp(:,i)/out%d(i)
                else
                    out%v(:,i)=0.0_dp
                end if
            end do
        end if
        out%info=0
    end function big_svd

    function svd_predict(x,fit,rows,cols) result(pc)
        type(fbm_real), intent(in) :: x
        type(big_svd_result), intent(in) :: fit
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), allocatable :: pc(:,:)
        pc=big_prod_mat(x,fit%v,rows,cols,fit%center,fit%scale)
    end function svd_predict

    subroutine sym_eigen(a,w,info)
        real(dp), intent(inout) :: a(:,:)
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: info
        real(dp), allocatable :: work(:)
        integer :: n,lwork
        n=size(a,1)
        lwork=max(1,3*n-1)
        allocate(work(lwork))
        call dsyev('V','U',n,a,n,w,work,lwork,info)
    end subroutine sym_eigen

    subroutine fbm_op_prod(self,x,y)
        class(fbm_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp), allocatable :: col(:)
        integer :: j
        allocate(col(size(self%rows)))
        y=0.0_dp
        do j=1,size(self%cols)
            call self%x%read_col(self%cols(j),col,self%rows)
            y=y+x(j)*col
        end do
    end subroutine fbm_op_prod

    subroutine fbm_op_tprod(self,x,y)
        class(fbm_operator), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp), allocatable :: col(:)
        integer :: j
        allocate(col(size(self%rows)))
        do j=1,size(self%cols)
            call self%x%read_col(self%cols(j),col,self%rows)
            y(j)=dot_product(col,x)
        end do
    end subroutine fbm_op_tprod

    subroutine resolve_indices(x,rows,cols,rr,cc)
        type(fbm_real), intent(in) :: x
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
    end subroutine resolve_indices

end module bigstatsr_svd
