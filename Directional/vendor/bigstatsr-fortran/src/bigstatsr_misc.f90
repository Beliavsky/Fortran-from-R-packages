! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr_misc
    use bigstatsr_kinds, only: dp
    use bigstatsr_utils, only: median_dp
    implicit none
    private
    public :: block_size, rows_along, cols_along, get_beta

contains

    integer function block_size(n,ncores,block_gb) result(m)
        integer, intent(in) :: n
        integer, intent(in), optional :: ncores
        real(dp), intent(in), optional :: block_gb
        integer :: nc
        real(dp) :: gb
        nc=1
        if (present(ncores)) nc=max(1,ncores)
        gb=1.0_dp
        if (present(block_gb)) gb=max(block_gb,tiny(1.0_dp))
        m=max(1,int(floor((gb/real(nc,dp))*1024.0_dp**3/(8.0_dp*real(max(1,n),dp)))))
    end function block_size

    function rows_along(n) result(ind)
        integer, intent(in) :: n
        integer, allocatable :: ind(:)
        integer :: i
        allocate(ind(n))
        ind=[(i,i=1,n)]
    end function rows_along

    function cols_along(n) result(ind)
        integer, intent(in) :: n
        integer, allocatable :: ind(:)
        integer :: i
        allocate(ind(n))
        ind=[(i,i=1,n)]
    end function cols_along

    function get_beta(betas,method,tol,maxiter) result(beta)
        real(dp), intent(in) :: betas(:,:)
        character(len=*), intent(in), optional :: method
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: old(:),new(:),dist(:)
        real(dp) :: eps,wden
        integer :: it,mi,j
        character(len=32) :: meth
        meth='geometric-median'
        if (present(method)) meth=trim(method)
        allocate(beta(size(betas,1)))
        select case(meth)
        case('mean-wise','mean')
            beta=sum(betas,dim=2)/real(size(betas,2),dp)
        case('median-wise','median')
            do j=1,size(betas,1)
                beta(j)=median_dp(betas(j,:))
            end do
        case default
            allocate(old(size(beta)),new(size(beta)),dist(size(betas,2)))
            old=sum(betas,dim=2)/real(size(betas,2),dp)
            eps=1.0e-10_dp
            if (present(tol)) eps=tol
            mi=10000
            if (present(maxiter)) mi=maxiter
            do it=1,mi
                do j=1,size(betas,2)
                    dist(j)=sqrt(sum((betas(:,j)-old)**2))
                end do
                if (any(dist<=sqrt(epsilon(1.0_dp)))) then
                    beta=betas(:,minloc(dist,dim=1))
                    return
                end if
                new=0.0_dp
                wden=sum(1.0_dp/dist)
                do j=1,size(betas,2)
                    new=new+betas(:,j)/dist(j)
                end do
                new=new/wden
                if (maxval(abs(new-old))<eps) exit
                old=new
            end do
            beta=new
        end select
    end function get_beta

end module bigstatsr_misc
