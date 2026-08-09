! Upstream license declaration: GPL (version unspecified)
module limsolve_linear
    use limsolve_kinds, only: dp
    use limsolve_types, only: LS_SUCCESS, LS_INVALID, LS_SINGULAR
    use limsolve_linalg, only: dense_solve
    implicit none
    private
    public :: solve_tridiag, solve_banded, solve_block

contains

    subroutine solve_tridiag(diam1,dia,diap1,b,x,status)
        real(dp), intent(in) :: diam1(:),dia(:),diap1(:),b(:,:)
        real(dp), intent(out) :: x(:,:)
        integer, intent(out) :: status
        real(dp), allocatable :: a(:,:), rhs(:), sol(:)
        integer :: n,nrhs,i,j,st
        n=size(dia); nrhs=size(b,2)
        status=LS_INVALID; x=0.0_dp
        if (size(diam1)/=max(0,n-1) .or. size(diap1)/=max(0,n-1) .or. &
            size(b,1)/=n .or. size(x,1)/=n .or. size(x,2)/=nrhs) return
        allocate(a(n,n),rhs(n),sol(n)); a=0.0_dp
        do i=1,n
            a(i,i)=dia(i)
        end do
        do i=1,n-1
            a(i+1,i)=diam1(i); a(i,i+1)=diap1(i)
        end do
        do j=1,nrhs
            rhs=b(:,j); call dense_solve(a,rhs,sol,st)
            if (st/=0) then; status=LS_SINGULAR; return; end if
            x(:,j)=sol
        end do
        status=LS_SUCCESS
    end subroutine solve_tridiag

    subroutine solve_banded(abd,nup,nlow,b,x,status,full)
        real(dp), intent(in) :: abd(:,:),b(:,:)
        integer, intent(in) :: nup,nlow
        real(dp), intent(out) :: x(:,:)
        integer, intent(out) :: status
        logical, intent(in), optional :: full
        real(dp), allocatable :: a(:,:),rhs(:),sol(:)
        integer :: n,nrhs,i,j,k,st
        logical :: isfull
        n=size(abd,2); nrhs=size(b,2); isfull=.false.; if(present(full)) isfull=full
        status=LS_INVALID; x=0.0_dp
        if(size(b,1)/=n .or. size(x,1)/=n .or. size(x,2)/=nrhs) return
        allocate(a(n,n),rhs(n),sol(n)); a=0.0_dp
        if(isfull) then
            if(size(abd,1)/=n) return
            a=abd
        else
            if(size(abd,1)/=nup+nlow+1) return
            do j=1,n
                do k=1,nup+nlow+1
                    i=j+k-nup-1
                    if(i>=1 .and. i<=n) a(i,j)=abd(k,j)
                end do
            end do
        end if
        do j=1,nrhs
            rhs=b(:,j); call dense_solve(a,rhs,sol,st)
            if(st/=0) then; status=LS_SINGULAR; return; end if
            x(:,j)=sol
        end do
        status=LS_SUCCESS
    end subroutine solve_banded

    subroutine solve_block(top,blocks,bot,b,overlap,x,status)
        real(dp), intent(in) :: top(:,:),blocks(:,:,:),bot(:,:),b(:,:)
        integer, intent(in) :: overlap
        real(dp), intent(out) :: x(:,:)
        integer, intent(out) :: status
        real(dp), allocatable :: a(:,:),rhs(:),sol(:)
        integer :: nrwt,nrwb,nrowb,ncolb,nb,n,nrhs,k,r0,c0,j,st
        nrwt=size(top,1); nrwb=size(bot,1); nrowb=size(blocks,1)
        ncolb=size(blocks,2); nb=size(blocks,3)
        n=nb*nrowb+overlap; nrhs=size(b,2)
        status=LS_INVALID; x=0.0_dp
        if(size(top,2)/=overlap .or. size(bot,2)/=overlap .or. &
            ncolb/=nrowb+overlap .or. overlap/=nrwt+nrwb .or. &
            size(b,1)/=n .or. size(x,1)/=n .or. size(x,2)/=nrhs) return
        allocate(a(n,n),rhs(n),sol(n)); a=0.0_dp
        if(nrwt>0) a(1:nrwt,1:overlap)=top
        do k=1,nb
            r0=nrwt+(k-1)*nrowb+1
            c0=(k-1)*nrowb+1
            a(r0:r0+nrowb-1,c0:c0+ncolb-1)=blocks(:,:,k)
        end do
        if(nrwb>0) a(n-nrwb+1:n,n-overlap+1:n)=bot
        do j=1,nrhs
            rhs=b(:,j); call dense_solve(a,rhs,sol,st)
            if(st/=0) then; status=LS_SINGULAR; return; end if
            x(:,j)=sol
        end do
        status=LS_SUCCESS
    end subroutine solve_block

end module limsolve_linear
