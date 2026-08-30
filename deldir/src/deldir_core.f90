module deldir_core
use deldir_kinds, only: dp
use deldir_kernel, only: master, binsrt
use deldir_types, only: deldir_result
implicit none
private
public :: deldir_compute
contains

subroutine deldir_compute(x, y, result, rw, eps, sort_points, status)
    real(dp), intent(in) :: x(:), y(:)
    type(deldir_result), intent(out) :: result
    real(dp), intent(in), optional :: rw(4), eps
    logical, intent(in), optional :: sort_points
    integer, intent(out), optional :: status
    real(dp), allocatable :: xu(:), yu(:), xs(:), ys(:), tx(:), ty(:)
    real(dp), allocatable :: xk(:), yk(:), delsgs(:,:), dirsgs(:,:), delsum(:,:), dirsum(:,:)
    integer, allocatable :: orig(:), ind(:), rind(:), ilst(:), nadj(:,:)
    real(dp) :: win(4), tol, dx, dy, total
    integer :: n, nu, i, j, k, ntot, madj, capdel, capdir, ndel, ndir, incadj, incseg, srow, i1, i2
    integer :: attempt
    logical :: do_sort, duplicate

    if (present(status)) status = 0
    if (size(x) /= size(y)) then
        call fail(1, 'deldir_compute: x and y lengths differ', status); return
    end if
    n = size(x)
    if (n < 2) then
        call fail(2, 'deldir_compute: at least two points are required', status); return
    end if
    tol = 1.0e-9_dp
    if (present(eps)) tol = eps
    if (tol <= 0.0_dp) then
        call fail(3, 'deldir_compute: eps must be positive', status); return
    end if
    do_sort = .true.; if (present(sort_points)) do_sort = sort_points

    if (present(rw)) then
        win = rw
        if (win(1) >= win(2) .or. win(3) >= win(4)) then
            call fail(4, 'deldir_compute: invalid rectangular window', status); return
        end if
    else
        dx = maxval(x)-minval(x); dy = maxval(y)-minval(y)
        if (dx <= 0.0_dp .or. dy <= 0.0_dp) then
            call fail(5, 'deldir_compute: zero coordinate range requires an explicit window', status); return
        end if
        win = [minval(x)-0.1_dp*dx, maxval(x)+0.1_dp*dx, minval(y)-0.1_dp*dy, maxval(y)+0.1_dp*dy]
    end if

    allocate(xu(n), yu(n), orig(n))
    nu = 0
    do i=1,n
        if (x(i) < win(1) .or. x(i) > win(2) .or. y(i) < win(3) .or. y(i) > win(4)) cycle
        duplicate = .false.
        do j=1,nu
            if (x(i) <= xu(j) .and. x(i) >= xu(j) .and. y(i) <= yu(j) .and. y(i) >= yu(j)) then
                duplicate = .true.; exit
            end if
        end do
        if (.not. duplicate) then
            nu=nu+1; xu(nu)=x(i); yu(nu)=y(i); orig(nu)=i
        end if
    end do
    if (nu < 2) then
        call fail(6, 'deldir_compute: fewer than two distinct points lie inside the window', status); return
    end if
    xu=xu(:nu); yu=yu(:nu); orig=orig(:nu)
    allocate(xs(nu),ys(nu),tx(nu),ty(nu),ind(nu),rind(nu),ilst(nu))
    xs=xu; ys=yu
    if (do_sort) then
        call binsrt(xs,ys,win,nu,ind,rind,tx,ty,ilst)
    else
        do i=1,nu; ind(i)=i; rind(i)=i; end do
    end if

    ntot=nu+4
    madj=max(20, ceiling(3.0_dp*sqrt(real(ntot,dp))))
    capdel=madj*(madj+1)/2
    capdir=capdel
    do attempt=1,30
        allocate(xk(-3:ntot),yk(-3:ntot),nadj(-3:ntot,0:madj))
        allocate(delsgs(6,capdel),dirsgs(10,capdir),delsum(nu,4),dirsum(nu,3))
        xk=0.0_dp; yk=0.0_dp; xk(1:nu)=xs; yk(1:nu)=ys
        delsgs=0.0_dp; dirsgs=0.0_dp; delsum=0.0_dp; dirsum=0.0_dp
        ndel=capdel; ndir=capdir; incadj=0; incseg=0
        call master(xk,yk,win,nu,ntot,nadj,madj,tol,delsgs,ndel,delsum,dirsgs,ndir,dirsum,incadj,incseg)
        if (incadj==0 .and. incseg==0) exit
        deallocate(xk,yk,nadj,delsgs,dirsgs,delsum,dirsum)
        if (incadj==1) madj=max(madj+1,ceiling(1.2_dp*real(madj,dp)))
        if (incseg==1) then
            capdel=max(capdel+1,ceiling(1.2_dp*real(capdel,dp)))
            capdir=max(capdir+1,ceiling(1.2_dp*real(capdir,dp)))
        end if
        capdel=max(capdel,madj*(madj+1)/2); capdir=max(capdir,capdel)
    end do
    if (incadj/=0 .or. incseg/=0) then
        call fail(7,'deldir_compute: internal storage could not be enlarged sufficiently',status); return
    end if

    result%n_data=nu; result%rw=win
    allocate(result%delsgs(ndel), result%dirsgs(ndir), result%summary(nu), result%ind_orig(nu))
    result%ind_orig=orig
    do k=1,ndel
        i1=rind(nint(delsgs(5,k))); i2=rind(nint(delsgs(6,k)))
        result%delsgs(k)%x1=delsgs(1,k); result%delsgs(k)%y1=delsgs(2,k)
        result%delsgs(k)%x2=delsgs(3,k); result%delsgs(k)%y2=delsgs(4,k)
        result%delsgs(k)%ind1=i1; result%delsgs(k)%ind2=i2
    end do
    do k=1,ndir
        i1=rind(nint(dirsgs(5,k))); i2=rind(nint(dirsgs(6,k)))
        result%dirsgs(k)%x1=dirsgs(1,k); result%dirsgs(k)%y1=dirsgs(2,k)
        result%dirsgs(k)%x2=dirsgs(3,k); result%dirsgs(k)%y2=dirsgs(4,k)
        result%dirsgs(k)%ind1=i1; result%dirsgs(k)%ind2=i2
        result%dirsgs(k)%bp1=dirsgs(7,k)>0.5_dp; result%dirsgs(k)%bp2=dirsgs(8,k)>0.5_dp
        j=nint(dirsgs(9,k)); if (j>0) j=rind(j); result%dirsgs(k)%thirdv1=j
        j=nint(dirsgs(10,k)); if (j>0) j=rind(j); result%dirsgs(k)%thirdv2=j
    end do
    result%del_area=sum(delsum(:,4)); result%dir_area=sum(dirsum(:,3))
    do j=1,nu
        srow=ind(j)
        result%summary(j)%x=xu(j); result%summary(j)%y=yu(j); result%summary(j)%original_index=orig(j)
        result%summary(j)%n_tri=nint(delsum(srow,3)); result%summary(j)%del_area=delsum(srow,4)
        result%summary(j)%n_tside=nint(dirsum(srow,1)); result%summary(j)%nbpt=nint(dirsum(srow,2))
        result%summary(j)%dir_area=dirsum(srow,3)
    end do
    total=result%del_area
    if (total>0.0_dp) then
        do j=1,nu; result%summary(j)%del_wt=result%summary(j)%del_area/total; end do
    end if
    total=result%dir_area
    if (total>0.0_dp) then
        do j=1,nu; result%summary(j)%dir_wt=result%summary(j)%dir_area/total; end do
    end if
end subroutine deldir_compute

subroutine fail(code,message,status)
    integer,intent(in)::code
    character(len=*),intent(in)::message
    integer,intent(out),optional::status
    if (present(status)) then
        status=code
    else
        error stop message
    end if
end subroutine fail
end module deldir_core
