module actuar_mde_v03
    use actuar_kinds, only: dp
    use actuar_grouped_v02, only: elev_grouped
    implicit none
    private
    public :: mde_result_t, mde_cvm, mde_grouped_cvm, mde_grouped_chisq, mde_grouped_las

    abstract interface
        function mde_model_fun(x, par) result(v)
            import dp
            real(dp), intent(in) :: x
            real(dp), intent(in) :: par(:)
            real(dp) :: v
        end function mde_model_fun
    end interface

    type :: mde_result_t
        real(dp), allocatable :: estimate(:)
        real(dp) :: distance = huge(1.0_dp)
        logical :: converged = .false.
        integer :: iterations = 0
    end type mde_result_t

    procedure(mde_model_fun), pointer :: active_fun => null()
    real(dp), allocatable :: ctx_x(:), ctx_target(:), ctx_weights(:)
    integer :: ctx_mode = 0
    real(dp) :: ctx_total = 0.0_dp
contains
    function mde_cvm(x, fun, start, weights, lower, upper, tol, maxit) result(res)
        real(dp), intent(in) :: x(:), start(:)
        procedure(mde_model_fun) :: fun
        real(dp), intent(in), optional :: weights(:), lower(:), upper(:), tol
        integer, intent(in), optional :: maxit
        type(mde_result_t) :: res
        real(dp), allocatable :: knots(:)
        integer :: i, n, m
        n=size(x); if(n<1 .or. size(start)<1) return
        allocate(knots(n));knots=x;call sort_real(knots)
        m=1
        do i=2,n
            if(knots(i)/=knots(m)) then;m=m+1;knots(m)=knots(i);end if
        end do
        call set_context(fun,knots(:m),1,m)
        do i=1,m;ctx_target(i)=real(count(x<=ctx_x(i)),dp)/real(n,dp);end do
        if(.not.set_weights(weights,m)) return
        call nelder_mead(mde_objective,start,res,lower,upper,tol,maxit)
        call clear_context()
    end function mde_cvm

    function mde_grouped_cvm(boundaries,counts,fun,start,weights,lower,upper,tol,maxit) result(res)
        real(dp),intent(in)::boundaries(:),counts(:),start(:)
        procedure(mde_model_fun)::fun
        real(dp),intent(in),optional::weights(:),lower(:),upper(:),tol
        integer,intent(in),optional::maxit
        type(mde_result_t)::res
        integer::m,i
        m=size(boundaries)
        if(size(counts)/=m-1 .or. any(counts<0.0_dp) .or. sum(counts)<=0.0_dp) return
        call set_context(fun,boundaries,1,m)
        ctx_target(1)=0.0_dp
        do i=2,m;ctx_target(i)=sum(counts(:i-1))/sum(counts);end do
        if(.not.set_weights(weights,m)) return
        call nelder_mead(mde_objective,start,res,lower,upper,tol,maxit)
        call clear_context()
    end function mde_grouped_cvm

    function mde_grouped_chisq(boundaries,counts,fun,start,weights,lower,upper,tol,maxit) result(res)
        real(dp),intent(in)::boundaries(:),counts(:),start(:)
        procedure(mde_model_fun)::fun
        real(dp),intent(in),optional::weights(:),lower(:),upper(:),tol
        integer,intent(in),optional::maxit
        type(mde_result_t)::res
        integer::r
        r=size(counts)
        if(size(boundaries)/=r+1 .or. any(counts<=0.0_dp)) return
        call set_context(fun,boundaries,2,r)
        ctx_target=counts;ctx_total=sum(counts)
        ctx_weights=1.0_dp/counts
        if(present(weights)) then
            if(size(weights)==1) then;ctx_weights=weights(1)
            else if(size(weights)==r) then;ctx_weights=weights
            else;call clear_context();return;end if
        end if
        call nelder_mead(mde_objective,start,res,lower,upper,tol,maxit)
        call clear_context()
    end function mde_grouped_chisq

    function mde_grouped_las(boundaries,counts,levfun,start,weights,lower,upper,tol,maxit) result(res)
        real(dp),intent(in)::boundaries(:),counts(:),start(:)
        procedure(mde_model_fun)::levfun
        real(dp),intent(in),optional::weights(:),lower(:),upper(:),tol
        integer,intent(in),optional::maxit
        type(mde_result_t)::res
        real(dp),allocatable::emp(:)
        integer::r,k
        r=size(counts)
        if(size(boundaries)/=r+1 .or. any(counts<0.0_dp)) return
        allocate(emp(r+1));do k=1,r+1;emp(k)=elev_grouped(boundaries(k),boundaries,counts);end do
        call set_context(levfun,boundaries,3,r)
        do k=1,r;ctx_target(k)=emp(k+1)-emp(k);end do
        if(.not.set_weights(weights,r)) return
        call nelder_mead(mde_objective,start,res,lower,upper,tol,maxit)
        call clear_context()
    end function mde_grouped_las

    subroutine set_context(fun,x,mode,ntarget)
        procedure(mde_model_fun)::fun
        real(dp),intent(in)::x(:)
        integer,intent(in)::mode,ntarget
        active_fun=>fun;ctx_mode=mode;ctx_total=0.0_dp
        if(allocated(ctx_x))deallocate(ctx_x,ctx_target,ctx_weights)
        allocate(ctx_x(size(x)),ctx_target(ntarget),ctx_weights(ntarget))
        ctx_x=x;ctx_target=0.0_dp;ctx_weights=1.0_dp
    end subroutine set_context

    logical function set_weights(weights,n) result(ok)
        real(dp),intent(in),optional::weights(:)
        integer,intent(in)::n
        ok=.true.;ctx_weights=1.0_dp
        if(present(weights)) then
            if(size(weights)==1) then;ctx_weights=weights(1)
            else if(size(weights)==n) then;ctx_weights=weights
            else;ok=.false.;call clear_context();end if
        end if
    end function set_weights

    subroutine clear_context()
        nullify(active_fun);ctx_mode=0;ctx_total=0.0_dp
        if(allocated(ctx_x))deallocate(ctx_x,ctx_target,ctx_weights)
    end subroutine clear_context

    function mde_objective(par) result(v)
        real(dp),intent(in)::par(:)
        real(dp)::v,d,model
        integer::k
        v=0.0_dp
        if(.not.associated(active_fun)) then;v=huge(1.0_dp);return;end if
        select case(ctx_mode)
        case(1)
            do k=1,size(ctx_target)
                d=active_fun(ctx_x(k),par)-ctx_target(k)
                if(.not.finite(d)) then;v=huge(1.0_dp);return;end if
                v=v+ctx_weights(k)*d*d
            end do
        case(2)
            do k=1,size(ctx_target)
                model=ctx_total*(active_fun(ctx_x(k+1),par)-active_fun(ctx_x(k),par))
                d=model-ctx_target(k)
                if(.not.finite(d)) then;v=huge(1.0_dp);return;end if
                v=v+ctx_weights(k)*d*d
            end do
        case(3)
            do k=1,size(ctx_target)
                model=active_fun(ctx_x(k+1),par)-active_fun(ctx_x(k),par)
                d=model-ctx_target(k)
                if(.not.finite(d)) then;v=huge(1.0_dp);return;end if
                v=v+ctx_weights(k)*d*d
            end do
        case default
            v=huge(1.0_dp)
        end select
    end function mde_objective

    subroutine nelder_mead(fn,start,res,lower,upper,tol,maxit)
        interface
            function fn(par) result(v)
                import dp
                real(dp),intent(in)::par(:)
                real(dp)::v
            end function fn
        end interface
        real(dp),intent(in)::start(:)
        type(mde_result_t),intent(out)::res
        real(dp),intent(in),optional::lower(:),upper(:),tol
        integer,intent(in),optional::maxit
        integer::n,it,mi,j,ilo,ihi,inhi
        real(dp)::eps,fr,fe,fc,spread,scale,xspread
        real(dp),allocatable::simp(:,:),f(:),cent(:),xr(:),xe(:),xc(:)
        n=size(start);mi=1000;if(present(maxit))mi=maxit
        eps=1.0e-9_dp;if(present(tol))eps=tol
        allocate(simp(n,n+1),f(n+1),cent(n),xr(n),xe(n),xc(n),res%estimate(n))
        simp(:,1)=start;call enforce_bounds(simp(:,1),lower,upper)
        do j=1,n
            simp(:,j+1)=simp(:,1);scale=0.05_dp*max(abs(start(j)),1.0_dp)
            simp(j,j+1)=simp(j,j+1)+scale;call enforce_bounds(simp(:,j+1),lower,upper)
        end do
        do j=1,n+1;f(j)=fn(simp(:,j));end do
        do it=1,mi
            ilo=minloc(f,dim=1);ihi=maxloc(f,dim=1);inhi=0
            do j=1,n+1
                if(j/=ihi) then
                    if(inhi==0)then
                        inhi=j
                    else if(f(j)>f(inhi))then
                        inhi=j
                    end if
                end if
            end do
            spread=maxval(abs(f-f(ilo)))/max(1.0_dp,abs(f(ilo)))
            xspread=0.0_dp
            do j=1,n+1;xspread=max(xspread,maxval(abs(simp(:,j)-simp(:,ilo))));end do
            if(spread<eps .and. xspread<sqrt(eps)*max(1.0_dp,maxval(abs(simp(:,ilo))))) exit
            cent=0.0_dp
            do j=1,n+1;if(j/=ihi)cent=cent+simp(:,j);end do
            cent=cent/real(n,dp)
            xr=cent+(cent-simp(:,ihi));call enforce_bounds(xr,lower,upper);fr=fn(xr)
            if(fr<f(ilo)) then
                xe=cent+2.0_dp*(xr-cent);call enforce_bounds(xe,lower,upper);fe=fn(xe)
                if(fe<fr) then;simp(:,ihi)=xe;f(ihi)=fe
                else;simp(:,ihi)=xr;f(ihi)=fr;end if
            else if(fr<f(inhi)) then
                simp(:,ihi)=xr;f(ihi)=fr
            else
                if(fr<f(ihi)) then;xc=cent+0.5_dp*(xr-cent)
                else;xc=cent+0.5_dp*(simp(:,ihi)-cent);end if
                call enforce_bounds(xc,lower,upper);fc=fn(xc)
                if(fc<min(fr,f(ihi))) then
                    simp(:,ihi)=xc;f(ihi)=fc
                else
                    do j=1,n+1
                        if(j/=ilo) then
                            simp(:,j)=simp(:,ilo)+0.5_dp*(simp(:,j)-simp(:,ilo))
                            call enforce_bounds(simp(:,j),lower,upper);f(j)=fn(simp(:,j))
                        end if
                    end do
                end if
            end if
        end do
        ilo=minloc(f,dim=1);res%estimate=simp(:,ilo);res%distance=f(ilo)
        res%iterations=min(it,mi);res%converged=(it<=mi .and. finite(res%distance))
    end subroutine nelder_mead

    subroutine enforce_bounds(x,lower,upper)
        real(dp),intent(inout)::x(:)
        real(dp),intent(in),optional::lower(:),upper(:)
        if(present(lower)) then
            if(size(lower)==1) then;x=max(x,lower(1))
            else if(size(lower)==size(x)) then;x=max(x,lower);end if
        end if
        if(present(upper)) then
            if(size(upper)==1) then;x=min(x,upper(1))
            else if(size(upper)==size(x)) then;x=min(x,upper);end if
        end if
    end subroutine enforce_bounds

    pure logical function finite(x)
        real(dp),intent(in)::x
        finite=(x==x .and. abs(x)<huge(x))
    end function finite

    subroutine sort_real(x)
        real(dp),intent(inout)::x(:)
        integer::i,j
        real(dp)::v
        do i=2,size(x)
            v=x(i);j=i-1
            do while(j>=1)
                if(x(j)<=v)exit
                x(j+1)=x(j);j=j-1
            end do
            x(j+1)=v
        end do
    end subroutine sort_real
end module actuar_mde_v03
