! Upstream license declaration: GPL (version unspecified)
module limsolve_ranges
    use limsolve_kinds, only: dp
    use limsolve_types
    use lpsolve, only: lp_result, solve_lp, LP_MIN, LP_MAX, LP_EQ, LP_GE, &
        LP_OPTIMAL, LP_SUBOPTIMAL, LP_INFEASIBLE, LP_UNBOUNDED
    implicit none
    private
    public :: linp, xranges, varranges, varsample

contains

    subroutine linp(e,f,g,h,cost,result,ispos,int_vec,lower,upper)
        real(dp), intent(in) :: e(:,:),f(:),g(:,:),h(:),cost(:)
        type(solve_result), intent(out) :: result
        logical, intent(in), optional :: ispos
        integer, intent(in), optional :: int_vec(:)
        real(dp), intent(in), optional :: lower(:),upper(:)
        real(dp), allocatable :: gg(:,:),hh(:),con(:,:),con2(:,:),rhs(:),obj(:)
        integer, allocatable :: sense(:),ints(:)
        type(lp_result) :: lpres
        logical :: positive
        integer :: n,me,mg,m

        n=size(cost); me=size(e,1); mg=size(g,1); positive=.true.
        if(present(ispos)) positive=ispos
        call init_result(result,n)
        if(size(e,2)/=n .or. size(g,2)/=n .or. size(f)/=me .or. size(h)/=mg) return
        call augment_bounds_local(g,h,n,gg,hh,lower,upper)
        m=me+size(gg,1)
        allocate(con(m,n),rhs(m),sense(m))
        if(me>0) then
            con(1:me,:)=e; rhs(1:me)=f; sense(1:me)=LP_EQ
        end if
        if(size(gg,1)>0) then
            con(me+1:m,:)=gg; rhs(me+1:m)=hh; sense(me+1:m)=LP_GE
        end if
        if(positive) then
            obj=cost
            if(present(int_vec)) then
                call solve_lp(LP_MIN,obj,con,sense,rhs,lpres,integer_variables=int_vec)
            else
                call solve_lp(LP_MIN,obj,con,sense,rhs,lpres)
            end if
            if(allocated(lpres%solution)) result%x=lpres%solution
        else
            allocate(obj(2*n)); obj(1:n)=cost; obj(n+1:2*n)=-cost
            allocate(con2(m,2*n))
            con2(:,1:n)=con
            con2(:,n+1:2*n)=-con
            call move_alloc(con2,con)
            if(present(int_vec)) then
                allocate(ints(2*size(int_vec)))
                ints(1:size(int_vec))=int_vec
                ints(size(int_vec)+1:)=int_vec+n
                call solve_lp(LP_MIN,obj,con,sense,rhs,lpres,integer_variables=ints)
            else
                call solve_lp(LP_MIN,obj,con,sense,rhs,lpres)
            end if
            if(allocated(lpres%solution)) result%x=lpres%solution(1:n)-lpres%solution(n+1:2*n)
        end if
        result%solution_norm=lpres%objective
        result%residual_norm=0.0_dp
        if(me>0) result%residual_norm=result%residual_norm+sum(abs(matmul(e,result%x)-f))
        if(size(gg,1)>0) result%residual_norm=result%residual_norm-sum(min(matmul(gg,result%x)-hh,0.0_dp))
        result%numiter=lpres%simplex_iterations
        select case(lpres%status)
        case(LP_OPTIMAL,LP_SUBOPTIMAL)
            result%status=LS_SUCCESS; result%is_error=.false.
        case(LP_INFEASIBLE)
            result%status=LS_INFEASIBLE; result%is_error=.true.
        case(LP_UNBOUNDED)
            result%status=LS_UNBOUNDED; result%is_error=.true.
        case default
            result%status=LS_NUMERICAL; result%is_error=.true.
        end select
    end subroutine linp

    function xranges(e,f,g,h,ispos,central,full,lower,upper) result(res)
        real(dp), intent(in) :: e(:,:),f(:),g(:,:),h(:)
        logical, intent(in), optional :: ispos,central,full
        real(dp), intent(in), optional :: lower(:),upper(:)
        type(range_result) :: res
        type(solve_result) :: smin,smax
        real(dp), allocatable :: c(:),solutions(:,:)
        logical :: pos,wantc,wantf
        integer :: n,i,ns
        n=max(size(e,2),size(g,2)); pos=.false.; if(present(ispos)) pos=ispos
        wantc=.false.; if(present(central)) wantc=central
        wantf=.false.; if(present(full)) wantf=full
        allocate(res%range(n,2)); res%range=0.0_dp
        if(wantc) then; allocate(res%central(n)); res%central=0.0_dp; end if
        if(wantf) then; allocate(solutions(n,2*n)); solutions=0.0_dp; end if
        allocate(c(n)); ns=0; res%status=LS_SUCCESS
        do i=1,n
            c=0.0_dp; c(i)=1.0_dp
            call linp_bounds_dispatch(e,f,g,h,c,smin,pos,lower,upper)
            c(i)=-1.0_dp
            call linp_bounds_dispatch(e,f,g,h,c,smax,pos,lower,upper)
            if(smin%status==LS_SUCCESS) then
                res%range(i,1)=smin%x(i)
                if(wantc) res%central=res%central+smin%x
                if(wantf) then; ns=ns+1; solutions(:,ns)=smin%x; end if
            else if(smin%status==LS_UNBOUNDED) then
                res%range(i,1)=-huge(1.0_dp)
            else
                res%range(i,1)=huge(1.0_dp); res%status=smin%status
            end if
            if(smax%status==LS_SUCCESS) then
                res%range(i,2)=smax%x(i)
                if(wantc) res%central=res%central+smax%x
                if(wantf) then; ns=ns+1; solutions(:,ns)=smax%x; end if
            else if(smax%status==LS_UNBOUNDED) then
                res%range(i,2)=huge(1.0_dp)
            else
                res%range(i,2)=-huge(1.0_dp); res%status=smax%status
            end if
        end do
        if(wantc) res%central=res%central/real(max(1,2*n),dp)
        if(wantf) then
            allocate(res%all_x(n,ns)); if(ns>0) res%all_x=solutions(:,1:ns)
        end if
    end function xranges

    function varranges(e,f,g,h,eqa,eqb,ispos,lower,upper) result(res)
        real(dp), intent(in) :: e(:,:),f(:),g(:,:),h(:),eqa(:,:)
        real(dp), intent(in), optional :: eqb(:)
        logical, intent(in), optional :: ispos
        real(dp), intent(in), optional :: lower(:),upper(:)
        type(range_result) :: res
        type(solve_result) :: smin,smax
        real(dp), allocatable :: c(:)
        logical :: pos
        integer :: nv,n,i
        n=size(eqa,2); nv=size(eqa,1); pos=.false.; if(present(ispos)) pos=ispos
        allocate(res%range(nv,2),c(n)); res%range=0.0_dp; res%status=LS_SUCCESS
        do i=1,nv
            c=eqa(i,:)
            call linp_bounds_dispatch(e,f,g,h,c,smin,pos,lower,upper)
            c=-eqa(i,:)
            call linp_bounds_dispatch(e,f,g,h,c,smax,pos,lower,upper)
            if(smin%status==LS_SUCCESS) then
                res%range(i,1)=dot_product(eqa(i,:),smin%x)
            else if(smin%status==LS_UNBOUNDED) then
                res%range(i,1)=-huge(1.0_dp)
            else
                res%status=smin%status
            end if
            if(smax%status==LS_SUCCESS) then
                res%range(i,2)=dot_product(eqa(i,:),smax%x)
            else if(smax%status==LS_UNBOUNDED) then
                res%range(i,2)=huge(1.0_dp)
            else
                res%status=smax%status
            end if
        end do
        if(present(eqb)) then
            if(size(eqb)==nv) then
                res%range(:,1)=res%range(:,1)-eqb
                res%range(:,2)=res%range(:,2)-eqb
            end if
        end if
    end function varranges

    subroutine varsample(x,eqa,var,eqb,status)
        real(dp), intent(in) :: x(:,:),eqa(:,:)
        real(dp), intent(out) :: var(:,:)
        real(dp), intent(in), optional :: eqb(:)
        integer, intent(out), optional :: status
        integer :: i
        if(size(x,2)/=size(eqa,2) .or. size(var,1)/=size(x,1) .or. size(var,2)/=size(eqa,1)) then
            var=0.0_dp; if(present(status)) status=LS_INVALID; return
        end if
        var=matmul(x,transpose(eqa))
        if(present(eqb)) then
            if(size(eqb)==size(eqa,1)) then
                do i=1,size(var,1); var(i,:)=var(i,:)-eqb; end do
            end if
        end if
        if(present(status)) status=LS_SUCCESS
    end subroutine varsample

    subroutine linp_bounds_dispatch(e,f,g,h,c,res,pos,lower,upper)
        real(dp),intent(in)::e(:,:),f(:),g(:,:),h(:),c(:)
        type(solve_result),intent(out)::res
        logical,intent(in)::pos
        real(dp),intent(in),optional::lower(:),upper(:)
        if(present(lower).and.present(upper)) then
            call linp(e,f,g,h,c,res,pos,lower=lower,upper=upper)
        else if(present(lower)) then
            call linp(e,f,g,h,c,res,pos,lower=lower)
        else if(present(upper)) then
            call linp(e,f,g,h,c,res,pos,upper=upper)
        else
            call linp(e,f,g,h,c,res,pos)
        end if
    end subroutine linp_bounds_dispatch

    subroutine augment_bounds_local(g,h,n,gg,hh,lower,upper)
        real(dp),intent(in)::g(:,:),h(:); integer,intent(in)::n
        real(dp),allocatable,intent(out)::gg(:,:),hh(:)
        real(dp),intent(in),optional::lower(:),upper(:)
        integer::m,nl,nu,i,k; logical::scalar
        m=size(g,1); nl=0; nu=0
        if(present(lower)) then
            if(size(lower)==1) then; if(abs(lower(1))<huge(1.0_dp)/10.0_dp) nl=n
            else; nl=count(abs(lower)<huge(1.0_dp)/10.0_dp); end if
        end if
        if(present(upper)) then
            if(size(upper)==1) then; if(abs(upper(1))<huge(1.0_dp)/10.0_dp) nu=n
            else; nu=count(abs(upper)<huge(1.0_dp)/10.0_dp); end if
        end if
        allocate(gg(m+nl+nu,n),hh(m+nl+nu)); gg=0.0_dp; hh=0.0_dp
        if(m>0) then; gg(1:m,:)=g; hh(1:m)=h; end if; k=m
        if(present(lower)) then
            scalar=size(lower)==1
            do i=1,n
                if(scalar) then
                    if(abs(lower(1))>=huge(1.0_dp)/10.0_dp) cycle
                    k=k+1; gg(k,i)=1.0_dp; hh(k)=lower(1)
                else if(i<=size(lower)) then
                    if(abs(lower(i))>=huge(1.0_dp)/10.0_dp) cycle
                    k=k+1; gg(k,i)=1.0_dp; hh(k)=lower(i)
                end if
            end do
        end if
        if(present(upper)) then
            scalar=size(upper)==1
            do i=1,n
                if(scalar) then
                    if(abs(upper(1))>=huge(1.0_dp)/10.0_dp) cycle
                    k=k+1; gg(k,i)=-1.0_dp; hh(k)=-upper(1)
                else if(i<=size(upper)) then
                    if(abs(upper(i))>=huge(1.0_dp)/10.0_dp) cycle
                    k=k+1; gg(k,i)=-1.0_dp; hh(k)=-upper(i)
                end if
            end do
        end if
    end subroutine augment_bounds_local

    subroutine init_result(result,n)
        type(solve_result),intent(out)::result; integer,intent(in)::n
        allocate(result%x(max(0,n))); result%x=0.0_dp
        result%status=LS_INVALID; result%is_error=.true.; result%residual_norm=0.0_dp
        result%solution_norm=0.0_dp; result%numiter=0
    end subroutine init_result

end module limsolve_ranges
