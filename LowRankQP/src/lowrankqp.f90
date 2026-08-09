! SPDX-License-Identifier: GPL-2.0-or-later
module lowrankqp
    use lowrankqp_kinds, only : dp
    use lowrankqp_linalg, only : chol_factor, chol_solve, lu_factor, lu_solve
    implicit none
    private

    integer, parameter, public :: LRQP_LU = 1
    integer, parameter, public :: LRQP_CHOL = 2
    integer, parameter, public :: LRQP_SMW = 3
    integer, parameter, public :: LRQP_PFCF = 4

    real(dp), parameter :: EPS_INIT = 1.0e-2_dp
    real(dp), parameter :: EPS_IPM = 1.0e-2_dp
    real(dp), parameter :: EPS_PERT = 1.0e-14_dp

    type, public :: lowrankqp_options
        integer :: method = LRQP_PFCF
        integer :: max_iter = 200
        real(dp) :: tol = 1.0e-8_dp
        logical :: verbose = .false.
    end type lowrankqp_options

    type, public :: lowrankqp_result
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: xi(:)
        real(dp), allocatable :: zeta(:)
        integer :: iterations = 0
        integer :: status = 0
        logical :: converged = .false.
        real(dp) :: primal_feasibility = huge(1.0_dp)
        real(dp) :: dual_feasibility = huge(1.0_dp)
        real(dp) :: complementarity = huge(1.0_dp)
        real(dp) :: duality_gap = huge(1.0_dp)
        real(dp) :: termination = huge(1.0_dp)
    end type lowrankqp_result

    type :: factor_state
        integer :: method = 0
        real(dp), allocatable :: mfac(:,:)
        integer, allocatable :: piv(:)
        real(dp), allocatable :: p(:,:), beta(:,:), lambda(:)
    end type factor_state

    public :: solve_low_rank_qp, lowrankqp_objective, lowrankqp_method

contains

pure integer function lowrankqp_method(name) result(method)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: key
    key = trim(adjustl(name))
    select case(key)
    case('LU','lu')
        method = LRQP_LU
    case('CHOL','chol','Chol')
        method = LRQP_CHOL
    case('SMW','smw')
        method = LRQP_SMW
    case('PFCF','pfcf')
        method = LRQP_PFCF
    case default
        method = 0
    end select
end function lowrankqp_method

pure real(dp) function lowrankqp_objective(v, d, alpha) result(value)
    real(dp), intent(in) :: v(:,:), d(:), alpha(:)
    real(dp), allocatable :: w(:)
    integer :: n, m
    n = size(v,1)
    m = size(v,2)
    if (m == n) then
        value = dot_product(d,alpha) + 0.5_dp*dot_product(alpha,matmul(v,alpha))
    else
        w = matmul(transpose(v),alpha)
        value = dot_product(d,alpha) + 0.5_dp*dot_product(w,w)
    end if
end function lowrankqp_objective

subroutine solve_low_rank_qp(v, d, a, b, u, result, options)
    real(dp), intent(in) :: v(:,:), d(:), a(:,:), b(:), u(:)
    type(lowrankqp_result), intent(out) :: result
    type(lowrankqp_options), intent(in), optional :: options
    type(lowrankqp_options) :: opt
    integer :: n, m, p, iter, info
    real(dp) :: mult, t
    real(dp), allocatable :: alpha(:), beta(:), xi(:), zeta(:)
    real(dp), allocatable :: da(:), db(:), dxi(:), dzeta(:)
    real(dp), allocatable :: uma(:), xiu(:), za(:), w(:), r1(:), r2(:), diagd(:)
    type(factor_state) :: fac

    opt = lowrankqp_options()
    if (present(options)) opt = options
    n = size(v,1); m = size(v,2); p = size(a,1)
    if (size(d) /= n .or. size(u) /= n .or. size(a,2) /= n .or. size(b) /= p) then
        result%status = -1
        return
    end if
    if (n < 1 .or. m < 1 .or. any(u <= 0.0_dp) .or. opt%max_iter < 1 .or. opt%tol <= 0.0_dp) then
        result%status = -2
        return
    end if
    if (opt%method < LRQP_LU .or. opt%method > LRQP_PFCF) then
        result%status = -3
        return
    end if

    allocate(alpha(n), xi(n), zeta(n), da(n), dxi(n), dzeta(n), uma(n), xiu(n), za(n), w(m), r1(n), diagd(n))
    allocate(beta(p), db(p), r2(p))
    beta = 0.0_dp; db = 0.0_dp; r2 = 0.0_dp
    call init_point(v,d,u,alpha,xi,zeta,w,r1)
    mult = 0.0_dp
    info = 0

    do iter = 1, opt%max_iter
        call calc_stats(v,d,a,b,u,alpha,beta,xi,zeta,mult,uma,xiu,za,w,r1,r2,diagd, &
            result%primal_feasibility,result%dual_feasibility,result%complementarity, &
            result%duality_gap,result%termination,t)
        if (opt%verbose) then
            write(*,'(i4,5(1x,es14.6))') iter, result%primal_feasibility, result%dual_feasibility, &
                result%complementarity, result%duality_gap, result%termination
        end if
        result%iterations = iter
        if (result%termination < opt%tol) then
            result%converged = .true.
            exit
        end if
        call factorize_system(v,diagd,opt%method,fac,info)
        if (info /= 0) then
            result%status = 2
            exit
        end if
        call calc_direction(v,a,alpha,xi,zeta,uma,xiu,za,r1,r2,diagd,fac,t,.false., &
            da,db,dxi,dzeta,info)
        if (info /= 0) then
            result%status = 3
            exit
        end if
        call calc_direction(v,a,alpha,xi,zeta,uma,xiu,za,r1,r2,diagd,fac,t,.true., &
            da,db,dxi,dzeta,info)
        if (info /= 0) then
            result%status = 3
            exit
        end if
        call take_step(alpha,beta,xi,zeta,da,db,dxi,dzeta,uma,mult)
    end do
    if (result%converged) result%status = 0
    if (.not. result%converged .and. result%status == 0) result%status = 1

    result%alpha = alpha
    result%beta = beta
    result%xi = xi
    result%zeta = zeta
end subroutine solve_low_rank_qp

subroutine init_point(v,d,u,alpha,xi,zeta,w,temp)
    real(dp), intent(in) :: v(:,:), d(:), u(:)
    real(dp), intent(out) :: alpha(:), xi(:), zeta(:), w(:), temp(:)
    integer :: n, m
    n=size(v,1); m=size(v,2)
    alpha = min(EPS_INIT,u*EPS_INIT)
    w = matmul(transpose(v),alpha)
    if (n == m) then
        temp = -w-d
    else
        temp = -matmul(v,w)-d
    end if
    xi = max(EPS_INIT,temp)
    zeta = max(EPS_INIT,xi-temp)
end subroutine init_point

subroutine calc_stats(v,d,a,b,u,alpha,beta,xi,zeta,mult,uma,xiu,za,w,r1,r2,diagd, &
    prim,dual,comp,gap,term,t)
    real(dp), intent(in) :: v(:,:), d(:), a(:,:), b(:), u(:), alpha(:), beta(:), xi(:), zeta(:), mult
    real(dp), intent(out) :: uma(:), xiu(:), za(:), w(:), r1(:), r2(:), diagd(:)
    real(dp), intent(out) :: prim,dual,comp,gap,term,t
    real(dp) :: quad, ca, tmp
    integer :: n,m,p
    n=size(v,1); m=size(v,2); p=size(a,1)
    w=matmul(transpose(v),alpha)
    uma=u-alpha
    xiu=xi/uma
    za=zeta/alpha
    if (n == m) then
        r1=-d-w-xi+zeta
        quad=dot_product(alpha,w)
    else
        r1=-d-matmul(v,w)-xi+zeta
        quad=dot_product(w,w)
    end if
    if (p > 0) then
        r1=r1-matmul(transpose(a),beta)
        r2=b-matmul(a,alpha)
        dual=sum(abs(r2))
    else
        dual=0.0_dp
    end if
    prim=sum(abs(r1))
    comp=dot_product(alpha,zeta)+dot_product(uma,xi)
    ca=dot_product(d,alpha)
    gap=abs(quad+ca+dot_product(u,xi))
    if (p > 0) gap=abs(quad+ca+dot_product(u,xi)+dot_product(b,beta))
    term=comp/(abs(0.5_dp*quad+ca)+1.0_dp)
    tmp=(1.0_dp-mult+EPS_IPM)/(10.0_dp+mult)
    t=comp*tmp*tmp/(2.0_dp*real(n,dp))
    diagd=xiu+za+EPS_PERT
end subroutine calc_stats

subroutine factorize_system(v,dg,method,fac,info)
    real(dp), intent(in) :: v(:,:), dg(:)
    integer, intent(in) :: method
    type(factor_state), intent(out) :: fac
    integer, intent(out) :: info
    integer :: n,m,i
    real(dp), allocatable :: dinvq(:,:)
    n=size(v,1); m=size(v,2)
    fac%method=method; info=0
    select case(method)
    case(LRQP_LU,LRQP_CHOL)
        allocate(fac%mfac(n,n))
        if (n == m) then
            fac%mfac=v
        else
            fac%mfac=matmul(v,transpose(v))
        end if
        do i=1,n
            fac%mfac(i,i)=fac%mfac(i,i)+dg(i)
        end do
        if (method == LRQP_LU) then
            allocate(fac%piv(n))
            call lu_factor(fac%mfac,fac%piv,info)
        else
            call chol_factor(fac%mfac,info)
        end if
    case(LRQP_SMW)
        allocate(dinvq(n,m),fac%mfac(m,m))
        do i=1,n
            dinvq(i,:)=v(i,:)/dg(i)
        end do
        fac%mfac=matmul(transpose(v),dinvq)
        do i=1,m
            fac%mfac(i,i)=fac%mfac(i,i)+1.0_dp
        end do
        call chol_factor(fac%mfac,info)
    case(LRQP_PFCF)
        call pfcf_factorize(v,dg,fac%p,fac%beta,fac%lambda)
    end select
end subroutine factorize_system

subroutine solve_system(v,dg,fac,rhs,sol,info)
    real(dp), intent(in) :: v(:,:), dg(:), rhs(:,:)
    type(factor_state), intent(in) :: fac
    real(dp), intent(out) :: sol(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: tmp(:,:)
    integer :: i,j,n,m
    n=size(v,1); m=size(v,2); info=0
    sol=rhs
    select case(fac%method)
    case(LRQP_LU)
        call lu_solve(fac%mfac,fac%piv,sol,info)
    case(LRQP_CHOL)
        call chol_solve(fac%mfac,sol,info)
    case(LRQP_SMW)
        do i=1,n
            sol(i,:)=sol(i,:)/dg(i)
        end do
        tmp=matmul(transpose(v),sol)
        call chol_solve(fac%mfac,tmp,info)
        if (info /= 0) return
        sol=rhs-matmul(v,tmp)
        do i=1,n
            sol(i,:)=sol(i,:)/dg(i)
        end do
    case(LRQP_PFCF)
        do j=1,size(sol,2)
            do i=1,m
                call pfcf_solve(fac%p(:,i),fac%beta(:,i),sol(:,j),.false.)
            end do
            sol(:,j)=sol(:,j)/fac%lambda
            do i=m,1,-1
                call pfcf_solve(fac%p(:,i),fac%beta(:,i),sol(:,j),.true.)
            end do
        end do
    end select
end subroutine solve_system

subroutine calc_direction(v,a,alpha,xi,zeta,uma,xiu,za,r1,r2,dg,fac,t,corrector, &
    da,db,dxi,dzeta,info)
    real(dp), intent(in) :: v(:,:), a(:,:), alpha(:), xi(:), zeta(:), uma(:), xiu(:), za(:), r1(:), r2(:), dg(:), t
    type(factor_state), intent(in) :: fac
    logical, intent(in) :: corrector
    real(dp), intent(inout) :: da(:), db(:), dxi(:), dzeta(:)
    integer, intent(out) :: info
    real(dp), allocatable :: r3(:),r4(:),r5(:),rr(:,:),rhs1(:,:),rmat(:,:),schur(:,:),brhs(:,:)
    integer :: n,p
    n=size(alpha); p=size(a,1); info=0
    allocate(r3(n),r4(n),r5(n))
    r3=-zeta; r4=-xi
    if (corrector) then
        r3=r3+(t-da*dzeta)/alpha
        r4=r4+(t+da*dxi)/uma
    end if
    r5=r1+r3-r4
    if (p > 0) then
        allocate(rmat(n,p),rr(n,1),rhs1(n,1),schur(p,p),brhs(p,1))
        call solve_system(v,dg,fac,transpose(a),rmat,info)
        if (info /= 0) return
        rhs1(:,1)=r5
        call solve_system(v,dg,fac,rhs1,rr,info)
        if (info /= 0) return
        schur=matmul(a,rmat)
        brhs(:,1)=matmul(a,rr(:,1))-r2
        call chol_factor(schur,info)
        if (info /= 0) return
        call chol_solve(schur,brhs,info)
        if (info /= 0) return
        db=brhs(:,1)
        da=rr(:,1)-matmul(rmat,db)
    else
        allocate(rr(n,1),rhs1(n,1)); rhs1(:,1)=r5
        call solve_system(v,dg,fac,rhs1,rr,info)
        if (info /= 0) return
        da=rr(:,1)
    end if
    dzeta=r3-za*da
    dxi=r4+xiu*da
end subroutine calc_direction

subroutine take_step(alpha,beta,xi,zeta,da,db,dxi,dzeta,uma,mult)
    real(dp), intent(inout) :: alpha(:),beta(:),xi(:),zeta(:)
    real(dp), intent(in) :: da(:),db(:),dxi(:),dzeta(:),uma(:)
    real(dp), intent(out) :: mult
    integer :: i
    mult=1.0_dp
    do i=1,size(alpha)
        if (da(i) < 0.0_dp) mult=min(mult,-alpha(i)/da(i))
        if (da(i) > 0.0_dp) mult=min(mult,uma(i)/da(i))
        if (dxi(i) < 0.0_dp) mult=min(mult,-xi(i)/dxi(i))
        if (dzeta(i) < 0.0_dp) mult=min(mult,-zeta(i)/dzeta(i))
    end do
    mult=0.99_dp*mult
    alpha=alpha+mult*da
    if (size(beta)>0) beta=beta+mult*db
    xi=xi+mult*dxi
    zeta=zeta+mult*dzeta
end subroutine take_step

subroutine pfcf_solve(p,beta,r,trans)
    real(dp), intent(in) :: p(:),beta(:)
    real(dp), intent(inout) :: r(:)
    logical, intent(in) :: trans
    integer :: j,n
    real(dp) :: sigma
    n=size(r)
    if (n == 1) return
    if (trans) then
        sigma=r(n)*p(n)
        do j=n-1,2,-1
            r(j)=r(j)-sigma*beta(j)
            sigma=sigma+r(j)*p(j)
        end do
        r(1)=r(1)-sigma*beta(1)
    else
        sigma=r(1)*beta(1)
        do j=2,n-1
            r(j)=r(j)-sigma*p(j)
            sigma=sigma+r(j)*beta(j)
        end do
        r(n)=r(n)-sigma*p(n)
    end if
end subroutine pfcf_solve

subroutine pfcf_factorize(q,d,p,beta,lambda)
    real(dp), intent(in) :: q(:,:),d(:)
    real(dp), allocatable, intent(out) :: p(:,:),beta(:,:),lambda(:)
    real(dp), allocatable :: lt(:),t(:)
    real(dp) :: tmp
    integer :: i,j,n,m
    real(dp), parameter :: infinity=1.0e12_dp, epsilon=1.0e-12_dp
    n=size(q,1); m=size(q,2)
    allocate(p(n,m),beta(n,m),lambda(n),lt(n),t(n+1))
    p=q; beta=0.0_dp; lt=d; lambda=d
    do i=1,m
        do j=1,i-1
            call pfcf_solve(p(:,j),beta(:,j),p(:,i),.false.)
        end do
        t(1)=1.0_dp
        do j=1,n
            tmp=p(j,i)
            if (abs(t(j)) < infinity) then
                if (abs(lt(j)) > epsilon) then
                    t(j+1)=t(j)+tmp*tmp/lt(j)
                    lambda(j)=lt(j)*t(j+1)/t(j)
                    beta(j,i)=tmp/(lt(j)*t(j+1))
                else if (abs(tmp) > epsilon) then
                    t(j+1)=infinity
                    lambda(j)=tmp*tmp/t(j)
                    beta(j,i)=1.0_dp/tmp
                else
                    t(j+1)=t(j); lambda(j)=0.0_dp; beta(j,i)=0.0_dp
                end if
            else
                t(j+1)=infinity; lambda(j)=lt(j); beta(j,i)=0.0_dp
            end if
        end do
        lt=lambda
    end do
end subroutine pfcf_factorize

end module lowrankqp
