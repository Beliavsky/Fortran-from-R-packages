! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_solver
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use ecos_types, only : dp, ecos_problem, ecos_result, ecos_settings, ecos_sparse_cache, &
        ECOS_OPTIMAL, ECOS_PINF, ECOS_DINF, ECOS_INACC_OFFSET, ECOS_MAXIT, ECOS_NUMERICS
    use ecos_linalg, only : solve_kkt, least_norm_equalities, vecnorm2
    use ecos_cones, only : cone_scalar_eval, cone_slack, cone_dual_from_scalar, cone_violation
    use ecos_sparse_solver, only : solve_pd_sparse, densify_problem
    use ecos_sparse, only : sparse_triplet_builder, triplet_to_csc, csc_to_csr
    implicit none
    private
    public :: ecos_solve_continuous

contains

    recursive subroutine ecos_solve_continuous(prob, result, settings, diagnose, cache)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        type(ecos_settings), intent(in), optional :: settings
        logical, intent(in), optional :: diagnose
        type(ecos_sparse_cache), intent(inout), optional :: cache
        type(ecos_settings) :: stg
        logical :: do_diag
        integer :: raw_status
        type(ecos_problem) :: dprob

        stg = ecos_settings()
        if (present(settings)) stg = settings
        do_diag = .true.
        if (present(diagnose)) do_diag = diagnose

        if (.not. prob%valid()) then
            call allocate_result(prob,result)
            result%exitflag = ECOS_NUMERICS
            result%numerr = 1
            result%infostring = 'Invalid problem data'
            return
        end if

        if (prob%sparse_storage) then
            if (stg%sparse_kkt .and. prob%dims%scalar_inequalities() > 0) then
                if (present(cache)) then
                    call solve_pd_sparse(prob,result,stg,raw_status,cache)
                else
                    call solve_pd_sparse(prob,result,stg,raw_status)
                end if
            else
                call densify_problem(prob,dprob)
                call solve_pd(dprob,result,stg,raw_status)
            end if
        else
            call solve_pd(prob,result,stg,raw_status)
        end if
        if (raw_status == ECOS_OPTIMAL .or. raw_status == ECOS_OPTIMAL+ECOS_INACC_OFFSET) return

        if (do_diag) then
            if (prob%sparse_storage) then
                if (sparse_unbounded_certificate(prob,stg,result)) then
                    result%exitflag = ECOS_DINF
                    result%infostring = 'Dual infeasible / primal unbounded (sparse ray certificate)'
                    result%dinf = 1.0_dp
                    result%pinf = 0.0_dp
                    call sanitize_result(result)
                    return
                end if
                if (sparse_primal_infeasible_certificate(prob,stg,result)) then
                    result%exitflag = ECOS_PINF
                    result%infostring = 'Primal infeasible (sparse dual-ray certificate)'
                    result%pinf = 1.0_dp
                    result%dinf = 0.0_dp
                    call sanitize_result(result)
                    return
                end if
                if (prob%nvar()+prob%ncone()+prob%neq() <= stg%dense_diagnostic_limit) then
                    call densify_problem(prob,dprob)
                    if (is_unbounded_direction(dprob,stg)) then
                        result%exitflag = ECOS_DINF
                        result%infostring = 'Dual infeasible / primal unbounded'
                        result%dinf = 1.0_dp
                        call sanitize_result(result)
                        return
                    end if
                    if (.not. phase1_feasible(dprob,stg)) then
                        result%exitflag = ECOS_PINF
                        result%infostring = 'Primal infeasible'
                        result%pinf = 1.0_dp
                        call sanitize_result(result)
                        return
                    end if
                end if
                return
            end if
            if (is_unbounded_direction(prob,stg)) then
                result%exitflag = ECOS_DINF
                result%infostring = 'Dual infeasible / primal unbounded'
                result%dinf = 1.0_dp
                result%pinf = 0.0_dp
                call sanitize_result(result)
                return
            end if
            if (.not. phase1_feasible(prob,stg)) then
                result%exitflag = ECOS_PINF
                result%infostring = 'Primal infeasible'
                result%pinf = 1.0_dp
                result%dinf = 0.0_dp
                call sanitize_result(result)
                return
            end if
        end if
    end subroutine ecos_solve_continuous


    subroutine initialize_sparse_problem_storage(prob, n, m, p)
        type(ecos_problem), intent(inout) :: prob
        integer, intent(in) :: n, m, p
        allocate(prob%c(n), prob%h(m), prob%b(p))
        prob%c = 0.0_dp; prob%h = 0.0_dp; prob%b = 0.0_dp
        prob%sparse_storage = .true.
    end subroutine initialize_sparse_problem_storage

    logical function sparse_unbounded_certificate(prob, stg, result) result(found)
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in) :: stg
        type(ecos_result), intent(inout) :: result
        type(ecos_problem) :: cp
        type(ecos_result) :: rr
        type(ecos_settings) :: cstg
        type(sparse_triplet_builder) :: tb
        integer :: n, m, p, i, j, k
        real(dp) :: cdot
        found = .false.
        n=prob%nvar(); m=prob%ncone(); p=prob%neq()
        if (m <= 0) return
        call initialize_sparse_problem_storage(cp,n,m,p+1)
        cp%dims=prob%dims
        cp%h=0.0_dp
        cp%b=0.0_dp; cp%b(p+1)=-1.0_dp
        cp%g_csc=prob%g_csc; call csc_to_csr(cp%g_csc,cp%g_csr)
        call tb%init(p+1,n,max(16,size(prob%a_csc%values)+n))
        do j=1,prob%a_csc%ncol
            do k=prob%a_csc%colptr(j),prob%a_csc%colptr(j+1)-1
                i=prob%a_csc%rowind(k)
                call tb%add(i,j,prob%a_csc%values(k))
            end do
        end do
        do j=1,n
            if (abs(prob%c(j))>0.0_dp) call tb%add(p+1,j,prob%c(j))
        end do
        call triplet_to_csc(tb,cp%a_csc); call csc_to_csr(cp%a_csc,cp%a_csr)
        cstg=stg; cstg%maxit=min(max(20,stg%certificate_maxit),200)
        call ecos_solve_continuous(cp,rr,cstg,.false.)
        if (rr%exitflag /= ECOS_OPTIMAL) return
        cdot=dot_product(prob%c,rr%x)
        if (cdot < -0.5_dp .and. cone_violation(cp,rr%x) <= 1.0e-5_dp) then
            if (allocated(result%dual_certificate)) deallocate(result%dual_certificate)
            allocate(result%dual_certificate(n)); result%dual_certificate=rr%x
            result%dual_certificate_valid=.true.; result%dinfres=max(rr%pres,rr%dres)
            found=.true.
        end if
    end function sparse_unbounded_certificate

    logical function sparse_primal_infeasible_certificate(prob, stg, result) result(found)
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in) :: stg
        type(ecos_result), intent(inout) :: result
        type(ecos_problem) :: cp
        type(ecos_result) :: rr
        type(ecos_settings) :: cstg
        type(sparse_triplet_builder) :: tg, ta
        integer :: n, m, p, nv, i, j, k, row, iq, qn, zcol
        real(dp) :: certdot
        real(dp), parameter :: euler = exp(1.0_dp)
        found=.false.
        n=prob%nvar(); m=prob%ncone(); p=prob%neq(); nv=p+m
        if (m <= 0 .or. nv <= 0) return
        call initialize_sparse_problem_storage(cp,nv,m,n+1)
        cp%dims=prob%dims; cp%h=0.0_dp; cp%b=0.0_dp; cp%b(n+1)=-1.0_dp

        ! Impose z in K*: linear and SOC blocks are self-dual.  For an
        ! exponential block (u,v,w), (-v,-u,e*w) lies in the primal cone.
        call tg%init(m,nv,max(16,m))
        row=1
        do i=1,prob%dims%l
            call tg%add(row,p+row,-1.0_dp); row=row+1
        end do
        if (allocated(prob%dims%q)) then
            do iq=1,size(prob%dims%q)
                qn=prob%dims%q(iq)
                do i=0,qn-1
                    call tg%add(row+i,p+row+i,-1.0_dp)
                end do
                row=row+qn
            end do
        end if
        do iq=1,prob%dims%e
            zcol=p+row
            call tg%add(row,zcol+1, 1.0_dp)
            call tg%add(row+1,zcol, 1.0_dp)
            call tg%add(row+2,zcol+2,-euler)
            row=row+3
        end do
        call triplet_to_csc(tg,cp%g_csc); call csc_to_csr(cp%g_csc,cp%g_csr)

        ! A^T y + G^T z = 0, then normalize b^T y + h^T z = -1.
        call ta%init(n+1,nv,max(32,size(prob%a_csc%values)+size(prob%g_csc%values)+nv))
        do j=1,prob%a_csc%ncol
            do k=prob%a_csc%colptr(j),prob%a_csc%colptr(j+1)-1
                i=prob%a_csc%rowind(k)
                call ta%add(j,i,prob%a_csc%values(k))
            end do
        end do
        do j=1,prob%g_csc%ncol
            do k=prob%g_csc%colptr(j),prob%g_csc%colptr(j+1)-1
                i=prob%g_csc%rowind(k)
                call ta%add(j,p+i,prob%g_csc%values(k))
            end do
        end do
        do i=1,p
            if (abs(prob%b(i))>0.0_dp) call ta%add(n+1,i,prob%b(i))
        end do
        do i=1,m
            if (abs(prob%h(i))>0.0_dp) call ta%add(n+1,p+i,prob%h(i))
        end do
        call triplet_to_csc(ta,cp%a_csc); call csc_to_csr(cp%a_csc,cp%a_csr)
        cstg=stg; cstg%maxit=min(max(20,stg%certificate_maxit),200)
        call ecos_solve_continuous(cp,rr,cstg,.false.)
        if (rr%exitflag /= ECOS_OPTIMAL) return
        certdot=0.0_dp
        if(p>0) certdot=certdot+dot_product(prob%b,rr%x(1:p))
        certdot=certdot+dot_product(prob%h,rr%x(p+1:p+m))
        if (certdot < -0.5_dp .and. rr%pres <= 1.0e-5_dp) then
            if (allocated(result%primal_certificate)) deallocate(result%primal_certificate)
            allocate(result%primal_certificate(nv)); result%primal_certificate=rr%x
            result%primal_certificate_valid=.true.; result%pinfres=max(rr%pres,rr%dres)
            found=.true.
        end if
    end function sparse_primal_infeasible_certificate

    subroutine sanitize_result(result)
        type(ecos_result), intent(inout) :: result
        integer :: i
        do i=1,size(result%x)
            if (.not.ieee_is_finite(result%x(i))) result%x(i)=0.0_dp
        end do
        do i=1,size(result%y)
            if (.not.ieee_is_finite(result%y(i))) result%y(i)=0.0_dp
        end do
        do i=1,size(result%s)
            if (.not.ieee_is_finite(result%s(i))) result%s(i)=0.0_dp
        end do
        do i=1,size(result%z)
            if (.not.ieee_is_finite(result%z(i))) result%z(i)=0.0_dp
        end do
    end subroutine sanitize_result

    subroutine allocate_result(prob,result)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        integer :: n,m,p
        n = prob%nvar()
        m = prob%ncone()
        p = prob%neq()
        allocate(result%x(n), result%y(p), result%s(m), result%z(m))
        result%x = 0.0_dp
        result%y = 0.0_dp
        result%s = 0.0_dp
        result%z = 0.0_dp
        result%r0 = 0.0_dp
    end subroutine allocate_result

    subroutine solve_pd(prob,result,stg,status)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        type(ecos_settings), intent(in) :: stg
        integer, intent(out) :: status
        real(dp), allocatable :: x(:), y(:), sl(:), lam(:)
        real(dp), allocatable :: g(:), jac(:,:), hlag(:,:), hbar(:,:)
        real(dp), allocatable :: req(:), rineq(:), rdual(:), rcent(:)
        real(dp), allocatable :: dx(:),dy(:),ds(:),dl(:)
        real(dp), allocatable :: dxa(:),dya(:),dsa(:),dla(:)
        real(dp), allocatable :: rhsx(:), w(:), tmp(:)
        real(dp) :: eqinit, pres, dres, mu, mua, sigma
        real(dp) :: ap, ad, apa, ada, scale_p, scale_d, reg, gap_est
        integer :: n,p,ni,i,it,info
        logical :: converged

        n = prob%nvar()
        p = prob%neq()
        ni = prob%dims%scalar_inequalities()
        call allocate_result(prob,result)
        allocate(x(n), y(p))
        y = 0.0_dp
        if (p > 0) then
            call least_norm_equalities(prob%amat,prob%b,x,eqinit,info)
            if (info /= 0) x = 0.0_dp
        else
            x = 0.0_dp
            eqinit = 0.0_dp
        end if
        if (prob%dims%e > 0) call repair_exp_log_arguments(prob,x)

        if (ni == 0) then
            call solve_no_inequalities(prob,result,stg,x)
            status = result%exitflag
            return
        end if

        allocate(sl(ni),lam(ni),g(ni),jac(ni,n),hlag(n,n),hbar(n,n))
        allocate(req(p),rineq(ni),rdual(n),rcent(ni))
        allocate(dx(n),dy(p),ds(ni),dl(ni),dxa(n),dya(p),dsa(ni),dla(ni))
        allocate(rhsx(n),w(ni),tmp(ni))
        lam = 1.0_dp
        call cone_scalar_eval(prob,x,lam,g,jac,hlag)
        sl = max(1.0_dp, -g + 1.0_dp)
        reg = max(stg%regularization,1.0e-12_dp)
        converged = .false.
        status = ECOS_MAXIT

        do it = 0, stg%maxit
            call cone_scalar_eval(prob,x,lam,g,jac,hlag)
            if (p > 0) then
                req = matmul(prob%amat,x)-prob%b
            end if
            rineq = g + sl
            rdual = prob%c + matmul(transpose(jac),lam)
            if (p > 0) rdual = rdual + matmul(transpose(prob%amat),y)
            mu = dot_product(sl,lam)/real(ni,dp)
            scale_p = 1.0_dp + max(vecnorm2(prob%h),merge(vecnorm2(prob%b),0.0_dp,p>0))
            scale_d = 1.0_dp + vecnorm2(prob%c)
            pres = vecnorm2(rineq)/scale_p
            if (p > 0) pres = max(pres,vecnorm2(req)/scale_p)
            dres = vecnorm2(rdual)/scale_d
            gap_est = real(ni,dp)*mu
            if (pres <= stg%feastol .and. dres <= stg%feastol .and. &
                (gap_est <= stg%abstol .or. gap_est <= stg%reltol*(1.0_dp+abs(dot_product(prob%c,x))))) then
                converged = .true.
                status = ECOS_OPTIMAL
                exit
            end if
            if (it >= stg%maxit) exit

            ! Affine predictor: sigma = 0.
            rcent = sl*lam
            call newton_direction(prob,jac,hlag,sl,lam,req,rineq,rdual,rcent,reg, &
                                  dxa,dya,dsa,dla,info)
            if (info /= 0) then
                reg = min(1.0e-2_dp,max(100.0_dp*reg,1.0e-8_dp))
                call newton_direction(prob,jac,hlag,sl,lam,req,rineq,rdual,rcent,reg, &
                                      dxa,dya,dsa,dla,info)
                if (info /= 0) then
                    status = ECOS_NUMERICS
                    exit
                end if
            end if
            apa = max_step(sl,dsa)
            ada = max_step(lam,dla)
            mua = dot_product(sl+apa*dsa,lam+ada*dla)/real(ni,dp)
            if (mu > tiny(1.0_dp)) then
                sigma = min(1.0_dp,max(0.0_dp,(mua/mu)**3))
            else
                sigma = 0.0_dp
            end if

            rcent = sl*lam + dsa*dla - sigma*mu
            call newton_direction(prob,jac,hlag,sl,lam,req,rineq,rdual,rcent,reg, &
                                  dx,dy,ds,dl,info)
            if (info /= 0) then
                status = ECOS_NUMERICS
                exit
            end if
            ap = min(1.0_dp,0.995_dp*max_step(sl,ds))
            ad = min(1.0_dp,0.995_dp*max_step(lam,dl))
            ! Protect the exponential cone's log arguments during the step.
            call limit_exp_step(prob,x,dx,ap)
            x = x + ap*dx
            sl = sl + ap*ds
            y = y + ad*dy
            lam = lam + ad*dl
            do i = 1, ni
                sl(i) = max(sl(i),1.0e-14_dp)
                lam(i) = max(lam(i),1.0e-14_dp)
            end do
        end do

        result%x = x
        result%y = y
        if (prob%ncone() > 0) then
            call cone_slack(prob,x,result%s)
            call cone_dual_from_scalar(prob,x,lam,result%z)
        end if
        result%iter = min(it,stg%maxit)
        result%exitflag = status
        if (status == ECOS_OPTIMAL) then
            result%infostring = 'Optimal solution found'
        else if (status == ECOS_MAXIT) then
            result%infostring = 'Maximum number of iterations reached'
        else
            result%infostring = 'Numerical problems encountered'
            result%numerr = 1
        end if
        call fill_summary(prob,result)
    end subroutine solve_pd


    subroutine repair_exp_log_arguments(prob,x)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(inout) :: x(:)
        real(dp), allocatable :: s(:), grad(:), proj(:), gram(:,:), rhs(:), yy(:)
        real(dp) :: delta, alpha, ng
        integer :: n,p,row,k,it,info,i
        n=size(x); p=prob%neq(); delta=1.0_dp
        allocate(s(prob%ncone()),grad(n),proj(n))
        if (p>0) then
            allocate(gram(p,p),rhs(p),yy(p))
            gram=matmul(prob%amat,transpose(prob%amat))
            do i=1,p
                gram(i,i)=gram(i,i)+1.0e-12_dp*max(1.0_dp,maxval(abs(gram)))
            end do
        end if
        row=prob%dims%l
        if (allocated(prob%dims%q)) row=row+sum(prob%dims%q)
        do it=1,40
            call cone_slack(prob,x,s)
            grad=0.0_dp
            row=prob%dims%l
            if (allocated(prob%dims%q)) row=row+sum(prob%dims%q)
            do k=1,prob%dims%e
                if (s(row+2)<delta) grad=grad+(delta-s(row+2))*prob%gmat(row+2,:)
                if (s(row+3)<delta) grad=grad+(delta-s(row+3))*prob%gmat(row+3,:)
                row=row+3
            end do
            ng=vecnorm2(grad)
            if (ng<=1.0e-10_dp) exit
            proj=grad
            if (p>0) then
                rhs=matmul(prob%amat,grad)
                call local_solve(gram,rhs,yy,info)
                if (info==0) proj=grad-matmul(transpose(prob%amat),yy)
            end if
            ng=vecnorm2(proj)
            if (ng<=1.0e-12_dp) exit
            alpha=min(1.0_dp,1.0_dp/(ng+1.0e-12_dp))
            x=x-alpha*proj
        end do
    end subroutine repair_exp_log_arguments

    subroutine newton_direction(prob,jac,hlag,sl,lam,req,rineq,rdual,rcent,reg, &
                                dx,dy,ds,dl,info)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: jac(:,:), hlag(:,:), sl(:), lam(:)
        real(dp), intent(in) :: req(:), rineq(:), rdual(:), rcent(:), reg
        real(dp), intent(out) :: dx(:),dy(:),ds(:),dl(:)
        integer, intent(out) :: info
        real(dp), allocatable :: hbar(:,:), rhsx(:), rhs_eq(:), ratio(:), v(:)
        integer :: n,p,i,ni
        n = size(dx)
        p = size(dy)
        ni = size(sl)
        allocate(hbar(n,n),rhsx(n),rhs_eq(p),ratio(ni),v(ni))
        ratio = lam/sl
        hbar = hlag
        do i = 1, ni
            hbar = hbar + ratio(i)*outer(jac(i,:),jac(i,:))
        end do
        v = (rcent - lam*rineq)/sl
        rhsx = -rdual + matmul(transpose(jac),v)
        if (p > 0) rhs_eq = -req
        call solve_kkt(hbar,prob%amat,rhsx,rhs_eq,dx,dy,reg,info)
        if (info /= 0) then
            ds = 0.0_dp
            dl = 0.0_dp
            return
        end if
        ds = -rineq - matmul(jac,dx)
        dl = (-rcent-lam*ds)/sl
    end subroutine newton_direction

    pure function outer(a,b) result(c)
        real(dp), intent(in) :: a(:),b(:)
        real(dp) :: c(size(a),size(b))
        c = spread(a,2,size(b))*spread(b,1,size(a))
    end function outer

    pure real(dp) function max_step(v,dv) result(a)
        real(dp), intent(in) :: v(:),dv(:)
        integer :: i
        a = 1.0_dp
        do i = 1, size(v)
            if (dv(i) < 0.0_dp) a = min(a,-v(i)/dv(i))
        end do
        a = max(0.0_dp,min(1.0_dp,a))
    end function max_step

    subroutine limit_exp_step(prob,x,dx,alpha)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:),dx(:)
        real(dp), intent(inout) :: alpha
        real(dp), allocatable :: s(:), dsx(:)
        real(dp) :: aa
        integer :: row, k
        if (prob%dims%e <= 0) return
        allocate(s(prob%ncone()),dsx(prob%ncone()))
        s = prob%h - matmul(prob%gmat,x)
        dsx = -matmul(prob%gmat,dx)
        row = prob%dims%l
        if (allocated(prob%dims%q)) row = row + sum(prob%dims%q)
        aa = alpha
        do k = 1, prob%dims%e
            if (dsx(row+2) < 0.0_dp) aa = min(aa,0.99_dp*max(s(row+2),1.0e-14_dp)/(-dsx(row+2)))
            if (dsx(row+3) < 0.0_dp) aa = min(aa,0.99_dp*max(s(row+3),1.0e-14_dp)/(-dsx(row+3)))
            row = row + 3
        end do
        alpha = max(1.0e-8_dp,min(alpha,aa))
    end subroutine limit_exp_step

    subroutine solve_no_inequalities(prob,result,stg,x0)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(inout) :: result
        type(ecos_settings), intent(in) :: stg
        real(dp), intent(in) :: x0(:)
        real(dp), allocatable :: proj(:), rhs(:), gram(:,:), yy(:)
        real(dp) :: nr
        integer :: p,n,info,i
        n = prob%nvar()
        p = prob%neq()
        result%x = x0
        if (p == 0) then
            if (vecnorm2(prob%c) > stg%feastol) then
                result%exitflag = ECOS_DINF
                result%infostring = 'Dual infeasible / primal unbounded'
                result%dinf = 1.0_dp
            else
                result%exitflag = ECOS_OPTIMAL
                result%infostring = 'Optimal solution found'
            end if
            call fill_summary(prob,result)
            return
        end if
        allocate(gram(p,p),yy(p),rhs(p),proj(n))
        gram = matmul(prob%amat,transpose(prob%amat))
        do i = 1,p
            gram(i,i) = gram(i,i) + 1.0e-12_dp*max(1.0_dp,maxval(abs(gram)))
        end do
        rhs = matmul(prob%amat,prob%c)
        call local_solve(gram,rhs,yy,info)
        if (info /= 0) then
            result%exitflag = ECOS_NUMERICS
            result%infostring = 'Numerical problems encountered'
            result%numerr = 1
            return
        end if
        proj = prob%c - matmul(transpose(prob%amat),yy)
        nr = vecnorm2(proj)
        if (nr > stg%feastol*(1.0_dp+vecnorm2(prob%c))) then
            result%exitflag = ECOS_DINF
            result%infostring = 'Dual infeasible / primal unbounded'
            result%dinf = 1.0_dp
        else
            result%y = -yy
            result%exitflag = ECOS_OPTIMAL
            result%infostring = 'Optimal solution found'
        end if
        call fill_summary(prob,result)
    end subroutine solve_no_inequalities

    subroutine local_solve(a,b,x,info)
        use ecos_linalg, only : solve_linear
        real(dp), intent(in) :: a(:,:),b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        call solve_linear(a,b,x,info)
    end subroutine local_solve

    subroutine fill_summary(prob,result)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(inout) :: result
        real(dp), allocatable :: rd(:), re(:)
        real(dp) :: cv
        integer :: n,p
        n = prob%nvar()
        p = prob%neq()
        result%pcost = dot_product(prob%c,result%x)
        result%dcost = -dot_product(prob%h,result%z)
        if (p > 0) result%dcost = result%dcost - dot_product(prob%b,result%y)
        allocate(rd(n))
        rd = prob%c
        if (prob%ncone() > 0) rd = rd + matmul(transpose(prob%gmat),result%z)
        if (p > 0) rd = rd + matmul(transpose(prob%amat),result%y)
        result%dres = vecnorm2(rd)/(1.0_dp+vecnorm2(prob%c))
        result%pres = cone_violation(prob,result%x)/(1.0_dp+merge(vecnorm2(prob%h),0.0_dp,prob%ncone()>0))
        if (p > 0) then
            allocate(re(p))
            re = matmul(prob%amat,result%x)-prob%b
            result%pres = max(result%pres,vecnorm2(re)/(1.0_dp+vecnorm2(prob%b)))
        end if
        result%gap = max(0.0_dp,result%pcost-result%dcost)
        result%relgap = result%gap/(1.0_dp+max(abs(result%pcost),abs(result%dcost)))
        result%pinfres = result%pres
        result%dinfres = result%dres
        result%r0 = 0.0_dp
        cv = cone_violation(prob,result%x)
        if (cv > 1.0e2_dp) result%pinfres = cv
    end subroutine fill_summary

    logical function is_unbounded_direction(prob,stg) result(unb)
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in) :: stg
        type(ecos_problem) :: dprob
        type(ecos_result) :: dr
        type(ecos_settings) :: dstg
        integer :: n,p,m
        unb = .false.
        n = prob%nvar()
        if (vecnorm2(prob%c) <= stg%feastol) return
        p = prob%neq()
        m = prob%ncone()
        allocate(dprob%c(n),dprob%gmat(m,n),dprob%h(m),dprob%amat(p+1,n),dprob%b(p+1))
        dprob%c = 0.0_dp
        if (m > 0) then
            dprob%gmat = prob%gmat
            dprob%h = 0.0_dp
        end if
        if (p > 0) dprob%amat(1:p,:) = prob%amat
        dprob%amat(p+1,:) = prob%c
        if (p > 0) dprob%b(1:p) = 0.0_dp
        dprob%b(p+1) = -1.0_dp
        dprob%dims = prob%dims
        dstg = stg
        dstg%maxit = max(50,min(150,stg%maxit+30))
        call ecos_solve_continuous(dprob,dr,dstg,.false.)
        unb = dr%exitflag == ECOS_OPTIMAL .and. dr%pres < max(1.0e-6_dp,100.0_dp*stg%feastol)
    end function is_unbounded_direction

    logical function phase1_feasible(prob,stg) result(feas)
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in) :: stg
        type(ecos_problem) :: pp
        type(ecos_result) :: rr
        type(ecos_settings) :: pstg
        real(dp), allocatable :: xeq(:), g(:), jac(:,:), hh(:,:), lam(:)
        real(dp) :: r0, eqr
        integer :: n,m,p,ni,info
        ! A conservative diagnosis: equality consistency plus a relaxed cone solve.
        n = prob%nvar()
        m = prob%ncone()
        p = prob%neq()
        ni = prob%dims%scalar_inequalities()
        allocate(xeq(n))
        if (p > 0) then
            call least_norm_equalities(prob%amat,prob%b,xeq,eqr,info)
            if (info /= 0 .or. eqr > 1.0e-6_dp*(1.0_dp+vecnorm2(prob%b))) then
                feas = .false.
                return
            end if
        else
            xeq = 0.0_dp
        end if
        if (ni == 0) then
            feas = .true.
            return
        end if
        allocate(g(ni),jac(ni,n),hh(n,n),lam(ni))
        lam = 0.0_dp
        call cone_scalar_eval(prob,xeq,lam,g,jac,hh)
        r0 = max(1.0_dp,maxval(g)+1.0_dp)
        ! If the scalar violations are modest, the primal-dual solver often has a feasible point nearby.
        if (r0 < 1.0e8_dp) then
            ! Try minimizing a common linear-cone relaxation by augmenting x with r.
            call build_relaxed_problem(prob,r0,pp)
            pstg = stg
            pstg%maxit = max(80,min(250,stg%maxit+80))
            call ecos_solve_continuous(pp,rr,pstg,.false.)
            if (rr%exitflag == ECOS_OPTIMAL) then
                feas = rr%x(n+1) <= max(1.0e-6_dp,100.0_dp*stg%feastol)
                return
            end if
        end if
        feas = .false.
    end function phase1_feasible



    subroutine build_relaxed_problem(prob,r0,pp)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: r0
        type(ecos_problem), intent(out) :: pp
        integer :: n,m,p,l,row_new,row_old,iq,qd,k
        n = prob%nvar()
        m = prob%ncone()
        p = prob%neq()
        l = prob%dims%l
        allocate(pp%c(n+1),pp%gmat(m+1,n+1),pp%h(m+1),pp%amat(p,n+1),pp%b(p))
        pp%c = 0.0_dp
        pp%c(n+1) = 1.0_dp
        pp%gmat = 0.0_dp
        pp%h = 0.0_dp

        ! Preserve ECOS cone row order: linear, SOC, exponential.
        if (l > 0) then
            pp%gmat(1:l,1:n) = prob%gmat(1:l,:)
            pp%gmat(1:l,n+1) = -1.0_dp
            pp%h(1:l) = prob%h(1:l)
        end if
        ! Insert r >= 0 as the final positive-orthant row.
        pp%gmat(l+1,n+1) = -1.0_dp
        pp%h(l+1) = 0.0_dp

        row_old = l
        row_new = l+1
        if (allocated(prob%dims%q)) then
            do iq=1,size(prob%dims%q)
                qd=prob%dims%q(iq)
                if (qd>0) then
                    pp%gmat(row_new+1:row_new+qd,1:n) = prob%gmat(row_old+1:row_old+qd,:)
                    pp%h(row_new+1:row_new+qd) = prob%h(row_old+1:row_old+qd)
                    ! Shift only the SOC head: (t+r,u) remains in Q for sufficiently large r.
                    pp%gmat(row_new+1,n+1) = -1.0_dp
                end if
                row_old=row_old+qd
                row_new=row_new+qd
            end do
        end if
        do k=1,prob%dims%e
            pp%gmat(row_new+1:row_new+3,1:n) = prob%gmat(row_old+1:row_old+3,:)
            pp%h(row_new+1:row_new+3) = prob%h(row_old+1:row_old+3)
            ! Move b and c positively; for sufficiently large r this gives a usable cone relaxation.
            pp%gmat(row_new+2,n+1) = -1.0_dp
            pp%gmat(row_new+3,n+1) = -1.0_dp
            row_old=row_old+3
            row_new=row_new+3
        end do

        if (p > 0) then
            pp%amat(:,1:n) = prob%amat
            pp%amat(:,n+1) = 0.0_dp
            pp%b = prob%b
        end if
        pp%dims%l = l+1
        if (allocated(prob%dims%q)) pp%dims%q = prob%dims%q
        pp%dims%e = prob%dims%e
        if (r0 < 0.0_dp) pp%c(n+1)=1.0_dp
    end subroutine build_relaxed_problem


end module ecos_solver
