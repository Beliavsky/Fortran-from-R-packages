! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_sparse_solver
    use ecos_types, only : dp, ecos_problem, ecos_result, ecos_settings, ecos_sparse_cache, &
        ecos_csc_matrix, ecos_csr_matrix, ECOS_OPTIMAL, ECOS_INACC_OFFSET, ECOS_MAXIT, ECOS_NUMERICS
    use ecos_sparse, only : sparse_triplet_builder, sparse_ldl_factor, triplet_to_csc, &
        csr_matvec, csr_tmatvec, sparse_structure_equal, sparse_pattern_hash
    use ecos_sparse_cones, only : sparse_cone_linearize, sparse_cone_slack, &
        sparse_cone_values, sparse_limit_exp_step
    use ecos_cones, only : cone_dual_from_scalar, cone_violation
    use ecos_linalg, only : vecnorm2
    use ecos_equilibration, only : ecos_scaling, equilibrate_problem_sparse, unscale_result
    implicit none
    private
    public :: solve_pd_sparse, densify_problem

contains

    subroutine densify_problem(prob,dprob)
        type(ecos_problem), intent(in) :: prob
        type(ecos_problem), intent(out) :: dprob
        integer :: n,m,p
        n=prob%nvar(); m=prob%ncone(); p=prob%neq()
        allocate(dprob%c(n),dprob%h(m),dprob%b(p))
        dprob%c=prob%c; dprob%h=prob%h; dprob%b=prob%b; dprob%dims=prob%dims
        if(m>0) then
            call prob%g_csc%to_dense(dprob%gmat)
        else
            allocate(dprob%gmat(0,n))
        end if
        if(p>0) then
            call prob%a_csc%to_dense(dprob%amat)
        else
            allocate(dprob%amat(0,n))
        end if
        dprob%sparse_storage=.false.
        if(allocated(prob%bool_vars)) then
            allocate(dprob%bool_vars(size(prob%bool_vars))); dprob%bool_vars=prob%bool_vars
        end if
        if(allocated(prob%int_vars)) then
            allocate(dprob%int_vars(size(prob%int_vars))); dprob%int_vars=prob%int_vars
        end if
    end subroutine densify_problem

    subroutine allocate_sparse_result(prob,result)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        allocate(result%x(prob%nvar()),result%y(prob%neq()),result%s(prob%ncone()),result%z(prob%ncone()))
        result%x=0.0_dp; result%y=0.0_dp; result%s=0.0_dp; result%z=0.0_dp
        result%sparse_backend_used=.true.
    end subroutine allocate_sparse_result

    subroutine build_eq_kkt(a,n,reg,kkt)
        type(ecos_csr_matrix), intent(in) :: a
        integer, intent(in) :: n
        real(dp), intent(in) :: reg
        type(ecos_csc_matrix), intent(out) :: kkt
        type(sparse_triplet_builder) :: tb
        integer :: p,i,k,j
        p=a%nrow
        call tb%init(n+p,n+p,max(64,n+p+size(a%values)))
        do i=1,n; call tb%add(i,i,1.0_dp); end do
        do i=1,p
            call tb%add(n+i,n+i,-max(reg,1.0e-12_dp))
            do k=a%rowptr(i),a%rowptr(i+1)-1
                j=a%colind(k); call tb%add(j,n+i,a%values(k))
            end do
        end do
        call triplet_to_csc(tb,kkt,.true.)
    end subroutine build_eq_kkt

    subroutine sparse_least_norm_equalities(prob,b,x,resnorm,reg,stg,info)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: b(:),reg
        real(dp), intent(out) :: x(:),resnorm
        type(ecos_settings), intent(in) :: stg
        integer, intent(out) :: info
        type(ecos_csc_matrix) :: kkt
        type(sparse_ldl_factor) :: fac
        real(dp), allocatable :: rhs(:),sol(:),sgn(:),ax(:)
        integer :: n,p,nref
        n=prob%nvar(); p=prob%neq(); info=0
        if(p==0) then; x=0.0_dp; resnorm=0.0_dp; return; end if
        call build_eq_kkt(prob%a_csr,n,reg,kkt)
        allocate(rhs(n+p),sol(n+p),sgn(n+p),ax(p))
        rhs=0.0_dp; rhs(n+1:n+p)=b
        sgn(1:n)=1.0_dp; sgn(n+1:n+p)=-1.0_dp
        call fac%analyze(kkt,stg%sparse_rcm,info,stg%sparse_amd); if(info/=0) return
        call fac%factorize(kkt,sgn,max(reg,1.0e-10_dp),1.0e-14_dp,info); if(info/=0) return
        call fac%solve_refined(kkt,rhs,sol,stg%iterative_refinement,stg%refinement_tol,nref,info)
        if(info/=0) return
        x=sol(1:n)
        call csr_matvec(prob%a_csr,x,ax)
        resnorm=vecnorm2(ax-b)
    end subroutine sparse_least_norm_equalities

    subroutine sparse_project_null(prob,v,proj,reg,stg,info)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: v(:),reg
        real(dp), intent(out) :: proj(:)
        type(ecos_settings), intent(in) :: stg
        integer, intent(out) :: info
        type(ecos_csc_matrix) :: kkt
        type(sparse_ldl_factor) :: fac
        real(dp), allocatable :: rhs(:),sol(:),sgn(:)
        integer :: n,p,nref
        n=prob%nvar(); p=prob%neq(); info=0
        if(p==0) then; proj=v; return; end if
        call build_eq_kkt(prob%a_csr,n,reg,kkt)
        allocate(rhs(n+p),sol(n+p),sgn(n+p))
        rhs=0.0_dp; rhs(1:n)=v
        sgn(1:n)=1.0_dp; sgn(n+1:n+p)=-1.0_dp
        call fac%analyze(kkt,stg%sparse_rcm,info,stg%sparse_amd); if(info/=0) return
        call fac%factorize(kkt,sgn,max(reg,1.0e-10_dp),1.0e-14_dp,info); if(info/=0) return
        call fac%solve_refined(kkt,rhs,sol,stg%iterative_refinement,stg%refinement_tol,nref,info)
        if(info==0) proj=sol(1:n)
    end subroutine sparse_project_null

    subroutine repair_exp_sparse(prob,x,stg,reg)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(inout) :: x(:)
        type(ecos_settings), intent(in) :: stg
        real(dp), intent(in) :: reg
        real(dp), allocatable :: s(:),grad(:),proj(:)
        real(dp) :: ng,alpha,delta
        integer :: row,k,it,j,kk,info
        if(prob%dims%e<=0) return
        allocate(s(prob%ncone()),grad(prob%nvar()),proj(prob%nvar()))
        delta=1.0_dp
        do it=1,30
            call sparse_cone_slack(prob,x,s); grad=0.0_dp
            row=prob%dims%l
            if(allocated(prob%dims%q)) row=row+sum(prob%dims%q)
            do k=1,prob%dims%e
                if(s(row+2)<delta) then
                    do kk=prob%g_csr%rowptr(row+2),prob%g_csr%rowptr(row+3)-1
                        j=prob%g_csr%colind(kk)
                        grad(j)=grad(j)+(delta-s(row+2))*prob%g_csr%values(kk)
                    end do
                end if
                if(s(row+3)<delta) then
                    do kk=prob%g_csr%rowptr(row+3),prob%g_csr%rowptr(row+4)-1
                        j=prob%g_csr%colind(kk)
                        grad(j)=grad(j)+(delta-s(row+3))*prob%g_csr%values(kk)
                    end do
                end if
                row=row+3
            end do
            ng=vecnorm2(grad); if(ng<=1.0e-10_dp) exit
            call sparse_project_null(prob,grad,proj,reg,stg,info)
            if(info/=0) proj=grad
            ng=vecnorm2(proj); if(ng<=1.0e-12_dp) exit
            alpha=min(1.0_dp,1.0_dp/(ng+1.0e-12_dp)); x=x-alpha*proj
        end do
    end subroutine repair_exp_sparse

    subroutine build_sparse_kkt(prob,jac,hlag,sl,lam,reg,kkt)
        type(ecos_problem), intent(in) :: prob
        type(ecos_csr_matrix), intent(in) :: jac
        type(ecos_csc_matrix), intent(in) :: hlag
        real(dp), intent(in) :: sl(:),lam(:),reg
        type(ecos_csc_matrix), intent(out) :: kkt
        type(sparse_triplet_builder) :: tb
        integer :: n,p,ni,ntot,j,k,i,row,idx
        n=prob%nvar(); p=prob%neq(); ni=size(sl); ntot=n+p+ni
        call tb%init(ntot,ntot,max(64,n+p+ni+size(hlag%values)+size(jac%values)+size(prob%a_csr%values)))
        do j=1,n
            call tb%add(j,j,reg)
            do k=hlag%colptr(j),hlag%colptr(j+1)-1
                i=hlag%rowind(k); call tb%add(i,j,hlag%values(k))
            end do
        end do
        do row=1,p
            idx=n+row; call tb%add(idx,idx,-reg)
            do k=prob%a_csr%rowptr(row),prob%a_csr%rowptr(row+1)-1
                call tb%add(prob%a_csr%colind(k),idx,prob%a_csr%values(k))
            end do
        end do
        do row=1,ni
            idx=n+p+row
            call tb%add(idx,idx,-sl(row)/max(lam(row),1.0e-18_dp)-reg)
            do k=jac%rowptr(row),jac%rowptr(row+1)-1
                call tb%add(jac%colind(k),idx,jac%values(k))
            end do
        end do
        call triplet_to_csc(tb,kkt,.true.)
    end subroutine build_sparse_kkt

    subroutine solve_direction_sparse(prob,jac,sl,lam,req,rineq,rdual,rcent,kkt,fac,stg, &
                                      dx,dy,ds,dl,nref,info)
        type(ecos_problem), intent(in) :: prob
        type(ecos_csr_matrix), intent(in) :: jac
        real(dp), intent(in) :: sl(:),lam(:),req(:),rineq(:),rdual(:),rcent(:)
        type(ecos_csc_matrix), intent(in) :: kkt
        type(sparse_ldl_factor), intent(in) :: fac
        type(ecos_settings), intent(in) :: stg
        real(dp), intent(out) :: dx(:),dy(:),ds(:),dl(:)
        integer, intent(out) :: nref,info
        real(dp), allocatable :: rhs(:),sol(:),jdx(:)
        integer :: n,p,ni
        n=prob%nvar(); p=prob%neq(); ni=size(sl)
        allocate(rhs(n+p+ni),sol(n+p+ni),jdx(ni)); rhs=0.0_dp
        rhs(1:n)=-rdual
        if(p>0) rhs(n+1:n+p)=-req
        rhs(n+p+1:n+p+ni)=rcent/max(lam,1.0e-18_dp)-rineq
        call fac%solve_refined(kkt,rhs,sol,stg%iterative_refinement,stg%refinement_tol,nref,info)
        if(info/=0) then; dx=0.0_dp; dy=0.0_dp; ds=0.0_dp; dl=0.0_dp; return; end if
        dx=sol(1:n)
        if(p>0) dy=sol(n+1:n+p)
        dl=sol(n+p+1:n+p+ni)
        call csr_matvec(jac,dx,jdx)
        ds=-rineq-jdx
    end subroutine solve_direction_sparse

    pure real(dp) function max_step(v,dv) result(a)
        real(dp), intent(in) :: v(:),dv(:)
        integer :: i
        a=1.0_dp
        do i=1,size(v)
            if(dv(i)<0.0_dp) a=min(a,-v(i)/dv(i))
        end do
        a=max(0.0_dp,min(1.0_dp,a))
    end function max_step

    subroutine fill_summary_sparse(prob,result)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(inout) :: result
        real(dp), allocatable :: rd(:),tmp(:),re(:)
        real(dp) :: cv
        integer :: n,p
        n=prob%nvar(); p=prob%neq()
        result%pcost=dot_product(prob%c,result%x)
        result%dcost=-dot_product(prob%h,result%z)
        if(p>0) result%dcost=result%dcost-dot_product(prob%b,result%y)
        allocate(rd(n),tmp(n)); rd=prob%c
        if(prob%ncone()>0) then
            call csr_tmatvec(prob%g_csr,result%z,tmp); rd=rd+tmp
        end if
        if(p>0) then
            call csr_tmatvec(prob%a_csr,result%y,tmp); rd=rd+tmp
        end if
        result%dres=vecnorm2(rd)/(1.0_dp+vecnorm2(prob%c))
        result%pres=cone_violation(prob,result%x)/(1.0_dp+merge(vecnorm2(prob%h),0.0_dp,prob%ncone()>0))
        if(p>0) then
            allocate(re(p)); call csr_matvec(prob%a_csr,result%x,re); re=re-prob%b
            result%pres=max(result%pres,vecnorm2(re)/(1.0_dp+vecnorm2(prob%b)))
        end if
        result%gap=max(0.0_dp,result%pcost-result%dcost)
        result%relgap=result%gap/(1.0_dp+max(abs(result%pcost),abs(result%dcost)))
        result%pinfres=result%pres; result%dinfres=result%dres; result%r0=0.0_dp
        cv=cone_violation(prob,result%x); if(cv>1.0e2_dp) result%pinfres=cv
    end subroutine fill_summary_sparse

    subroutine solve_pd_sparse(prob,result,stg,status,cache)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        type(ecos_settings), intent(in) :: stg
        integer, intent(out) :: status
        type(ecos_sparse_cache), intent(inout), optional :: cache
        type(ecos_problem) :: scaled
        type(ecos_result) :: scaled_result
        type(ecos_scaling) :: scale
        integer :: info

        if (stg%equilibrate) then
            call equilibrate_problem_sparse(prob,stg,scaled,scale,info)
            if (info /= 0) then
                status = ECOS_NUMERICS
                call allocate_sparse_result(prob,result)
                result%exitflag = status
                result%infostring = 'Sparse equilibration failed'
                return
            end if
            if (present(cache)) then
                call solve_pd_sparse_core(scaled,scaled_result,stg,status,cache)
            else
                call solve_pd_sparse_core(scaled,scaled_result,stg,status)
            end if
            result = scaled_result
            call unscale_result(scale,result)
            result%min_col_scale = minval(scale%x)
            result%max_col_scale = maxval(scale%x)
            if (size(scale%eq) > 0 .or. size(scale%cone) > 0) then
                result%min_row_scale = huge(1.0_dp)
                result%max_row_scale = 0.0_dp
                if (size(scale%eq) > 0) then
                    result%min_row_scale = min(result%min_row_scale,minval(scale%eq))
                    result%max_row_scale = max(result%max_row_scale,maxval(scale%eq))
                end if
                if (size(scale%cone) > 0) then
                    result%min_row_scale = min(result%min_row_scale,minval(scale%cone))
                    result%max_row_scale = max(result%max_row_scale,maxval(scale%cone))
                end if
            end if
            call fill_summary_sparse(prob,result)
        else
            if (present(cache)) then
                call solve_pd_sparse_core(prob,result,stg,status,cache)
            else
                call solve_pd_sparse_core(prob,result,stg,status)
            end if
        end if
    end subroutine solve_pd_sparse

    subroutine solve_pd_sparse_core(prob,result,stg,status,cache)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        type(ecos_settings), intent(in) :: stg
        integer, intent(out) :: status
        type(ecos_sparse_cache), intent(inout), optional :: cache
        type(ecos_csr_matrix) :: jac
        type(ecos_csc_matrix) :: hlag,kkt,kkt_prev
        type(sparse_ldl_factor) :: fac
        real(dp), allocatable :: x(:),y(:),sl(:),lam(:),g(:),req(:),rineq(:),rdual(:),rcent(:)
        real(dp), allocatable :: dx(:),dy(:),ds(:),dl(:),dxa(:),dya(:),dsa(:),dla(:),tmpn(:),sgn(:)
        real(dp), allocatable :: bestx(:),besty(:),bestsl(:),bestlam(:)
        real(dp) :: eqinit,pres,dres,mu,mua,sigma,ap,ad,apa,ada,scale_p,scale_d,reg,gap_est
        real(dp) :: merit,best_merit,best_pres,best_dres,best_gap,t0,t1,t2,t3
        integer :: n,p,ni,it,info,i,nref,nref2,total_ref,pattern_hash,reg_try
        logical :: converged,need_analyze,warm_used

        n=prob%nvar(); p=prob%neq(); ni=prob%dims%scalar_inequalities()
        call allocate_sparse_result(prob,result)
        if(ni==0) then
            status=ECOS_NUMERICS; result%exitflag=status
            result%infostring='Sparse backend requires at least one cone inequality'; return
        end if
        allocate(x(n),y(p),sl(ni),lam(ni),g(ni),req(p),rineq(ni),rdual(n),rcent(ni))
        allocate(dx(n),dy(p),ds(ni),dl(ni),dxa(n),dya(p),dsa(ni),dla(ni),tmpn(n))
        allocate(bestx(n),besty(p),bestsl(ni),bestlam(ni))
        reg=max(stg%regularization,1.0e-12_dp); y=0.0_dp
        warm_used=.false.
        if (present(cache)) then
            if (cache%warm_valid .and. allocated(cache%warm_x) .and. allocated(cache%warm_y) .and. &
                allocated(cache%warm_sl) .and. allocated(cache%warm_lam)) then
                if (size(cache%warm_x)==n .and. size(cache%warm_y)==p .and. &
                    size(cache%warm_sl)==ni .and. size(cache%warm_lam)==ni) then
                    x=cache%warm_x; y=cache%warm_y
                    sl=max(cache%warm_sl,1.0e-12_dp); lam=max(cache%warm_lam,1.0e-12_dp)
                    warm_used=.true.; result%cached_warm_starts=1
                end if
            end if
        end if
        if (.not.warm_used) then
            if(p>0) then
                call sparse_least_norm_equalities(prob,prob%b,x,eqinit,reg,stg,info)
                if(info/=0) x=0.0_dp
            else
                x=0.0_dp; eqinit=0.0_dp
            end if
            if(prob%dims%e>0) call repair_exp_sparse(prob,x,stg,reg)
            lam=1.0_dp
            call sparse_cone_linearize(prob,x,lam,g,jac,hlag)
            sl=max(1.0_dp,-g+1.0_dp)
        end if
        status=ECOS_MAXIT; converged=.false.; total_ref=0; need_analyze=.true.
        best_merit=huge(1.0_dp); best_pres=huge(1.0_dp)
        best_dres=huge(1.0_dp); best_gap=huge(1.0_dp)
        bestx=x; besty=y; bestsl=sl; bestlam=lam

        do it=0,stg%maxit
            call sparse_cone_linearize(prob,x,lam,g,jac,hlag)
            if(p>0) then
                call csr_matvec(prob%a_csr,x,req); req=req-prob%b
            end if
            rineq=g+sl
            call csr_tmatvec(jac,lam,rdual); rdual=prob%c+rdual
            if(p>0) then
                call csr_tmatvec(prob%a_csr,y,tmpn); rdual=rdual+tmpn
            end if
            mu=dot_product(sl,lam)/real(ni,dp)
            scale_p=1.0_dp+max(vecnorm2(prob%h),merge(vecnorm2(prob%b),0.0_dp,p>0))
            scale_d=1.0_dp+vecnorm2(prob%c)
            pres=vecnorm2(rineq)/scale_p
            if(p>0) pres=max(pres,vecnorm2(req)/scale_p)
            dres=vecnorm2(rdual)/scale_d; gap_est=real(ni,dp)*mu
            merit=pres+dres+gap_est/(1.0_dp+abs(dot_product(prob%c,x)))
            if (merit < best_merit .and. all(abs(x) < huge(1.0_dp))) then
                best_merit=merit; best_pres=pres; best_dres=dres; best_gap=gap_est
                bestx=x; besty=y; bestsl=sl; bestlam=lam
            end if
            if(pres<=stg%feastol .and. dres<=stg%feastol .and. &
               (gap_est<=stg%abstol .or. gap_est<=stg%reltol*(1.0_dp+abs(dot_product(prob%c,x))))) then
                converged=.true.; status=ECOS_OPTIMAL; exit
            end if
            if(it>=stg%maxit) exit

            call build_sparse_kkt(prob,jac,hlag,sl,lam,reg,kkt)
            pattern_hash=sparse_pattern_hash(kkt)
            if(need_analyze .or. .not.fac%analyzed) then
                info=1
                if (present(cache)) then
                    if (cache%symbolic_valid .and. cache%n==kkt%ncol .and. &
                        cache%structure_hash==pattern_hash) then
                        call fac%load_symbolic(cache,info)
                        if(info==0) then
                            result%cached_symbolic_reuses=result%cached_symbolic_reuses+1
                        end if
                    end if
                end if
                if(info/=0) then
                    call cpu_time(t0)
                    call fac%analyze(kkt,stg%sparse_rcm,info,stg%sparse_amd)
                    call cpu_time(t1); result%time_ordering=result%time_ordering+(t1-t0)
                    result%symbolic_analyses=result%symbolic_analyses+1
                    if(info/=0) then; status=ECOS_NUMERICS; exit; end if
                    if (present(cache)) call fac%save_symbolic(cache,pattern_hash)
                end if
                kkt_prev=kkt; need_analyze=.false.
            else if(.not.sparse_structure_equal(kkt,kkt_prev)) then
                call cpu_time(t0)
                call fac%analyze(kkt,stg%sparse_rcm,info,stg%sparse_amd)
                call cpu_time(t1); result%time_ordering=result%time_ordering+(t1-t0)
                result%symbolic_analyses=result%symbolic_analyses+1
                if(info/=0) then; status=ECOS_NUMERICS; exit; end if
                if (present(cache)) call fac%save_symbolic(cache,pattern_hash)
                kkt_prev=kkt
            end if
            allocate(sgn(n+p+ni)); sgn(1:n)=1.0_dp
            if(p+ni>0) sgn(n+1:n+p+ni)=-1.0_dp
            info=1
            do reg_try=0,max(0,stg%max_regularization_updates)
                call cpu_time(t0)
                call fac%factorize(kkt,sgn,max(reg,1.0e-12_dp),1.0e-14_dp,info)
                call cpu_time(t1); result%time_factorization=result%time_factorization+(t1-t0)
                result%numeric_factorizations=result%numeric_factorizations+1
                if(info==0) exit
                if(.not.stg%dynamic_regularization) exit
                reg=min(1.0e-2_dp,max(10.0_dp*reg,1.0e-10_dp))
                result%regularization_updates=result%regularization_updates+1
                call build_sparse_kkt(prob,jac,hlag,sl,lam,reg,kkt)
            end do
            deallocate(sgn)
            if(info/=0) then; status=ECOS_NUMERICS; exit; end if
            result%kkt_nnz=max(result%kkt_nnz,size(kkt%values))
            result%ldl_nnz=max(result%ldl_nnz,fac%symbolic_nnz)
            if(result%kkt_nnz>0) result%factor_fill_ratio= &
                real(result%ldl_nnz,dp)/real(result%kkt_nnz,dp)

            rcent=sl*lam
            call cpu_time(t2)
            call solve_direction_sparse(prob,jac,sl,lam,req,rineq,rdual,rcent,kkt,fac,stg, &
                                        dxa,dya,dsa,dla,nref,info)
            call cpu_time(t3); result%time_refinement=result%time_refinement+(t3-t2)
            total_ref=total_ref+nref
            if(info/=0) then; status=ECOS_NUMERICS; exit; end if
            apa=max_step(sl,dsa); ada=max_step(lam,dla)
            mua=dot_product(sl+apa*dsa,lam+ada*dla)/real(ni,dp)
            if(mu>tiny(1.0_dp)) then
                sigma=min(1.0_dp,max(0.0_dp,(mua/mu)**3))
            else
                sigma=0.0_dp
            end if
            rcent=sl*lam+dsa*dla-sigma*mu
            call cpu_time(t2)
            call solve_direction_sparse(prob,jac,sl,lam,req,rineq,rdual,rcent,kkt,fac,stg, &
                                        dx,dy,ds,dl,nref2,info)
            call cpu_time(t3); result%time_refinement=result%time_refinement+(t3-t2)
            total_ref=total_ref+nref2
            if(info/=0) then; status=ECOS_NUMERICS; exit; end if
            ap=min(1.0_dp,0.995_dp*max_step(sl,ds)); ad=min(1.0_dp,0.995_dp*max_step(lam,dl))
            call sparse_limit_exp_step(prob,x,dx,ap)
            x=x+ap*dx; sl=sl+ap*ds; y=y+ad*dy; lam=lam+ad*dl
            do i=1,ni
                sl(i)=max(sl(i),1.0e-14_dp); lam(i)=max(lam(i),1.0e-14_dp)
            end do
        end do

        if (.not.converged .and. best_merit < huge(1.0_dp)) then
            x=bestx; y=besty; sl=bestsl; lam=bestlam
            if (best_pres<=stg%feastol_inacc .and. best_dres<=stg%feastol_inacc .and. &
                (best_gap<=stg%abstol_inacc .or. &
                 best_gap<=stg%reltol_inacc*(1.0_dp+abs(dot_product(prob%c,x))))) then
                status=ECOS_OPTIMAL+ECOS_INACC_OFFSET
            end if
        end if
        if (present(cache)) then
            if (allocated(cache%warm_x)) deallocate(cache%warm_x)
            if (allocated(cache%warm_y)) deallocate(cache%warm_y)
            if (allocated(cache%warm_sl)) deallocate(cache%warm_sl)
            if (allocated(cache%warm_lam)) deallocate(cache%warm_lam)
            allocate(cache%warm_x(n),cache%warm_y(p),cache%warm_sl(ni),cache%warm_lam(ni))
            cache%warm_x=x; cache%warm_y=y; cache%warm_sl=sl; cache%warm_lam=lam
            cache%warm_valid=.true.
        end if
        result%x=x; result%y=y
        call sparse_cone_slack(prob,x,result%s)
        call cone_dual_from_scalar(prob,x,lam,result%z)
        result%iter=min(it,stg%maxit); result%exitflag=status; result%iterative_refinements=total_ref
        if(status==ECOS_OPTIMAL) then
            result%infostring='Optimal solution found (sparse LDL backend)'
        else if(status==ECOS_OPTIMAL+ECOS_INACC_OFFSET) then
            result%infostring='Inaccurate optimal solution found (sparse LDL backend)'
        else if(status==ECOS_MAXIT) then
            result%infostring='Maximum number of iterations reached (sparse LDL backend)'
        else
            result%infostring='Numerical problems encountered in sparse LDL backend'; result%numerr=1
        end if
        call fill_summary_sparse(prob,result)
        if(.not.converged .and. status==ECOS_OPTIMAL) status=ECOS_MAXIT
    end subroutine solve_pd_sparse_core

end module ecos_sparse_solver
