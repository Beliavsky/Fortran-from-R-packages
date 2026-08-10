! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_bb
    use ecos_types, only : dp, ecos_problem, ecos_result, ecos_settings, ecos_sparse_cache, &
        ECOS_OPTIMAL, ECOS_PINF, ECOS_DINF, ECOS_MAXIT
    use ecos_solver, only : ecos_solve_continuous
    use ecos_sparse, only : sparse_triplet_builder, triplet_to_csc, csc_to_csr
    implicit none
    private
    public :: ecos_solve_mixed_integer

contains

    subroutine ecos_solve_mixed_integer(prob,result,settings)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        type(ecos_settings), intent(in) :: settings
        real(dp), allocatable :: lb(:),ub(:),bestx(:)
        integer, allocatable :: ivars(:)
        real(dp) :: bestobj
        integer :: n,i,nodes,bb_reuses
        logical :: have_best
        type(ecos_sparse_cache) :: cache

        n = prob%nvar()
        allocate(lb(n),ub(n),bestx(n))
        lb = -huge(1.0_dp)
        ub = huge(1.0_dp)
        call collect_integer_vars(prob,ivars)
        if (allocated(prob%bool_vars)) then
            do i=1,size(prob%bool_vars)
                if (prob%bool_vars(i)>=1 .and. prob%bool_vars(i)<=n) then
                    lb(prob%bool_vars(i))=0.0_dp
                    ub(prob%bool_vars(i))=1.0_dp
                end if
            end do
        end if
        bestobj = huge(1.0_dp)
        bestx = 0.0_dp
        nodes = 0; bb_reuses=0
        have_best = .false.
        call branch(prob,settings,ivars,lb,ub,bestobj,bestx,have_best,nodes,cache,bb_reuses)

        call initialize_result_shape(prob,result)
        result%mi_iter = nodes
        result%bb_symbolic_reuses = bb_reuses
        result%sparse_backend_used = prob%sparse_storage
        if (have_best) then
            result%x = bestx
            result%pcost = bestobj
            result%exitflag = ECOS_OPTIMAL
            result%infostring = 'Optimal branch and bound solution found'
            ! Re-solve the fixed continuous problem to get useful dual/slack output.
            call solve_fixed(prob,settings,bestx,ivars,result)
            result%x = bestx
            result%pcost = bestobj
            result%mi_iter = nodes
            result%bb_symbolic_reuses = bb_reuses
            result%exitflag = ECOS_OPTIMAL
            result%infostring = 'Optimal branch and bound solution found'
        else if (nodes >= settings%mi_max_iters) then
            result%exitflag = ECOS_MAXIT
            result%infostring = 'Maximum mixed-integer iterations reached with no feasible solution'
        else
            result%exitflag = ECOS_PINF
            result%infostring = 'Mixed-integer problem is infeasible'
        end if
    end subroutine ecos_solve_mixed_integer

    recursive subroutine branch(prob,stg,ivars,lb,ub,bestobj,bestx,have_best,nodes,cache,bb_reuses)
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in) :: stg
        integer, intent(in) :: ivars(:)
        real(dp), intent(in) :: lb(:),ub(:)
        real(dp), intent(inout) :: bestobj,bestx(:)
        logical, intent(inout) :: have_best
        integer, intent(inout) :: nodes,bb_reuses
        type(ecos_sparse_cache), intent(inout) :: cache
        type(ecos_problem) :: nodeprob
        type(ecos_result) :: rr
        real(dp), allocatable :: lb2(:),ub2(:)
        real(dp) :: obj,frac,bfrac,v,flo,cei,cut
        integer :: k,bk,j

        if (nodes >= stg%mi_max_iters) return
        if (any(lb > ub + stg%mi_int_tol)) return
        nodes = nodes + 1
        if (prob%sparse_storage) then
            call add_box_constraints_sparse(prob,ivars,lb,ub,nodeprob)
            call ecos_solve_continuous(nodeprob,rr,stg,cache=cache)
            bb_reuses=bb_reuses+rr%cached_symbolic_reuses
        else
            call add_box_constraints(prob,lb,ub,nodeprob)
            call ecos_solve_continuous(nodeprob,rr,stg)
        end if
        if (rr%exitflag == ECOS_PINF) return
        if (rr%exitflag == ECOS_DINF) return
        if (rr%exitflag /= ECOS_OPTIMAL .and. rr%exitflag /= ECOS_MAXIT) return
        obj = dot_product(prob%c,rr%x)
        if (have_best) then
            cut = bestobj - max(stg%mi_abs_eps,stg%mi_rel_eps*abs(bestobj))
            if (obj >= cut) return
        end if
        bk = 0
        bfrac = 0.0_dp
        do k=1,size(ivars)
            j=ivars(k)
            v=rr%x(j)
            frac=abs(v-anint(v))
            if (frac > stg%mi_int_tol .and. frac > bfrac) then
                bfrac=frac
                bk=j
            end if
        end do
        if (bk == 0) then
            if (.not.have_best .or. obj < bestobj) then
                bestobj=obj
                bestx=rr%x
                have_best=.true.
            end if
            return
        end if
        v=rr%x(bk)
        flo=floor(v)
        cei=ceiling(v)
        allocate(lb2(size(lb)),ub2(size(ub)))
        ! Explore the closer branch first.
        if (v-flo <= cei-v) then
            lb2=lb; ub2=ub; ub2(bk)=min(ub2(bk),flo)
            call branch(prob,stg,ivars,lb2,ub2,bestobj,bestx,have_best,nodes,cache,bb_reuses)
            if (nodes < stg%mi_max_iters) then
                lb2=lb; ub2=ub; lb2(bk)=max(lb2(bk),cei)
                call branch(prob,stg,ivars,lb2,ub2,bestobj,bestx,have_best,nodes,cache,bb_reuses)
            end if
        else
            lb2=lb; ub2=ub; lb2(bk)=max(lb2(bk),cei)
            call branch(prob,stg,ivars,lb2,ub2,bestobj,bestx,have_best,nodes,cache,bb_reuses)
            if (nodes < stg%mi_max_iters) then
                lb2=lb; ub2=ub; ub2(bk)=min(ub2(bk),flo)
                call branch(prob,stg,ivars,lb2,ub2,bestobj,bestx,have_best,nodes,cache,bb_reuses)
            end if
        end if
    end subroutine branch

    subroutine collect_integer_vars(prob,ivars)
        type(ecos_problem), intent(in) :: prob
        integer, allocatable, intent(out) :: ivars(:)
        integer, allocatable :: tmp(:)
        integer :: n,k,i
        n=0
        if (allocated(prob%bool_vars)) n=n+size(prob%bool_vars)
        if (allocated(prob%int_vars)) n=n+size(prob%int_vars)
        allocate(tmp(max(1,n)))
        k=0
        if (allocated(prob%bool_vars)) then
            do i=1,size(prob%bool_vars)
                if (.not.any(tmp(1:k)==prob%bool_vars(i))) then
                    k=k+1; tmp(k)=prob%bool_vars(i)
                end if
            end do
        end if
        if (allocated(prob%int_vars)) then
            do i=1,size(prob%int_vars)
                if (.not.any(tmp(1:k)==prob%int_vars(i))) then
                    k=k+1; tmp(k)=prob%int_vars(i)
                end if
            end do
        end if
        allocate(ivars(k))
        if (k>0) ivars=tmp(1:k)
    end subroutine collect_integer_vars

    subroutine add_box_constraints_sparse(prob,ivars,lb,ub,node)
        type(ecos_problem), intent(in) :: prob
        integer, intent(in) :: ivars(:)
        real(dp), intent(in) :: lb(:),ub(:)
        type(ecos_problem), intent(out) :: node
        type(sparse_triplet_builder) :: tg
        integer :: n,m,p,nadd,j,k,r,newr,ii,row
        logical :: haslb,hasub
        n=prob%nvar(); m=prob%ncone(); p=prob%neq(); nadd=2*size(ivars)
        allocate(node%c(n),node%h(m+nadd),node%b(p))
        node%c=prob%c; node%b=prob%b; node%sparse_storage=.true.
        node%a_csc=prob%a_csc; node%a_csr=prob%a_csr
        node%dims=prob%dims; node%dims%l=prob%dims%l+nadd
        call tg%init(m+nadd,n,max(16,size(prob%g_csc%values)+nadd))
        node%h=0.0_dp
        do j=1,prob%g_csc%ncol
            do k=prob%g_csc%colptr(j),prob%g_csc%colptr(j+1)-1
                r=prob%g_csc%rowind(k)
                if(r<=prob%dims%l) then
                    newr=r
                else
                    newr=r+nadd
                end if
                call tg%add(newr,j,prob%g_csc%values(k))
            end do
        end do
        if(prob%dims%l>0) node%h(1:prob%dims%l)=prob%h(1:prob%dims%l)
        row=prob%dims%l
        do ii=1,size(ivars)
            j=ivars(ii)
            haslb=lb(j)>-0.5_dp*huge(1.0_dp)
            hasub=ub(j)< 0.5_dp*huge(1.0_dp)
            row=row+1
            if(haslb) then
                call tg%add(row,j,-1.0_dp); node%h(row)=-lb(j)
            else
                call tg%add(row,j,0.0_dp); node%h(row)=1.0_dp
            end if
            row=row+1
            if(hasub) then
                call tg%add(row,j,1.0_dp); node%h(row)=ub(j)
            else
                call tg%add(row,j,0.0_dp); node%h(row)=1.0_dp
            end if
        end do
        if(m>prob%dims%l) node%h(row+1:m+nadd)=prob%h(prob%dims%l+1:m)
        call triplet_to_csc(tg,node%g_csc); call csc_to_csr(node%g_csc,node%g_csr)
    end subroutine add_box_constraints_sparse

    subroutine add_box_constraints(prob,lb,ub,node)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: lb(:),ub(:)
        type(ecos_problem), intent(out) :: node
        integer :: n,m,p,nadd,i,row,rest
        logical :: haslb,hasub
        n=prob%nvar(); m=prob%ncone(); p=prob%neq()
        nadd=0
        do i=1,n
            if (lb(i)>-0.5_dp*huge(1.0_dp)) nadd=nadd+1
            if (ub(i)< 0.5_dp*huge(1.0_dp)) nadd=nadd+1
        end do
        allocate(node%c(n),node%gmat(m+nadd,n),node%h(m+nadd),node%amat(p,n),node%b(p))
        node%c=prob%c
        if (p>0) then
            node%amat=prob%amat; node%b=prob%b
        end if
        row=0
        if (prob%dims%l>0) then
            node%gmat(1:prob%dims%l,:)=prob%gmat(1:prob%dims%l,:)
            node%h(1:prob%dims%l)=prob%h(1:prob%dims%l)
            row=prob%dims%l
        end if
        do i=1,n
            haslb=lb(i)>-0.5_dp*huge(1.0_dp)
            hasub=ub(i)< 0.5_dp*huge(1.0_dp)
            if (haslb) then
                row=row+1; node%gmat(row,:)=0.0_dp; node%gmat(row,i)=-1.0_dp; node%h(row)=-lb(i)
            end if
            if (hasub) then
                row=row+1; node%gmat(row,:)=0.0_dp; node%gmat(row,i)=1.0_dp; node%h(row)=ub(i)
            end if
        end do
        rest=m-prob%dims%l
        if (rest>0) then
            node%gmat(row+1:row+rest,:)=prob%gmat(prob%dims%l+1:m,:)
            node%h(row+1:row+rest)=prob%h(prob%dims%l+1:m)
        end if
        node%dims%l=prob%dims%l+nadd
        if (allocated(prob%dims%q)) node%dims%q=prob%dims%q
        node%dims%e=prob%dims%e
    end subroutine add_box_constraints

    subroutine initialize_result_shape(prob,result)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        allocate(result%x(prob%nvar()),result%y(prob%neq()),result%s(prob%ncone()),result%z(prob%ncone()))
        result%x=0.0_dp; result%y=0.0_dp; result%s=0.0_dp; result%z=0.0_dp
    end subroutine initialize_result_shape

    subroutine solve_fixed(prob,stg,xfix,ivars,result)
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in) :: stg
        real(dp), intent(in) :: xfix(:)
        integer, intent(in) :: ivars(:)
        type(ecos_result), intent(inout) :: result
        type(ecos_problem) :: fp
        type(ecos_result) :: rr
        type(ecos_sparse_cache) :: cache
        real(dp), allocatable :: lb(:),ub(:)
        integer :: n,m,p,k,j,l,nadd
        n=prob%nvar(); m=prob%ncone(); p=prob%neq()
        if(prob%sparse_storage) then
            allocate(lb(n),ub(n)); lb=-huge(1.0_dp); ub=huge(1.0_dp)
            do k=1,size(ivars)
                j=ivars(k); lb(j)=xfix(j); ub(j)=xfix(j)
            end do
            call add_box_constraints_sparse(prob,ivars,lb,ub,fp)
            call ecos_solve_continuous(fp,rr,stg,.false.,cache)
            if(rr%exitflag==ECOS_OPTIMAL) then
                result%y=rr%y
                l=prob%dims%l; nadd=2*size(ivars)
                if(l>0) then
                    result%s(1:l)=rr%s(1:l); result%z(1:l)=rr%z(1:l)
                end if
                if(m>l) then
                    result%s(l+1:m)=rr%s(l+nadd+1:m+nadd)
                    result%z(l+1:m)=rr%z(l+nadd+1:m+nadd)
                end if
                result%dcost=rr%dcost; result%pres=rr%pres; result%dres=rr%dres
                result%gap=rr%gap; result%relgap=rr%relgap
            end if
            return
        end if
        allocate(fp%c(n),fp%gmat(m,n),fp%h(m),fp%amat(p+size(ivars),n),fp%b(p+size(ivars)))
        fp%c=prob%c
        if (m>0) then; fp%gmat=prob%gmat; fp%h=prob%h; end if
        if (p>0) then; fp%amat(1:p,:)=prob%amat; fp%b(1:p)=prob%b; end if
        do k=1,size(ivars)
            j=ivars(k); fp%amat(p+k,:)=0.0_dp; fp%amat(p+k,j)=1.0_dp; fp%b(p+k)=xfix(j)
        end do
        fp%dims=prob%dims
        call ecos_solve_continuous(fp,rr,stg,.false.)
        if (rr%exitflag==ECOS_OPTIMAL) then
            result%y=rr%y(1:p)
            result%s=rr%s; result%z=rr%z
            result%dcost=rr%dcost; result%pres=rr%pres; result%dres=rr%dres
            result%gap=rr%gap; result%relgap=rr%relgap
        end if
    end subroutine solve_fixed

end module ecos_bb
