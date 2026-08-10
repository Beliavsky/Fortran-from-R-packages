! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_api
    use ecos_types
    use ecos_solver, only : ecos_solve_continuous
    use ecos_bb, only : ecos_solve_mixed_integer
    use ecos_sparse, only : csc_to_csr, sparse_structure_equal
    implicit none
    private
    public :: ecos_csolve, ecos_setup, ecos_solve, ecos_update, ecos_cleanup
    public :: ecos_control, setup_problem_dense, setup_problem_csc
    public :: ecos_problem, ecos_result, ecos_settings, ecos_dims, ecos_workspace
    public :: ecos_csc_matrix, ecos_csr_matrix
    public :: dp, make_csc_matrix, csc_from_zero_based
    public :: ECOS_OPTIMAL, ECOS_PINF, ECOS_DINF, ECOS_INACC_OFFSET, ECOS_MAXIT
    public :: ECOS_NUMERICS, ECOS_OUTCONE, ECOS_SIGINT, ECOS_FATAL

    interface ecos_csolve
        module procedure ecos_csolve_problem
        module procedure ecos_csolve_dense
        module procedure ecos_csolve_csc
    end interface ecos_csolve

contains

    logical function has_integer_variables(prob) result(has_int)
        type(ecos_problem), intent(in) :: prob
        has_int = .false.
        if (allocated(prob%bool_vars)) then
            if (size(prob%bool_vars) > 0) has_int = .true.
        end if
        if (allocated(prob%int_vars)) then
            if (size(prob%int_vars) > 0) has_int = .true.
        end if
    end function has_integer_variables

    function ecos_control(maxit,feastol,reltol,abstol,feastol_inacc,abstol_inacc, &
                          reltol_inacc,verbose,mi_max_iters,mi_int_tol,mi_abs_eps,mi_rel_eps, &
                          sparse_kkt,sparse_rcm,iterative_refinement,refinement_tol,sparse_amd, &
                          equilibrate,equilibration_iters,dynamic_regularization, &
                          max_regularization_updates,certificate_maxit) result(ctrl)
        integer, intent(in), optional :: maxit,mi_max_iters,iterative_refinement
        integer, intent(in), optional :: equilibration_iters,max_regularization_updates,certificate_maxit
        real(dp), intent(in), optional :: feastol,reltol,abstol,feastol_inacc,abstol_inacc,reltol_inacc
        real(dp), intent(in), optional :: mi_int_tol,mi_abs_eps,mi_rel_eps
        logical, intent(in), optional :: verbose,sparse_kkt,sparse_rcm,sparse_amd,equilibrate
        logical, intent(in), optional :: dynamic_regularization
        real(dp), intent(in), optional :: refinement_tol
        type(ecos_settings) :: ctrl
        ctrl=ecos_settings()
        if (present(maxit)) ctrl%maxit=maxit
        if (present(feastol)) ctrl%feastol=feastol
        if (present(reltol)) ctrl%reltol=reltol
        if (present(abstol)) ctrl%abstol=abstol
        if (present(feastol_inacc)) ctrl%feastol_inacc=feastol_inacc
        if (present(abstol_inacc)) ctrl%abstol_inacc=abstol_inacc
        if (present(reltol_inacc)) ctrl%reltol_inacc=reltol_inacc
        if (present(verbose)) ctrl%verbose=verbose
        if (present(mi_max_iters)) ctrl%mi_max_iters=mi_max_iters
        if (present(mi_int_tol)) ctrl%mi_int_tol=mi_int_tol
        if (present(mi_abs_eps)) ctrl%mi_abs_eps=mi_abs_eps
        if (present(mi_rel_eps)) ctrl%mi_rel_eps=mi_rel_eps
        if (present(sparse_kkt)) ctrl%sparse_kkt=sparse_kkt
        if (present(sparse_rcm)) ctrl%sparse_rcm=sparse_rcm
        if (present(sparse_amd)) ctrl%sparse_amd=sparse_amd
        if (present(equilibrate)) ctrl%equilibrate=equilibrate
        if (present(equilibration_iters)) ctrl%equilibration_iters=max(0,equilibration_iters)
        if (present(dynamic_regularization)) ctrl%dynamic_regularization=dynamic_regularization
        if (present(max_regularization_updates)) then
            ctrl%max_regularization_updates=max(0,max_regularization_updates)
        end if
        if (present(certificate_maxit)) ctrl%certificate_maxit=max(10,certificate_maxit)
        if (present(iterative_refinement)) ctrl%iterative_refinement=max(0,iterative_refinement)
        if (present(refinement_tol)) ctrl%refinement_tol=max(0.0_dp,refinement_tol)
    end function ecos_control

    subroutine setup_problem_dense(prob,c,g,h,dims,a,b,bool_vars,int_vars,ierr)
        type(ecos_problem), intent(out) :: prob
        real(dp), intent(in) :: c(:),g(:,:),h(:)
        type(ecos_dims), intent(in) :: dims
        real(dp), intent(in), optional :: a(:,:),b(:)
        integer, intent(in), optional :: bool_vars(:),int_vars(:)
        integer, intent(out), optional :: ierr
        integer :: n,p
        if (present(ierr)) ierr=0
        n=size(c)
        prob%sparse_storage=.false.
        allocate(prob%c(n)); prob%c=c
        allocate(prob%gmat(size(g,1),size(g,2))); prob%gmat=g
        allocate(prob%h(size(h))); prob%h=h
        prob%dims=dims
        if (present(a)) then
            p=size(a,1); allocate(prob%amat(p,n)); prob%amat=a
            allocate(prob%b(p))
            if (present(b)) then; prob%b=b; else; prob%b=0.0_dp; end if
        else
            allocate(prob%amat(0,n),prob%b(0))
        end if
        if (present(bool_vars)) then; allocate(prob%bool_vars(size(bool_vars))); prob%bool_vars=bool_vars; end if
        if (present(int_vars)) then; allocate(prob%int_vars(size(int_vars))); prob%int_vars=int_vars; end if
        if (.not.prob%valid()) then
            if (present(ierr)) ierr=1
        end if
    end subroutine setup_problem_dense

    subroutine setup_problem_csc(prob,c,g,h,dims,a,b,bool_vars,int_vars,ierr)
        type(ecos_problem), intent(out) :: prob
        real(dp), intent(in) :: c(:),h(:)
        type(ecos_csc_matrix), intent(in) :: g
        type(ecos_dims), intent(in) :: dims
        type(ecos_csc_matrix), intent(in), optional :: a
        real(dp), intent(in), optional :: b(:)
        integer, intent(in), optional :: bool_vars(:),int_vars(:)
        integer, intent(out), optional :: ierr
        integer :: n,p
        if(present(ierr)) ierr=0
        n=size(c); prob%sparse_storage=.true.
        allocate(prob%c(n),prob%h(size(h))); prob%c=c; prob%h=h; prob%dims=dims
        prob%g_csc=g
        call csc_to_csr(prob%g_csc,prob%g_csr)
        if(present(a)) then
            p=a%nrow; prob%a_csc=a; call csc_to_csr(prob%a_csc,prob%a_csr)
            allocate(prob%b(p))
            if(present(b)) then; prob%b=b; else; prob%b=0.0_dp; end if
        else
            p=0; allocate(prob%b(0))
            prob%a_csc%nrow=0; prob%a_csc%ncol=n
            allocate(prob%a_csc%colptr(n+1),prob%a_csc%rowind(0),prob%a_csc%values(0))
            prob%a_csc%colptr=1
            call csc_to_csr(prob%a_csc,prob%a_csr)
        end if
        if(present(bool_vars)) then; allocate(prob%bool_vars(size(bool_vars))); prob%bool_vars=bool_vars; end if
        if(present(int_vars)) then; allocate(prob%int_vars(size(int_vars))); prob%int_vars=int_vars; end if
        if(.not.prob%valid()) then
            if(present(ierr)) ierr=1
        end if
    end subroutine setup_problem_csc

    subroutine ecos_csolve_problem(prob,result,control)
        type(ecos_problem), intent(in) :: prob
        type(ecos_result), intent(out) :: result
        type(ecos_settings), intent(in), optional :: control
        type(ecos_settings) :: stg
        stg=ecos_settings(); if (present(control)) stg=control
        if (has_integer_variables(prob)) then
            call ecos_solve_mixed_integer(prob,result,stg)
        else
            call ecos_solve_continuous(prob,result,stg)
        end if
    end subroutine ecos_csolve_problem

    subroutine ecos_csolve_dense(c,g,h,dims,result,a,b,bool_vars,int_vars,control,ierr)
        real(dp), intent(in) :: c(:),g(:,:),h(:)
        type(ecos_dims), intent(in) :: dims
        type(ecos_result), intent(out) :: result
        real(dp), intent(in), optional :: a(:,:),b(:)
        integer, intent(in), optional :: bool_vars(:),int_vars(:)
        type(ecos_settings), intent(in), optional :: control
        integer, intent(out), optional :: ierr
        type(ecos_problem) :: prob
        integer :: ie
        call setup_problem_dense(prob,c,g,h,dims,a,b,bool_vars,int_vars,ie)
        if (present(ierr)) ierr=ie
        if (ie/=0) then
            allocate(result%x(size(c)),result%y(0),result%s(size(h)),result%z(size(h)))
            result%exitflag=ECOS_FATAL; result%infostring='Invalid problem data'; return
        end if
        call ecos_csolve_problem(prob,result,control)
    end subroutine ecos_csolve_dense

    subroutine ecos_csolve_csc(c,g,h,dims,result,a,b,bool_vars,int_vars,control,ierr)
        real(dp), intent(in) :: c(:),h(:)
        type(ecos_csc_matrix), intent(in) :: g
        type(ecos_dims), intent(in) :: dims
        type(ecos_result), intent(out) :: result
        type(ecos_csc_matrix), intent(in), optional :: a
        real(dp), intent(in), optional :: b(:)
        integer, intent(in), optional :: bool_vars(:),int_vars(:)
        type(ecos_settings), intent(in), optional :: control
        integer, intent(out), optional :: ierr
        type(ecos_problem) :: prob
        integer :: ie
        call setup_problem_csc(prob,c,g,h,dims,a,b,bool_vars,int_vars,ie)
        if (present(ierr)) ierr=ie
        if (ie/=0) then
            allocate(result%x(size(c)),result%y(0),result%s(size(h)),result%z(size(h)))
            result%exitflag=ECOS_FATAL; result%infostring='Invalid problem data'; return
        end if
        call ecos_csolve_problem(prob,result,control)
    end subroutine ecos_csolve_csc

    subroutine ecos_setup(ws,prob,control)
        type(ecos_workspace), intent(out) :: ws
        type(ecos_problem), intent(in) :: prob
        type(ecos_settings), intent(in), optional :: control
        ws%problem=prob
        ws%settings=ecos_settings(); if (present(control)) ws%settings=control
        ws%initialized=prob%valid()
    end subroutine ecos_setup

    subroutine ecos_solve(ws,result,control,ierr)
        type(ecos_workspace), intent(inout) :: ws
        type(ecos_result), intent(out) :: result
        type(ecos_settings), intent(in), optional :: control
        integer, intent(out), optional :: ierr
        type(ecos_settings) :: stg
        if (.not.ws%initialized) then
            if (present(ierr)) ierr=1
            allocate(result%x(0),result%y(0),result%s(0),result%z(0))
            result%exitflag=ECOS_FATAL; result%infostring='Workspace is cleaned up or invalid'; return
        end if
        stg=ws%settings
        if (present(control)) then
            if (control%sparse_amd .neqv. ws%settings%sparse_amd .or. &
                control%sparse_rcm .neqv. ws%settings%sparse_rcm) then
                ws%sparse_cache%symbolic_valid=.false.
            end if
            if (control%equilibrate .neqv. ws%settings%equilibrate .or. &
                control%equilibration_iters /= ws%settings%equilibration_iters) then
                ws%sparse_cache%warm_valid=.false.
            end if
            stg=control
            ws%settings=control
        end if
        if (has_integer_variables(ws%problem)) then
            call ecos_csolve_problem(ws%problem,result,stg)
        else
            call ecos_solve_continuous(ws%problem,result,stg,cache=ws%sparse_cache)
        end if
        if (present(ierr)) ierr=0
    end subroutine ecos_solve

    subroutine ecos_update(ws,c,h,b,gmat,amat,g_csc,a_csc,ierr)
        type(ecos_workspace), intent(inout) :: ws
        real(dp), intent(in), optional :: c(:),h(:),b(:),gmat(:,:),amat(:,:)
        type(ecos_csc_matrix), intent(in), optional :: g_csc,a_csc
        integer, intent(out), optional :: ierr
        type(ecos_csc_matrix) :: tmpc
        integer :: ie
        logical :: matrix_values_changed,matrix_structure_changed,same_structure
        ie=0; matrix_values_changed=.false.; matrix_structure_changed=.false.
        if(.not.ws%initialized) then; ie=1; goto 100; end if
        if(present(c)) then
            if(size(c)/=size(ws%problem%c)) then; ie=2; goto 100; end if
            ws%problem%c=c
        end if
        if(present(h)) then
            if(size(h)/=size(ws%problem%h)) then; ie=3; goto 100; end if
            ws%problem%h=h
        end if
        if(present(b)) then
            if(size(b)/=size(ws%problem%b)) then; ie=4; goto 100; end if
            ws%problem%b=b
        end if
        if(present(gmat)) then
            if(size(gmat,1)/=ws%problem%ncone() .or. size(gmat,2)/=ws%problem%nvar()) then
                ie=5; goto 100
            end if
            if(ws%problem%sparse_storage) then
                call make_csc_matrix(gmat,tmpc)
                same_structure=sparse_structure_equal(ws%problem%g_csc,tmpc)
                ws%problem%g_csc=tmpc; matrix_values_changed=.true.
                if(.not.same_structure) matrix_structure_changed=.true.
                call csc_to_csr(ws%problem%g_csc,ws%problem%g_csr)
            else
                ws%problem%gmat=gmat
            end if
        end if
        if(present(amat)) then
            if(size(amat,1)/=ws%problem%neq() .or. size(amat,2)/=ws%problem%nvar()) then
                ie=6; goto 100
            end if
            if(ws%problem%sparse_storage) then
                call make_csc_matrix(amat,tmpc)
                same_structure=sparse_structure_equal(ws%problem%a_csc,tmpc)
                ws%problem%a_csc=tmpc; matrix_values_changed=.true.
                if(.not.same_structure) matrix_structure_changed=.true.
                call csc_to_csr(ws%problem%a_csc,ws%problem%a_csr)
            else
                ws%problem%amat=amat
            end if
        end if
        if(present(g_csc)) then
            if(g_csc%nrow/=ws%problem%ncone() .or. g_csc%ncol/=ws%problem%nvar()) then
                ie=7; goto 100
            end if
            if(ws%problem%sparse_storage) then
                same_structure=sparse_structure_equal(ws%problem%g_csc,g_csc)
                ws%problem%g_csc=g_csc; call csc_to_csr(ws%problem%g_csc,ws%problem%g_csr)
                matrix_values_changed=.true.; if(.not.same_structure) matrix_structure_changed=.true.
            else
                call g_csc%to_dense(ws%problem%gmat)
            end if
        end if
        if(present(a_csc)) then
            if(a_csc%nrow/=ws%problem%neq() .or. a_csc%ncol/=ws%problem%nvar()) then
                ie=8; goto 100
            end if
            if(ws%problem%sparse_storage) then
                same_structure=sparse_structure_equal(ws%problem%a_csc,a_csc)
                ws%problem%a_csc=a_csc; call csc_to_csr(ws%problem%a_csc,ws%problem%a_csr)
                matrix_values_changed=.true.; if(.not.same_structure) matrix_structure_changed=.true.
            else
                call a_csc%to_dense(ws%problem%amat)
            end if
        end if
100     if (ie==0 .and. matrix_values_changed) ws%sparse_cache%warm_valid=.false.
        if (ie==0 .and. matrix_structure_changed) ws%sparse_cache%symbolic_valid=.false.
        if(present(ierr)) ierr=ie
    end subroutine ecos_update

    subroutine ecos_cleanup(ws)
        type(ecos_workspace), intent(inout) :: ws
        ws%initialized=.false.
        ws%sparse_cache%symbolic_valid=.false.
        ws%sparse_cache%warm_valid=.false.
    end subroutine ecos_cleanup

end module ecos_api

module ecos
    use ecos_api
    implicit none
    public
end module ecos
