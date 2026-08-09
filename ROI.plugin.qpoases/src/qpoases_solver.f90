! SPDX-License-Identifier: LGPL-2.1-or-later
module qpoases_solver
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use qpoases_kinds, only : dp
    use qpoases_linalg, only : solve_linear, norm_inf
    use qpoases_types
    use qpoases_active_set, only : find_feasible_point, solve_working_set
    implicit none
    private
    public :: solve_general_qp

contains

    pure logical function finite_value(x)
        real(dp), intent(in) :: x
        finite_value = ieee_is_finite(x)
    end function finite_value

    subroutine build_standard_constraints(n, m, a, lb, ub, lba, uba, tol, &
                                          e_mat, e_rhs, e_id, c_mat, c_rhs, c_id, status)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: a(m,n), lb(n), ub(n), lba(m), uba(m), tol
        real(dp), allocatable, intent(out) :: e_mat(:,:), e_rhs(:)
        integer, allocatable, intent(out) :: e_id(:)
        real(dp), allocatable, intent(out) :: c_mat(:,:), c_rhs(:)
        integer, allocatable, intent(out) :: c_id(:)
        integer, intent(out) :: status

        real(dp), allocatable :: etmp(:,:), ctmp(:,:), ertmp(:), crtmp(:)
        integer, allocatable :: eitmp(:), citmp(:)
        integer :: me, mi, i
        real(dp) :: scale

        allocate(etmp(n,n+m), ertmp(n+m), eitmp(n+m))
        allocate(ctmp(n,2*n+2*m), crtmp(2*n+2*m), citmp(2*n+2*m))
        me = 0
        mi = 0
        status = successful_return

        do i = 1, n
            if (finite_value(lb(i)) .and. finite_value(ub(i))) then
                if (lb(i) > ub(i) + tol * max(1.0_dp,abs(lb(i)),abs(ub(i)))) then
                    status = ret_invalid_arguments
                    exit
                end if
                scale = max(1.0_dp,abs(lb(i)),abs(ub(i)))
                if (abs(ub(i)-lb(i)) <= tol * scale) then
                    me = me + 1
                    etmp(:,me) = 0.0_dp
                    etmp(i,me) = 1.0_dp
                    ertmp(me) = 0.5_dp * (lb(i)+ub(i))
                    eitmp(me) = i
                    cycle
                end if
            end if
            if (finite_value(lb(i))) then
                mi = mi + 1
                ctmp(:,mi) = 0.0_dp
                ctmp(i,mi) = 1.0_dp
                crtmp(mi) = lb(i)
                citmp(mi) = i
            end if
            if (finite_value(ub(i))) then
                mi = mi + 1
                ctmp(:,mi) = 0.0_dp
                ctmp(i,mi) = -1.0_dp
                crtmp(mi) = -ub(i)
                citmp(mi) = -i
            end if
        end do
        if (status /= successful_return) then
            allocate(e_mat(n,0),e_rhs(0),e_id(0),c_mat(n,0),c_rhs(0),c_id(0))
            return
        end if

        do i = 1, m
            if (finite_value(lba(i)) .and. finite_value(uba(i))) then
                if (lba(i) > uba(i) + tol * max(1.0_dp,abs(lba(i)),abs(uba(i)))) then
                    status = ret_invalid_arguments
                    exit
                end if
                scale = max(1.0_dp,abs(lba(i)),abs(uba(i)))
                if (abs(uba(i)-lba(i)) <= tol * scale) then
                    me = me + 1
                    etmp(:,me) = a(i,:)
                    ertmp(me) = 0.5_dp * (lba(i)+uba(i))
                    eitmp(me) = n + i
                    cycle
                end if
            end if
            if (finite_value(lba(i))) then
                mi = mi + 1
                ctmp(:,mi) = a(i,:)
                crtmp(mi) = lba(i)
                citmp(mi) = n + i
            end if
            if (finite_value(uba(i))) then
                mi = mi + 1
                ctmp(:,mi) = -a(i,:)
                crtmp(mi) = -uba(i)
                citmp(mi) = -(n+i)
            end if
        end do
        if (status /= successful_return) then
            allocate(e_mat(n,0),e_rhs(0),e_id(0),c_mat(n,0),c_rhs(0),c_id(0))
            return
        end if

        allocate(e_mat(n,me),e_rhs(me),e_id(me))
        allocate(c_mat(n,mi),c_rhs(mi),c_id(mi))
        if (me > 0) then
            e_mat = etmp(:,1:me)
            e_rhs = ertmp(1:me)
            e_id = eitmp(1:me)
        end if
        if (mi > 0) then
            c_mat = ctmp(:,1:mi)
            c_rhs = crtmp(1:mi)
            c_id = citmp(1:mi)
        end if
    end subroutine build_standard_constraints

    subroutine prepare_hessian(h_in, hessian_type, h)
        real(dp), intent(in) :: h_in(:,:)
        integer, intent(in) :: hessian_type
        real(dp), intent(out) :: h(size(h_in,1),size(h_in,2))
        integer :: i, n

        n = size(h,1)
        select case (hessian_type)
        case (hst_zero)
            h = 0.0_dp
        case (hst_identity)
            h = 0.0_dp
            do i = 1, n
                h(i,i) = 1.0_dp
            end do
        case default
            h = 0.5_dp * (h_in + transpose(h_in))
        end select
    end subroutine prepare_hessian

    subroutine initial_reference(h, g, hessian_type, xref)
        real(dp), intent(in) :: h(:,:), g(:)
        integer, intent(in) :: hessian_type
        real(dp), intent(out) :: xref(:)
        integer :: n, info
        real(dp), allocatable :: rhs(:)

        n = size(g)
        xref = 0.0_dp
        if (hessian_type == hst_zero) return
        if (hessian_type == hst_identity) then
            xref = -g
            return
        end if
        allocate(rhs(n))
        rhs = -g
        call solve_linear(h, rhs, xref, info)
        if (info /= 0) xref = 0.0_dp
    end subroutine initial_reference

    subroutine map_duals(n, m, e_id, c_id, active, le, li, y)
        integer, intent(in) :: n, m
        integer, intent(in) :: e_id(:), c_id(:)
        logical, intent(in) :: active(:)
        real(dp), intent(in) :: le(:), li(:)
        real(dp), intent(out) :: y(n+m)
        integer :: i, id

        y = 0.0_dp
        do i = 1, size(e_id)
            id = e_id(i)
            if (id > 0 .and. id <= n+m) y(id) = y(id) + le(i)
        end do
        do i = 1, size(c_id)
            if (.not. active(i)) cycle
            id = c_id(i)
            if (id > 0) then
                y(id) = y(id) + li(i)
            else
                y(-id) = y(-id) - li(i)
            end if
        end do
    end subroutine map_duals


    subroutine zero_hessian_unbounded(g, e_mat, c_mat, options, max_nwsr, unbounded)
        real(dp), intent(in) :: g(:), e_mat(:,:), c_mat(:,:)
        type(qpoases_options), intent(in) :: options
        integer, intent(in) :: max_nwsr
        logical, intent(out) :: unbounded
        real(dp), allocatable :: h(:,:), cdir(:,:), rhsdir(:), d(:), le(:), li(:)
        logical, allocatable :: active(:)
        integer :: n, mi, i, status, nwsr
        real(dp) :: scale

        n = size(g)
        mi = size(c_mat,2)
        unbounded = .false.
        scale = max(1.0_dp,maxval(abs(g)))
        if (maxval(abs(g)) <= sqrt(epsilon(1.0_dp))*scale) return

        allocate(h(n,n),cdir(n,mi+2*n),rhsdir(mi+2*n),d(n))
        h = 0.0_dp
        do i = 1, n
            h(i,i) = 1.0e-12_dp
        end do
        if (mi > 0) cdir(:,1:mi) = c_mat
        rhsdir = 0.0_dp
        do i = 1, n
            cdir(:,mi+2*i-1) = 0.0_dp
            cdir(i,mi+2*i-1) = 1.0_dp
            rhsdir(mi+2*i-1) = -1.0_dp
            cdir(:,mi+2*i) = 0.0_dp
            cdir(i,mi+2*i) = -1.0_dp
            rhsdir(mi+2*i) = -1.0_dp
        end do
        d = 0.0_dp
        block
            real(dp), allocatable :: erhs0(:)
            allocate(erhs0(size(e_mat,2)))
            erhs0 = 0.0_dp
            call solve_working_set(h,g,e_mat,erhs0,cdir,rhsdir,d,options, &
            max(50,min(max_nwsr,1000)),active=active,lambda_eq=le,lambda_ineq=li, &
                status=status,nwsr=nwsr)
        end block
        if (status == ret_qp_solved) then
            if (dot_product(g,d) < -sqrt(epsilon(1.0_dp))*scale .and. &
                maxval(abs(d)) > sqrt(epsilon(1.0_dp))) unbounded = .true.
        end if
    end subroutine zero_hessian_unbounded

    subroutine solve_general_qp(h_in, g, a, lb, ub, lba, uba, hessian_type, options, &
                                max_nwsr, result, x_start, warm_ids, active_ids)
        real(dp), intent(in) :: h_in(:,:), g(:), a(:,:), lb(:), ub(:), lba(:), uba(:)
        integer, intent(in) :: hessian_type
        type(qpoases_options), intent(in) :: options
        integer, intent(in) :: max_nwsr
        type(qpoases_result), intent(out) :: result
        real(dp), intent(in), optional :: x_start(:)
        integer, intent(in), optional :: warm_ids(:)
        integer, allocatable, intent(out), optional :: active_ids(:)

        real(dp), allocatable :: h(:,:), e_mat(:,:), e_rhs(:), c_mat(:,:), c_rhs(:)
        integer, allocatable :: e_id(:), c_id(:), ids(:)
        real(dp), allocatable :: x(:), xref(:), le(:), li(:)
        logical, allocatable :: active(:), warm_active(:)
        integer :: n, m, status, phase_nwsr, solve_nwsr, i, k
        real(dp) :: tol
        logical :: unb

        n = size(g)
        m = size(lba)
        allocate(result%x(n), result%y(n+m))
        result%x = 0.0_dp
        result%y = 0.0_dp
        result%status = ret_invalid_arguments

        if (n <= 0 .or. size(h_in,1) /= n .or. size(h_in,2) /= n) return
        if (size(a,1) /= m .or. size(a,2) /= n) return
        if (size(lb) /= n .or. size(ub) /= n .or. size(uba) /= m) return
        if (hessian_type < hst_zero .or. hessian_type > hst_unknown) return
        if (max_nwsr <= 0) return

        allocate(h(n,n), x(n), xref(n))
        call prepare_hessian(h_in,hessian_type,h)
        tol = max(options%bound_tolerance,100.0_dp*epsilon(1.0_dp))
        call build_standard_constraints(n,m,a,lb,ub,lba,uba,tol,e_mat,e_rhs,e_id, &
                                        c_mat,c_rhs,c_id,status)
        if (status /= successful_return) then
            result%status = status
            return
        end if

        if (present(x_start)) then
            if (size(x_start) == n) then
                xref = x_start
            else
                call initial_reference(h,g,hessian_type,xref)
            end if
        else
            call initial_reference(h,g,hessian_type,xref)
        end if

        call find_feasible_point(xref,e_mat,e_rhs,c_mat,c_rhs,options,max_nwsr,x, &
                                 status,phase_nwsr)
        if (status /= ret_qp_solved) then
            result%status = ret_qp_infeasible
            result%infeasible = .true.
            result%initialised = .true.
            result%nwsr = phase_nwsr
            result%x = x
            return
        end if

        if (hessian_type == hst_zero) then
            call zero_hessian_unbounded(g,e_mat,c_mat,options,max_nwsr,unb)
            if (unb) then
                result%status = ret_qp_unbounded
                result%unbounded = .true.
                result%initialised = .true.
                result%x = x
                result%objval = dot_product(g,x)
                result%nwsr = phase_nwsr
                return
            end if
        end if

        allocate(warm_active(size(c_rhs)))
        warm_active = .false.
        if (present(warm_ids)) then
            do i = 1, size(c_id)
                if (any(warm_ids == c_id(i))) warm_active(i) = .true.
            end do
        end if

        call solve_working_set(h,g,e_mat,e_rhs,c_mat,c_rhs,x,options, &
            max(1,max_nwsr-phase_nwsr),warm_active,active,le,li,status,solve_nwsr)

        result%x = x
        result%objval = 0.5_dp * dot_product(x,matmul(h,x)) + dot_product(g,x)
        result%status = status
        result%nwsr = phase_nwsr + solve_nwsr
        result%initialised = .true.
        result%solved = status == ret_qp_solved
        result%infeasible = status == ret_qp_infeasible
        result%unbounded = status == ret_qp_unbounded
        if (allocated(active)) then
            call map_duals(n,m,e_id,c_id,active,le,li,result%y)
            result%n_fixed = 0
            result%n_active_constraints = 0
            do i = 1, size(c_id)
                if (.not. active(i)) cycle
                if (abs(c_id(i)) <= n) then
                    result%n_fixed = result%n_fixed + 1
                else
                    result%n_active_constraints = result%n_active_constraints + 1
                end if
            end do
            do i = 1, size(e_id)
                if (e_id(i) <= n) then
                    result%n_fixed = result%n_fixed + 1
                else
                    result%n_active_constraints = result%n_active_constraints + 1
                end if
            end do
            result%n_fixed = min(n,result%n_fixed)
            result%n_free = n - result%n_fixed
            result%n_inactive_constraints = max(0,m-result%n_active_constraints)

            if (present(active_ids)) then
                k = count(active)
                allocate(ids(k))
                k = 0
                do i = 1, size(c_id)
                    if (active(i)) then
                        k = k + 1
                        ids(k) = c_id(i)
                    end if
                end do
                call move_alloc(ids,active_ids)
            end if
        else
            result%n_free = n
            result%n_fixed = 0
            result%n_active_constraints = count(abs(e_id) > n)
            result%n_inactive_constraints = max(0,m-result%n_active_constraints)
            if (present(active_ids)) allocate(active_ids(0))
        end if
    end subroutine solve_general_qp
end module qpoases_solver
