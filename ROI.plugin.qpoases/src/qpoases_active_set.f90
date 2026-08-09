! SPDX-License-Identifier: LGPL-2.1-or-later
module qpoases_active_set
    use qpoases_kinds, only : dp
    use qpoases_linalg, only : solve_linear, norm_inf, norm2_vec
    use qpoases_types, only : qpoases_options, ret_qp_solved, ret_qp_infeasible, &
        ret_max_nwsr_reached, ret_hessian_not_spd
    implicit none
    private
    public :: solve_working_set, find_feasible_point

contains

    subroutine project_equalities(xref, e_mat, e_rhs, x, ok)
        real(dp), intent(in) :: xref(:)
        real(dp), intent(in) :: e_mat(:,:), e_rhs(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: kkt(:,:), rhs(:), sol(:)
        integer :: n, me, info, i

        n = size(xref)
        me = size(e_rhs)
        if (me == 0) then
            x = xref
            ok = .true.
            return
        end if

        allocate(kkt(n+me,n+me), rhs(n+me), sol(n+me))
        kkt = 0.0_dp
        do i = 1, n
            kkt(i,i) = 1.0_dp
        end do
        kkt(1:n,n+1:n+me) = -e_mat
        kkt(n+1:n+me,1:n) = transpose(e_mat)
        rhs(1:n) = xref
        rhs(n+1:n+me) = e_rhs
        call solve_linear(kkt, rhs, sol, info)
        ok = info == 0
        if (ok) then
            x = sol(1:n)
            ok = norm_inf(matmul(transpose(e_mat),x) - e_rhs) <= &
                1000.0_dp * sqrt(epsilon(1.0_dp)) * max(1.0_dp, norm_inf(e_rhs))
        end if
    end subroutine project_equalities

    subroutine solve_working_set(h, g, e_mat, e_rhs, c_mat, c_rhs, x, options, &
                                 max_nwsr, warm_active, active, lambda_eq, &
                                 lambda_ineq, status, nwsr)
        real(dp), intent(in) :: h(:,:), g(:)
        real(dp), intent(in) :: e_mat(:,:), e_rhs(:)
        real(dp), intent(in) :: c_mat(:,:), c_rhs(:)
        real(dp), intent(inout) :: x(:)
        type(qpoases_options), intent(in) :: options
        integer, intent(in) :: max_nwsr
        logical, intent(in), optional :: warm_active(:)
        logical, allocatable, intent(out) :: active(:)
        real(dp), allocatable, intent(out) :: lambda_eq(:), lambda_ineq(:)
        integer, intent(out) :: status, nwsr

        real(dp), allocatable :: grad(:), p(:), aw(:,:), kkt(:,:), rhs(:), sol(:)
        real(dp), allocatable :: lamb(:), hwork(:,:)
        real(dp) :: tol, btol, alpha, ratio, cp, slack, reg, reg0
        integer :: n, me, mi, nw, i, j, info, blocker, drop_idx, attempt
        integer, allocatable :: widx(:)
        logical :: solved_linear

        n = size(g)
        me = size(e_rhs)
        mi = size(c_rhs)
        allocate(active(mi), lambda_eq(me), lambda_ineq(mi))
        active = .false.
        lambda_eq = 0.0_dp
        lambda_ineq = 0.0_dp
        if (present(warm_active)) then
            if (size(warm_active) == mi) active = warm_active
        end if

        tol = max(options%termination_tolerance, 100.0_dp * epsilon(1.0_dp))
        btol = max(options%bound_tolerance, 100.0_dp * epsilon(1.0_dp))

        ! Add constraints already tight at the feasible starting point.
        do j = 1, mi
            if (dot_product(c_mat(:,j),x) - c_rhs(j) <= btol) active(j) = .true.
        end do

        nwsr = 0
        status = ret_max_nwsr_reached
        allocate(grad(n), p(n), hwork(n,n))

        do while (nwsr < max_nwsr)
            nw = count(active)
            allocate(widx(nw))
            if (nw > 0) then
                i = 0
                do j = 1, mi
                    if (active(j)) then
                        i = i + 1
                        widx(i) = j
                    end if
                end do
            end if

            allocate(aw(n,me+nw), kkt(n+me+nw,n+me+nw), &
                     rhs(n+me+nw), sol(n+me+nw), lamb(me+nw))
            if (me > 0) aw(:,1:me) = e_mat
            if (nw > 0) then
                do i = 1, nw
                    aw(:,me+i) = c_mat(:,widx(i))
                end do
            end if

            grad = matmul(h,x) + g
            solved_linear = .false.
            reg0 = max(options%eps_regularisation, 1000.0_dp * epsilon(1.0_dp))
            reg = 0.0_dp
            do attempt = 0, max(2, options%num_regularisation_steps + 2)
                hwork = h
                if (attempt > 0) then
                    reg = reg0 * 10.0_dp**real(attempt-1,dp)
                    do i = 1, n
                        hwork(i,i) = hwork(i,i) + reg
                    end do
                end if
                kkt = 0.0_dp
                kkt(1:n,1:n) = hwork
                if (me + nw > 0) then
                    kkt(1:n,n+1:n+me+nw) = -aw
                    kkt(n+1:n+me+nw,1:n) = transpose(aw)
                end if
                rhs = 0.0_dp
                rhs(1:n) = -grad
                call solve_linear(kkt, rhs, sol, info)
                if (info == 0) then
                    solved_linear = .true.
                    exit
                end if
            end do

            if (.not. solved_linear) then
                ! A newly added inequality can make the working set dependent.
                if (nw > 0) then
                    active(widx(nw)) = .false.
                    nwsr = nwsr + 1
                    deallocate(widx, aw, kkt, rhs, sol, lamb)
                    cycle
                end if
                status = ret_hessian_not_spd
                deallocate(widx, aw, kkt, rhs, sol, lamb)
                return
            end if

            p = sol(1:n)
            if (me + nw > 0) lamb = sol(n+1:n+me+nw)

            if (norm2_vec(p) <= tol * max(1.0_dp,norm2_vec(x))) then
                if (me > 0) lambda_eq = lamb(1:me)
                lambda_ineq = 0.0_dp
                if (nw > 0) then
                    do i = 1, nw
                        lambda_ineq(widx(i)) = lamb(me+i)
                    end do
                end if

                drop_idx = 0
                if (nw > 0) then
                    do i = 1, nw
                        if (lamb(me+i) < -btol) then
                            if (drop_idx == 0) then
                                drop_idx = i
                            else if (lamb(me+i) < lamb(me+drop_idx)) then
                                drop_idx = i
                            end if
                        end if
                    end do
                end if

                if (drop_idx == 0) then
                    status = ret_qp_solved
                    deallocate(widx, aw, kkt, rhs, sol, lamb)
                    return
                end if

                active(widx(drop_idx)) = .false.
                nwsr = nwsr + 1
                deallocate(widx, aw, kkt, rhs, sol, lamb)
                cycle
            end if

            alpha = 1.0_dp
            blocker = 0
            do j = 1, mi
                if (active(j)) cycle
                cp = dot_product(c_mat(:,j), p)
                if (cp < -btol) then
                    slack = dot_product(c_mat(:,j),x) - c_rhs(j)
                    ratio = max(0.0_dp, slack) / (-cp)
                    if (ratio < alpha) then
                        alpha = ratio
                        blocker = j
                    end if
                end if
            end do

            x = x + alpha * p
            ! Remove tiny numerical infeasibility in equalities only through the KKT step.
            if (blocker > 0 .and. alpha < 1.0_dp - btol) then
                active(blocker) = .true.
                nwsr = nwsr + 1
            else
                nwsr = nwsr + 1
            end if

            deallocate(widx, aw, kkt, rhs, sol, lamb)
        end do
    end subroutine solve_working_set

    subroutine find_feasible_point(xref, e_mat, e_rhs, c_mat, c_rhs, options, &
                                   max_nwsr, x, status, nwsr)
        real(dp), intent(in) :: xref(:)
        real(dp), intent(in) :: e_mat(:,:), e_rhs(:)
        real(dp), intent(in) :: c_mat(:,:), c_rhs(:)
        type(qpoases_options), intent(in) :: options
        integer, intent(in) :: max_nwsr
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: status, nwsr

        real(dp), allocatable :: xeq(:), z(:), hz(:,:), gz(:), ez(:,:), cz(:,:)
        real(dp), allocatable :: le(:), li(:)
        real(dp) :: t0, viol, epsx, ftol
        logical, allocatable :: active(:)
        logical :: ok
        integer :: n, me, mi, j, inner_status, inner_nwsr

        n = size(xref)
        me = size(e_rhs)
        mi = size(c_rhs)
        allocate(xeq(n))
        call project_equalities(xref, e_mat, e_rhs, xeq, ok)
        if (.not. ok) then
            x = xref
            status = ret_qp_infeasible
            nwsr = 0
            return
        end if

        ftol = max(10.0_dp * options%bound_tolerance, 1.0e-10_dp)
        if (mi == 0) then
            x = xeq
            status = ret_qp_solved
            nwsr = 0
            return
        end if

        viol = maxval(c_rhs - matmul(transpose(c_mat),xeq))
        if (viol <= ftol) then
            x = xeq
            status = ret_qp_solved
            nwsr = 0
            return
        end if

        allocate(z(n+1), hz(n+1,n+1), gz(n+1), ez(n+1,me), cz(n+1,mi+1))
        t0 = max(0.0_dp, viol) + 1.0_dp
        z(1:n) = xeq
        z(n+1) = t0

        epsx = 1.0e-12_dp
        hz = 0.0_dp
        do j = 1, n
            hz(j,j) = epsx
        end do
        hz(n+1,n+1) = 1.0_dp
        gz = 0.0_dp
        gz(1:n) = -epsx * xeq

        if (me > 0) then
            ez = 0.0_dp
            ez(1:n,:) = e_mat
        end if
        cz = 0.0_dp
        cz(1:n,1:mi) = c_mat
        cz(n+1,1:mi) = 1.0_dp
        cz(n+1,mi+1) = 1.0_dp

        call solve_working_set(hz, gz, ez, e_rhs, cz, [c_rhs,0.0_dp], z, &
            options, max_nwsr, active=active, lambda_eq=le, lambda_ineq=li, &
            status=inner_status, nwsr=inner_nwsr)
        nwsr = inner_nwsr
        x = z(1:n)

        if (inner_status /= ret_qp_solved) then
            status = ret_qp_infeasible
            return
        end if
        if (z(n+1) > ftol) then
            status = ret_qp_infeasible
            return
        end if
        if (maxval(c_rhs - matmul(transpose(c_mat),x)) > 100.0_dp * ftol) then
            status = ret_qp_infeasible
            return
        end if
        status = ret_qp_solved
    end subroutine find_feasible_point
end module qpoases_active_set
