! SPDX-License-Identifier: LGPL-2.0-only
module lpsolve_simplex
    use lpsolve_types, only: dp, LP_LE, LP_GE, LP_EQ, LP_OPTIMAL, LP_INFEASIBLE, &
        LP_UNBOUNDED, LP_NUMFAILURE, LP_TIMEOUT, lp_control
    implicit none
    private

    public :: simplex_relaxation

contains

    subroutine simplex_relaxation(c, a, sense, rhs, control, x, obj, status, &
        iterations, row_dual, reduced_cost)
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: sense(:)
        real(dp), intent(in) :: rhs(:)
        type(lp_control), intent(in) :: control
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: obj
        integer, intent(out) :: status
        integer, intent(out) :: iterations
        real(dp), intent(out), optional :: row_dual(:)
        real(dp), intent(out), optional :: reduced_cost(:)

        integer :: n, m, nt, i, j, nart
        integer, allocatable :: basis(:), slack_col(:), surplus_col(:), art_col(:)
        integer, allocatable :: nsense(:)
        real(dp), allocatable :: tab(:,:), row_factor(:)
        logical, allocatable :: artificial(:), allow(:)
        real(dp) :: scale, q
        integer :: phase_iter, st

        n = size(c)
        m = size(rhs)
        x = 0.0_dp
        obj = 0.0_dp
        iterations = 0
        status = LP_NUMFAILURE

        if (size(a,1) /= m .or. size(a,2) /= n .or. size(sense) /= m) return
        if (size(x) /= n) return

        allocate(nsense(m), row_factor(m), basis(m), slack_col(m), surplus_col(m), art_col(m))
        nsense = sense
        row_factor = 1.0_dp
        slack_col = 0
        surplus_col = 0
        art_col = 0

        nart = 0
        nt = n
        do i = 1, m
            if (rhs(i) < 0.0_dp) then
                row_factor(i) = -1.0_dp
                select case (nsense(i))
                case (LP_LE)
                    nsense(i) = LP_GE
                case (LP_GE)
                    nsense(i) = LP_LE
                end select
            end if
            select case (nsense(i))
            case (LP_LE)
                nt = nt + 1
                slack_col(i) = nt
            case (LP_GE)
                nt = nt + 1
                surplus_col(i) = nt
                nt = nt + 1
                art_col(i) = nt
                nart = nart + 1
            case (LP_EQ)
                nt = nt + 1
                art_col(i) = nt
                nart = nart + 1
            case default
                status = LP_NUMFAILURE
                return
            end select
        end do

        allocate(tab(0:m,1:nt+1), artificial(nt), allow(nt))
        tab = 0.0_dp
        artificial = .false.
        allow = .true.

        do i = 1, m
            scale = 1.0_dp
            if (control%scale_rows) then
                scale = max(1.0_dp, abs(rhs(i)))
                if (n > 0) scale = max(scale, maxval(abs(a(i,:))))
            end if
            row_factor(i) = row_factor(i) / scale
            tab(i,1:n) = row_factor(i) * a(i,:)
            tab(i,nt+1) = row_factor(i) * rhs(i)
            if (slack_col(i) > 0) then
                tab(i,slack_col(i)) = 1.0_dp
                basis(i) = slack_col(i)
            else if (surplus_col(i) > 0) then
                tab(i,surplus_col(i)) = -1.0_dp
                tab(i,art_col(i)) = 1.0_dp
                artificial(art_col(i)) = .true.
                basis(i) = art_col(i)
            else
                tab(i,art_col(i)) = 1.0_dp
                artificial(art_col(i)) = .true.
                basis(i) = art_col(i)
            end if
        end do

        ! Phase I: maximize minus the sum of artificial variables.
        tab(0,:) = 0.0_dp
        do j = 1, nt
            if (artificial(j)) tab(0,j) = 1.0_dp
        end do
        do i = 1, m
            if (artificial(basis(i))) tab(0,:) = tab(0,:) - tab(i,:)
        end do

        allow = .not. artificial
        call simplex_iterate(tab, basis, allow, control, phase_iter, st)
        iterations = iterations + phase_iter
        if (st == LP_UNBOUNDED) then
            status = LP_NUMFAILURE
            return
        else if (st /= LP_OPTIMAL) then
            status = st
            return
        end if

        if (tab(0,nt+1) < -control%feasibility_tol) then
            status = LP_INFEASIBLE
            return
        end if

        ! Remove zero-valued artificial variables from the basis whenever possible.
        do i = 1, m
            if (.not. artificial(basis(i))) cycle
            if (tab(i,nt+1) > control%feasibility_tol) then
                status = LP_INFEASIBLE
                return
            end if
            do j = 1, nt
                if (artificial(j)) cycle
                if (abs(tab(i,j)) > control%feasibility_tol) then
                    call pivot_tableau(tab, basis, i, j, control%optimality_tol)
                    exit
                end if
            end do
        end do

        ! Phase II objective.
        tab(0,:) = 0.0_dp
        tab(0,1:n) = -c
        do i = 1, m
            q = tab(0,basis(i))
            if (abs(q) > control%optimality_tol) tab(0,:) = tab(0,:) - q * tab(i,:)
        end do
        allow = .not. artificial

        call simplex_iterate(tab, basis, allow, control, phase_iter, st)
        iterations = iterations + phase_iter
        status = st
        if (status /= LP_OPTIMAL) return

        x = 0.0_dp
        do i = 1, m
            if (basis(i) >= 1 .and. basis(i) <= n) then
                x(basis(i)) = tab(i,nt+1)
            end if
        end do
        where (abs(x) <= control%feasibility_tol) x = 0.0_dp
        obj = tab(0,nt+1)

        if (present(reduced_cost)) then
            if (size(reduced_cost) == n) reduced_cost = -tab(0,1:n)
        end if

        if (present(row_dual)) then
            if (size(row_dual) == m) then
                row_dual = 0.0_dp
                do i = 1, m
                    if (slack_col(i) > 0) then
                        row_dual(i) = tab(0,slack_col(i))
                    else if (surplus_col(i) > 0) then
                        row_dual(i) = -tab(0,surplus_col(i))
                    else if (art_col(i) > 0) then
                        row_dual(i) = tab(0,art_col(i))
                    end if
                    row_dual(i) = row_factor(i) * row_dual(i)
                end do
            end if
        end if

    end subroutine simplex_relaxation


    subroutine pivot_tableau(tab, basis, leave, enter, tol)
        real(dp), intent(inout) :: tab(0:,1:)
        integer, intent(inout) :: basis(:)
        integer, intent(in) :: leave, enter
        real(dp), intent(in) :: tol
        integer :: i, m
        real(dp) :: piv, coeff

        m = ubound(tab,1)
        piv = tab(leave,enter)
        tab(leave,:) = tab(leave,:) / piv
        do i = 0, m
            if (i == leave) cycle
            coeff = tab(i,enter)
            if (abs(coeff) > tol) tab(i,:) = tab(i,:) - coeff * tab(leave,:)
        end do
        basis(leave) = enter
    end subroutine pivot_tableau


    subroutine simplex_iterate(tab, basis, allow, control, iterations, status)
        real(dp), intent(inout) :: tab(0:,1:)
        integer, intent(inout) :: basis(:)
        logical, intent(in) :: allow(:)
        type(lp_control), intent(in) :: control
        integer, intent(out) :: iterations
        integer, intent(out) :: status

        integer :: m, nt, enter, leave, i, j
        real(dp) :: ratio, best_ratio, piv, coeff
        real(dp) :: tol, rt, start_time, now

        m = ubound(tab,1)
        nt = ubound(tab,2) - 1
        tol = control%optimality_tol
        iterations = 0
        status = LP_OPTIMAL
        call cpu_time(start_time)

        do while (iterations < control%max_simplex_iter)
            if (control%timeout_seconds > 0.0_dp) then
                if (mod(iterations,100) == 0) then
                    call cpu_time(now)
                    if (now - start_time >= control%timeout_seconds) then
                        status = LP_TIMEOUT
                        return
                    end if
                end if
            end if
            enter = 0
            if (control%bland_rule) then
                do j = 1, nt
                    if (allow(j) .and. tab(0,j) < -tol) then
                        enter = j
                        exit
                    end if
                end do
            else
                rt = -tol
                do j = 1, nt
                    if (allow(j) .and. tab(0,j) < rt) then
                        rt = tab(0,j)
                        enter = j
                    end if
                end do
            end if

            if (enter == 0) then
                status = LP_OPTIMAL
                return
            end if

            leave = 0
            best_ratio = huge(1.0_dp)
            do i = 1, m
                coeff = tab(i,enter)
                if (coeff > control%feasibility_tol) then
                    ratio = tab(i,nt+1) / coeff
                    if (ratio < best_ratio - control%feasibility_tol) then
                        best_ratio = ratio
                        leave = i
                    else if (abs(ratio - best_ratio) <= control%feasibility_tol) then
                        if (leave == 0 .or. basis(i) < basis(leave)) leave = i
                    end if
                end if
            end do

            if (leave == 0) then
                status = LP_UNBOUNDED
                return
            end if

            piv = tab(leave,enter)
            if (abs(piv) <= control%feasibility_tol) then
                status = LP_NUMFAILURE
                return
            end if

            call pivot_tableau(tab, basis, leave, enter, tol)
            iterations = iterations + 1
        end do

        status = LP_NUMFAILURE
    end subroutine simplex_iterate

end module lpsolve_simplex
