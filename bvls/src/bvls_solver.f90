! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern free-form translation of Stark & Parker's BVLS active-set
! algorithm. The original fixed-form implementation is preserved under
! original/bvls-master/src/bvls.f.
module bvls_solver
    use bvls_kinds, only : dp
    use bvls_qr, only : householder_lsq
    use bvls_types
    implicit none
    private
    public :: bvls_fit, bvls_core

contains

    subroutine bvls_fit(a, b, lower, upper, result, key, istate)
        real(dp), intent(in) :: a(:,:), b(:), lower(:), upper(:)
        type(bvls_result), intent(out) :: result
        integer, intent(in), optional :: key
        integer, intent(in), optional :: istate(:)

        real(dp), allocatable :: x(:), w(:)
        integer, allocatable :: state(:)
        integer :: n, m, key_use, loop_a, status

        m = size(a, 1)
        n = size(a, 2)
        allocate(result%x(n), result%fitted(m), result%residuals(m), &
                 result%gradient(n), result%istate(n+1), x(n), w(n), state(n+1))
        result%x = 0.0_dp
        result%fitted = 0.0_dp
        result%residuals = 0.0_dp
        result%gradient = 0.0_dp
        result%istate = 0
        x = 0.0_dp
        w = 0.0_dp
        state = 0
        key_use = 0
        if (present(key)) key_use = key
        if (present(istate)) then
            if (size(istate) == n + 1) state = istate
        end if

        call bvls_core(key_use, a, b, lower, upper, x, w, state, loop_a, status)

        result%x = x
        if (m > 0 .and. n > 0) result%fitted = matmul(a, x)
        if (m > 0) result%residuals = b - result%fitted
        if (n > 0 .and. m > 0) result%gradient = matmul(transpose(a), result%residuals)
        result%deviance = sum(result%residuals**2)
        result%residual_norm = sqrt(max(result%deviance, 0.0_dp))
        result%istate = state
        result%iterations = loop_a
        result%status = status
    end subroutine bvls_fit

    subroutine bvls_core(key, a, b, lower, upper, x, w, istate, loop_a, status)
        integer, intent(in) :: key
        real(dp), intent(in) :: a(:,:), b(:), lower(:), upper(:)
        real(dp), intent(out) :: x(:), w(:)
        integer, intent(inout) :: istate(:)
        integer, intent(out) :: loop_a, status

        real(dp), allocatable :: residual_work(:), active_a(:,:), zz(:)
        real(dp) :: bdiff, bnorm, obj, worst, bad, bound, resq
        real(dp) :: alpha, alf
        integer :: m, n, mm, nbound, nact, jj, ifrom5, iact, it
        integer :: j, k, k1, noldb, sj, qr_info, temp_state
        logical :: feasible
        real(dp), parameter :: eps = 1.0e-11_dp

        m = size(a, 1)
        n = size(a, 2)
        x = 0.0_dp
        w = 0.0_dp
        loop_a = 0
        status = bvls_success

        if (m <= 0 .or. n <= 0 .or. size(b) /= m .or. size(lower) /= n .or. &
            size(upper) /= n .or. size(x) /= n .or. size(w) /= n .or. &
            size(istate) /= n + 1) then
            status = bvls_invalid_dimensions
            return
        end if

        mm = min(m, n)
        allocate(residual_work(m), active_a(m,mm), zz(mm))
        jj = 0
        ifrom5 = 0

        bdiff = 0.0_dp
        do j = 1, n
            if (lower(j) > upper(j)) then
                status = bvls_inconsistent_bounds
                return
            end if
            bdiff = max(bdiff, upper(j) - lower(j))
        end do
        if (.not. (bdiff > 0.0_dp .or. bdiff < 0.0_dp)) then
            x = lower
            istate(1:n) = [( -j, j=1,n )]
            istate(n+1) = n
            status = bvls_no_free_variables
            return
        end if

        if (key == 0) then
            nbound = n
            nact = 0
            do j = 1, n
                istate(j) = -j
                x(j) = lower(j)
            end do
        else
            nbound = istate(n+1)
            if (nbound < 0 .or. nbound > n) then
                status = bvls_invalid_state
                return
            end if
            nact = n - nbound
            if (nact > mm .or. .not. valid_state(istate, n, nbound)) then
                status = bvls_invalid_state
                return
            end if
            do k = 1, nbound
                j = abs(istate(k))
                if (istate(k) < 0) x(j) = lower(j)
                if (istate(k) > 0) x(j) = upper(j)
            end do
            do k = nbound + 1, n
                j = istate(k)
                x(j) = 0.5_dp * (upper(j) + lower(j))
            end do
        end if

        bnorm = sqrt(sum(b**2))

        do loop_a = 1, 3*n
            obj = 0.0_dp
            w = 0.0_dp
            residual_work = b - matmul(a, x)
            obj = sum(residual_work**2)
            w = matmul(transpose(a), residual_work)

            if (sqrt(obj) <= bnorm*eps .or. (loop_a > 1 .and. nbound == 0)) then
                istate(n+1) = nbound
                if (n > 0) w(1) = sqrt(obj)
                status = bvls_success
                return
            end if

            do k = nbound + 1, n
                j = istate(k)
                residual_work = residual_work + a(:,j) * x(j)
            end do

            if (loop_a == 1 .and. key /= 0) then
                call solve_active_set()
                if (status /= bvls_success) return
                cycle
            end if

            do
                worst = 0.0_dp
                it = 1
                iact = 0
                do j = 1, nbound
                    k = abs(istate(j))
                    bad = w(k) * real(sign(1, istate(j)), dp)
                    if (bad < worst) then
                        it = j
                        worst = bad
                        iact = k
                    end if
                end do

                if (worst >= 0.0_dp) then
                    istate(n+1) = nbound
                    if (n > 0) w(1) = sqrt(obj)
                    status = bvls_success
                    return
                end if

                if (iact == jj) then
                    w(jj) = 0.0_dp
                    cycle
                end if

                if (istate(it) > 0) bound = upper(iact)
                if (istate(it) < 0) bound = lower(iact)
                residual_work = residual_work + bound * a(:,iact)
                ifrom5 = istate(it)

                istate(it) = istate(nbound)
                nbound = nbound - 1
                nact = nact + 1
                istate(nbound+1) = iact
                if (nact > mm) then
                    status = bvls_rank_failure
                    return
                end if

                call solve_active_set()
                if (status /= bvls_success) return
                exit
            end do
        end do

        loop_a = 3*n
        istate(n+1) = nbound
        status = bvls_max_iterations

    contains

        subroutine solve_active_set()
            do
                if (nact <= 0) return
                do k = 1, nact
                    j = istate(nbound + nact + 1 - k)
                    active_a(:,k) = a(:,j)
                end do

                call householder_lsq(active_a(:,1:nact), residual_work, zz(1:nact), resq, qr_info)
                iact = istate(nbound + 1)

                if (qr_info /= 0 .or. resq < 0.0_dp .or. &
                    (ifrom5 > 0 .and. zz(nact) > upper(iact)) .or. &
                    (ifrom5 < 0 .and. zz(nact) < lower(iact))) then
                    nbound = nbound + 1
                    if (x(iact) < upper(iact)) then
                        istate(nbound) = -abs(istate(nbound))
                    else
                        istate(nbound) = abs(istate(nbound))
                    end if
                    nact = nact - 1
                    residual_work = residual_work - x(iact) * a(:,iact)
                    ifrom5 = 0
                    w(iact) = 0.0_dp
                    return
                end if

                if (ifrom5 /= 0) jj = 0
                ifrom5 = 0

                feasible = .true.
                k1 = 1
                do k = 1, nact
                    j = istate(k + nbound)
                    if (zz(nact+1-k) < lower(j) .or. zz(nact+1-k) > upper(j)) then
                        feasible = .false.
                        k1 = k
                        exit
                    end if
                end do

                if (feasible) then
                    do k = 1, nact
                        j = istate(k + nbound)
                        x(j) = zz(nact+1-k)
                    end do
                    return
                end if

                alpha = 2.0_dp
                sj = 0
                do k = k1, nact
                    j = istate(k + nbound)
                    alf = alpha
                    if (zz(nact+1-k) > upper(j)) &
                        alf = (upper(j)-x(j)) / (zz(nact+1-k)-x(j))
                    if (zz(nact+1-k) < lower(j)) &
                        alf = (lower(j)-x(j)) / (zz(nact+1-k)-x(j))
                    if (alf < alpha) then
                        alpha = alf
                        jj = j
                        if (zz(nact+1-k) - lower(j) > 0.0_dp) then
                            sj = 1
                        else
                            sj = -1
                        end if
                    end if
                end do

                do k = 1, nact
                    j = istate(k + nbound)
                    x(j) = x(j) + alpha * (zz(nact+1-k) - x(j))
                end do

                noldb = nbound
                do k = 1, nact
                    j = istate(k + noldb)
                    if ((upper(j)-x(j) <= 0.0_dp) .or. (j == jj .and. sj > 0)) then
                        x(j) = upper(j)
                        temp_state = istate(nbound+1)
                        istate(k+noldb) = temp_state
                        istate(nbound+1) = j
                        nbound = nbound + 1
                        residual_work = residual_work - upper(j) * a(:,j)
                    else if ((x(j)-lower(j) <= 0.0_dp) .or. (j == jj .and. sj < 0)) then
                        x(j) = lower(j)
                        temp_state = istate(nbound+1)
                        istate(k+noldb) = temp_state
                        istate(nbound+1) = -j
                        nbound = nbound + 1
                        residual_work = residual_work - lower(j) * a(:,j)
                    end if
                end do
                nact = n - nbound
                if (nact <= 0) return
            end do
        end subroutine solve_active_set

        logical function valid_state(state, nvar, nbd) result(ok)
            integer, intent(in) :: state(:), nvar, nbd
            logical, allocatable :: seen(:)
            integer :: q, idx

            allocate(seen(nvar))
            seen = .false.
            ok = .true.
            do q = 1, nvar
                if (q <= nbd) then
                    idx = abs(state(q))
                else
                    idx = state(q)
                end if
                if (idx < 1 .or. idx > nvar) then
                    ok = .false.
                    return
                end if
                if (seen(idx)) then
                    ok = .false.
                    return
                end if
                seen(idx) = .true.
            end do
        end function valid_state

    end subroutine bvls_core

end module bvls_solver
