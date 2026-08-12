! SPDX-License-Identifier: LGPL-2.0-only
module lpsolve_core
    use lpsolve_types
    use lpsolve_simplex, only: simplex_relaxation
    implicit none
    private

    public :: solve_lp, solve_lp_sparse

contains

    subroutine solve_lp(direction, objective, a, sense, rhs, result, control, &
        integer_variables, binary_variables)
        integer, intent(in) :: direction
        real(dp), intent(in) :: objective(:)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: sense(:)
        real(dp), intent(in) :: rhs(:)
        type(lp_result), intent(out) :: result
        type(lp_control), intent(in), optional :: control
        integer, intent(in), optional :: integer_variables(:)
        integer, intent(in), optional :: binary_variables(:)

        type(lp_control) :: ctl
        logical, allocatable :: int_mask(:), bin_mask(:)
        integer :: n, k, want, found, total_nodes, total_iters
        real(dp), allocatable :: aw(:,:), bw(:)
        integer, allocatable :: sw(:)
        type(lp_result) :: one
        real(dp), allocatable :: pool(:,:), pobj(:)

        ctl = lp_control()
        if (present(control)) ctl = control
        call clear_result(result)

        n = size(objective)
        if (direction /= LP_MIN .and. direction /= LP_MAX) then
            result%status = LP_NUMFAILURE
            return
        end if
        if (size(a,2) /= n .or. size(a,1) /= size(rhs) .or. size(sense) /= size(rhs)) then
            result%status = LP_NUMFAILURE
            return
        end if

        allocate(int_mask(n), bin_mask(n))
        int_mask = .false.
        bin_mask = .false.
        result%status = LP_OPTIMAL
        if (present(integer_variables)) call set_index_mask(integer_variables, int_mask, result%status)
        if (result%status == LP_NUMFAILURE) return
        if (present(binary_variables)) call set_index_mask(binary_variables, bin_mask, result%status)
        if (result%status == LP_NUMFAILURE) return
        int_mask = int_mask .or. bin_mask

        if (.not. any(int_mask)) then
            call solve_continuous(direction, objective, a, sense, rhs, ctl, result)
            return
        end if

        want = 1
        if (all(bin_mask) .and. ctl%num_binary_solutions > 1) want = ctl%num_binary_solutions
        allocate(pool(n,want), pobj(want))
        pool = 0.0_dp
        pobj = 0.0_dp
        found = 0
        total_nodes = 0
        total_iters = 0

        do k = 1, want
            call append_exclusions(a, sense, rhs, pool, k-1, aw, sw, bw)
            call solve_mip_single(direction, objective, aw, sw, bw, int_mask, bin_mask, ctl, one)
            if (one%status /= LP_OPTIMAL .and. one%status /= LP_SUBOPTIMAL) then
                if (k == 1) result = one
                exit
            end if
            pool(:,k) = one%solution
            pobj(k) = one%objective
            found = k
            total_nodes = total_nodes + one%nodes
            total_iters = total_iters + one%simplex_iterations
            if (k == 1) result = one
        end do

        result%nodes = total_nodes
        result%simplex_iterations = total_iters
        result%solution_count = found
        if (found > 0) then
            allocate(result%solutions(n,found))
            allocate(result%solution_objectives(found))
            result%solutions = pool(:,1:found)
            result%solution_objectives = pobj(1:found)
        end if
    end subroutine solve_lp


    subroutine solve_lp_sparse(direction, objective, sparse, sense, rhs, result, &
        control, integer_variables, binary_variables)
        integer, intent(in) :: direction
        real(dp), intent(in) :: objective(:)
        type(sparse_constraints), intent(in) :: sparse
        integer, intent(in) :: sense(:)
        real(dp), intent(in) :: rhs(:)
        type(lp_result), intent(out) :: result
        type(lp_control), intent(in), optional :: control
        integer, intent(in), optional :: integer_variables(:)
        integer, intent(in), optional :: binary_variables(:)

        real(dp), allocatable :: a(:,:)
        integer :: k

        call clear_result(result)
        if (sparse%ncol /= size(objective) .or. sparse%nrow /= size(rhs)) then
            result%status = LP_NUMFAILURE
            return
        end if
        if (.not. allocated(sparse%row) .or. .not. allocated(sparse%col) .or. &
            .not. allocated(sparse%val)) then
            result%status = LP_NUMFAILURE
            return
        end if
        if (size(sparse%row) /= size(sparse%col) .or. size(sparse%row) /= size(sparse%val)) then
            result%status = LP_NUMFAILURE
            return
        end if

        allocate(a(sparse%nrow,sparse%ncol))
        a = 0.0_dp
        do k = 1, size(sparse%val)
            if (sparse%row(k) < 1 .or. sparse%row(k) > sparse%nrow .or. &
                sparse%col(k) < 1 .or. sparse%col(k) > sparse%ncol) then
                result%status = LP_NUMFAILURE
                return
            end if
            a(sparse%row(k),sparse%col(k)) = a(sparse%row(k),sparse%col(k)) + sparse%val(k)
        end do

        if (present(control)) then
            if (present(integer_variables) .and. present(binary_variables)) then
                call solve_lp(direction, objective, a, sense, rhs, result, control, &
                    integer_variables, binary_variables)
            else if (present(integer_variables)) then
                call solve_lp(direction, objective, a, sense, rhs, result, control, &
                    integer_variables=integer_variables)
            else if (present(binary_variables)) then
                call solve_lp(direction, objective, a, sense, rhs, result, control, &
                    binary_variables=binary_variables)
            else
                call solve_lp(direction, objective, a, sense, rhs, result, control)
            end if
        else
            if (present(integer_variables) .and. present(binary_variables)) then
                call solve_lp(direction, objective, a, sense, rhs, result, &
                    integer_variables=integer_variables, binary_variables=binary_variables)
            else if (present(integer_variables)) then
                call solve_lp(direction, objective, a, sense, rhs, result, &
                    integer_variables=integer_variables)
            else if (present(binary_variables)) then
                call solve_lp(direction, objective, a, sense, rhs, result, &
                    binary_variables=binary_variables)
            else
                call solve_lp(direction, objective, a, sense, rhs, result)
            end if
        end if
    end subroutine solve_lp_sparse


    subroutine solve_continuous(direction, objective, a, sense, rhs, control, result, &
        lower, upper)
        integer, intent(in) :: direction
        real(dp), intent(in) :: objective(:)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: sense(:)
        real(dp), intent(in) :: rhs(:)
        type(lp_control), intent(in) :: control
        type(lp_result), intent(out) :: result
        real(dp), intent(in), optional :: lower(:), upper(:)

        integer :: n, m, me, i, r, st, it
        real(dp) :: sign_dir, objmax
        real(dp), allocatable :: ae(:,:), be(:), cmax(:), x(:), y(:), rc(:)
        integer, allocatable :: se(:)
        real(dp), allocatable :: lo(:), up(:)

        call clear_result(result)
        n = size(objective)
        m = size(rhs)
        sign_dir = merge(1.0_dp, -1.0_dp, direction == LP_MAX)

        allocate(lo(n), up(n))
        lo = 0.0_dp
        up = LP_INFINITY
        if (present(lower)) lo = lower
        if (present(upper)) up = upper
        if (size(lo) /= n .or. size(up) /= n) then
            result%status = LP_NUMFAILURE
            return
        end if
        if (any(lo < -control%feasibility_tol) .or. &
            any(up < lo - control%feasibility_tol)) then
            result%status = LP_INFEASIBLE
            return
        end if

        me = m
        do i = 1, n
            if (lo(i) > control%feasibility_tol) me = me + 1
            if (up(i) < LP_INFINITY) me = me + 1
        end do
        allocate(ae(me,n), be(me), se(me))
        if (m > 0) then
            ae(1:m,:) = a
            be(1:m) = rhs
            se(1:m) = sense
        end if
        r = m
        do i = 1, n
            if (lo(i) > control%feasibility_tol) then
                r = r + 1
                ae(r,:) = 0.0_dp
                ae(r,i) = 1.0_dp
                be(r) = lo(i)
                se(r) = LP_GE
            end if
            if (up(i) < LP_INFINITY) then
                r = r + 1
                ae(r,:) = 0.0_dp
                ae(r,i) = 1.0_dp
                be(r) = up(i)
                se(r) = LP_LE
            end if
        end do

        allocate(cmax(n), x(n), y(me), rc(n))
        cmax = sign_dir * objective
        call simplex_relaxation(cmax, ae, se, be, control, x, objmax, st, it, y, rc)

        result%status = st
        result%simplex_iterations = it
        allocate(result%solution(n), result%duals(m), result%reduced_costs(n))
        result%solution = x
        result%objective = sign_dir * objmax
        if (m > 0) result%duals = sign_dir * y(1:m)
        result%reduced_costs = sign_dir * rc
        if (st == LP_OPTIMAL) result%solution_count = 1
    end subroutine solve_continuous


    subroutine solve_mip_single(direction, objective, a, sense, rhs, int_mask, &
        bin_mask, control, result)
        integer, intent(in) :: direction
        real(dp), intent(in) :: objective(:)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: sense(:)
        real(dp), intent(in) :: rhs(:)
        logical, intent(in) :: int_mask(:), bin_mask(:)
        type(lp_control), intent(in) :: control
        type(lp_result), intent(out) :: result

        integer :: n, global_status, nodes, total_iterations
        real(dp), allocatable :: lower(:), upper(:), best_x(:)
        real(dp) :: best_score, sign_dir, start_time
        logical :: have_best, stop_search
        type(lp_result) :: sens_result
        real(dp), allocatable :: fixlo(:), fixup(:)
        integer :: i

        call clear_result(result)
        n = size(objective)
        allocate(lower(n), upper(n), best_x(n))
        lower = 0.0_dp
        upper = LP_INFINITY
        where (bin_mask) upper = 1.0_dp
        best_x = 0.0_dp
        sign_dir = merge(1.0_dp, -1.0_dp, direction == LP_MAX)
        best_score = -huge(1.0_dp)
        have_best = .false.
        stop_search = .false.
        global_status = LP_OPTIMAL
        nodes = 0
        total_iterations = 0
        call cpu_time(start_time)

        call branch_node(lower, upper, 0)

        result%nodes = nodes
        result%simplex_iterations = total_iterations
        if (have_best) then
            allocate(result%solution(n))
            result%solution = best_x
            result%objective = sign_dir * best_score
            result%solution_count = 1
            if (global_status == LP_TIMEOUT .or. global_status == LP_SUBOPTIMAL) then
                result%status = global_status
            else
                result%status = LP_OPTIMAL
            end if

            ! Re-solve with integer variables fixed to provide useful LP duals/reduced costs.
            allocate(fixlo(n), fixup(n))
            fixlo = 0.0_dp
            fixup = LP_INFINITY
            where (bin_mask) fixup = 1.0_dp
            do i = 1, n
                if (int_mask(i)) then
                    fixlo(i) = nint(best_x(i))
                    fixup(i) = nint(best_x(i))
                end if
            end do
            call solve_continuous(direction, objective, a, sense, rhs, control, sens_result, &
                fixlo, fixup)
            if (allocated(sens_result%duals)) then
                allocate(result%duals(size(sens_result%duals)))
                result%duals = sens_result%duals
            end if
            if (allocated(sens_result%reduced_costs)) then
                allocate(result%reduced_costs(size(sens_result%reduced_costs)))
                result%reduced_costs = sens_result%reduced_costs
            end if
        else
            result%status = global_status
            if (result%status == LP_OPTIMAL) result%status = LP_INFEASIBLE
            allocate(result%solution(n))
            result%solution = 0.0_dp
        end if

    contains

        recursive subroutine branch_node(lo, up, depth)
            real(dp), intent(in) :: lo(:), up(:)
            integer, intent(in) :: depth
            type(lp_result) :: relax
            integer :: j, branch_var
            real(dp) :: frac, best_frac, score, now, fl, ce
            real(dp), allocatable :: lo2(:), up2(:)

            if (stop_search) return
            nodes = nodes + 1
            if (nodes > control%max_nodes) then
                global_status = LP_SUBOPTIMAL
                stop_search = .true.
                return
            end if
            if (control%timeout_seconds > 0.0_dp) then
                call cpu_time(now)
                if (now - start_time >= control%timeout_seconds) then
                    global_status = LP_TIMEOUT
                    stop_search = .true.
                    return
                end if
            end if

            call solve_continuous(direction, objective, a, sense, rhs, control, relax, lo, up)
            total_iterations = total_iterations + relax%simplex_iterations
            if (relax%status == LP_INFEASIBLE) return
            if (relax%status == LP_UNBOUNDED) then
                global_status = LP_UNBOUNDED
                stop_search = .true.
                return
            end if
            if (relax%status == LP_TIMEOUT) then
                global_status = LP_TIMEOUT
                stop_search = .true.
                return
            end if
            if (relax%status /= LP_OPTIMAL) then
                global_status = LP_SUBOPTIMAL
                return
            end if

            score = sign_dir * relax%objective
            if (have_best) then
                if (score <= best_score + control%optimality_tol) return
            end if

            branch_var = 0
            best_frac = control%integrality_tol
            do j = 1, n
                if (.not. int_mask(j)) cycle
                frac = abs(relax%solution(j) - nint(relax%solution(j)))
                if (frac > best_frac) then
                    best_frac = frac
                    branch_var = j
                end if
            end do

            if (branch_var == 0) then
                if (.not. have_best .or. score > best_score + control%optimality_tol) then
                    best_score = score
                    best_x = relax%solution
                    where (int_mask) best_x = real(nint(best_x),dp)
                    have_best = .true.
                end if
                return
            end if

            fl = floor(relax%solution(branch_var))
            ce = ceiling(relax%solution(branch_var))
            allocate(lo2(n), up2(n))

            ! Explore the closer side first; this often gets an incumbent quickly.
            if (relax%solution(branch_var) - fl <= ce - relax%solution(branch_var)) then
                lo2 = lo
                up2 = up
                up2(branch_var) = min(up2(branch_var), fl)
                if (up2(branch_var) >= lo2(branch_var) - control%feasibility_tol) then
                    call branch_node(lo2, up2, depth+1)
                end if
                if (stop_search) return
                lo2 = lo
                up2 = up
                lo2(branch_var) = max(lo2(branch_var), ce)
                if (up2(branch_var) >= lo2(branch_var) - control%feasibility_tol) then
                    call branch_node(lo2, up2, depth+1)
                end if
            else
                lo2 = lo
                up2 = up
                lo2(branch_var) = max(lo2(branch_var), ce)
                if (up2(branch_var) >= lo2(branch_var) - control%feasibility_tol) then
                    call branch_node(lo2, up2, depth+1)
                end if
                if (stop_search) return
                lo2 = lo
                up2 = up
                up2(branch_var) = min(up2(branch_var), fl)
                if (up2(branch_var) >= lo2(branch_var) - control%feasibility_tol) then
                    call branch_node(lo2, up2, depth+1)
                end if
            end if
        end subroutine branch_node

    end subroutine solve_mip_single


    subroutine append_exclusions(a, sense, rhs, previous, count, ao, so, bo)
        real(dp), intent(in) :: a(:,:), rhs(:), previous(:,:)
        integer, intent(in) :: sense(:), count
        real(dp), allocatable, intent(out) :: ao(:,:), bo(:)
        integer, allocatable, intent(out) :: so(:)
        integer :: m, n, k, j, ones

        m = size(rhs)
        n = size(a,2)
        allocate(ao(m+count,n), bo(m+count), so(m+count))
        if (m > 0) then
            ao(1:m,:) = a
            bo(1:m) = rhs
            so(1:m) = sense
        end if
        do k = 1, count
            ao(m+k,:) = 0.0_dp
            ones = 0
            do j = 1, n
                if (previous(j,k) > 0.5_dp) then
                    ao(m+k,j) = -1.0_dp
                    ones = ones + 1
                else
                    ao(m+k,j) = 1.0_dp
                end if
            end do
            bo(m+k) = 1.0_dp - real(ones,dp)
            so(m+k) = LP_GE
        end do
    end subroutine append_exclusions


    subroutine set_index_mask(indices, mask, status)
        integer, intent(in) :: indices(:)
        logical, intent(inout) :: mask(:)
        integer, intent(inout) :: status
        integer :: k

        status = LP_OPTIMAL
        do k = 1, size(indices)
            if (indices(k) < 1 .or. indices(k) > size(mask)) then
                status = LP_NUMFAILURE
                return
            end if
            mask(indices(k)) = .true.
        end do
    end subroutine set_index_mask


    subroutine clear_result(result)
        type(lp_result), intent(out) :: result
        result%status = LP_NUMFAILURE
        result%objective = 0.0_dp
        result%solution_count = 0
        result%simplex_iterations = 0
        result%nodes = 0
        result%sensitivity_ranges_available = .false.
    end subroutine clear_result

end module lpsolve_core
