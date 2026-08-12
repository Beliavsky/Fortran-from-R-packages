! Modern Fortran translation of computational code from TSP 1.2.7.
! Original Copyright (C) Michael Hahsler and Kurt Hornik.
! SPDX-License-Identifier: GPL-3.0-only
! See LICENSE, COPYING, and UPSTREAM.md for provenance and licensing.

module tsp_heuristics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use tsp_kinds, only : dp
    use tsp_types, only : tsp_control, tsp_tour, tsp_identity, tsp_random, &
        tsp_nearest_insertion, tsp_farthest_insertion, tsp_cheapest_insertion, &
        tsp_arbitrary_insertion, tsp_nn, tsp_repetitive_nn, tsp_two_opt_method, &
        tsp_sa_method, sa_reversal, sa_swap, sa_mixed
    use tsp_core, only : tour_length, insertion_cost, replace_infinite, is_square_matrix
    implicit none
    private

    public :: random_permutation, nearest_neighbor, repetitive_nearest_neighbor
    public :: insertion_heuristic, arbitrary_insertion
    public :: two_opt, two_opt_symmetric, simulated_annealing, solve_tsp, method_name

contains

    subroutine random_permutation(n, order)
        integer, intent(in) :: n
        integer, allocatable, intent(out) :: order(:)
        integer :: i, j, tmp
        real(dp) :: u

        allocate(order(n))
        order = [(i, i=1,n)]
        do i = n, 2, -1
            call random_number(u)
            j = 1 + int(u * real(i, dp))
            if (j > i) j = i
            tmp = order(i)
            order(i) = order(j)
            order(j) = tmp
        end do
    end subroutine random_permutation

    subroutine nearest_neighbor(cost, order, start)
        real(dp), intent(in) :: cost(:,:)
        integer, allocatable, intent(out) :: order(:)
        integer, intent(in), optional :: start
        logical, allocatable :: placed(:)
        integer :: n, cur, next, i, k, first
        real(dp) :: best, u
        integer, allocatable :: ties(:)
        integer :: ntie

        n = size(cost,1)
        allocate(order(n), placed(n), ties(n))
        placed = .false.
        if (n == 0) return

        if (present(start)) then
            first = start
        else
            call random_number(u)
            first = 1 + int(u * real(n,dp))
            if (first > n) first = n
        end if
        if (first < 1 .or. first > n) first = 1

        cur = first
        order(1) = cur
        placed(cur) = .true.
        do i = 2, n
            best = huge(1.0_dp)
            ntie = 0
            do k = 1, n
                if (placed(k)) cycle
                if (ntie == 0 .or. cost(cur,k) < best) then
                    best = cost(cur,k)
                    ntie = 1
                    ties(1) = k
                else if (same_value(cost(cur,k), best)) then
                    ntie = ntie + 1
                    ties(ntie) = k
                end if
            end do
            if (ntie <= 1) then
                next = ties(1)
            else
                call random_number(u)
                k = 1 + int(u * real(ntie,dp))
                if (k > ntie) k = ntie
                next = ties(k)
            end if
            order(i) = next
            placed(next) = .true.
            cur = next
        end do
    end subroutine nearest_neighbor

    subroutine repetitive_nearest_neighbor(cost, order)
        real(dp), intent(in) :: cost(:,:)
        integer, allocatable, intent(out) :: order(:)
        integer, allocatable :: candidate(:)
        integer :: n, s
        real(dp) :: best_len, len

        n = size(cost,1)
        allocate(order(n))
        if (n == 0) return
        best_len = huge(1.0_dp)
        do s = 1, n
            call nearest_neighbor(cost, candidate, s)
            len = tour_length(cost, candidate)
            if (len < best_len) then
                best_len = len
                order = candidate
            end if
        end do
    end subroutine repetitive_nearest_neighbor

    subroutine insertion_heuristic(cost, kind, order, start)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: kind
        integer, allocatable, intent(out) :: order(:)
        integer, intent(in), optional :: start
        logical, allocatable :: placed(:)
        integer, allocatable :: partial(:), candidates(:), ties(:)
        real(dp), allocatable :: delta(:)
        integer :: n, m, first, k, pos, i, j, nc, ntie, winner
        real(dp) :: u, score, best, d

        n = size(cost,1)
        allocate(order(n))
        if (n == 0) return
        allocate(placed(n), partial(n), candidates(n), ties(n))
        placed = .false.

        if (present(start)) then
            first = start
        else
            call random_number(u)
            first = 1 + int(u * real(n,dp))
            if (first > n) first = n
        end if
        if (first < 1 .or. first > n) first = 1
        partial(1) = first
        placed(first) = .true.
        m = 1

        do while (m < n)
            nc = 0
            do i = 1, n
                if (.not. placed(i)) then
                    nc = nc + 1
                    candidates(nc) = i
                end if
            end do

            select case (kind)
            case (tsp_nearest_insertion, tsp_farthest_insertion)
                if (kind == tsp_nearest_insertion) then
                    best = huge(1.0_dp)
                else
                    best = -huge(1.0_dp)
                end if
                ntie = 0
                do i = 1, nc
                    k = candidates(i)
                    score = huge(1.0_dp)
                    do j = 1, m
                        d = min(cost(k,partial(j)), cost(partial(j),k))
                        score = min(score, d)
                    end do
                    if (ntie == 0 .or. &
                        (kind == tsp_nearest_insertion .and. score < best) .or. &
                        (kind == tsp_farthest_insertion .and. score > best)) then
                        best = score
                        ntie = 1
                        ties(1) = k
                    else if (same_value(score, best)) then
                        ntie = ntie + 1
                        ties(ntie) = k
                    end if
                end do
                winner = pick_one(ties(:ntie))
            case (tsp_cheapest_insertion)
                best = huge(1.0_dp)
                ntie = 0
                do i = 1, nc
                    k = candidates(i)
                    call insertion_cost(cost, partial(:m), k, delta)
                    score = minval(delta)
                    if (ntie == 0 .or. score < best) then
                        best = score
                        ntie = 1
                        ties(1) = k
                    else if (same_value(score, best)) then
                        ntie = ntie + 1
                        ties(ntie) = k
                    end if
                end do
                winner = pick_one(ties(:ntie))
            case (tsp_arbitrary_insertion)
                winner = pick_one(candidates(:nc))
            case default
                winner = candidates(1)
            end select

            placed(winner) = .true.
            if (m == 1) then
                m = 2
                partial(2) = winner
            else
                call insertion_cost(cost, partial(:m), winner, delta)
                best = minval(delta)
                ntie = 0
                do i = 1, m
                    if (same_value(delta(i), best)) then
                        ntie = ntie + 1
                        ties(ntie) = i
                    end if
                end do
                pos = pick_one(ties(:ntie))
                do i = m, pos + 1, -1
                    partial(i+1) = partial(i)
                end do
                partial(pos+1) = winner
                m = m + 1
            end if
        end do
        order = partial(:n)
    end subroutine insertion_heuristic

    subroutine arbitrary_insertion(cost, order)
        real(dp), intent(in) :: cost(:,:)
        integer, allocatable, intent(out) :: order(:)
        integer, allocatable :: rorder(:), partial(:)
        real(dp), allocatable :: shuffled(:,:), delta(:)
        integer :: n, i, j, pos, k

        n = size(cost,1)
        if (n <= 2) then
            call random_permutation(n, order)
            return
        end if
        call random_permutation(n, rorder)
        allocate(shuffled(n,n))
        do j = 1, n
            do i = 1, n
                shuffled(i,j) = cost(rorder(i), rorder(j))
            end do
        end do
        allocate(partial(n))
        partial(1) = 1
        partial(2) = 2
        do k = 3, n
            call insertion_cost(shuffled, partial(:k-1), k, delta)
            pos = minloc(delta, dim=1)
            do i = k - 1, pos + 1, -1
                partial(i+1) = partial(i)
            end do
            partial(pos+1) = k
        end do
        allocate(order(n))
        do i = 1, n
            order(i) = rorder(partial(i))
        end do
    end subroutine arbitrary_insertion

    pure logical function same_value(a, b) result(equal)
        real(dp), intent(in) :: a, b
        equal = (a <= b .and. a >= b)
    end function same_value

    integer function pick_one(values) result(v)
        integer, intent(in) :: values(:)
        real(dp) :: u
        integer :: k
        if (size(values) <= 1) then
            v = values(1)
            return
        end if
        call random_number(u)
        k = 1 + int(u * real(size(values),dp))
        if (k > size(values)) k = size(values)
        v = values(k)
    end function pick_one

    subroutine two_opt(cost, initial, order)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: initial(:)
        integer, allocatable, intent(out) :: order(:)
        integer :: i, j, n, swap1, swap2, k, tmp
        integer :: swaps
        real(dp) :: imp, imp_tmp, imp_best
        real(dp), parameter :: epsilon = 1.0e-7_dp

        n = size(initial)
        allocate(order, source=initial)
        if (n <= 2) return

        do
            swaps = 0
            swap1 = 0
            swap2 = 0
            imp_best = 0.0_dp
            do i = 2, n - 1
                imp = cost(order(i-1), order(i)) + cost(order(i), order(i+1))
                do j = i + 1, n - 1
                    imp = imp + cost(order(j), order(j+1))
                    imp = imp - cost(order(j), order(j-1))
                    imp_tmp = imp - cost(order(i-1), order(j)) - cost(order(i), order(j+1))
                    if (imp_tmp > epsilon) then
                        swaps = swaps + 1
                        if (imp_tmp > imp_best) then
                            imp_best = imp_tmp
                            swap1 = i
                            swap2 = j
                        end if
                    end if
                end do
                j = n
                imp_tmp = imp + cost(order(j), order(1)) - cost(order(j), order(j-1)) &
                    - cost(order(i-1), order(j)) - cost(order(i), order(1))
                if (imp_tmp > epsilon) then
                    swaps = swaps + 1
                    if (imp_tmp > imp_best) then
                        imp_best = imp_tmp
                        swap1 = i
                        swap2 = j
                    end if
                end if
            end do
            if (swaps == 0) exit
            do k = 0, (swap2 - swap1) / 2
                if (swap1 + k >= swap2 - k) exit
                tmp = order(swap1+k)
                order(swap1+k) = order(swap2-k)
                order(swap2-k) = tmp
            end do
        end do
    end subroutine two_opt

    subroutine two_opt_symmetric(cost, initial, order)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: initial(:)
        integer, allocatable, intent(out) :: order(:)
        integer :: i, j, n, swaps, swap1, swap2, k, tmp
        real(dp) :: e1, e2, e1_swap, e2_swap, improvement, current_improvement

        n = size(initial)
        allocate(order, source=initial)
        if (n <= 2) return
        do
            swaps = 0
            swap1 = 0
            swap2 = 0
            improvement = 0.0_dp
            do i = 1, n - 2
                e1 = cost(order(i),order(i+1))
                do j = i + 1, n - 1
                    e2 = cost(order(j),order(j+1))
                    e1_swap = cost(order(i),order(j))
                    e2_swap = cost(order(i+1),order(j+1))
                    current_improvement = (e1 + e2) - (e1_swap + e2_swap)
                    if (current_improvement > 0.0_dp) then
                        swaps = swaps + 1
                        if (current_improvement > improvement) then
                            improvement = current_improvement
                            swap1 = i + 1
                            swap2 = j
                        end if
                    end if
                end do
                e2 = cost(order(n),order(1))
                e1_swap = cost(order(i),order(n))
                e2_swap = cost(order(i+1),order(1))
                current_improvement = (e1 + e2) - (e1_swap + e2_swap)
                if (current_improvement > 0.0_dp) then
                    swaps = swaps + 1
                    if (current_improvement > improvement) then
                        improvement = current_improvement
                        swap1 = i + 1
                        swap2 = n
                    end if
                end if
            end do
            if (swaps == 0) exit
            do k = 0, (swap2 - swap1) / 2
                if (swap1 + k >= swap2 - k) exit
                tmp = order(swap1+k)
                order(swap1+k) = order(swap2-k)
                order(swap2-k) = tmp
            end do
        end do
    end subroutine two_opt_symmetric

    subroutine simulated_annealing(cost, initial, order, temp, tmax, maxit, move_kind)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: initial(:)
        integer, allocatable, intent(out) :: order(:)
        real(dp), intent(in), optional :: temp
        integer, intent(in), optional :: tmax, maxit, move_kind
        integer, allocatable :: current(:), proposal(:), best_order(:)
        integer :: i, ntmax, nmax, kind
        real(dp) :: initial_temp, temperature, cur_cost, new_cost, best_cost, u

        allocate(current, source=initial)
        allocate(best_order, source=initial)
        ntmax = 10
        nmax = 10000
        kind = sa_reversal
        if (present(tmax)) ntmax = max(1, tmax)
        if (present(maxit)) nmax = max(1, maxit)
        if (present(move_kind)) kind = move_kind
        initial_temp = tour_length(cost, current) / real(max(1,size(current)),dp)
        if (present(temp)) then
            if (temp > 0.0_dp) initial_temp = temp
        end if
        if (.not. ieee_is_finite(initial_temp) .or. initial_temp <= 0.0_dp) initial_temp = 1.0_dp

        cur_cost = tour_length(cost, current)
        best_cost = cur_cost
        do i = 1, nmax
            temperature = initial_temp / log(real(((i - 1) / ntmax) * ntmax, dp) + exp(1.0_dp))
            call sa_local_move(current, kind, proposal)
            new_cost = tour_length(cost, proposal)
            if (new_cost <= cur_cost) then
                current = proposal
                cur_cost = new_cost
            else if (temperature > 0.0_dp .and. ieee_is_finite(new_cost) .and. ieee_is_finite(cur_cost)) then
                call random_number(u)
                if (u < exp((cur_cost - new_cost) / temperature)) then
                    current = proposal
                    cur_cost = new_cost
                end if
            end if
            if (cur_cost < best_cost) then
                best_cost = cur_cost
                best_order = current
            end if
        end do
        allocate(order, source=best_order)
    end subroutine simulated_annealing

    subroutine sa_local_move(tour, kind, new_tour)
        integer, intent(in) :: tour(:)
        integer, intent(in) :: kind
        integer, allocatable, intent(out) :: new_tour(:)
        integer :: i, j, k, tmp, actual_kind, n
        real(dp) :: u

        n = size(tour)
        allocate(new_tour, source=tour)
        if (n < 2) return
        actual_kind = kind
        if (actual_kind == sa_mixed) then
            call random_number(u)
            if (u > 0.5_dp) then
                actual_kind = sa_reversal
            else
                actual_kind = sa_swap
            end if
        end if
        call random_number(u)
        i = 1 + int(u * real(n,dp))
        if (i > n) i = n
        do
            call random_number(u)
            j = 1 + int(u * real(n,dp))
            if (j > n) j = n
            if (j /= i) exit
        end do
        if (actual_kind == sa_swap) then
            tmp = new_tour(i)
            new_tour(i) = new_tour(j)
            new_tour(j) = tmp
        else
            if (i > j) then
                k = i
                i = j
                j = k
            end if
            do k = 0, (j - i) / 2
                if (i + k >= j - k) exit
                tmp = new_tour(i+k)
                new_tour(i+k) = new_tour(j-k)
                new_tour(j-k) = tmp
            end do
        end if
    end subroutine sa_local_move

    subroutine two_opt_with_restarts(cost, ctl, order)
        real(dp), intent(in) :: cost(:,:)
        type(tsp_control), intent(in) :: ctl
        integer, allocatable, intent(out) :: order(:)
        integer, allocatable :: initial(:), candidate(:), best(:)
        integer :: r, nr
        real(dp) :: len, best_len

        if (allocated(ctl%tour)) then
            call two_opt(cost, ctl%tour, order)
            return
        end if
        nr = max(1, ctl%two_opt_repetitions)
        best_len = huge(1.0_dp)
        do r = 1, nr
            call random_permutation(size(cost,1), initial)
            call two_opt(cost, initial, candidate)
            len = tour_length(cost, candidate)
            if (r == 1 .or. len < best_len) then
                best_len = len
                best = candidate
            end if
        end do
        allocate(order, source=best)
    end subroutine two_opt_with_restarts

    function solve_tsp(cost, method, control) result(solution)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in), optional :: method
        type(tsp_control), intent(in), optional :: control
        type(tsp_tour) :: solution
        type(tsp_control) :: ctl
        real(dp), allocatable :: work(:,:)
        integer, allocatable :: order(:), candidate(:), initial(:), refined(:)
        integer :: m, r, reps, ierr, start, i
        real(dp) :: best_len, len

        if (.not. is_square_matrix(cost) .or. any(ieee_is_nan(cost))) then
            allocate(solution%order(0))
            solution%length = huge(1.0_dp)
            solution%method = "invalid"
            return
        end if
        ctl = tsp_control()
        if (present(control)) ctl = control
        m = tsp_arbitrary_insertion
        if (present(method)) m = method
        if (.not. present(method)) ctl%two_opt = .true.
        call replace_infinite(cost, work, ierr=ierr)
        if (ierr /= 0) work = cost

        reps = max(1, ctl%rep)
        if (m == tsp_repetitive_nn) reps = 1
        best_len = huge(1.0_dp)
        do r = 1, reps
            start = ctl%start
            select case (m)
            case (tsp_identity)
                allocate(candidate(size(cost,1)))
                candidate = [(i, i=1,size(cost,1))]
            case (tsp_random)
                call random_permutation(size(cost,1), candidate)
            case (tsp_nearest_insertion, tsp_farthest_insertion, tsp_cheapest_insertion)
                if (start > 0) then
                    call insertion_heuristic(work, m, candidate, start)
                else
                    call insertion_heuristic(work, m, candidate)
                end if
            case (tsp_arbitrary_insertion)
                call arbitrary_insertion(work, candidate)
            case (tsp_nn)
                if (start > 0) then
                    call nearest_neighbor(work, candidate, start)
                else
                    call nearest_neighbor(work, candidate)
                end if
            case (tsp_repetitive_nn)
                call repetitive_nearest_neighbor(work, candidate)
            case (tsp_two_opt_method)
                call two_opt_with_restarts(work, ctl, candidate)
            case (tsp_sa_method)
                if (allocated(ctl%tour)) then
                    initial = ctl%tour
                else
                    call random_permutation(size(cost,1), initial)
                end if
                if (ctl%temp > 0.0_dp) then
                    call simulated_annealing(work, initial, candidate, ctl%temp, ctl%tmax, ctl%maxit, ctl%sa_move)
                else
                    call simulated_annealing(work, initial, candidate, tmax=ctl%tmax, maxit=ctl%maxit, move_kind=ctl%sa_move)
                end if
            case default
                call arbitrary_insertion(work, candidate)
            end select

            if (ctl%two_opt .and. m /= tsp_two_opt_method) then
                call two_opt(work, candidate, refined)
                candidate = refined
            end if
            len = tour_length(cost, candidate)
            if (r == 1 .or. len < best_len) then
                best_len = len
                order = candidate
            end if
        end do
        solution%order = order
        solution%length = tour_length(cost, order)
        solution%method = method_name(m)
        if (ctl%two_opt .and. m /= tsp_two_opt_method) solution%method = trim(solution%method)//"+two_opt"
    end function solve_tsp

    pure function method_name(method) result(name)
        integer, intent(in) :: method
        character(len=48) :: name
        select case(method)
        case(tsp_identity); name = "identity"
        case(tsp_random); name = "random"
        case(tsp_nearest_insertion); name = "nearest_insertion"
        case(tsp_farthest_insertion); name = "farthest_insertion"
        case(tsp_cheapest_insertion); name = "cheapest_insertion"
        case(tsp_arbitrary_insertion); name = "arbitrary_insertion"
        case(tsp_nn); name = "nn"
        case(tsp_repetitive_nn); name = "repetitive_nn"
        case(tsp_two_opt_method); name = "two_opt"
        case(tsp_sa_method); name = "sa"
        case default; name = "unknown"
        end select
    end function method_name

end module tsp_heuristics
