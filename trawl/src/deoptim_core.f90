! Computational translation of DEoptim 2.2-8 src/de4_0.c.
! Original package license: GPL (>= 2). See LICENSE-GPL-2 and original/.
module deoptim_core
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use deoptim_kinds, only : dp, i8
    use deoptim_rng, only : de_rng
    use deoptim_types, only : de_control, de_result, de_objective, de_map, &
        de_success, de_invalid_input, de_unsupported, de_objective_nan, de_map_error
    implicit none
    private
    public :: deoptim_solve

contains

    subroutine deoptim_solve(objective, lower, upper, result, control, initialpop, map)
        procedure(de_objective) :: objective
        real(dp), intent(in) :: lower(:), upper(:)
        type(de_result), intent(out) :: result
        type(de_control), intent(in), optional :: control
        real(dp), intent(in), optional :: initialpop(:,:)
        procedure(de_map), optional :: map

        type(de_control) :: ctrl
        type(de_rng) :: rng
        integer :: d, np, i, j, k, iter, iter_tol, p_np, pbest
        integer :: r1, r2, r3, nstoremax, store_count
        integer :: urn(4)
        integer, allocatable :: work_urn(:), sort_index(:)
        real(dp), allocatable :: oldp(:,:), trialp(:,:)
        real(dp), allocatable :: oldc(:), trialc(:)
        real(dp), allocatable :: best(:), best_before(:), bestmem_tmp(:,:), bestval_tmp(:)
        real(dp), allocatable :: store_tmp(:,:,:), temp_cost(:)
        real(dp) :: bestval, f_now, cr_now, dither, jitter, cross_draw
        real(dp) :: mean_cr, mean_f, good_cr, good_f, good_f2
        integer :: good_np
        logical :: map_ok

        call init_empty_result(result)
        ctrl = de_control()
        if (present(control)) ctrl = control

        d = size(lower)
        if (d <= 0 .or. size(upper) /= d) then
            call fail(result, de_invalid_input, "lower and upper must have the same positive length")
            return
        end if
        if (any(lower > upper)) then
            call fail(result, de_invalid_input, "every lower bound must be <= its upper bound")
            return
        end if

        if (present(initialpop)) then
            if (ctrl%np <= 0 .or. ctrl%np /= size(initialpop,1)) ctrl%np = size(initialpop,1)
        end if
        call normalize_control(ctrl, d)
        np = ctrl%np
        if (np < 4) then
            call fail(result, de_invalid_input, "NP must be at least 4")
            return
        end if
        if (ctrl%bs) then
            call fail(result, de_unsupported, "bs=.true. is not supported by DEoptim 2.2-8's active C engine")
            return
        end if
        if (present(initialpop)) then
            if (size(initialpop,1) /= np .or. size(initialpop,2) /= d) then
                call fail(result, de_invalid_input, "initialpop must have shape (NP,npar)")
                return
            end if
        else
            if (any(.not. finite_interval(lower, upper))) then
                call fail(result, de_invalid_input, "finite bounds are required when initialpop is absent")
                return
            end if
        end if

        allocate(oldp(np,d), trialp(np,d))
        allocate(oldc(np), trialc(np))
        allocate(best(d), best_before(d), sort_index(np), temp_cost(np), work_urn(np))
        allocate(bestmem_tmp(d, ctrl%itermax), bestval_tmp(ctrl%itermax))
        bestmem_tmp = 0.0_dp
        bestval_tmp = 0.0_dp

        if (ctrl%storepopfrom <= ctrl%itermax) then
            nstoremax = 1 + max(0, (ctrl%itermax - ctrl%storepopfrom) / ctrl%storepopfreq)
        else
            nstoremax = 0
        end if
        allocate(store_tmp(np,d,nstoremax))
        if (nstoremax > 0) store_tmp = 0.0_dp

        call rng%seed(ctrl%seed)
        if (present(initialpop)) then
            oldp = initialpop
        else
            do j = 1, d
                do i = 1, np
                    oldp(i,j) = lower(j) + rng%uniform() * (upper(j) - lower(j))
                end do
            end do
        end if

        if (present(map)) then
            call map_population(oldp, lower, upper, map, rng, result%map_duplicates_remaining, map_ok)
            if (.not. map_ok) then
                call fail(result, de_map_error, "mapping callback produced NaN values")
                return
            end if
        end if
        call evaluate_population(oldp, oldc, objective, result%nfeval, map_ok)
        if (.not. map_ok) then
            call fail(result, de_objective_nan, "objective callback returned NaN")
            return
        end if

        bestval = huge(1.0_dp)
        best = oldp(1,:)
        do i = 1, np
            if (oldc(i) <= bestval) then
                bestval = oldc(i)
                best = oldp(i,:)
            end if
        end do

        p_np = nint(ctrl%p * real(np,dp))
        p_np = max(2, min(np, p_np))
        cr_now = ctrl%cr
        f_now = ctrl%f
        mean_cr = ctrl%cr
        mean_f = ctrl%f
        good_cr = 0.0_dp
        good_f = 0.0_dp
        good_f2 = 0.0_dp
        good_np = 0
        iter = 0
        iter_tol = 0
        store_count = 0

        do while (iter < ctrl%itermax .and. bestval > ctrl%vtr .and. iter_tol <= ctrl%steptol)
            ! The C code stores the current population before creating generation iter+1.
            if (mod(iter, ctrl%storepopfreq) == 0 .and. iter >= ctrl%storepopfrom - 1) then
                if (store_count < nstoremax) then
                    store_count = store_count + 1
                    store_tmp(:,:,store_count) = oldp
                end if
            end if

            iter = iter + 1
            bestmem_tmp(:,iter) = best
            bestval_tmp(iter) = bestval
            best_before = best

            if (ctrl%strategy == 5) dither = f_now + rng%uniform() * (1.0_dp - f_now)

            if (ctrl%strategy == 6) then
                temp_cost = oldc
                call sort_with_index(temp_cost, sort_index)
            end if

            do i = 1, np
                trialp(i,:) = oldp(i,:)
                call permute_urn(rng, urn, np, i, work_urn)
                r1 = urn(2)
                r2 = urn(3)
                r3 = urn(4)

                if (ctrl%c <= 0.0_dp) then
                    cr_now = ctrl%cr
                    f_now = ctrl%f
                else
                    cr_now = min(1.0_dp, max(0.0_dp, rng%normal(mean_cr, 0.1_dp)))
                    do
                        f_now = min(1.0_dp, rng%cauchy(mean_f, 0.1_dp))
                        if (f_now > 0.0_dp) exit
                    end do
                end if

                j = rng%randint(d)
                k = 0
                do
                    select case (ctrl%strategy)
                    case (1)
                        trialp(i,j) = oldp(r1,j) + f_now * (oldp(r2,j) - oldp(r3,j))
                    case (2)
                        trialp(i,j) = trialp(i,j) + f_now * (best_before(j) - trialp(i,j)) + &
                            f_now * (oldp(r2,j) - oldp(r3,j))
                    case (3)
                        jitter = 0.0001_dp * rng%uniform() + f_now
                        trialp(i,j) = best_before(j) + jitter * (oldp(r1,j) - oldp(r2,j))
                    case (4)
                        trialp(i,j) = oldp(r1,j) + (f_now + rng%uniform() * (1.0_dp - f_now)) * &
                            (oldp(r2,j) - oldp(r3,j))
                    case (5)
                        trialp(i,j) = oldp(r1,j) + dither * (oldp(r2,j) - oldp(r3,j))
                    case (6)
                        pbest = sort_index(rng%randint(p_np))
                        trialp(i,j) = oldp(i,j) + f_now * (oldp(pbest,j) - oldp(i,j)) + &
                            f_now * (oldp(r1,j) - oldp(r2,j))
                    end select
                    j = modulo(j, d) + 1
                    k = k + 1
                    cross_draw = rng%uniform()
                    if (.not. (cross_draw < cr_now .and. k < d)) exit
                end do

                do j = 1, d
                    if (trialp(i,j) < lower(j)) then
                        trialp(i,j) = lower(j) + rng%uniform() * (upper(j) - lower(j))
                    end if
                    if (trialp(i,j) > upper(j)) then
                        trialp(i,j) = upper(j) - rng%uniform() * (upper(j) - lower(j))
                    end if
                end do

            end do

            if (present(map)) then
                call map_population(trialp, lower, upper, map, rng, result%map_duplicates_remaining, map_ok)
                if (.not. map_ok) then
                    call fail(result, de_map_error, "mapping callback produced NaN values")
                    return
                end if
            end if

            call select_generation(trialp, oldp, oldc, objective, trialc, best, bestval, &
                result%nfeval, ctrl%c, cr_now, f_now, mean_cr, mean_f, &
                good_cr, good_f, good_f2, good_np, map_ok)
            if (.not. map_ok) then
                call fail(result, de_objective_nan, "objective callback returned NaN")
                return
            end if
            oldp = trialp
            oldc = trialc

            if (ctrl%trace > 0) then
                if (mod(iter, ctrl%trace) == 0) call print_progress(iter, bestval, best)
            end if

            if (bestval_tmp(iter) - bestval < ctrl%reltol * (abs(bestval_tmp(iter)) + ctrl%reltol)) then
                iter_tol = iter_tol + 1
            else
                iter_tol = 0
            end if
        end do

        result%bestval = bestval
        result%iter = iter
        allocate(result%bestmem(d)); result%bestmem = best
        allocate(result%pop(np,d)); result%pop = oldp
        allocate(result%bestmemit(d,iter)); result%bestmemit = bestmem_tmp(:,1:iter)
        allocate(result%bestvalit(iter)); result%bestvalit = bestval_tmp(1:iter)
        result%nstore = store_count
        allocate(result%storepop(np,d,store_count))
        if (store_count > 0) result%storepop = store_tmp(:,:,1:store_count)
        result%status = de_success
        result%message = "success"
    end subroutine deoptim_solve

    subroutine select_generation(trialp, oldp, oldc, objective, outc, best, bestval, &
            nfeval, adapt_c, cr_used, f_used, mean_cr, mean_f, good_cr, good_f, good_f2, good_np, ok)
        real(dp), intent(inout) :: trialp(:,:)
        real(dp), intent(in) :: oldp(:,:), oldc(:)
        procedure(de_objective) :: objective
        real(dp), intent(inout) :: outc(:)
        real(dp), intent(inout) :: best(:), bestval
        integer(i8), intent(inout) :: nfeval
        real(dp), intent(in) :: adapt_c, cr_used, f_used
        real(dp), intent(inout) :: mean_cr, mean_f, good_cr, good_f, good_f2
        integer, intent(inout) :: good_np
        logical, intent(out) :: ok
        real(dp), allocatable :: values(:)
        integer :: i

        allocate(values(size(oldc)))
        call evaluate_population(trialp, values, objective, nfeval, ok)
        if (.not. ok) return

        do i = 1, size(oldc)
            if (values(i) <= oldc(i)) then
                outc(i) = values(i)
                if (values(i) <= bestval) then
                    best = trialp(i,:)
                    bestval = values(i)
                end if
                if (adapt_c > 0.0_dp) then
                    good_np = good_np + 1
                    good_cr = good_cr + cr_used / real(good_np,dp)
                    good_f = good_f + f_used
                    good_f2 = good_f2 + f_used * f_used
                end if
            else
                trialp(i,:) = oldp(i,:)
                outc(i) = oldc(i)
            end if
        end do

        if (adapt_c > 0.0_dp .and. good_f > 0.0_dp) then
            mean_cr = (1.0_dp - adapt_c) * mean_cr + adapt_c * good_cr
            mean_f = (1.0_dp - adapt_c) * mean_f + adapt_c * good_f2 / good_f
        end if
    end subroutine select_generation

    subroutine evaluate_population(pop, values, objective, nfeval, ok)
        real(dp), intent(in) :: pop(:,:)
        real(dp), intent(out) :: values(:)
        procedure(de_objective) :: objective
        integer(i8), intent(inout) :: nfeval
        logical, intent(out) :: ok
        integer :: i

        ok = .true.
        do i = 1, size(pop,1)
            values(i) = objective(pop(i,:))
            nfeval = nfeval + 1_i8
            if (ieee_is_nan(values(i))) then
                ok = .false.
                return
            end if
        end do
    end subroutine evaluate_population

    subroutine map_population(pop, lower, upper, map, rng, duplicates_remaining, ok)
        real(dp), intent(inout) :: pop(:,:)
        real(dp), intent(in) :: lower(:), upper(:)
        procedure(de_map) :: map
        type(de_rng), intent(inout) :: rng
        integer, intent(out) :: duplicates_remaining
        logical, intent(out) :: ok
        logical, allocatable :: dup(:)
        integer :: i, j, tries

        ok = .true.
        do i = 1, size(pop,1)
            call map(pop(i,:))
            if (any(ieee_is_nan(pop(i,:)))) then
                ok = .false.
                return
            end if
        end do

        allocate(dup(size(pop,1)))
        call mark_duplicates(pop, dup)
        tries = 0
        do while (tries < 5 .and. any(dup))
            do i = 1, size(pop,1)
                if (dup(i)) then
                    do j = 1, size(pop,2)
                        pop(i,j) = lower(j) + rng%uniform() * (upper(j) - lower(j))
                    end do
                    call map(pop(i,:))
                    if (any(ieee_is_nan(pop(i,:)))) then
                        ok = .false.
                        return
                    end if
                end if
            end do
            call mark_duplicates(pop, dup)
            tries = tries + 1
        end do
        duplicates_remaining = count(dup)
    end subroutine map_population

    subroutine mark_duplicates(pop, dup)
        real(dp), intent(in) :: pop(:,:)
        logical, intent(out) :: dup(:)
        integer :: i, j

        dup = .false.
        do i = 2, size(pop,1)
            do j = 1, i - 1
                if (all(abs(pop(i,:) - pop(j,:)) <= 0.0_dp)) then
                    dup(i) = .true.
                    exit
                end if
            end do
        end do
    end subroutine mark_duplicates

    subroutine permute_urn(rng, urn2, np, avoid, urn1)
        type(de_rng), intent(inout) :: rng
        integer, intent(out) :: urn2(4)
        integer, intent(in) :: np, avoid
        integer, intent(inout) :: urn1(np)
        integer :: k, pos, idx, i

        do i = 1, np
            urn1(i) = i
        end do
        k = np
        pos = 1
        idx = avoid
        do while (k > np - 4)
            urn2(pos) = urn1(idx)
            urn1(idx) = urn1(k)
            k = k - 1
            pos = pos + 1
            if (pos <= 4) idx = rng%randint(k)
        end do
    end subroutine permute_urn

    subroutine sort_with_index(values, index)
        real(dp), intent(inout) :: values(:)
        integer, intent(out) :: index(:)
        integer :: i, j, idx
        real(dp) :: key

        do i = 1, size(values)
            index(i) = i
        end do
        do i = 2, size(values)
            key = values(i)
            idx = index(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j+1) = values(j)
                index(j+1) = index(j)
                j = j - 1
            end do
            values(j+1) = key
            index(j+1) = idx
        end do
    end subroutine sort_with_index

    elemental logical function finite_interval(lo, hi)
        real(dp), intent(in) :: lo, hi
        finite_interval = abs(lo) <= huge(lo) .and. abs(hi) <= huge(hi) .and. &
            .not. ieee_is_nan(lo) .and. .not. ieee_is_nan(hi)
    end function finite_interval

    subroutine normalize_control(ctrl, d)
        type(de_control), intent(inout) :: ctrl
        integer, intent(in) :: d

        if (ctrl%itermax <= 0) ctrl%itermax = 200
        if (ctrl%f < 0.0_dp .or. ctrl%f > 2.0_dp) ctrl%f = 0.8_dp
        if (ctrl%cr < 0.0_dp .or. ctrl%cr > 1.0_dp) ctrl%cr = 0.5_dp
        if (ctrl%strategy < 1 .or. ctrl%strategy > 6) ctrl%strategy = 2
        if (ctrl%np <= 0 .or. ctrl%np < 4) ctrl%np = 10 * d
        ctrl%storepopfreq = max(1, ctrl%storepopfreq)
        if (ctrl%storepopfreq > ctrl%itermax) ctrl%storepopfreq = 1
        if (ctrl%storepopfrom <= 0) ctrl%storepopfrom = ctrl%itermax + 1
        if (ctrl%p <= 0.0_dp .or. ctrl%p > 1.0_dp) ctrl%p = 0.2_dp
        if (ctrl%c < 0.0_dp .or. ctrl%c > 1.0_dp) ctrl%c = 0.0_dp
        if (ctrl%reltol < 0.0_dp) ctrl%reltol = sqrt(epsilon(1.0_dp))
        if (ctrl%steptol <= 0) ctrl%steptol = ctrl%itermax
        ctrl%trace = max(0, ctrl%trace)
    end subroutine normalize_control

    subroutine print_progress(iter, bestval, best)
        integer, intent(in) :: iter
        real(dp), intent(in) :: bestval, best(:)
        integer :: j
        write(*,'(a,i0,a,es16.8,a)',advance='no') "Iteration: ", iter, " bestvalit: ", bestval, " bestmemit:"
        do j = 1, size(best)
            write(*,'(1x,es14.6)',advance='no') best(j)
        end do
        write(*,*)
    end subroutine print_progress

    subroutine init_empty_result(result)
        type(de_result), intent(out) :: result
        result%bestval = huge(1.0_dp)
        result%nfeval = 0_i8
        result%iter = 0
        result%nstore = 0
        result%status = de_success
        result%message = ""
        result%map_duplicates_remaining = 0
    end subroutine init_empty_result

    subroutine fail(result, status, message)
        type(de_result), intent(inout) :: result
        integer, intent(in) :: status
        character(len=*), intent(in) :: message
        result%status = status
        result%message = message
    end subroutine fail

end module deoptim_core
