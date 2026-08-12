! genalg-fortran -- translation of rbga() and rbga.bin() from genalg 0.2.1.
! License: GPL-2.0-only. See LICENSE.
module genalg_core
    use iso_fortran_env, only : int64
    use genalg_kinds, only : dp
    use genalg_rng, only : rng_state
    implicit none
    private

    type, public :: rbga_control
        integer :: pop_size = 200
        integer :: iters = 100
        real(dp) :: mutation_chance = -1.0_dp
        integer :: elitism = -1
        integer(int64) :: seed = 123456789_int64
        logical :: verbose = .false.
    end type rbga_control

    type, public :: rbga_bin_control
        integer :: pop_size = 200
        integer :: iters = 100
        real(dp) :: mutation_chance = -1.0_dp
        integer :: elitism = -1
        real(dp) :: zero_to_one_ratio = 10.0_dp
        integer(int64) :: seed = 123456789_int64
        logical :: verbose = .false.
        ! genalg 0.2.1 does not invalidate cached fitness after binary mutation.
        logical :: legacy_binary_eval_cache = .true.
    end type rbga_bin_control

    type, public :: rbga_result
        real(dp), allocatable :: population(:,:)
        real(dp), allocatable :: evaluations(:)
        real(dp), allocatable :: best(:)
        real(dp), allocatable :: mean(:)
        real(dp), allocatable :: best_chromosome(:)
        real(dp) :: best_value = huge(1.0_dp)
        integer :: nfe = 0
        integer :: mutation_count = 0
    end type rbga_result

    type, public :: rbga_bin_result
        integer, allocatable :: population(:,:)
        real(dp), allocatable :: evaluations(:)
        real(dp), allocatable :: best(:)
        real(dp), allocatable :: mean(:)
        integer, allocatable :: best_chromosome(:)
        real(dp) :: best_value = huge(1.0_dp)
        integer :: nfe = 0
        integer :: mutation_count = 0
    end type rbga_bin_result

    abstract interface
        function real_eval_func(x) result(value)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function real_eval_func

        function binary_eval_func(x) result(value)
            import dp
            integer, intent(in) :: x(:)
            real(dp) :: value
        end function binary_eval_func

        subroutine real_monitor_func(iter, population, evaluations, best, mean)
            import dp
            integer, intent(in) :: iter
            real(dp), intent(in) :: population(:,:), evaluations(:)
            real(dp), intent(in) :: best(:), mean(:)
        end subroutine real_monitor_func

        subroutine binary_monitor_func(iter, population, evaluations, best, mean)
            import dp
            integer, intent(in) :: iter
            integer, intent(in) :: population(:,:)
            real(dp), intent(in) :: evaluations(:)
            real(dp), intent(in) :: best(:), mean(:)
        end subroutine binary_monitor_func
    end interface

    public :: rbga, rbga_bin
    public :: real_eval_func, binary_eval_func
    public :: real_monitor_func, binary_monitor_func

contains

    subroutine rbga(string_min, string_max, eval_func, result, control, suggestions, monitor_func)
        real(dp), intent(in) :: string_min(:), string_max(:)
        procedure(real_eval_func) :: eval_func
        type(rbga_result), intent(out) :: result
        type(rbga_control), intent(in), optional :: control
        real(dp), intent(in), optional :: suggestions(:,:)
        procedure(real_monitor_func), optional :: monitor_func

        type(rbga_control) :: ctl
        type(rng_state) :: rng
        integer :: nvar, pop_size, iters, elitism
        real(dp) :: mutation_chance
        real(dp), allocatable :: population(:,:), new_population(:,:)
        real(dp), allocatable :: evals(:), new_evals(:), sorted_evals(:)
        real(dp), allocatable :: parent_weights(:)
        logical, allocatable :: valid(:), new_valid(:)
        integer, allocatable :: order(:), sample_ids(:)
        integer :: i, j, iter, child, p1, p2, cross, nsuggest, idx
        integer :: direction, mutation_count
        real(dp) :: damp, mutation_val, mutation

        ctl = rbga_control()
        if (present(control)) ctl = control
        nvar = size(string_min)
        if (size(string_max) /= nvar) error stop "rbga: string_min and string_max sizes differ"
        if (nvar < 1) error stop "rbga: chromosome size must be positive"
        if (any(string_max < string_min)) error stop "rbga: upper bound below lower bound"
        pop_size = ctl%pop_size
        iters = ctl%iters
        if (pop_size < 5) error stop "rbga: pop_size must be at least 5"
        if (iters < 1) error stop "rbga: iters must be at least 1"
        elitism = ctl%elitism
        if (elitism < 0) elitism = pop_size / 5
        if (elitism < 0 .or. elitism >= pop_size) error stop "rbga: invalid elitism"
        mutation_chance = ctl%mutation_chance
        if (mutation_chance < 0.0_dp) mutation_chance = 1.0_dp / real(nvar + 1, dp)
        if (mutation_chance < 0.0_dp .or. mutation_chance > 1.0_dp) then
            error stop "rbga: mutation_chance must be in [0,1]"
        end if

        call rng%seed(ctl%seed)
        allocate(population(pop_size,nvar), evals(pop_size), valid(pop_size))
        allocate(result%best(iters), result%mean(iters))
        evals = 0.0_dp
        valid = .false.

        nsuggest = 0
        if (present(suggestions)) then
            if (size(suggestions,2) /= nvar) error stop "rbga: suggestion width mismatch"
            nsuggest = size(suggestions,1)
            if (nsuggest > pop_size) error stop "rbga: too many suggestions"
            if (nsuggest > 0) population(1:nsuggest,:) = suggestions
        end if
        do i = nsuggest + 1, pop_size
            do j = 1, nvar
                population(i,j) = string_min(j) + rng%uniform() * (string_max(j) - string_min(j))
            end do
        end do

        allocate(order(pop_size), sample_ids(pop_size), sorted_evals(pop_size), parent_weights(pop_size))
        do i = 1, pop_size
            parent_weights(i) = exp(-0.5_dp * (real(i,dp) / (real(pop_size,dp)/3.0_dp))**2)
        end do
        mutation_count = 0
        result%nfe = 0

        do iter = 1, iters
            do i = 1, pop_size
                if (.not. valid(i)) then
                    evals(i) = eval_func(population(i,:))
                    valid(i) = .true.
                    result%nfe = result%nfe + 1
                end if
            end do
            result%best(iter) = minval(evals)
            result%mean(iter) = sum(evals) / real(pop_size, dp)

            if (present(monitor_func)) then
                call monitor_func(iter, population, evals, result%best(1:iter), result%mean(1:iter))
            end if
            if (ctl%verbose) then
                write(*,'(a,i0,a,es14.6)') "iteration ", iter, ": best = ", result%best(iter)
            end if
            if (iter == iters) exit

            call sort_index(evals, order)
            do i = 1, pop_size
                sorted_evals(i) = evals(order(i))
            end do
            allocate(new_population(pop_size,nvar), new_evals(pop_size), new_valid(pop_size))
            new_evals = 0.0_dp
            new_valid = .false.
            if (elitism > 0) then
                do i = 1, elitism
                    new_population(i,:) = population(order(i),:)
                    new_evals(i) = sorted_evals(i)
                    new_valid(i) = .true.
                end do
            end if

            if (nvar > 1) then
                do child = elitism + 1, pop_size
                    call weighted_pair_without_replacement(parent_weights, rng, p1, p2)
                    cross = rng%integer(0, nvar)
                    if (cross == 0) then
                        new_population(child,:) = population(order(p2),:)
                        new_evals(child) = sorted_evals(p2)
                        new_valid(child) = .true.
                    else if (cross == nvar) then
                        new_population(child,:) = population(order(p1),:)
                        new_evals(child) = sorted_evals(p1)
                        new_valid(child) = .true.
                    else
                        new_population(child,1:cross) = population(order(p1),1:cross)
                        new_population(child,cross+1:nvar) = population(order(p2),cross+1:nvar)
                    end if
                end do
            else
                call random_subset(pop_size, pop_size-elitism, rng, sample_ids)
                do child = elitism + 1, pop_size
                    idx = order(sample_ids(child-elitism))
                    new_population(child,:) = population(idx,:)
                end do
            end if

            call move_alloc(new_population, population)
            call move_alloc(new_evals, evals)
            call move_alloc(new_valid, valid)

            if (mutation_chance > 0.0_dp) then
                damp = real(iters - iter, dp) / real(iters, dp)
                do i = elitism + 1, pop_size
                    do j = 1, nvar
                        if (rng%uniform() < mutation_chance) then
                            if (rng%integer(0,1) == 0) then
                                direction = -1
                            else
                                direction = 1
                            end if
                            ! Preserve the exact genalg expression and operator precedence.
                            mutation_val = string_max(j) - string_min(j) * 0.67_dp
                            mutation = population(i,j) + real(direction,dp) * mutation_val * damp
                            if (mutation < string_min(j) .or. mutation > string_max(j)) then
                                mutation = string_min(j) + rng%uniform() * (string_max(j) - string_min(j))
                            end if
                            population(i,j) = mutation
                            valid(i) = .false.
                            mutation_count = mutation_count + 1
                        end if
                    end do
                end do
            end if
        end do

        allocate(result%population(pop_size,nvar), result%evaluations(pop_size))
        result%population = population
        result%evaluations = evals
        idx = minloc(evals, dim=1)
        allocate(result%best_chromosome(nvar))
        result%best_chromosome = population(idx,:)
        result%best_value = evals(idx)
        result%mutation_count = mutation_count
    end subroutine rbga

    subroutine rbga_bin(size_chromosome, eval_func, result, control, suggestions, monitor_func)
        integer, intent(in) :: size_chromosome
        procedure(binary_eval_func) :: eval_func
        type(rbga_bin_result), intent(out) :: result
        type(rbga_bin_control), intent(in), optional :: control
        integer, intent(in), optional :: suggestions(:,:)
        procedure(binary_monitor_func), optional :: monitor_func

        type(rbga_bin_control) :: ctl
        type(rng_state) :: rng
        integer :: nvar, pop_size, iters, elitism, nzero
        real(dp) :: mutation_chance
        integer, allocatable :: population(:,:), new_population(:,:)
        real(dp), allocatable :: evals(:), new_evals(:), sorted_evals(:), parent_weights(:)
        logical, allocatable :: valid(:), new_valid(:)
        integer, allocatable :: order(:), sample_ids(:)
        integer :: i, j, iter, child, p1, p2, cross, nsuggest, idx, mutation_count

        ctl = rbga_bin_control()
        if (present(control)) ctl = control
        nvar = size_chromosome
        if (nvar < 1) error stop "rbga_bin: size must be positive"
        pop_size = ctl%pop_size
        iters = ctl%iters
        if (pop_size < 5) error stop "rbga_bin: pop_size must be at least 5"
        if (iters < 1) error stop "rbga_bin: iters must be at least 1"
        elitism = ctl%elitism
        if (elitism < 0) elitism = pop_size / 5
        if (elitism < 0 .or. elitism >= pop_size) error stop "rbga_bin: invalid elitism"
        mutation_chance = ctl%mutation_chance
        if (mutation_chance < 0.0_dp) mutation_chance = 1.0_dp / real(nvar + 1, dp)
        if (mutation_chance < 0.0_dp .or. mutation_chance > 1.0_dp) then
            error stop "rbga_bin: mutation_chance must be in [0,1]"
        end if
        ! R's rep(0, times) coerces the count to an integer. This mirrors that intent.
        nzero = max(0, int(ctl%zero_to_one_ratio))

        call rng%seed(ctl%seed)
        allocate(population(pop_size,nvar), evals(pop_size), valid(pop_size))
        allocate(result%best(iters), result%mean(iters))
        evals = 0.0_dp
        valid = .false.

        nsuggest = 0
        if (present(suggestions)) then
            if (size(suggestions,2) /= nvar) error stop "rbga_bin: suggestion width mismatch"
            nsuggest = size(suggestions,1)
            if (nsuggest > pop_size) error stop "rbga_bin: too many suggestions"
            if (nsuggest > 0) population(1:nsuggest,:) = suggestions
        end if
        do i = nsuggest + 1, pop_size
            call random_nonzero_binary(population(i,:), nzero, rng)
        end do

        allocate(order(pop_size), sample_ids(pop_size), sorted_evals(pop_size), parent_weights(pop_size))
        do i = 1, pop_size
            parent_weights(i) = exp(-0.5_dp * (real(i,dp) / (real(pop_size,dp)/3.0_dp))**2)
        end do
        mutation_count = 0
        result%nfe = 0

        do iter = 1, iters
            do i = 1, pop_size
                if (.not. valid(i)) then
                    evals(i) = eval_func(population(i,:))
                    valid(i) = .true.
                    result%nfe = result%nfe + 1
                end if
            end do
            result%best(iter) = minval(evals)
            result%mean(iter) = sum(evals) / real(pop_size, dp)

            if (present(monitor_func)) then
                call monitor_func(iter, population, evals, result%best(1:iter), result%mean(1:iter))
            end if
            if (ctl%verbose) then
                write(*,'(a,i0,a,es14.6)') "iteration ", iter, ": best = ", result%best(iter)
            end if
            if (iter == iters) exit

            call sort_index(evals, order)
            do i = 1, pop_size
                sorted_evals(i) = evals(order(i))
            end do
            allocate(new_population(pop_size,nvar), new_evals(pop_size), new_valid(pop_size))
            new_evals = 0.0_dp
            new_valid = .false.
            if (elitism > 0) then
                do i = 1, elitism
                    new_population(i,:) = population(order(i),:)
                    new_evals(i) = sorted_evals(i)
                    new_valid(i) = .true.
                end do
            end if

            if (nvar > 1) then
                do child = elitism + 1, pop_size
                    call weighted_pair_without_replacement(parent_weights, rng, p1, p2)
                    cross = rng%integer(0, nvar)
                    if (cross == 0) then
                        new_population(child,:) = population(order(p2),:)
                        new_evals(child) = sorted_evals(p2)
                        new_valid(child) = .true.
                    else if (cross == nvar) then
                        new_population(child,:) = population(order(p1),:)
                        new_evals(child) = sorted_evals(p1)
                        new_valid(child) = .true.
                    else
                        new_population(child,1:cross) = population(order(p1),1:cross)
                        new_population(child,cross+1:nvar) = population(order(p2),cross+1:nvar)
                        if (sum(new_population(child,:)) == 0) then
                            call random_nonzero_binary(new_population(child,:), nzero, rng)
                        end if
                    end if
                end do
            else
                call random_subset(pop_size, pop_size-elitism, rng, sample_ids)
                do child = elitism + 1, pop_size
                    idx = order(sample_ids(child-elitism))
                    new_population(child,:) = population(idx,:)
                end do
            end if

            call move_alloc(new_population, population)
            call move_alloc(new_evals, evals)
            call move_alloc(new_valid, valid)

            if (mutation_chance > 0.0_dp) then
                do i = elitism + 1, pop_size
                    do j = 1, nvar
                        if (rng%uniform() < mutation_chance) then
                            population(i,j) = random_binary_gene(nzero, rng)
                            mutation_count = mutation_count + 1
                            if (.not. ctl%legacy_binary_eval_cache) valid(i) = .false.
                        end if
                    end do
                end do
            end if
        end do

        allocate(result%population(pop_size,nvar), result%evaluations(pop_size))
        result%population = population
        result%evaluations = evals
        idx = minloc(evals, dim=1)
        allocate(result%best_chromosome(nvar))
        result%best_chromosome = population(idx,:)
        result%best_value = evals(idx)
        result%mutation_count = mutation_count
    end subroutine rbga_bin

    subroutine sort_index(values, order)
        real(dp), intent(in) :: values(:)
        integer, intent(out) :: order(size(values))
        integer :: i, j, key

        do i = 1, size(values)
            order(i) = i
        end do
        do i = 2, size(values)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (values(order(j)) <= values(key)) exit
                order(j+1) = order(j)
                j = j - 1
            end do
            order(j+1) = key
        end do
    end subroutine sort_index

    subroutine weighted_pair_without_replacement(weights, rng, first, second)
        real(dp), intent(in) :: weights(:)
        type(rng_state), intent(inout) :: rng
        integer, intent(out) :: first, second
        real(dp) :: total, target, accum
        integer :: i

        total = sum(weights)
        target = rng%uniform() * total
        accum = 0.0_dp
        first = size(weights)
        do i = 1, size(weights)
            accum = accum + weights(i)
            if (target <= accum) then
                first = i
                exit
            end if
        end do

        total = sum(weights) - weights(first)
        target = rng%uniform() * total
        accum = 0.0_dp
        second = 0
        do i = 1, size(weights)
            if (i == first) cycle
            accum = accum + weights(i)
            if (target <= accum) then
                second = i
                exit
            end if
        end do
        if (second == 0) then
            do i = size(weights), 1, -1
                if (i /= first) then
                    second = i
                    exit
                end if
            end do
        end if
    end subroutine weighted_pair_without_replacement

    subroutine random_subset(n, k, rng, out)
        integer, intent(in) :: n, k
        type(rng_state), intent(inout) :: rng
        integer, intent(out) :: out(:)
        integer, allocatable :: pool(:)
        integer :: i, j, tmp

        if (k < 0 .or. k > n) error stop "random_subset: invalid k"
        if (size(out) < k) error stop "random_subset: output too small"
        allocate(pool(n))
        do i = 1, n
            pool(i) = i
        end do
        do i = 1, k
            j = rng%integer(i, n)
            tmp = pool(i)
            pool(i) = pool(j)
            pool(j) = tmp
            out(i) = pool(i)
        end do
    end subroutine random_subset

    function random_binary_gene(nzero, rng) result(bit)
        integer, intent(in) :: nzero
        type(rng_state), intent(inout) :: rng
        integer :: bit

        if (rng%integer(1, nzero + 1) == nzero + 1) then
            bit = 1
        else
            bit = 0
        end if
    end function random_binary_gene

    subroutine random_nonzero_binary(x, nzero, rng)
        integer, intent(out) :: x(:)
        integer, intent(in) :: nzero
        type(rng_state), intent(inout) :: rng
        integer :: j

        do
            do j = 1, size(x)
                x(j) = random_binary_gene(nzero, rng)
            end do
            if (sum(x) > 0) exit
        end do
    end subroutine random_nonzero_binary

end module genalg_core
