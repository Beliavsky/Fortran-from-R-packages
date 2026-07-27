! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmc_pricing
    use lsmc_european, only : eu_put_bs
    use lsmc_kinds, only : dp
    use lsmc_linear_algebra, only : least_squares
    use lsmc_math, only : mean_value, sample_standard_error
    use lsmc_random, only : seed_random_number
    use lsmc_simulation, only : simulate_antithetic_correlated_gbm_paths
    use lsmc_simulation, only : simulate_antithetic_gbm_paths
    use lsmc_simulation, only : simulate_correlated_gbm_paths
    use lsmc_simulation, only : simulate_gbm_paths
    use lsmc_types, only : option_result
    implicit none
    private

    interface american_put_lsmc
        module procedure amer_put_lsm
    end interface american_put_lsmc

    interface american_put_lsmc_antithetic
        module procedure amer_put_lsm_av
    end interface american_put_lsmc_antithetic

    interface american_put_lsmc_control
        module procedure amer_put_lsm_cv
    end interface american_put_lsmc_control

    interface asian_american_put_lsmc
        module procedure asian_amer_put_lsm
    end interface asian_american_put_lsmc

    interface quanto_american_put_lsmc
        module procedure quanto_amer_put_lsm
    end interface quanto_american_put_lsmc

    interface quanto_american_put_lsmc_antithetic
        module procedure quanto_amer_put_lsm_av
    end interface quanto_american_put_lsmc_antithetic

    public :: amer_put_lsm
    public :: amer_put_lsm_av
    public :: amer_put_lsm_cv
    public :: asian_amer_put_lsm
    public :: quanto_amer_put_lsm
    public :: quanto_amer_put_lsm_av
    public :: american_put_lsmc
    public :: american_put_lsmc_antithetic
    public :: american_put_lsmc_control
    public :: asian_american_put_lsmc
    public :: quanto_american_put_lsmc
    public :: quanto_american_put_lsmc_antithetic

contains

    function amer_put_lsm(spot, sigma, n, m, strike, rate, dividend, maturity, seed) result(result)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: strike
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        integer, intent(in), optional :: seed
        type(option_result) :: result
        real(dp), allocatable :: paths(:, :)
        real(dp), allocatable :: present_values(:)
        real(dp) :: d
        real(dp) :: k
        real(dp) :: q
        real(dp) :: s
        real(dp) :: t
        real(dp) :: v
        integer :: np
        integer :: nt

        call resolve_basic_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, s, v, np, nt, k, d, q, t)
        if (present(seed)) call seed_random_number(seed)
        allocate(paths(np, nt), present_values(np))
        call simulate_gbm_paths(s, v, np, nt, d, q, t, paths)
        call plain_lsm_values(paths, k, d, t, present_values)
        call set_result(result, present_values, s, v, np, np, nt, k, d, q, t, &
            "American put", "Least-squares Monte Carlo")
    end function amer_put_lsm

    function amer_put_lsm_av(spot, sigma, n, m, strike, rate, dividend, maturity, seed) result(result)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: strike
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        integer, intent(in), optional :: seed
        type(option_result) :: result
        real(dp), allocatable :: pair_values(:)
        real(dp), allocatable :: paths(:, :)
        real(dp), allocatable :: present_values(:)
        real(dp) :: d
        real(dp) :: k
        real(dp) :: q
        real(dp) :: s
        real(dp) :: t
        real(dp) :: v
        integer :: np
        integer :: nt

        call resolve_basic_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, s, v, np, nt, k, d, q, t)
        if (present(seed)) call seed_random_number(seed)
        allocate(paths(2 * np, nt), present_values(2 * np), pair_values(np))
        call simulate_antithetic_gbm_paths(s, v, np, nt, d, q, t, paths)
        call plain_lsm_values(paths, k, d, t, present_values)
        pair_values = 0.5_dp * (present_values(1:np) + present_values(np + 1:2 * np))
        call set_result(result, pair_values, s, v, np, 2 * np, nt, k, d, q, t, &
            "American put", "LSMC with antithetic variates")
    end function amer_put_lsm_av

    function amer_put_lsm_cv(spot, sigma, n, m, strike, rate, dividend, maturity, seed) result(result)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: strike
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        integer, intent(in), optional :: seed
        type(option_result) :: result
        real(dp), allocatable :: adjusted_values(:)
        real(dp), allocatable :: paths(:, :)
        real(dp), allocatable :: present_values(:)
        real(dp) :: bs_put
        real(dp) :: d
        real(dp) :: k
        real(dp) :: q
        real(dp) :: s
        real(dp) :: t
        real(dp) :: v
        integer :: np
        integer :: nt

        call resolve_basic_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, s, v, np, nt, k, d, q, t)
        if (present(seed)) call seed_random_number(seed)
        allocate(paths(np, nt), present_values(np), adjusted_values(np))
        call simulate_gbm_paths(s, v, np, nt, d, q, t, paths)
        call plain_lsm_values(paths, k, d, t, present_values)
        bs_put = eu_put_bs(s, v, k, d, q, t)
        adjusted_values = present_values - exp(-d * t) * max(k - paths(:, nt), 0.0_dp) + bs_put
        call set_result(result, adjusted_values, s, v, np, np, nt, k, d, q, t, &
            "American put", "LSMC with European-put control variate")
    end function amer_put_lsm_cv

    function asian_amer_put_lsm(spot, sigma, n, m, strike, rate, dividend, maturity, seed) result(result)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: strike
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        integer, intent(in), optional :: seed
        type(option_result) :: result
        real(dp), allocatable :: paths(:, :)
        real(dp), allocatable :: present_values(:)
        real(dp) :: d
        real(dp) :: k
        real(dp) :: q
        real(dp) :: s
        real(dp) :: t
        real(dp) :: v
        integer :: np
        integer :: nt

        call resolve_basic_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, s, v, np, nt, k, d, q, t)
        if (present(seed)) call seed_random_number(seed)
        allocate(paths(np, nt), present_values(np))
        call simulate_gbm_paths(s, v, np, nt, d, q, t, paths)
        call asian_lsm_values(paths, k, d, t, present_values)
        call set_result(result, present_values, s, v, np, np, nt, k, d, q, t, &
            "Asian American put", "Least-squares Monte Carlo")
    end function asian_amer_put_lsm

    function quanto_amer_put_lsm(spot, sigma, n, m, strike, rate, dividend, maturity, spot2, sigma2, rate2, &
            dividend2, rho, seed) result(result)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: strike
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        real(dp), intent(in), optional :: spot2
        real(dp), intent(in), optional :: sigma2
        real(dp), intent(in), optional :: rate2
        real(dp), intent(in), optional :: dividend2
        real(dp), intent(in), optional :: rho
        integer, intent(in), optional :: seed
        type(option_result) :: result
        real(dp), allocatable :: paths1(:, :)
        real(dp), allocatable :: paths2(:, :)
        real(dp), allocatable :: present_values(:)
        real(dp) :: d
        real(dp) :: d2
        real(dp) :: k
        real(dp) :: q
        real(dp) :: q2
        real(dp) :: correlation
        real(dp) :: s
        real(dp) :: s2
        real(dp) :: t
        real(dp) :: v
        real(dp) :: v2
        integer :: np
        integer :: nt

        call resolve_quanto_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, spot2, sigma2, rate2, &
            dividend2, rho, s, v, np, nt, k, d, q, t, s2, v2, d2, q2, correlation)
        if (present(seed)) call seed_random_number(seed)
        allocate(paths1(np, nt), paths2(np, nt), present_values(np))
        call simulate_correlated_gbm_paths(s, v, d, q, s2, v2, d2, q2, correlation, np, nt, t, paths1, paths2)
        call quanto_lsm_values(paths1, paths2, k, d, t, present_values)
        call set_result(result, present_values, s, v, np, np, nt, k, d, q, t, &
            "Quanto American put", "Least-squares Monte Carlo")
        call set_quanto_result(result, s2, v2, d2, q2, correlation)
    end function quanto_amer_put_lsm

    function quanto_amer_put_lsm_av(spot, sigma, n, m, strike, rate, dividend, maturity, spot2, sigma2, rate2, &
            dividend2, rho, seed) result(result)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: strike
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        real(dp), intent(in), optional :: spot2
        real(dp), intent(in), optional :: sigma2
        real(dp), intent(in), optional :: rate2
        real(dp), intent(in), optional :: dividend2
        real(dp), intent(in), optional :: rho
        integer, intent(in), optional :: seed
        type(option_result) :: result
        real(dp), allocatable :: pair_values(:)
        real(dp), allocatable :: paths1(:, :)
        real(dp), allocatable :: paths2(:, :)
        real(dp), allocatable :: present_values(:)
        real(dp) :: d
        real(dp) :: d2
        real(dp) :: k
        real(dp) :: q
        real(dp) :: q2
        real(dp) :: correlation
        real(dp) :: s
        real(dp) :: s2
        real(dp) :: t
        real(dp) :: v
        real(dp) :: v2
        integer :: np
        integer :: nt

        call resolve_quanto_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, spot2, sigma2, rate2, &
            dividend2, rho, s, v, np, nt, k, d, q, t, s2, v2, d2, q2, correlation)
        if (present(seed)) call seed_random_number(seed)
        allocate(paths1(2 * np, nt), paths2(2 * np, nt), present_values(2 * np), pair_values(np))
        call simulate_antithetic_correlated_gbm_paths(s, v, d, q, s2, v2, d2, q2, correlation, np, nt, t, paths1, paths2)
        call quanto_lsm_values(paths1, paths2, k, d, t, present_values)
        pair_values = 0.5_dp * (present_values(1:np) + present_values(np + 1:2 * np))
        call set_result(result, pair_values, s, v, np, 2 * np, nt, k, d, q, t, &
            "Quanto American put", "LSMC with antithetic variates")
        call set_quanto_result(result, s2, v2, d2, q2, correlation)
    end function quanto_amer_put_lsm_av

    subroutine plain_lsm_values(paths, strike, rate, maturity, present_values)
        real(dp), intent(in) :: paths(:, :)
        real(dp), intent(in) :: strike
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: maturity
        real(dp), intent(out) :: present_values(:)
        real(dp), allocatable :: basis(:, :)
        real(dp), allocatable :: continuation_targets(:)
        real(dp), allocatable :: exercise_values(:)
        real(dp), allocatable :: cashflows(:)
        real(dp), allocatable :: continuation(:)
        real(dp) :: beta(3)
        real(dp) :: dt
        integer, allocatable :: exercise_step(:)
        integer, allocatable :: indices(:)
        integer :: i
        integer :: j
        integer :: n_in
        integer :: n_paths
        integer :: n_steps
        logical :: ok

        n_paths = size(paths, 1)
        n_steps = size(paths, 2)
        if (size(present_values) /= n_paths) error stop "plain_lsm_values: size mismatch"
        dt = maturity / real(n_steps, dp)
        allocate(cashflows(n_paths), exercise_step(n_paths), exercise_values(n_paths), indices(n_paths))
        cashflows = max(strike - paths(:, n_steps), 0.0_dp)
        exercise_step = n_steps

        do i = n_steps - 1, 1, -1
            exercise_values = max(strike - paths(:, i), 0.0_dp)
            n_in = pack_indices(exercise_values > 0.0_dp, indices)
            if (n_in < 3) cycle
            allocate(basis(n_in, 3), continuation_targets(n_in), continuation(n_in))
            do j = 1, n_in
                basis(j, :) = [1.0_dp, paths(indices(j), i), paths(indices(j), i)**2]
                continuation_targets(j) = cashflows(indices(j)) * &
                    exp(-rate * dt * real(exercise_step(indices(j)) - i, dp))
            end do
            call least_squares(basis, continuation_targets, beta, ok)
            if (ok) then
                continuation = matmul(basis, beta)
                do j = 1, n_in
                    if (exercise_values(indices(j)) > continuation(j)) then
                        cashflows(indices(j)) = exercise_values(indices(j))
                        exercise_step(indices(j)) = i
                    end if
                end do
            end if
            deallocate(basis, continuation_targets, continuation)
        end do
        present_values = cashflows * exp(-rate * dt * real(exercise_step, dp))
    end subroutine plain_lsm_values

    subroutine asian_lsm_values(paths, strike, rate, maturity, present_values)
        real(dp), intent(in) :: paths(:, :)
        real(dp), intent(in) :: strike
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: maturity
        real(dp), intent(out) :: present_values(:)
        real(dp), allocatable :: averages(:, :)
        real(dp), allocatable :: basis(:, :)
        real(dp), allocatable :: continuation_targets(:)
        real(dp), allocatable :: exercise_values(:)
        real(dp), allocatable :: cashflows(:)
        real(dp), allocatable :: continuation(:)
        real(dp) :: beta(3)
        real(dp) :: dt
        integer, allocatable :: exercise_step(:)
        integer, allocatable :: indices(:)
        integer :: i
        integer :: j
        integer :: n_in
        integer :: n_paths
        integer :: n_steps
        logical :: ok

        n_paths = size(paths, 1)
        n_steps = size(paths, 2)
        if (size(present_values) /= n_paths) error stop "asian_lsm_values: size mismatch"
        dt = maturity / real(n_steps, dp)
        allocate(averages(n_paths, n_steps), cashflows(n_paths), exercise_step(n_paths), &
            exercise_values(n_paths), indices(n_paths))
        averages(:, 1) = paths(:, 1)
        do i = 2, n_steps
            averages(:, i) = (real(i - 1, dp) * averages(:, i - 1) + paths(:, i)) / real(i, dp)
        end do
        cashflows = max(strike - averages(:, n_steps), 0.0_dp)
        exercise_step = n_steps

        do i = n_steps - 1, 1, -1
            exercise_values = max(strike - averages(:, i), 0.0_dp)
            n_in = pack_indices(exercise_values > 0.0_dp, indices)
            if (n_in < 3) cycle
            allocate(basis(n_in, 3), continuation_targets(n_in), continuation(n_in))
            do j = 1, n_in
                basis(j, :) = [1.0_dp, paths(indices(j), i), paths(indices(j), i)**2]
                continuation_targets(j) = cashflows(indices(j)) * &
                    exp(-rate * dt * real(exercise_step(indices(j)) - i, dp))
            end do
            call least_squares(basis, continuation_targets, beta, ok)
            if (ok) then
                continuation = matmul(basis, beta)
                do j = 1, n_in
                    if (exercise_values(indices(j)) > continuation(j)) then
                        cashflows(indices(j)) = exercise_values(indices(j))
                        exercise_step(indices(j)) = i
                    end if
                end do
            end if
            deallocate(basis, continuation_targets, continuation)
        end do
        present_values = cashflows * exp(-rate * dt * real(exercise_step, dp))
    end subroutine asian_lsm_values

    subroutine quanto_lsm_values(paths1, paths2, strike, rate, maturity, present_values)
        real(dp), intent(in) :: paths1(:, :)
        real(dp), intent(in) :: paths2(:, :)
        real(dp), intent(in) :: strike
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: maturity
        real(dp), intent(out) :: present_values(:)
        real(dp), allocatable :: basis(:, :)
        real(dp), allocatable :: continuation_targets(:)
        real(dp), allocatable :: exercise_values(:)
        real(dp), allocatable :: cashflows(:)
        real(dp), allocatable :: continuation(:)
        real(dp) :: beta(6)
        real(dp) :: dt
        integer, allocatable :: exercise_step(:)
        integer, allocatable :: indices(:)
        integer :: i
        integer :: j
        integer :: n_in
        integer :: n_paths
        integer :: n_steps
        logical :: ok

        if (any(shape(paths1) /= shape(paths2))) error stop "quanto_lsm_values: path shape mismatch"
        n_paths = size(paths1, 1)
        n_steps = size(paths1, 2)
        if (size(present_values) /= n_paths) error stop "quanto_lsm_values: size mismatch"
        dt = maturity / real(n_steps, dp)
        allocate(cashflows(n_paths), exercise_step(n_paths), exercise_values(n_paths), indices(n_paths))
        cashflows = max(strike - paths1(:, n_steps), 0.0_dp) * paths2(:, n_steps)
        exercise_step = n_steps

        do i = n_steps - 1, 1, -1
            exercise_values = max(strike - paths1(:, i), 0.0_dp) * paths2(:, i)
            n_in = pack_indices(paths1(:, i) < strike, indices)
            if (n_in < 6) cycle
            allocate(basis(n_in, 6), continuation_targets(n_in), continuation(n_in))
            do j = 1, n_in
                basis(j, :) = [1.0_dp, paths1(indices(j), i), paths1(indices(j), i)**2, &
                    paths2(indices(j), i), paths2(indices(j), i)**2, &
                    paths1(indices(j), i) * paths2(indices(j), i)]
                continuation_targets(j) = cashflows(indices(j)) * &
                    exp(-rate * dt * real(exercise_step(indices(j)) - i, dp))
            end do
            call least_squares(basis, continuation_targets, beta, ok)
            if (ok) then
                continuation = matmul(basis, beta)
                do j = 1, n_in
                    if (exercise_values(indices(j)) > continuation(j)) then
                        cashflows(indices(j)) = exercise_values(indices(j))
                        exercise_step(indices(j)) = i
                    end if
                end do
            end if
            deallocate(basis, continuation_targets, continuation)
        end do
        present_values = cashflows * exp(-rate * dt * real(exercise_step, dp))
    end subroutine quanto_lsm_values

    integer function pack_indices(mask, indices) result(count)
        logical, intent(in) :: mask(:)
        integer, intent(out) :: indices(:)
        integer :: i

        if (size(indices) < size(mask)) error stop "pack_indices: output too small"
        count = 0
        do i = 1, size(mask)
            if (mask(i)) then
                count = count + 1
                indices(count) = i
            end if
        end do
    end function pack_indices

    subroutine resolve_basic_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, s, v, np, nt, k, d, q, t)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: strike
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        real(dp), intent(out) :: s
        real(dp), intent(out) :: v
        integer, intent(out) :: np
        integer, intent(out) :: nt
        real(dp), intent(out) :: k
        real(dp), intent(out) :: d
        real(dp), intent(out) :: q
        real(dp), intent(out) :: t

        s = 1.0_dp
        v = 0.2_dp
        np = 1000
        nt = 365
        k = 1.1_dp
        d = 0.06_dp
        q = 0.0_dp
        t = 1.0_dp
        if (present(spot)) s = spot
        if (present(sigma)) v = sigma
        if (present(n)) np = n
        if (present(m)) nt = m
        if (present(strike)) k = strike
        if (present(rate)) d = rate
        if (present(dividend)) q = dividend
        if (present(maturity)) t = maturity
        if (s <= 0.0_dp .or. k <= 0.0_dp) error stop "LSMC: spot and strike must be positive"
        if (v < 0.0_dp) error stop "LSMC: sigma must be nonnegative"
        if (np <= 0 .or. nt <= 1) error stop "LSMC: n must be positive and m must exceed one"
        if (t <= 0.0_dp) error stop "LSMC: maturity must be positive"
    end subroutine resolve_basic_inputs

    subroutine resolve_quanto_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, spot2, sigma2, rate2, &
            dividend2, rho, s, v, np, nt, k, d, q, t, s2, v2, d2, q2, correlation)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: strike
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        real(dp), intent(in), optional :: spot2
        real(dp), intent(in), optional :: sigma2
        real(dp), intent(in), optional :: rate2
        real(dp), intent(in), optional :: dividend2
        real(dp), intent(in), optional :: rho
        real(dp), intent(out) :: s
        real(dp), intent(out) :: v
        integer, intent(out) :: np
        integer, intent(out) :: nt
        real(dp), intent(out) :: k
        real(dp), intent(out) :: d
        real(dp), intent(out) :: q
        real(dp), intent(out) :: t
        real(dp), intent(out) :: s2
        real(dp), intent(out) :: v2
        real(dp), intent(out) :: d2
        real(dp), intent(out) :: q2
        real(dp), intent(out) :: correlation

        call resolve_basic_inputs(spot, sigma, n, m, strike, rate, dividend, maturity, s, v, np, nt, k, d, q, t)
        s2 = 1.0_dp
        v2 = 0.2_dp
        d2 = 0.0_dp
        q2 = 0.0_dp
        correlation = 0.0_dp
        if (present(spot2)) s2 = spot2
        if (present(sigma2)) v2 = sigma2
        if (present(rate2)) d2 = rate2
        if (present(dividend2)) q2 = dividend2
        if (present(rho)) correlation = rho
        if (s2 <= 0.0_dp) error stop "Quanto LSMC: spot2 must be positive"
        if (v2 < 0.0_dp) error stop "Quanto LSMC: sigma2 must be nonnegative"
        if (abs(correlation) > 1.0_dp) error stop "Quanto LSMC: abs(rho) must not exceed one"
    end subroutine resolve_quanto_inputs

    subroutine set_result(result, values, spot, sigma, n_paths, effective_paths, n_steps, strike, rate, dividend, &
            maturity, option_type, method)
        type(option_result), intent(out) :: result
        real(dp), intent(in) :: values(:)
        real(dp), intent(in) :: spot
        real(dp), intent(in) :: sigma
        integer, intent(in) :: n_paths
        integer, intent(in) :: effective_paths
        integer, intent(in) :: n_steps
        real(dp), intent(in) :: strike
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: dividend
        real(dp), intent(in) :: maturity
        character(len=*), intent(in) :: option_type
        character(len=*), intent(in) :: method

        result%price = mean_value(values)
        result%standard_error = sample_standard_error(values)
        result%spot = spot
        result%sigma = sigma
        result%n_paths = n_paths
        result%effective_paths = effective_paths
        result%n_steps = n_steps
        result%strike = strike
        result%rate = rate
        result%dividend = dividend
        result%maturity = maturity
        result%option_type = option_type
        result%method = method
    end subroutine set_result

    subroutine set_quanto_result(result, spot2, sigma2, rate2, dividend2, rho)
        type(option_result), intent(inout) :: result
        real(dp), intent(in) :: spot2
        real(dp), intent(in) :: sigma2
        real(dp), intent(in) :: rate2
        real(dp), intent(in) :: dividend2
        real(dp), intent(in) :: rho

        result%spot2 = spot2
        result%sigma2 = sigma2
        result%rate2 = rate2
        result%dividend2 = dividend2
        result%rho = rho
    end subroutine set_quanto_result

end module lsmc_pricing
