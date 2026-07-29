! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of cvar 0.6 by Georgi N. Boshnakov.
module cvar_garch
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use cvar_kinds, only : dp
    use cvar_status, only : cvar_ok, cvar_invalid_model, cvar_allocation_failure
    use cvar_special, only : normal_quantile, std_student_t_quantile, ged_quantile
    use cvar_rng, only : rng_state, rng_seed, rng_normal, rng_std_student_t, rng_ged
    use cvar_risk, only : empirical_quantile
    implicit none
    private

    integer, parameter, public :: innovation_normal = 1
    integer, parameter, public :: innovation_std_t = 2
    integer, parameter, public :: innovation_ged = 3

    type, public :: garch11_model
        real(dp) :: omega = 0.0_dp
        real(dp) :: alpha = 0.0_dp
        real(dp) :: beta = 0.0_dp
        integer :: innovation = innovation_normal
        real(dp) :: shape = 5.0_dp
        real(dp) :: eps0 = 0.0_dp
        real(dp) :: h0 = 0.0_dp
        real(dp) :: eps0sq = 0.0_dp
        logical :: has_eps0 = .false.
        logical :: has_h0 = .false.
        logical :: has_eps0sq = .false.
    contains
        procedure :: valid => garch11_valid
        procedure :: unconditional_variance => garch11_unconditional_variance
    end type garch11_model

    type, public :: garch11_simulation
        real(dp), allocatable :: eps(:)
        real(dp), allocatable :: h(:)
        real(dp), allocatable :: eta(:)
        integer :: status = cvar_ok
    end type garch11_simulation

    type, public :: garch11_forecast
        real(dp), allocatable :: eps(:)
        real(dp), allocatable :: h(:)
        real(dp), allocatable :: plugin_interval(:, :)
        real(dp), allocatable :: simulation_interval(:, :)
        real(dp), allocatable :: simulated_eps(:, :)
        real(dp), allocatable :: simulated_h(:, :)
        integer :: status = cvar_ok
    end type garch11_forecast

    public :: make_garch11, simulate_garch11, forecast_garch11
    public :: innovation_quantile

contains

    function make_garch11(omega, alpha, beta, innovation, shape, eps0, h0, eps0sq) result(model)
        real(dp), intent(in) :: omega, alpha, beta
        integer, intent(in), optional :: innovation
        real(dp), intent(in), optional :: shape, eps0, h0, eps0sq
        type(garch11_model) :: model

        model%omega = omega
        model%alpha = alpha
        model%beta = beta
        if (present(innovation)) model%innovation = innovation
        if (present(shape)) model%shape = shape
        if (present(eps0)) then
            model%eps0 = eps0
            model%has_eps0 = .true.
        end if
        if (present(h0)) then
            model%h0 = h0
            model%has_h0 = .true.
        end if
        if (present(eps0sq)) then
            model%eps0sq = eps0sq
            model%has_eps0sq = .true.
        end if
    end function make_garch11

    pure function garch11_valid(self) result(is_valid)
        class(garch11_model), intent(in) :: self
        logical :: is_valid

        is_valid = ieee_is_finite(self%omega) .and. ieee_is_finite(self%alpha) .and. &
                   ieee_is_finite(self%beta) .and. self%omega > 0.0_dp .and. &
                   self%alpha >= 0.0_dp .and. self%beta >= 0.0_dp .and. &
                   self%alpha + self%beta < 1.0_dp
        if (.not. is_valid) return
        select case (self%innovation)
        case (innovation_normal)
            continue
        case (innovation_std_t)
            is_valid = ieee_is_finite(self%shape) .and. self%shape > 2.0_dp
        case (innovation_ged)
            is_valid = ieee_is_finite(self%shape) .and. self%shape > 0.0_dp
        case default
            is_valid = .false.
        end select
        if (self%has_h0) is_valid = is_valid .and. self%h0 > 0.0_dp
        if (self%has_eps0sq) is_valid = is_valid .and. self%eps0sq >= 0.0_dp
    end function garch11_valid

    pure function garch11_unconditional_variance(self) result(value)
        class(garch11_model), intent(in) :: self
        real(dp) :: value

        value = self%omega / (1.0_dp - self%alpha - self%beta)
    end function garch11_unconditional_variance

    subroutine simulate_garch11(model, n, result, burnin, seed)
        type(garch11_model), intent(in) :: model
        integer, intent(in) :: n
        type(garch11_simulation), intent(out) :: result
        integer, intent(in), optional :: burnin
        integer(int64), intent(in), optional :: seed
        type(rng_state) :: rng
        integer :: use_burnin
        integer(int64) :: use_seed

        use_burnin = 0
        if (present(burnin)) use_burnin = max(0, burnin)
        use_seed = 5489_int64
        if (present(seed)) use_seed = seed
        call rng_seed(rng, use_seed)
        call simulate_with_rng(model, n, use_burnin, rng, result)
    end subroutine simulate_garch11

    subroutine simulate_with_rng(model, n, burnin, rng, result)
        type(garch11_model), intent(in) :: model
        integer, intent(in) :: n, burnin
        type(rng_state), intent(inout) :: rng
        type(garch11_simulation), intent(out) :: result
        real(dp), allocatable :: all_eps(:), all_h(:), all_eta(:)
        real(dp) :: initial_h, initial_eps_sq
        integer :: total, i, alloc_status

        if (.not. model%valid() .or. n < 1 .or. burnin < 0) then
            result%status = cvar_invalid_model
            return
        end if
        total = n + burnin
        allocate(all_eps(total), all_h(total), all_eta(total), stat=alloc_status)
        if (alloc_status /= 0) then
            result%status = cvar_allocation_failure
            return
        end if

        initial_h = model%unconditional_variance()
        if (model%has_h0) initial_h = model%h0
        initial_eps_sq = model%unconditional_variance()
        if (model%has_eps0) initial_eps_sq = model%eps0 * model%eps0
        if (model%has_eps0sq) initial_eps_sq = model%eps0sq

        all_h(1) = model%omega + model%alpha * initial_eps_sq + model%beta * initial_h
        all_eta(1) = draw_innovation(model, rng)
        all_eps(1) = sqrt(all_h(1)) * all_eta(1)
        do i = 2, total
            all_h(i) = model%omega + model%alpha * all_eps(i - 1)**2 + &
                       model%beta * all_h(i - 1)
            all_eta(i) = draw_innovation(model, rng)
            all_eps(i) = sqrt(all_h(i)) * all_eta(i)
        end do

        allocate(result%eps(n), result%h(n), result%eta(n), stat=alloc_status)
        if (alloc_status /= 0) then
            result%status = cvar_allocation_failure
            return
        end if
        result%eps = all_eps(burnin + 1:total)
        result%h = all_h(burnin + 1:total)
        result%eta = all_eta(burnin + 1:total)
        result%status = cvar_ok
    end subroutine simulate_with_rng

    subroutine forecast_garch11(model, eps, sigmasq, n_ahead, result, nsim, seed, &
                                lower_probability, upper_probability)
        type(garch11_model), intent(in) :: model
        real(dp), intent(in) :: eps(:), sigmasq(:)
        integer, intent(in) :: n_ahead
        type(garch11_forecast), intent(out) :: result
        integer, intent(in), optional :: nsim
        integer(int64), intent(in), optional :: seed
        real(dp), intent(in), optional :: lower_probability, upper_probability
        type(garch11_model) :: future_model
        type(garch11_simulation) :: path
        type(rng_state) :: rng
        integer :: use_nsim, i, alloc_status
        integer(int64) :: use_seed
        real(dp) :: p_lower, p_upper, q_lower, q_upper

        if (.not. model%valid() .or. size(eps) < 1 .or. size(sigmasq) < 1 .or. &
            n_ahead < 1 .or. sigmasq(size(sigmasq)) <= 0.0_dp) then
            result%status = cvar_invalid_model
            return
        end if
        use_nsim = 1000
        if (present(nsim)) use_nsim = nsim
        if (use_nsim < 1) then
            result%status = cvar_invalid_model
            return
        end if
        p_lower = 0.025_dp
        p_upper = 0.975_dp
        if (present(lower_probability)) p_lower = lower_probability
        if (present(upper_probability)) p_upper = upper_probability
        if (p_lower <= 0.0_dp .or. p_upper >= 1.0_dp .or. p_lower >= p_upper) then
            result%status = cvar_invalid_model
            return
        end if

        allocate(result%eps(n_ahead), result%h(n_ahead), &
                 result%plugin_interval(n_ahead, 2), &
                 result%simulation_interval(n_ahead, 2), &
                 result%simulated_eps(n_ahead, use_nsim), &
                 result%simulated_h(n_ahead, use_nsim), stat=alloc_status)
        if (alloc_status /= 0) then
            result%status = cvar_allocation_failure
            return
        end if

        result%eps = 0.0_dp
        result%h(1) = model%omega + model%alpha * eps(size(eps))**2 + &
                      model%beta * sigmasq(size(sigmasq))
        do i = 2, n_ahead
            result%h(i) = model%omega + (model%alpha + model%beta) * result%h(i - 1)
        end do

        q_lower = innovation_quantile(model, p_lower)
        q_upper = innovation_quantile(model, p_upper)
        result%plugin_interval(:, 1) = q_lower * sqrt(result%h)
        result%plugin_interval(:, 2) = q_upper * sqrt(result%h)

        future_model = model
        future_model%eps0 = eps(size(eps))
        future_model%eps0sq = future_model%eps0**2
        future_model%h0 = sigmasq(size(sigmasq))
        future_model%has_eps0 = .true.
        future_model%has_eps0sq = .true.
        future_model%has_h0 = .true.

        use_seed = 5489_int64
        if (present(seed)) use_seed = seed
        call rng_seed(rng, use_seed)
        do i = 1, use_nsim
            call simulate_with_rng(future_model, n_ahead, 0, rng, path)
            if (path%status /= cvar_ok) then
                result%status = path%status
                return
            end if
            result%simulated_eps(:, i) = path%eps
            result%simulated_h(:, i) = path%h
        end do
        do i = 1, n_ahead
            result%simulation_interval(i, 1) = &
                empirical_quantile(result%simulated_eps(i, :), p_lower)
            result%simulation_interval(i, 2) = &
                empirical_quantile(result%simulated_eps(i, :), p_upper)
        end do
        result%status = cvar_ok
    end subroutine forecast_garch11

    function draw_innovation(model, rng) result(value)
        type(garch11_model), intent(in) :: model
        type(rng_state), intent(inout) :: rng
        real(dp) :: value

        select case (model%innovation)
        case (innovation_normal)
            value = rng_normal(rng)
        case (innovation_std_t)
            value = rng_std_student_t(rng, model%shape)
        case (innovation_ged)
            value = rng_ged(rng, model%shape)
        case default
            value = 0.0_dp
        end select
    end function draw_innovation

    pure function innovation_quantile(model, p) result(value)
        type(garch11_model), intent(in) :: model
        real(dp), intent(in) :: p
        real(dp) :: value

        select case (model%innovation)
        case (innovation_normal)
            value = normal_quantile(p)
        case (innovation_std_t)
            value = std_student_t_quantile(p, model%shape)
        case (innovation_ged)
            value = ged_quantile(p, model%shape)
        case default
            value = 0.0_dp
        end select
    end function innovation_quantile

end module cvar_garch
