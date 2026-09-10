! SPDX-License-Identifier: GPL-2.0-or-later
module dlm_random
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use r_linalg, only : thin_svd
    use dlm_types, only : dp, dlm_model, dlm_filter_result, dlm_gibbs_dig_result
    use dlm_types, only : dlm_success, dlm_invalid_shape, dlm_invalid_argument, dlm_linalg_failure
    use dlm_types, only : dlm_model_is_valid, dlm_model_is_time_varying
    use dlm_kalman, only : dlm_bsample_covariance, dlm_filter
    implicit none
    private

    public :: dlm_bsample
    public :: dlm_gibbs_dig
    public :: dlm_random_model
    public :: rwishart
    public :: seed_dlm_rng

contains

    subroutine seed_dlm_rng(seed)
        integer, intent(in) :: seed !! Deterministic scalar seed expanded to the compiler's intrinsic RNG seed vector.

        integer, allocatable :: put(:)
        integer :: i, n

        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(seed + 104729 * i + 7919 * i * i, huge(1) - 1)
            if (put(i) == 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine seed_dlm_rng

    subroutine normal_vector(z)
        real(dp), intent(out) :: z(:) !! Independent standard-normal variates generated with a Box-Muller transform.

        real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)
        real(dp) :: r, theta, u1, u2
        integer :: i

        i = 1
        do while (i <= size(z))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(1.0_dp))
            r = sqrt(-2.0_dp * log(u1))
            theta = two_pi * u2
            z(i) = r * cos(theta)
            if (i + 1 <= size(z)) z(i + 1) = r * sin(theta)
            i = i + 2
        end do
    end subroutine normal_vector

    subroutine covariance_draw(mean, covariance, draw, info)
        real(dp), intent(in) :: mean(:) !! Mean vector of the multivariate normal draw.
        real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-semidefinite covariance matrix.
        real(dp), intent(out) :: draw(:) !! Generated normal vector with the requested mean and covariance.
        integer, intent(out) :: info !! Zero on success or a linear-algebra/shape error code.

        real(dp), allocatable :: u(:, :), singular(:), vt(:, :), z(:)
        integer :: i, n

        n = size(mean)
        if (size(covariance, 1) /= n .or. size(covariance, 2) /= n .or. size(draw) /= n) then
            info = dlm_invalid_shape
            return
        end if
        call thin_svd(0.5_dp * (covariance + transpose(covariance)), u, singular, vt, info)
        if (info /= 0) then
            info = dlm_linalg_failure
            return
        end if
        allocate(z(n))
        call normal_vector(z)
        do i = 1, n
            z(i) = sqrt(max(0.0_dp, singular(i))) * z(i)
        end do
        draw = mean + matmul(u, z)
        info = dlm_success
    end subroutine covariance_draw

    subroutine dlm_bsample(filtered, theta, info, seed)
        type(dlm_filter_result), intent(in) :: filtered !! Kalman filter result defining the backward-sampling distribution.
        real(dp), allocatable, intent(out) :: theta(:, :) !! Sampled path by time row, including the initial state.
        integer, intent(out) :: info !! Zero on success or a `dlm_*`/linear-algebra status code.
        integer, intent(in), optional :: seed !! Optional deterministic RNG seed applied before drawing the path.

        real(dp), allocatable :: conditional_mean(:), conditional_covariance(:, :)
        integer :: i, n, p_state, status

        if (present(seed)) call seed_dlm_rng(seed)
        n = size(filtered%a, 1)
        p_state = size(filtered%m, 2)
        allocate(theta(n + 1, p_state), conditional_mean(p_state), conditional_covariance(p_state, p_state))
        call covariance_draw(filtered%m(n + 1, :), filtered%c(:, :, n + 1), theta(n + 1, :), status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        do i = n, 1, -1
            call dlm_bsample_covariance(filtered, i, conditional_mean, conditional_covariance, theta(i + 1, :), status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            call covariance_draw(conditional_mean, conditional_covariance, theta(i, :), status)
            if (status /= dlm_success) then
                info = status
                return
            end if
        end do
        info = dlm_success
    end subroutine dlm_bsample

    subroutine rwishart(df, p, wishart, info, sigma, sqrt_sigma, seed)
        integer, intent(in) :: df !! Integer Wishart degrees of freedom; must be at least `p`.
        integer, intent(in) :: p !! Matrix dimension and number of Bartlett factors; must be positive.
        real(dp), allocatable, intent(out) :: wishart(:, :) !! Generated `p` by `p` Wishart random matrix.
        integer, intent(out) :: info !! Zero on success or a `dlm_*`/linear-algebra status code.
        real(dp), intent(in), optional :: sigma(:, :) !! Scale covariance `Sigma`; mutually exclusive with `sqrt_sigma`.
        real(dp), intent(in), optional :: sqrt_sigma(:, :) !! Factor satisfying `transpose(sqrt_sigma)*sqrt_sigma = Sigma`.
        integer, intent(in), optional :: seed !! Optional deterministic RNG seed applied before this draw.

        real(dp), allocatable :: z(:, :), root(:, :), normals(:), u(:, :), singular(:), vt(:, :), work(:, :)
        integer :: i, j, status

        info = dlm_success
        if (p < 1 .or. df < p) then
            info = dlm_invalid_argument
            return
        end if
        if (present(sigma) .and. present(sqrt_sigma)) then
            info = dlm_invalid_argument
            return
        end if
        if (present(sigma)) then
            if (size(sigma, 1) /= p .or. size(sigma, 2) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(sqrt_sigma)) then
            if (size(sqrt_sigma, 1) /= p .or. size(sqrt_sigma, 2) /= p) then
                info = dlm_invalid_shape
                return
            end if
        end if
        if (present(seed)) call seed_dlm_rng(seed)
        allocate(z(p, p), root(p, p), normals(max(df, p)))
        z = 0.0_dp
        do i = 1, p
            call normal_vector(normals(1:df - i + 1))
            z(i, i) = sqrt(sum(normals(1:df - i + 1) ** 2))
            if (i < p) then
                call normal_vector(normals(1:p - i))
                do j = i + 1, p
                    z(j, i) = normals(j - i)
                end do
            end if
        end do

        root = 0.0_dp
        do i = 1, p
            root(i, i) = 1.0_dp
        end do
        if (present(sqrt_sigma)) root = sqrt_sigma
        if (present(sigma)) then
            call thin_svd(0.5_dp * (sigma + transpose(sigma)), u, singular, vt, status)
            if (status /= 0) then
                info = dlm_linalg_failure
                return
            end if
            root = 0.0_dp
            do i = 1, p
                root(i, :) = sqrt(max(0.0_dp, singular(i))) * u(:, i)
            end do
        end if
        work = matmul(z, root)
        wishart = matmul(transpose(work), work)
        wishart = 0.5_dp * (wishart + transpose(wishart))
    end subroutine rwishart

    recursive subroutine gamma_rate_draw(shape, rate, draw, info)
        real(dp), intent(in) :: shape !! Positive gamma shape parameter.
        real(dp), intent(in) :: rate !! Positive gamma rate parameter, inverse of the scale.
        real(dp), intent(out) :: draw !! Generated gamma variate.
        integer, intent(out) :: info !! Zero on success or `dlm_invalid_argument` for nonpositive parameters.

        real(dp) :: c, d, u, v, z(1), x
        integer :: iteration

        if (shape <= 0.0_dp .or. rate <= 0.0_dp) then
            draw = 0.0_dp
            info = dlm_invalid_argument
            return
        end if

        if (shape < 1.0_dp) then
            call gamma_rate_draw(shape + 1.0_dp, rate, draw, info)
            if (info /= dlm_success) return
            call random_number(u)
            u = max(u, tiny(1.0_dp))
            draw = draw * u ** (1.0_dp / shape)
            return
        end if

        d = shape - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do iteration = 1, 100000
            call normal_vector(z)
            x = 1.0_dp + c * z(1)
            if (x <= 0.0_dp) cycle
            v = x ** 3
            call random_number(u)
            if (u < 1.0_dp - 0.0331_dp * z(1) ** 4) then
                draw = d * v / rate
                info = dlm_success
                return
            end if
            u = max(u, tiny(1.0_dp))
            if (log(u) < 0.5_dp * z(1) ** 2 + d * (1.0_dp - v + log(v))) then
                draw = d * v / rate
                info = dlm_success
                return
            end if
        end do
        draw = 0.0_dp
        info = dlm_invalid_argument
    end subroutine gamma_rate_draw

    pure subroutine expand_positive_vector(input, n, output, info)
        real(dp), intent(in) :: input(:) !! Positive prior parameters of length one or requested output length.
        integer, intent(in) :: n !! Requested output length; must be positive.
        real(dp), allocatable, intent(out) :: output(:) !! Expanded prior parameter vector of length `n`.
        integer, intent(out) :: info !! Zero on success or `dlm_invalid_shape`/`dlm_invalid_argument`.

        if (n < 1 .or. any(input <= 0.0_dp)) then
            allocate(output(0))
            info = dlm_invalid_argument
            return
        end if
        if (size(input) == 1) then
            allocate(output(n))
            output = input(1)
        else if (size(input) == n) then
            output = input
        else
            allocate(output(0))
            info = dlm_invalid_shape
            return
        end if
        info = dlm_success
    end subroutine expand_positive_vector

    subroutine dlm_gibbs_dig(y, model, n_sample, result, info, thin, indices, save_states, &
                             a_y, b_y, a_theta, b_theta, shape_y, rate_y, shape_theta, &
                             rate_theta, seed)
        real(dp), intent(in) :: y(:) !! Univariate observation sequence; missing values are not supported by this sampler.
        type(dlm_model), intent(in) :: model !! Constant univariate DLM with diagonal system covariance `W`.
        integer, intent(in) :: n_sample !! Number of retained Gibbs draws; must be positive.
        type(dlm_gibbs_dig_result), intent(out) :: result !! Draws of `V`, selected diagonal `W`, and optional states.
        integer, intent(out) :: info !! Zero on success or a `dlm_*`/linear-algebra status code.
        integer, intent(in), optional :: thin !! Sweeps discarded between retained draws; default zero.
        integer, intent(in), optional :: indices(:) !! One-based state indices whose diagonal `W` entries are sampled.
        logical, intent(in), optional :: save_states !! Retain sampled state paths when true; default true.
        real(dp), intent(in), optional :: a_y !! Prior mean of observation precision; requires `b_y`.
        real(dp), intent(in), optional :: b_y !! Prior variance of observation precision; requires `a_y`.
        real(dp), intent(in), optional :: a_theta(:) !! Precision prior means; scalar or length `p`, recycled independently.
        real(dp), intent(in), optional :: b_theta(:) !! Precision prior variances; scalar or length `p`, recycled independently.
        real(dp), intent(in), optional :: shape_y !! Gamma prior shape for observation precision; requires `rate_y`.
        real(dp), intent(in), optional :: rate_y !! Gamma prior rate for observation precision; requires `shape_y`.
        real(dp), intent(in), optional :: shape_theta(:) !! Gamma prior shapes for system precisions; length one or selected count.
        real(dp), intent(in), optional :: rate_theta(:) !! Gamma prior rates for system precisions; length one or selected count.
        integer, intent(in), optional :: seed !! Optional deterministic RNG seed applied before Gibbs sampling.

        type(dlm_model) :: work_model
        type(dlm_filter_result) :: filtered
        real(dp), allocatable :: theta(:, :), y_matrix(:, :)
        real(dp), allocatable :: prior_shape_theta(:), prior_rate_theta(:)
        real(dp), allocatable :: theta_mean(:), theta_variance(:)
        real(dp) :: draw, prior_rate_y, prior_shape_y, residual, ss_theta, ss_y
        logical, allocatable :: selected(:)
        integer, allocatable :: idx(:)
        logical :: keep_states
        integer :: every, i, it, j, n, p, r, retained, status, sweeps

        info = dlm_success
        if (.not. dlm_model_is_valid(model) .or. dlm_model_is_time_varying(model)) then
            info = dlm_invalid_argument
            return
        end if
        if (size(model%ff, 1) /= 1 .or. n_sample < 1 .or. size(y) < 1) then
            info = dlm_invalid_argument
            return
        end if
        if (any(ieee_is_nan(y))) then
            info = dlm_invalid_argument
            return
        end if
        r = size(model%m0)
        do i = 1, r
            do j = 1, r
                if (i /= j .and. abs(model%w(i, j)) > 100.0_dp * epsilon(1.0_dp)) then
                    info = dlm_invalid_argument
                    return
                end if
            end do
        end do

        every = 1
        if (present(thin)) then
            if (thin < 0) then
                info = dlm_invalid_argument
                return
            end if
            every = thin + 1
        end if
        sweeps = n_sample * every
        keep_states = .true.
        if (present(save_states)) keep_states = save_states

        allocate(selected(r))
        selected = .false.
        if (present(indices)) then
            if (size(indices) < 1) then
                info = dlm_invalid_argument
                return
            end if
            do i = 1, size(indices)
                if (indices(i) < 1 .or. indices(i) > r) then
                    info = dlm_invalid_argument
                    return
                end if
                selected(indices(i)) = .true.
            end do
        else
            selected = .true.
        end if
        p = count(selected)
        allocate(idx(p))
        idx = pack([(i, i = 1, r)], selected)

        if (present(a_y) .neqv. present(b_y)) then
            info = dlm_invalid_argument
            return
        end if
        if (present(a_y)) then
            if (a_y <= 0.0_dp .or. b_y <= 0.0_dp) then
                info = dlm_invalid_argument
                return
            end if
            prior_shape_y = a_y ** 2 / b_y
            prior_rate_y = a_y / b_y
        else
            if (.not. present(shape_y) .or. .not. present(rate_y)) then
                info = dlm_invalid_argument
                return
            end if
            if (shape_y <= 0.0_dp .or. rate_y <= 0.0_dp) then
                info = dlm_invalid_argument
                return
            end if
            prior_shape_y = shape_y
            prior_rate_y = rate_y
        end if

        if (present(a_theta) .neqv. present(b_theta)) then
            info = dlm_invalid_argument
            return
        end if
        if (present(a_theta)) then
            call expand_positive_vector(a_theta, p, theta_mean, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            call expand_positive_vector(b_theta, p, theta_variance, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            prior_shape_theta = theta_mean ** 2 / theta_variance
            prior_rate_theta = theta_mean / theta_variance
        else
            if (.not. present(shape_theta) .or. .not. present(rate_theta)) then
                info = dlm_invalid_argument
                return
            end if
            call expand_positive_vector(shape_theta, p, prior_shape_theta, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            call expand_positive_vector(rate_theta, p, prior_rate_theta, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
        end if

        if (present(seed)) call seed_dlm_rng(seed)
        n = size(y)
        allocate(y_matrix(n, 1))
        y_matrix(:, 1) = y
        allocate(result%d_v(n_sample), result%d_w(n_sample, p))
        if (keep_states) allocate(result%theta(n + 1, r, n_sample))
        work_model = model
        retained = 0

        do it = 1, sweeps
            call dlm_filter(y_matrix, work_model, filtered, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            call dlm_bsample(filtered, theta, status)
            if (status /= dlm_success) then
                info = status
                return
            end if

            ss_y = 0.0_dp
            do i = 1, n
                residual = y(i) - dot_product(work_model%ff(1, :), theta(i + 1, :))
                ss_y = ss_y + residual ** 2
            end do
            call gamma_rate_draw(prior_shape_y + 0.5_dp * real(n, dp), &
                                 prior_rate_y + 0.5_dp * ss_y, draw, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            work_model%v(1, 1) = 1.0_dp / draw

            do j = 1, p
                ss_theta = 0.0_dp
                do i = 1, n
                    residual = theta(i + 1, idx(j)) - &
                               dot_product(work_model%gg(idx(j), :), theta(i, :))
                    ss_theta = ss_theta + residual ** 2
                end do
                call gamma_rate_draw(prior_shape_theta(j) + 0.5_dp * real(n, dp), &
                                     prior_rate_theta(j) + 0.5_dp * ss_theta, draw, status)
                if (status /= dlm_success) then
                    info = status
                    return
                end if
                work_model%w(idx(j), idx(j)) = 1.0_dp / draw
            end do

            if (modulo(it, every) == 0) then
                retained = retained + 1
                result%d_v(retained) = work_model%v(1, 1)
                do j = 1, p
                    result%d_w(retained, j) = work_model%w(idx(j), idx(j))
                end do
                if (keep_states) result%theta(:, :, retained) = theta
            end if
        end do
        info = dlm_success
    end subroutine dlm_gibbs_dig

    subroutine dlm_random_model(m_obs, p_state, model, info, seed)
        integer, intent(in) :: m_obs !! Observation dimension of the generated constant DLM; must be positive.
        integer, intent(in) :: p_state !! State dimension of the generated constant DLM; must be positive.
        type(dlm_model), intent(out) :: model !! Random constant DLM with stable transition scaling and positive-definite variances.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status code.
        integer, intent(in), optional :: seed !! Optional deterministic RNG seed applied before model generation.

        real(dp), allocatable :: z(:), vdraw(:, :), wdraw(:, :)
        real(dp) :: bound
        integer :: i, status

        if (m_obs < 1 .or. p_state < 1) then
            info = dlm_invalid_argument
            return
        end if
        if (present(seed)) call seed_dlm_rng(seed)
        allocate(model%m0(p_state), model%c0(p_state, p_state), model%ff(m_obs, p_state))
        allocate(model%v(m_obs, m_obs), model%gg(p_state, p_state), model%w(p_state, p_state))
        allocate(z(max(m_obs * p_state, p_state * p_state)))
        call normal_vector(z(1:m_obs * p_state))
        model%ff = reshape(z(1:m_obs * p_state), [m_obs, p_state])
        call normal_vector(z(1:p_state * p_state))
        model%gg = reshape(z(1:p_state * p_state), [p_state, p_state])
        bound = maxval(sum(abs(model%gg), dim=2))
        if (bound > 0.95_dp) model%gg = 0.95_dp * model%gg / bound
        call rwishart(2 * m_obs, m_obs, vdraw, status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        call rwishart(2 * p_state, p_state, wdraw, status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        model%v = vdraw
        model%w = wdraw
        model%m0 = 0.0_dp
        model%c0 = 0.0_dp
        do i = 1, p_state
            model%c0(i, i) = 100.0_dp
        end do
        info = dlm_success
    end subroutine dlm_random_model

end module dlm_random
