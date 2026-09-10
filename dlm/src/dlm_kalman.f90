! SPDX-License-Identifier: GPL-2.0-or-later
module dlm_kalman
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use r_linalg, only : solve_spd, thin_svd
    use dlm_types, only : dp, dlm_model, dlm_filter_result, dlm_smooth_result, dlm_forecast_result
    use dlm_types, only : dlm_success, dlm_invalid_shape, dlm_invalid_model, dlm_linalg_failure
    use dlm_types, only : dlm_model_is_valid, dlm_model_is_time_varying, dlm_matrices_at
    implicit none
    private

    public :: dlm_bsample_covariance
    public :: dlm_filter
    public :: dlm_forecast
    public :: dlm_ll
    public :: dlm_residuals
    public :: dlm_smooth
    public :: dlm_svd2var

contains

    subroutine pseudo_solve(a, b, x, info)
        real(dp), intent(in) :: a(:, :) !! Square matrix solved in Moore-Penrose sense using an SVD.
        real(dp), intent(in) :: b(:, :) !! Right-hand sides with shape `(size(a,1), nrhs)`.
        real(dp), intent(out) :: x(:, :) !! Pseudoinverse solution with shape `(size(a,2), nrhs)`.
        integer, intent(out) :: info !! Zero on success or the SVD status code.

        real(dp), allocatable :: u(:, :), singular(:), vt(:, :), temp(:, :)
        real(dp) :: tolerance
        integer :: i, n

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b, 1) /= n .or. size(x, 1) /= n .or. size(x, 2) /= size(b, 2)) then
            info = dlm_invalid_shape
            return
        end if
        call thin_svd(a, u, singular, vt, info)
        if (info /= 0) return
        allocate(temp(n, size(b, 2)))
        temp = matmul(transpose(u), b)
        tolerance = 0.0_dp
        if (n > 0) tolerance = real(n, dp) * epsilon(1.0_dp) * maxval(singular)
        do i = 1, n
            if (singular(i) > tolerance) then
                temp(i, :) = temp(i, :) / singular(i)
            else
                temp(i, :) = 0.0_dp
            end if
        end do
        x = matmul(transpose(vt), temp)
    end subroutine pseudo_solve

    subroutine stable_solve(a, b, x, info)
        real(dp), intent(in) :: a(:, :) !! Symmetric covariance matrix, preferably positive definite.
        real(dp), intent(in) :: b(:, :) !! Right-hand sides with leading dimension equal to `size(a,1)`.
        real(dp), intent(out) :: x(:, :) !! Solutions of `a*x=b`, using a pseudoinverse fallback for semidefinite `a`.
        integer, intent(out) :: info !! Zero on success or a linear-algebra status code.

        integer :: lapack_info

        call solve_spd(a, b, x, lapack_info)
        if (lapack_info == 0) then
            info = 0
        else
            call pseudo_solve(a, b, x, info)
        end if
    end subroutine stable_solve

    subroutine dlm_filter(y, model, result, info)
        real(dp), intent(in) :: y(:, :) !! Observation matrix with time in rows; IEEE NaNs denote missing components.
        type(dlm_model), intent(in) :: model !! Dynamic linear model, including optional time-varying index maps.
        type(dlm_filter_result), intent(out) :: result !! Prior/filtered means and covariance arrays plus one-step forecasts.
        integer, intent(out) :: info !! Zero on success or a `dlm_*`/linear-algebra status code.

        real(dp), allocatable :: ff(:, :), v(:, :), gg(:, :), w(:, :)
        real(dp), allocatable :: ff_good(:, :), v_good(:, :), q_good(:, :), rhs(:, :), solved(:, :)
        real(dp), allocatable :: gain(:, :), innovation(:), identity(:, :), temp(:, :)
        integer, allocatable :: good(:)
        integer :: i, j, k, m_obs, n, p_state, ngood, status

        info = dlm_success
        if (.not. dlm_model_is_valid(model)) then
            info = dlm_invalid_model
            return
        end if
        n = size(y, 1)
        m_obs = size(model%ff, 1)
        p_state = size(model%ff, 2)
        if (size(y, 2) /= m_obs) then
            info = dlm_invalid_shape
            return
        end if
        if (dlm_model_is_time_varying(model)) then
            if (.not. allocated(model%x)) then
                info = dlm_invalid_model
                return
            end if
            if (size(model%x, 1) < n) then
                info = dlm_invalid_shape
                return
            end if
        end if

        result%model = model
        result%y = y
        allocate(result%m(n + 1, p_state), result%c(p_state, p_state, n + 1))
        allocate(result%a(n, p_state), result%r(p_state, p_state, n))
        allocate(result%f(n, m_obs), result%q(m_obs, m_obs, n))
        result%m = 0.0_dp
        result%c = 0.0_dp
        result%a = 0.0_dp
        result%r = 0.0_dp
        result%f = 0.0_dp
        result%q = 0.0_dp
        result%m(1, :) = model%m0
        result%c(:, :, 1) = 0.5_dp * (model%c0 + transpose(model%c0))

        do i = 1, n
            call dlm_matrices_at(model, i, ff, v, gg, w, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            result%a(i, :) = matmul(gg, result%m(i, :))
            result%r(:, :, i) = matmul(gg, matmul(result%c(:, :, i), transpose(gg))) + w
            result%r(:, :, i) = 0.5_dp * (result%r(:, :, i) + transpose(result%r(:, :, i)))
            result%f(i, :) = matmul(ff, result%a(i, :))
            result%q(:, :, i) = matmul(ff, matmul(result%r(:, :, i), transpose(ff))) + v
            result%q(:, :, i) = 0.5_dp * (result%q(:, :, i) + transpose(result%q(:, :, i)))

            ngood = count(.not. ieee_is_nan(y(i, :)))
            if (ngood == 0) then
                result%m(i + 1, :) = result%a(i, :)
                result%c(:, :, i + 1) = result%r(:, :, i)
                cycle
            end if

            allocate(good(ngood))
            k = 0
            do j = 1, m_obs
                if (.not. ieee_is_nan(y(i, j))) then
                    k = k + 1
                    good(k) = j
                end if
            end do
            allocate(ff_good(ngood, p_state), v_good(ngood, ngood), q_good(ngood, ngood))
            allocate(rhs(ngood, p_state), solved(ngood, p_state), gain(p_state, ngood))
            allocate(innovation(ngood), identity(p_state, p_state), temp(p_state, p_state))
            do j = 1, ngood
                ff_good(j, :) = ff(good(j), :)
                innovation(j) = y(i, good(j)) - result%f(i, good(j))
                do k = 1, ngood
                    v_good(j, k) = v(good(j), good(k))
                end do
            end do
            q_good = matmul(ff_good, matmul(result%r(:, :, i), transpose(ff_good))) + v_good
            rhs = matmul(ff_good, result%r(:, :, i))
            call stable_solve(q_good, rhs, solved, status)
            if (status /= 0) then
                info = dlm_linalg_failure
                return
            end if
            gain = transpose(solved)
            result%m(i + 1, :) = result%a(i, :) + matmul(gain, innovation)

            identity = 0.0_dp
            do j = 1, p_state
                identity(j, j) = 1.0_dp
            end do
            temp = identity - matmul(gain, ff_good)
            result%c(:, :, i + 1) = matmul(temp, matmul(result%r(:, :, i), transpose(temp))) + &
                                      matmul(gain, matmul(v_good, transpose(gain)))
            result%c(:, :, i + 1) = 0.5_dp * (result%c(:, :, i + 1) + transpose(result%c(:, :, i + 1)))
            deallocate(good, ff_good, v_good, q_good, rhs, solved, gain, innovation, identity, temp)
        end do
    end subroutine dlm_filter

    subroutine dlm_ll(y, model, value, info)
        real(dp), intent(in) :: y(:, :) !! Observation matrix with IEEE NaNs marking missing components.
        type(dlm_model), intent(in) :: model !! DLM whose R-style negative log likelihood is evaluated.
        real(dp), intent(out) :: value !! Negative log likelihood without the constant `n*log(2*pi)/2`, matching R `dlmLL`.
        integer, intent(out) :: info !! Zero on success; singular observed forecast covariance returns a linear-algebra error.

        type(dlm_filter_result) :: filtered
        real(dp), allocatable :: q_good(:, :), rhs(:, :), solved(:, :), u(:, :), singular(:), vt(:, :)
        real(dp), allocatable :: innovation(:)
        real(dp) :: tolerance
        integer :: i, j, k, m_obs, ngood, status

        value = 0.0_dp
        call dlm_filter(y, model, filtered, status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        m_obs = size(y, 2)
        do i = 1, size(y, 1)
            ngood = count(.not. ieee_is_nan(y(i, :)))
            if (ngood == 0) cycle
            allocate(q_good(ngood, ngood), rhs(ngood, 1), solved(ngood, 1), innovation(ngood))
            k = 0
            do j = 1, m_obs
                if (.not. ieee_is_nan(y(i, j))) then
                    k = k + 1
                    innovation(k) = y(i, j) - filtered%f(i, j)
                end if
            end do
            k = 0
            do j = 1, m_obs
                if (ieee_is_nan(y(i, j))) cycle
                k = k + 1
                q_good(k, :) = pack(filtered%q(j, :, i), .not. ieee_is_nan(y(i, :)))
            end do
            call thin_svd(q_good, u, singular, vt, status)
            if (status /= 0) then
                info = dlm_linalg_failure
                return
            end if
            tolerance = real(ngood, dp) * epsilon(1.0_dp) * maxval(singular)
            if (any(singular <= tolerance)) then
                info = dlm_linalg_failure
                return
            end if
            rhs(:, 1) = innovation
            call solve_spd(q_good, rhs, solved, status)
            if (status /= 0) then
                info = dlm_linalg_failure
                return
            end if
            value = value + 0.5_dp * (sum(log(singular)) + dot_product(innovation, solved(:, 1)))
            deallocate(q_good, rhs, solved, innovation, u, singular, vt)
        end do
        info = dlm_success
    end subroutine dlm_ll

    subroutine dlm_smooth(filtered, result, info)
        type(dlm_filter_result), intent(in) :: filtered !! Filter output providing posterior/prior moments and original model.
        type(dlm_smooth_result), intent(out) :: result !! Rauch-Tung-Striebel smoothed state means and covariances.
        integer, intent(out) :: info !! Zero on success or a `dlm_*`/linear-algebra status code.

        real(dp), allocatable :: ff(:, :), v(:, :), gg(:, :), w(:, :), rhs(:, :), solved(:, :), gain(:, :)
        integer :: i, n, p_state, status

        n = size(filtered%a, 1)
        p_state = size(filtered%m, 2)
        if (size(filtered%m, 1) /= n + 1) then
            info = dlm_invalid_shape
            return
        end if
        allocate(result%s(n + 1, p_state), result%covariance(p_state, p_state, n + 1))
        result%s(n + 1, :) = filtered%m(n + 1, :)
        result%covariance(:, :, n + 1) = filtered%c(:, :, n + 1)
        do i = n, 1, -1
            call dlm_matrices_at(filtered%model, i, ff, v, gg, w, status)
            if (status /= dlm_success) then
                info = status
                return
            end if
            allocate(rhs(p_state, p_state), solved(p_state, p_state), gain(p_state, p_state))
            rhs = matmul(gg, filtered%c(:, :, i))
            call stable_solve(filtered%r(:, :, i), rhs, solved, status)
            if (status /= 0) then
                info = dlm_linalg_failure
                return
            end if
            gain = transpose(solved)
            result%s(i, :) = filtered%m(i, :) + matmul(gain, result%s(i + 1, :) - filtered%a(i, :))
            result%covariance(:, :, i) = filtered%c(:, :, i) + &
                matmul(gain, matmul(result%covariance(:, :, i + 1) - filtered%r(:, :, i), transpose(gain)))
            result%covariance(:, :, i) = 0.5_dp * &
                (result%covariance(:, :, i) + transpose(result%covariance(:, :, i)) )
            deallocate(rhs, solved, gain)
        end do
        info = dlm_success
    end subroutine dlm_smooth

    pure subroutine dlm_forecast(model, n_ahead, result, info)
        type(dlm_model), intent(in) :: model !! Constant DLM supplying current state prior and system matrices.
        integer, intent(in) :: n_ahead !! Number of future time steps; must be nonnegative.
        type(dlm_forecast_result), intent(out) :: result !! Forecast state/observation means and covariance arrays.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status code.

        integer :: i, m_obs, p_state
        real(dp), allocatable :: previous_mean(:), previous_covariance(:, :)

        if (.not. dlm_model_is_valid(model)) then
            info = dlm_invalid_model
            return
        end if
        if (dlm_model_is_time_varying(model)) then
            info = dlm_invalid_model
            return
        end if
        if (n_ahead < 0) then
            info = dlm_invalid_shape
            return
        end if
        p_state = size(model%m0)
        m_obs = size(model%ff, 1)
        allocate(result%a(n_ahead, p_state), result%r(p_state, p_state, n_ahead))
        allocate(result%f(n_ahead, m_obs), result%q(m_obs, m_obs, n_ahead))
        allocate(previous_mean(p_state), previous_covariance(p_state, p_state))
        previous_mean = model%m0
        previous_covariance = model%c0
        do i = 1, n_ahead
            result%a(i, :) = matmul(model%gg, previous_mean)
            result%r(:, :, i) = matmul(model%gg, matmul(previous_covariance, transpose(model%gg))) + model%w
            result%r(:, :, i) = 0.5_dp * (result%r(:, :, i) + transpose(result%r(:, :, i)))
            result%f(i, :) = matmul(model%ff, result%a(i, :))
            result%q(:, :, i) = matmul(model%ff, matmul(result%r(:, :, i), transpose(model%ff))) + model%v
            result%q(:, :, i) = 0.5_dp * (result%q(:, :, i) + transpose(result%q(:, :, i)))
            previous_mean = result%a(i, :)
            previous_covariance = result%r(:, :, i)
        end do
        info = dlm_success
    end subroutine dlm_forecast

    pure subroutine dlm_svd2var(u, d, covariance, info)
        real(dp), intent(in) :: u(:, :) !! Square orthogonal singular-vector matrix.
        real(dp), intent(in) :: d(:) !! Standard-deviation singular values whose squares are covariance eigenvalues.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Reconstructed covariance `u*diag(d**2)*u'`.
        integer, intent(out) :: info !! Zero on success or `dlm_invalid_shape`.

        real(dp), allocatable :: scaled(:, :)
        integer :: n

        n = size(d)
        if (size(u, 1) /= n .or. size(u, 2) /= n) then
            allocate(covariance(0, 0))
            info = dlm_invalid_shape
            return
        end if
        allocate(scaled(n, n), covariance(n, n))
        scaled = u * spread(d, 1, n)
        covariance = matmul(scaled, transpose(scaled))
        covariance = 0.5_dp * (covariance + transpose(covariance))
        info = dlm_success
    end subroutine dlm_svd2var

    pure subroutine dlm_residuals(filtered, residual, standard_deviation, standardized)
        type(dlm_filter_result), intent(in) :: filtered !! Filter result with observations and forecast moments.
        real(dp), allocatable, intent(out) :: residual(:, :) !! One-step errors; input NaNs remain NaN.
        real(dp), allocatable, intent(out) :: standard_deviation(:, :) !! Marginal one-step forecast standard deviations.
        logical, intent(in), optional :: standardized !! Standardize by forecast SD when true; default true.

        logical :: use_standardized
        integer :: i, j, m_obs, n

        use_standardized = .true.
        if (present(standardized)) use_standardized = standardized
        n = size(filtered%y, 1)
        m_obs = size(filtered%y, 2)
        allocate(residual(n, m_obs), standard_deviation(n, m_obs))
        do i = 1, n
            do j = 1, m_obs
                standard_deviation(i, j) = sqrt(max(0.0_dp, filtered%q(j, j, i)))
                if (ieee_is_nan(filtered%y(i, j))) then
                    residual(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
                else
                    residual(i, j) = filtered%y(i, j) - filtered%f(i, j)
                    if (use_standardized .and. standard_deviation(i, j) > 0.0_dp) then
                        residual(i, j) = residual(i, j) / standard_deviation(i, j)
                    end if
                end if
            end do
        end do
    end subroutine dlm_residuals

    subroutine dlm_bsample_covariance(filtered, time_index, conditional_mean, conditional_covariance, next_state, info)
        type(dlm_filter_result), intent(in) :: filtered !! Filter result defining the backward conditional distribution.
        integer, intent(in) :: time_index !! State index from 1 through the number of observations.
        real(dp), intent(out) :: conditional_mean(:) !! Mean of state `time_index` conditional on `next_state`.
        real(dp), intent(out) :: conditional_covariance(:, :) !! Conditional covariance of that backward-sampling state.
        real(dp), intent(in) :: next_state(:) !! State at the following time point used to condition the backward distribution.
        integer, intent(out) :: info !! Zero on success or a `dlm_*`/linear-algebra status code.

        real(dp), allocatable :: ff(:, :), v(:, :), gg(:, :), w(:, :), rhs(:, :), solved(:, :), gain(:, :)
        integer :: p_state, status

        p_state = size(filtered%m, 2)
        if (time_index < 1 .or. time_index > size(filtered%a, 1)) then
            info = dlm_invalid_shape
            return
        end if
        if (size(next_state) /= p_state .or. size(conditional_mean) /= p_state) then
            info = dlm_invalid_shape
            return
        end if
        if (size(conditional_covariance, 1) /= p_state .or. size(conditional_covariance, 2) /= p_state) then
            info = dlm_invalid_shape
            return
        end if
        call dlm_matrices_at(filtered%model, time_index, ff, v, gg, w, status)
        if (status /= dlm_success) then
            info = status
            return
        end if
        allocate(rhs(p_state, p_state), solved(p_state, p_state), gain(p_state, p_state))
        rhs = matmul(gg, filtered%c(:, :, time_index))
        call stable_solve(filtered%r(:, :, time_index), rhs, solved, status)
        if (status /= 0) then
            info = dlm_linalg_failure
            return
        end if
        gain = transpose(solved)
        conditional_mean = filtered%m(time_index, :) + &
                           matmul(gain, next_state - filtered%a(time_index, :))
        conditional_covariance = filtered%c(:, :, time_index) - &
                                 matmul(gain, matmul(filtered%r(:, :, time_index), transpose(gain)))
        conditional_covariance = 0.5_dp * (conditional_covariance + transpose(conditional_covariance))
        info = dlm_success
    end subroutine dlm_bsample_covariance

end module dlm_kalman
