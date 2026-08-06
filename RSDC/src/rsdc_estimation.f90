! SPDX-License-Identifier: GPL-3.0-only
module rsdc_estimation
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_model, rsdc_control, rsdc_filter_result
    use rsdc_types, only: rsdc_const, rsdc_nox, rsdc_tvtp
    use rsdc_parameters, only: lower_tri_size, partial_to_correlation, correlation_to_partial
    use rsdc_parameters, only: matrix_to_correlations, correlations_to_matrix, transition_at
    use rsdc_parameters, only: expected_parameter_count, unpack_natural_parameters
    use rsdc_likelihood, only: rsdc_negative_log_likelihood
    use rsdc_filter, only: rsdc_hamilton
    use rsdc_linalg, only: sample_correlation, inverse_matrix, mean_columns
    use deoptimr, only: jde_optimize, jde_control, de_result, seed_rng
    implicit none
    private
    public :: rsdc_estimate, rsdc_natural_from_theta, rsdc_theta_from_model
    public :: rsdc_theta_from_natural

    real(dp), allocatable, save :: context_y(:, :), context_x(:, :)
    integer, save :: context_method = rsdc_const
    integer, save :: context_n = 1
    integer, save :: context_k = 0
    integer, save :: context_p = 0

contains

    subroutine rsdc_estimate(method, residuals, n_regimes, model, x, control, out_of_sample, ok)
        integer, intent(in) :: method, n_regimes
        real(dp), intent(in) :: residuals(:, :)
        type(rsdc_model), intent(out) :: model
        real(dp), intent(in), optional :: x(:, :)
        type(rsdc_control), intent(in), optional :: control
        logical, intent(in), optional :: out_of_sample
        logical, intent(out), optional :: ok

        type(rsdc_control) :: cfg
        type(jde_control) :: de_cfg
        type(de_result) :: de_out
        real(dp), allocatable :: lower(:), upper(:), theta0(:), init(:, :), best_theta(:), natural(:)
        real(dp), allocatable :: yfit(:, :), xfit(:, :), yoos(:, :), xoos(:, :)
        real(dp) :: value, best_value
        integer :: d, start, nstarts, cut, n, k, p
        logical :: oos, good

        cfg = rsdc_control()
        if (present(control)) cfg = control
        oos = .false.; if (present(out_of_sample)) oos = out_of_sample
        k = size(residuals, 2)
        n = merge(1, n_regimes, method == rsdc_const)
        p = 0; if (present(x)) p = size(x, 2)
        good = size(residuals, 1) >= 2 .and. k >= 2
        good = good .and. (method == rsdc_const .or. n >= 2)
        if (method == rsdc_tvtp) good = good .and. present(x) .and. size(x, 1) == size(residuals, 1)
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if

        if (oos) then
            cut = nint(0.7_dp * real(size(residuals, 1), dp))
            cut = max(2, min(cut, size(residuals, 1) - 1))
            allocate(yfit(cut, k), yoos(size(residuals, 1) - cut, k))
            yfit = residuals(1:cut, :); yoos = residuals(cut + 1:, :)
            if (method == rsdc_tvtp) then
                allocate(xfit(cut, p), xoos(size(residuals, 1) - cut, p))
                xfit = x(1:cut, :); xoos = x(cut + 1:, :)
            else
                allocate(xfit(0, 0), xoos(0, 0))
            end if
        else
            allocate(yfit(size(residuals, 1), k)); yfit = residuals
            allocate(xfit(0, 0), yoos(0, k), xoos(0, 0))
            if (method == rsdc_tvtp) then
                deallocate(xfit); allocate(xfit(size(x, 1), p)); xfit = x
            end if
        end if

        call make_theta_bounds(method, n, k, p, lower, upper)
        d = size(lower)
        if (allocated(cfg%start)) then
            if (size(cfg%start) /= expected_parameter_count(method, n, k, p)) then
                if (present(ok)) ok = .false.
                return
            end if
            call rsdc_theta_from_natural(cfg%start, method, n, k, p, theta0, good)
            if (.not. good) then
                if (present(ok)) ok = .false.
                return
            end if
            lower = min(lower, theta0)
            upper = max(upper, theta0)
        else
            call initial_theta(method, n, yfit, xfit, theta0)
        end if
        allocate(init(d, 1), best_theta(d)); init(:, 1) = theta0
        best_value = huge(1.0_dp)
        nstarts = max(1, cfg%n_starts)

        context_method = method; context_n = n; context_k = k; context_p = p
        if (allocated(context_y)) deallocate(context_y)
        if (allocated(context_x)) deallocate(context_x)
        allocate(context_y(size(yfit, 1), k)); context_y = yfit
        if (method == rsdc_tvtp) then
            allocate(context_x(size(xfit, 1), p)); context_x = xfit
        else
            allocate(context_x(0, 0))
        end if

        if (allocated(cfg%start)) then
            best_theta = theta0
            call local_pattern_search(best_theta, lower, upper, cfg, best_value)
            de_out%iterations = 0
        else
            do start = 1, nstarts
                call seed_rng(cfg%seed + start - 1)
                de_cfg = jde_control()
                de_cfg%population_size = cfg%population_size
                de_cfg%max_iterations = cfg%max_global_iterations
                de_cfg%tolerance = max(cfg%tolerance * 0.01_dp, 1.0e-12_dp)
                de_cfg%trace = cfg%trace
                de_cfg%save_population = .false.
                if (start == 1) then
                    call jde_optimize(lower, upper, theta_objective, de_out, de_cfg, initial_population=init)
                else
                    call jde_optimize(lower, upper, theta_objective, de_out, de_cfg)
                end if
                call local_pattern_search(de_out%parameters, lower, upper, cfg, value)
                if (value < best_value) then
                    best_value = value
                    best_theta = de_out%parameters
                end if
            end do
        end if

        call rsdc_natural_from_theta(best_theta, method, n, k, p, natural)
        call build_model(natural, method, n, k, p, yfit, xfit, model, good)
        model%iterations = de_out%iterations
        model%convergence = merge(0, 1, good .and. best_value < 1.0e9_dp)
        if (oos .and. good) call compute_oos_likelihood(model, yfit, yoos, xfit, xoos)
        if (cfg%compute_vcov .and. good) call estimate_vcov(model, yfit, xfit)
        if (present(ok)) ok = good
    end subroutine rsdc_estimate

    real(dp) function theta_objective(theta) result(value)
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable :: natural(:)
        call rsdc_natural_from_theta(theta, context_method, context_n, context_k, context_p, natural)
        if (context_method == rsdc_tvtp) then
            value = rsdc_negative_log_likelihood(natural, context_y, context_method, context_n, context_x)
        else
            value = rsdc_negative_log_likelihood(natural, context_y, context_method, context_n)
        end if
    end function theta_objective

    subroutine make_theta_bounds(method, n, k, p, lower, upper)
        integer, intent(in) :: method, n, k, p
        real(dp), allocatable, intent(out) :: lower(:), upper(:)
        integer :: c, nh, d
        c = lower_tri_size(k)
        nh = 0
        d = 0
        select case (method)
        case (rsdc_const)
            nh = 0; d = c
        case (rsdc_nox)
            nh = n * (n - 1); d = nh + n * c
        case (rsdc_tvtp)
            nh = merge(n * p, n * (n - 1) * p, n == 2)
            d = nh + n * c
        end select
        allocate(lower(d), upper(d))
        if (nh > 0) then
            if (method == rsdc_tvtp) then
                lower(1:nh) = -10.0_dp; upper(1:nh) = 10.0_dp
            else
                lower(1:nh) = -4.6_dp; upper(1:nh) = 4.6_dp
            end if
        end if
        lower(nh + 1:) = -0.999_dp; upper(nh + 1:) = 0.999_dp
    end subroutine make_theta_bounds

    subroutine initial_theta(method, n, y, x, theta)
        integer, intent(in) :: method, n
        real(dp), intent(in) :: y(:, :), x(:, :)
        real(dp), allocatable, intent(out) :: theta(:)
        real(dp), allocatable :: cor(:, :), z(:), lower(:), upper(:)
        integer :: c, nh, s, lo, hi, i, j, p
        call make_theta_bounds(method, n, size(y, 2), size(x, 2), lower, upper)
        allocate(theta(size(lower))); theta = 0.0_dp
        c = lower_tri_size(size(y, 2)); p = size(x, 2); nh = size(theta) - n * c
        allocate(cor(size(y, 2), size(y, 2)), z(c))
        call sample_correlation(y, cor)
        call correlation_to_partial(cor, z)
        if (method == rsdc_nox) then
            if (n == 2) then
                theta(1:2) = log(0.95_dp / 0.05_dp)
            else
                do i = 1, n
                    lo = (i - 1) * (n - 1) + 1
                    hi = i * (n - 1)
                    theta(lo:hi) = -2.0_dp
                    if (i < n) theta(lo + i - 1) = 2.0_dp
                end do
            end if
        else if (method == rsdc_tvtp .and. p > 0) then
            if (n == 2) then
                do i = 1, n
                    theta((i - 1) * p + 1) = log(0.95_dp / 0.05_dp)
                end do
            else
                do i = 1, n
                    do j = 1, n - 1
                        if (i < n .and. i == j) then
                            theta((i - 1) * (n - 1) * p + (j - 1) * p + 1) = &
                                log(0.95_dp * real(n - 1, dp) / 0.05_dp)
                        else if (i == n) then
                            theta((i - 1) * (n - 1) * p + (j - 1) * p + 1) = &
                                log(0.05_dp / (0.95_dp * real(n - 1, dp)))
                        end if
                    end do
                end do
            end if
        end if
        do s = 1, n
            lo = nh + (s - 1) * c + 1
            hi = nh + s * c
            theta(lo:hi) = max(0.45_dp, 1.0_dp - 0.25_dp * real(s - 1, dp)) * z
        end do
        theta = min(max(theta, lower), upper)
    end subroutine initial_theta

    subroutine rsdc_natural_from_theta(theta, method, n, k, p, natural)
        real(dp), intent(in) :: theta(:)
        integer, intent(in) :: method, n, k, p
        real(dp), allocatable, intent(out) :: natural(:)
        real(dp), allocatable :: r(:, :), rho(:), logits(:), e(:)
        integer :: c, nh, s, lo, hi, i, j, q
        real(dp) :: cmax
        c = lower_tri_size(k)
        nh = size(theta) - n * c
        allocate(natural(expected_parameter_count(method, n, k, p)))
        if (nh > 0) then
            if (method == rsdc_tvtp) then
                natural(1:nh) = theta(1:nh)
            else if (n == 2) then
                do i = 1, 2
                    natural(i) = logistic_local(theta(i))
                end do
            else
                allocate(logits(n), e(n))
                q = 0
                do i = 1, n
                    logits = 0.0_dp
                    do j = 1, n - 1
                        q = q + 1
                        logits(j) = theta(q)
                    end do
                    cmax = maxval(logits); e = exp(logits - cmax); e = e / sum(e)
                    natural((i - 1) * (n - 1) + 1:i * (n - 1)) = e(1:n - 1)
                end do
            end if
        end if
        allocate(r(k, k), rho(c))
        do s = 1, n
            lo = nh + (s - 1) * c + 1
            hi = nh + s * c
            call partial_to_correlation(theta(lo:hi), k, r)
            call matrix_to_correlations(r, rho)
            natural(nh + (s - 1) * c + 1:nh + s * c) = rho
        end do
    end subroutine rsdc_natural_from_theta

    subroutine rsdc_theta_from_natural(parameters, method, n, k, p, theta, ok)
        real(dp), intent(in) :: parameters(:)
        integer, intent(in) :: method, n, k, p
        real(dp), allocatable, intent(out) :: theta(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: beta(:, :), rho(:, :), pmat(:, :), lower(:), upper(:)
        real(dp), allocatable :: r(:, :), z(:)
        integer :: c, nh, i, j, q, s, lo, hi
        call unpack_natural_parameters(parameters, method, n, k, p, beta, rho, pmat, ok)
        if (.not. ok) return
        call make_theta_bounds(method, n, k, p, lower, upper)
        allocate(theta(size(lower))); theta = 0.0_dp
        c = lower_tri_size(k); nh = size(theta) - n * c
        select case (method)
        case (rsdc_nox)
            if (n == 2) then
                theta(1) = logit_local(pmat(1, 1))
                theta(2) = logit_local(pmat(2, 2))
            else
                q = 0
                do i = 1, n
                    do j = 1, n - 1
                        q = q + 1
                        theta(q) = log(max(pmat(i, j), tiny(1.0_dp)) / &
                                       max(pmat(i, n), tiny(1.0_dp)))
                    end do
                end do
            end if
        case (rsdc_tvtp)
            theta(1:nh) = reshape(transpose(beta), [nh])
        end select
        allocate(r(k, k), z(c))
        do s = 1, n
            call correlations_to_matrix(rho(s, :), k, r)
            call correlation_to_partial(r, z)
            lo = nh + (s - 1) * c + 1
            hi = nh + s * c
            theta(lo:hi) = z
        end do
        ok = all(abs(theta(nh + 1:)) < 1.0_dp)
    end subroutine rsdc_theta_from_natural

    subroutine rsdc_theta_from_model(model, theta)
        type(rsdc_model), intent(in) :: model
        real(dp), allocatable, intent(out) :: theta(:)
        real(dp), allocatable :: lower(:), upper(:), r(:, :), z(:)
        integer :: c, nh, s, lo, hi, i, j, q, n, k, p
        n = model%n_regimes; k = model%n_series; p = model%n_covariates; c = lower_tri_size(k)
        call make_theta_bounds(model%method, n, k, p, lower, upper)
        allocate(theta(size(lower))); theta = 0.0_dp; nh = size(theta) - n * c
        if (model%method == rsdc_nox) then
            if (n == 2) then
                theta(1) = logit_local(model%transition_matrix(1, 1))
                theta(2) = logit_local(model%transition_matrix(2, 2))
            else
                q = 0
                do i = 1, n
                    do j = 1, n - 1
                        q = q + 1
                        theta(q) = log(max(model%transition_matrix(i, j), tiny(1.0_dp)) / &
                                           max(model%transition_matrix(i, n), tiny(1.0_dp)))
                    end do
                end do
            end if
        else if (model%method == rsdc_tvtp) then
            theta(1:nh) = reshape(transpose(model%beta), [nh])
        end if
        allocate(r(k, k), z(c))
        do s = 1, n
            call correlations_to_matrix(model%correlations(s, :), k, r)
            call correlation_to_partial(r, z)
            lo = nh + (s - 1) * c + 1; hi = nh + s * c
            theta(lo:hi) = z
        end do
    end subroutine rsdc_theta_from_model

    subroutine local_pattern_search(theta, lower, upper, cfg, value)
        real(dp), intent(inout) :: theta(:)
        real(dp), intent(in) :: lower(:), upper(:)
        type(rsdc_control), intent(in) :: cfg
        real(dp), intent(out) :: value
        real(dp), allocatable :: trial(:), step(:)
        real(dp) :: candidate
        integer :: iter, j, sign
        logical :: improved
        allocate(trial(size(theta)), step(size(theta)))
        step = cfg%initial_step * (upper - lower)
        value = theta_objective(theta)
        do iter = 1, cfg%max_local_iterations
            improved = .false.
            do j = 1, size(theta)
                do sign = -1, 1, 2
                    trial = theta
                    trial(j) = min(max(theta(j) + real(sign, dp) * step(j), lower(j)), upper(j))
                    candidate = theta_objective(trial)
                    if (candidate < value) then
                        theta = trial; value = candidate; improved = .true.
                    end if
                end do
            end do
            if (.not. improved) step = 0.5_dp * step
            if (maxval(step) < cfg%tolerance) exit
        end do
    end subroutine local_pattern_search

    subroutine build_model(parameters, method, n, k, p, y, x, model, ok)
        real(dp), intent(in) :: parameters(:), y(:, :), x(:, :)
        integer, intent(in) :: method, n, k, p
        type(rsdc_model), intent(out) :: model
        logical, intent(out) :: ok
        real(dp), allocatable :: beta(:, :), rho(:, :), pmat(:, :)
        integer :: s
        call unpack_natural_parameters(parameters, method, n, k, p, beta, rho, pmat, ok)
        if (.not. ok) return
        model%method = method; model%n_regimes = n; model%n_series = k; model%n_covariates = p
        allocate(model%parameters(size(parameters))); model%parameters = parameters
        allocate(model%correlations(n, lower_tri_size(k))); model%correlations = rho
        allocate(model%covariance(k, k, n))
        do s = 1, n
            call correlations_to_matrix(rho(s, :), k, model%covariance(:, :, s))
        end do
        if (method == rsdc_tvtp) then
            allocate(model%beta(size(beta, 1), size(beta, 2))); model%beta = beta
            allocate(model%average_x(p)); model%average_x = mean_columns(x)
            allocate(model%transition_matrix(n, n))
            call transition_at(beta, model%average_x, model%transition_matrix)
            model%log_likelihood = -rsdc_negative_log_likelihood(parameters, y, method, n, x)
        else
            allocate(model%beta(0, 0), model%average_x(0))
            allocate(model%transition_matrix(size(pmat, 1), size(pmat, 2))); model%transition_matrix = pmat
            model%log_likelihood = -rsdc_negative_log_likelihood(parameters, y, method, n)
        end if
        call order_regimes(model)
        ok = model%log_likelihood > -1.0e9_dp
    end subroutine build_model

    subroutine order_regimes(model)
        type(rsdc_model), intent(inout) :: model
        integer, allocatable :: ord(:)
        real(dp), allocatable :: score(:), rho(:, :), cov(:, :, :), pmat(:, :), beta(:, :), natural(:)
        integer :: i, j, tmp, n, c, p, ref, row, col, lo, hi
        n = model%n_regimes
        if (n <= 1) return
        c = size(model%correlations, 2)
        allocate(ord(n), score(n))
        do i = 1, n
            ord(i) = i; score(i) = sum(model%correlations(i, :)) / real(c, dp)
        end do
        do i = 1, n - 1
            do j = i + 1, n
                if (score(ord(j)) < score(ord(i))) then
                    tmp = ord(i); ord(i) = ord(j); ord(j) = tmp
                end if
            end do
        end do
        if (all(ord == [(i, i=1,n)])) return
        allocate(rho(n, c), cov(model%n_series, model%n_series, n), pmat(n, n))
        rho = model%correlations(ord, :); cov = model%covariance(:, :, ord)
        pmat = model%transition_matrix(ord, ord)
        model%correlations = rho; model%covariance = cov; model%transition_matrix = pmat
        if (model%method == rsdc_tvtp) then
            p = model%n_covariates
            allocate(beta(n, size(model%beta, 2)))
            if (n == 2) then
                beta = model%beta(ord, :)
            else
                ref = ord(n)
                do row = 1, n
                    do col = 1, n - 1
                        lo = (col - 1) * p + 1; hi = col * p
                        beta(row, lo:hi) = beta_column(model%beta, ord(row), ord(col), n, p) - &
                                           beta_column(model%beta, ord(row), ref, n, p)
                    end do
                end do
            end if
            model%beta = beta
        end if
        call pack_model(model, natural)
        model%parameters = natural
    end subroutine order_regimes

    function beta_column(beta, row, col, n, p) result(v)
        real(dp), intent(in) :: beta(:, :)
        integer, intent(in) :: row, col, n, p
        real(dp) :: v(p)
        if (col == n) then
            v = 0.0_dp
        else
            v = beta(row, (col - 1) * p + 1:col * p)
        end if
    end function beta_column

    subroutine pack_model(model, parameters)
        type(rsdc_model), intent(in) :: model
        real(dp), allocatable, intent(out) :: parameters(:)
        integer :: n, k, p, c, nh, i
        n = model%n_regimes; k = model%n_series; p = model%n_covariates; c = lower_tri_size(k)
        allocate(parameters(expected_parameter_count(model%method, n, k, p)))
        select case (model%method)
        case (rsdc_const)
            parameters = model%correlations(1, :)
        case (rsdc_nox)
            nh = n * (n - 1)
            if (n == 2) then
                parameters(1:2) = [model%transition_matrix(1, 1), model%transition_matrix(2, 2)]
            else
                do i = 1, n
                    parameters((i - 1) * (n - 1) + 1:i * (n - 1)) = model%transition_matrix(i, 1:n - 1)
                end do
            end if
            parameters(nh + 1:) = reshape(transpose(model%correlations), [n * c])
        case (rsdc_tvtp)
            nh = size(model%beta)
            parameters(1:nh) = reshape(transpose(model%beta), [nh])
            parameters(nh + 1:) = reshape(transpose(model%correlations), [n * c])
        end select
    end subroutine pack_model

    subroutine compute_oos_likelihood(model, yis, yoos, xis, xoos)
        type(rsdc_model), intent(inout) :: model
        real(dp), intent(in) :: yis(:, :), yoos(:, :), xis(:, :), xoos(:, :)
        type(rsdc_filter_result) :: a, b
        real(dp), allocatable :: xi(:)
        if (size(yoos, 1) == 0) return
        select case (model%method)
        case (rsdc_const)
            model%log_likelihood_oos = -rsdc_negative_log_likelihood(model%parameters, yoos, &
                rsdc_const, 1)
        case (rsdc_nox)
            call rsdc_hamilton(yis, model%correlations, a, pmat=model%transition_matrix)
            if (.not. a%ok) return
            xi = a%filtered(:, size(yis, 1))
            call rsdc_hamilton(yoos, model%correlations, b, pmat=model%transition_matrix, xi_init=xi)
            if (b%ok) model%log_likelihood_oos = b%log_likelihood
        case (rsdc_tvtp)
            call rsdc_hamilton(yis, model%correlations, a, x=xis, beta=model%beta)
            if (.not. a%ok) return
            xi = a%filtered(:, size(yis, 1))
            call rsdc_hamilton(yoos, model%correlations, b, x=xoos, beta=model%beta, xi_init=xi)
            if (b%ok) model%log_likelihood_oos = b%log_likelihood
        end select
    end subroutine compute_oos_likelihood

    subroutine estimate_vcov(model, y, x)
        type(rsdc_model), intent(inout) :: model
        real(dp), intent(in) :: y(:, :), x(:, :)
        real(dp), allocatable :: hess(:, :), inv(:, :), p0(:), pp(:), pm(:), ppp(:), ppm(:), pmp(:), pmm(:)
        real(dp), allocatable :: h(:)
        real(dp) :: f0, fp, fm
        integer :: i, j, np
        logical :: ok
        np = size(model%parameters)
        allocate(hess(np, np), inv(np, np), p0(np), pp(np), pm(np), ppp(np), ppm(np), pmp(np), pmm(np), h(np))
        p0 = model%parameters; h = 1.0e-4_dp * (1.0_dp + abs(p0))
        f0 = natural_objective(p0)
        hess = 0.0_dp
        do i = 1, np
            pp = p0; pm = p0; pp(i) = pp(i) + h(i); pm(i) = pm(i) - h(i)
            fp = natural_objective(pp); fm = natural_objective(pm)
            hess(i, i) = (fp - 2.0_dp * f0 + fm) / h(i) ** 2
            do j = i + 1, np
                ppp = p0; ppm = p0; pmp = p0; pmm = p0
                ppp(i)=ppp(i)+h(i); ppp(j)=ppp(j)+h(j)
                ppm(i)=ppm(i)+h(i); ppm(j)=ppm(j)-h(j)
                pmp(i)=pmp(i)-h(i); pmp(j)=pmp(j)+h(j)
                pmm(i)=pmm(i)-h(i); pmm(j)=pmm(j)-h(j)
                hess(i, j) = (natural_objective(ppp)-natural_objective(ppm)- &
                    natural_objective(pmp)+natural_objective(pmm))/(4.0_dp*h(i)*h(j))
                hess(j, i) = hess(i, j)
            end do
        end do
        call inverse_matrix(hess, inv, ok)
        if (ok) then
            allocate(model%vcov(np, np), model%standard_errors(np))
            model%vcov = inv
            do i = 1, np
                model%standard_errors(i) = sqrt(max(inv(i, i), 0.0_dp))
            end do
        end if
    contains
        real(dp) function natural_objective(par) result(v)
            real(dp), intent(in) :: par(:)
            if (model%method == rsdc_tvtp) then
                v = rsdc_negative_log_likelihood(par, y, model%method, model%n_regimes, x)
            else
                v = rsdc_negative_log_likelihood(par, y, model%method, model%n_regimes)
            end if
        end function natural_objective
    end subroutine estimate_vcov

    pure real(dp) function logistic_local(x) result(p)
        real(dp), intent(in) :: x
        if (x >= 0.0_dp) then
            p = 1.0_dp / (1.0_dp + exp(-x))
        else
            p = exp(x) / (1.0_dp + exp(x))
        end if
    end function logistic_local

    pure real(dp) function logit_local(p) result(x)
        real(dp), intent(in) :: p
        real(dp) :: q
        q = max(min(p, 1.0_dp - 1.0e-12_dp), 1.0e-12_dp)
        x = log(q / (1.0_dp - q))
    end function logit_local
end module rsdc_estimation
