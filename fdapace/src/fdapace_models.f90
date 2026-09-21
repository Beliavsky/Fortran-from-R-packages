module fdapace_models
    use fdapace_covariance, only : get_cr_cov_yz
    use fdapace_fpca, only : fpca_dense
    use fdapace_kinds, only : dp
    use fdapace_math, only : linear_interp, mean_value
    use fdapace_numerics, only : empirical_quantile, inverse_matrix, least_squares, logistic_regression, &
        sample_indices_with_replacement
    use fdapace_smoothing, only : lwls1d, trapz_rcpp
    use fdapace_types, only : fam_result, fcreg_result, flm_ci_result, flm_result, fopt_result, &
        fpc_quantile_result, fpca_result, multifam_result, sbf_result, tvam_result, vcam_result
    implicit none
    private

    public :: fam
    public :: fc_reg
    public :: flm_ci_scalar_dense
    public :: flm_scalar_dense
    public :: f_opt_des
    public :: fpc_quantile
    public :: multi_fam
    public :: sb_fitting
    public :: tvam
    public :: vcam

contains

    subroutine sb_fitting(y, x_eval, x_obs, result, h, support, max_iter, tolerance)
        real(dp), intent(in) :: y(:) !! Scalar responses for additive backfitting.
        real(dp), intent(in) :: x_eval(:,:) !! Evaluation points with one column per additive component.
        real(dp), intent(in) :: x_obs(:,:) !! Observed multivariate predictors with subjects in rows.
        type(sbf_result), intent(out) :: result !! Additive component estimates, marginal NW fits, mean response, and convergence metadata.
        real(dp), intent(in), optional :: h(:) !! Positive component bandwidths; defaults to 0.25*n^(-1/5) times each support range.
        real(dp), intent(in), optional :: support(:,:) !! Component supports as d-by-2 bounds; defaults to observed minima and maxima.
        integer, intent(in), optional :: max_iter !! Maximum backfitting iterations; defaults to 50.
        real(dp), intent(in), optional :: tolerance !! Maximum RMS component-change tolerance; defaults to 5e-5.
        real(dp), allocatable :: bandwidth(:)
        real(dp), allocatable :: current(:,:)
        real(dp), allocatable :: old(:,:)
        real(dp), allocatable :: residual(:)
        real(dp), allocatable :: supp(:,:)
        real(dp) :: change
        real(dp) :: tol
        integer :: d
        integer :: i
        integer :: iter
        integer :: itmax
        integer :: j
        integer :: l
        integer :: n
        integer :: ne

        n = size(x_obs, 1)
        d = size(x_obs, 2)
        ne = size(x_eval, 1)
        if (size(y) /= n .or. size(x_eval, 2) /= d .or. n < 2 .or. d < 1 .or. ne < 1) then
            error stop "sb_fitting: incompatible dimensions"
        end if
        allocate(supp(d, 2), bandwidth(d))
        if (present(support)) then
            if (size(support, 1) /= d .or. size(support, 2) /= 2) error stop "sb_fitting: support shape mismatch"
            supp = support
        else
            do j = 1, d
                supp(j, 1) = minval(x_obs(:, j))
                supp(j, 2) = maxval(x_obs(:, j))
            end do
        end if
        if (present(h)) then
            if (size(h) /= d .or. any(h <= 0.0_dp)) error stop "sb_fitting: h must contain d positive bandwidths"
            bandwidth = h
        else
            do j = 1, d
                bandwidth(j) = 0.25_dp * real(n, dp)**(-0.2_dp) * max(supp(j, 2) - supp(j, 1), epsilon(1.0_dp))
            end do
        end if
        itmax = 50
        if (present(max_iter)) itmax = max_iter
        tol = 5.0e-5_dp
        if (present(tolerance)) tol = tolerance
        if (itmax < 1 .or. tol <= 0.0_dp) error stop "sb_fitting: invalid convergence controls"
        result%mean_y = mean_value(y)
        allocate(current(ne, d), old(ne, d), result%nw(ne, d), residual(n))
        current = 0.0_dp
        do j = 1, d
            result%nw(:, j) = smooth_unordered(bandwidth(j), "epan", x_obs(:, j), y - result%mean_y, &
                x_eval(:, j), 0, 0)
        end do
        do iter = 1, itmax
            old = current
            do j = 1, d
                residual = y - result%mean_y
                do l = 1, d
                    if (l == j) cycle
                    do i = 1, n
                        residual(i) = residual(i) - linear_interp(x_eval(:, l), current(:, l), x_obs(i, l), .true.)
                    end do
                end do
                current(:, j) = smooth_unordered(bandwidth(j), "epan", x_obs(:, j), residual, x_eval(:, j), 0, 0)
                current(:, j) = current(:, j) - sum(current(:, j)) / real(ne, dp)
            end do
            change = 0.0_dp
            do j = 1, d
                change = max(change, sqrt(sum((current(:, j) - old(:, j))**2) / real(ne, dp)))
            end do
            if (change <= tol) exit
        end do
        result%fit = current
        result%iterations = iter
        result%error = change
        result%converged = change <= tol
    end subroutine sb_fitting

    subroutine fam(y, x, t, result, n_eval, alpha, support, fve_threshold)
        real(dp), intent(in) :: y(:) !! Scalar responses paired with rows of the functional predictor x.
        real(dp), intent(in) :: x(:,:) !! Dense regular functional predictor with subjects in rows.
        real(dp), intent(in) :: t(:) !! Common regular support for columns of x.
        type(fam_result), intent(out) :: result !! FPC-score additive component curves, score grids, bandwidths, and FPCA basis metadata.
        integer, intent(in), optional :: n_eval !! Number of evaluation points per component; defaults to 51.
        real(dp), intent(in), optional :: alpha !! Positive bandwidth shrinkage factor; defaults to 0.7.
        real(dp), intent(in), optional :: support(:) !! Two-element studentized-score support; defaults to [-2,2].
        real(dp), intent(in), optional :: fve_threshold !! FPCA FVE threshold; defaults to 0.99.
        type(fpca_result) :: fit
        real(dp), allocatable :: score_std(:,:)
        real(dp) :: a
        real(dp) :: lo
        real(dp) :: hi
        real(dp) :: fve
        real(dp) :: hstd
        integer :: i
        integer :: j
        integer :: ne
        integer :: nk

        if (size(y) /= size(x, 1)) error stop "fam: response length must match curves"
        ne = 51
        if (present(n_eval)) ne = n_eval
        if (ne < 2) error stop "fam: n_eval must be at least two"
        a = 0.7_dp
        if (present(alpha)) a = alpha
        if (a <= 0.0_dp) error stop "fam: alpha must be positive"
        lo = -2.0_dp
        hi = 2.0_dp
        if (present(support)) then
            if (size(support) /= 2 .or. support(2) <= support(1)) error stop "fam: support must contain increasing bounds"
            lo = support(1)
            hi = support(2)
        end if
        fve = 0.99_dp
        if (present(fve_threshold)) fve = fve_threshold
        call fpca_dense(x, t, fit, fve_threshold=fve, assume_error=.false.)
        nk = size(fit%lambda)
        allocate(score_std(size(x, 1), nk), result%fam(ne, nk), result%xi(ne, nk), result%bw(nk))
        score_std = fit%xi_est
        do j = 1, nk
            score_std(:, j) = score_std(:, j) / sqrt(max(fit%lambda(j), epsilon(1.0_dp)))
        end do
        result%mu = mean_value(y)
        hstd = a * 0.25_dp * real(size(y), dp)**(-0.2_dp) * (hi - lo)
        do j = 1, nk
            do i = 1, ne
                result%xi(i, j) = (lo + real(i - 1, dp) * (hi - lo) / real(ne - 1, dp)) * sqrt(fit%lambda(j))
            end do
            result%fam(:, j) = smooth_unordered(hstd, "epan", score_std(:, j), y - result%mu, &
                result%xi(:, j) / sqrt(fit%lambda(j)), 1, 0)
            result%fam(:, j) = result%fam(:, j) - sum(result%fam(:, j)) / real(ne, dp)
            result%bw(j) = hstd * sqrt(fit%lambda(j))
        end do
        result%lambda = fit%lambda
        result%phi = fit%phi
        result%work_grid = fit%work_grid
    end subroutine fam

    subroutine multi_fam(y, x, t, result, n_eval, alpha, support, fve_threshold)
        real(dp), intent(in) :: y(:) !! Scalar responses paired with rows of all functional predictors.
        real(dp), intent(in) :: x(:,:,:) !! Dense functional predictors shaped subjects-by-grid-by-predictor.
        real(dp), intent(in) :: t(:) !! Common regular support used by every functional predictor.
        type(multifam_result), intent(out) :: result !! Smooth-backfitting component curves on concatenated FPC-score grids and score metadata.
        integer, intent(in), optional :: n_eval !! Common number of evaluation points per FPC component; defaults to 51.
        real(dp), intent(in), optional :: alpha !! Positive bandwidth shrinkage factor; defaults to 0.7.
        real(dp), intent(in), optional :: support(:) !! Two-element studentized-score support; defaults to [-2,2].
        real(dp), intent(in), optional :: fve_threshold !! FVE threshold used for each predictor FPCA; defaults to 0.95.
        type(fpca_result), allocatable :: fits(:)
        type(sbf_result) :: sbf
        real(dp), allocatable :: h(:)
        real(dp), allocatable :: scores(:,:)
        real(dp), allocatable :: score_grid(:,:)
        real(dp) :: a
        real(dp) :: fve
        real(dp) :: hi
        real(dp) :: lo
        integer :: d
        integer :: i
        integer :: j
        integer :: k
        integer :: ne
        integer :: offset
        integer :: total_k

        if (size(y) /= size(x, 1) .or. size(x, 2) /= size(t)) error stop "multi_fam: incompatible dimensions"
        d = size(x, 3)
        if (d < 1) error stop "multi_fam: at least one functional predictor is required"
        ne = 51
        if (present(n_eval)) ne = n_eval
        if (ne < 2) error stop "multi_fam: n_eval must be at least two"
        a = 0.7_dp
        if (present(alpha)) a = alpha
        lo = -2.0_dp
        hi = 2.0_dp
        if (present(support)) then
            if (size(support) /= 2 .or. support(2) <= support(1)) error stop "multi_fam: invalid support"
            lo = support(1)
            hi = support(2)
        end if
        fve = 0.95_dp
        if (present(fve_threshold)) fve = fve_threshold
        allocate(fits(d))
        total_k = 0
        do j = 1, d
            call fpca_dense(x(:, :, j), t, fits(j), fve_threshold=fve, assume_error=.false.)
            total_k = total_k + size(fits(j)%lambda)
        end do
        if (total_k > size(y)) error stop "multi_fam: selected FPC score count exceeds sample size"
        allocate(scores(size(y), total_k), score_grid(ne, total_k), result%lambda(total_k), h(total_k))
        offset = 0
        do j = 1, d
            do k = 1, size(fits(j)%lambda)
                scores(:, offset + k) = fits(j)%xi_est(:, k) / sqrt(max(fits(j)%lambda(k), epsilon(1.0_dp)))
                result%lambda(offset + k) = fits(j)%lambda(k)
                do i = 1, ne
                    score_grid(i, offset + k) = lo + real(i - 1, dp) * (hi - lo) / real(ne - 1, dp)
                end do
                h(offset + k) = a * 0.25_dp * real(size(y), dp)**(-0.2_dp) * (hi - lo)
            end do
            offset = offset + size(fits(j)%lambda)
        end do
        call sb_fitting(y, score_grid, scores, sbf, h=h)
        result%mu = sbf%mean_y
        result%sbfit = sbf%fit
        allocate(result%xi(ne, total_k), result%bw(total_k))
        do k = 1, total_k
            result%xi(:, k) = score_grid(:, k) * sqrt(result%lambda(k))
            result%bw(k) = h(k) * sqrt(result%lambda(k))
        end do
    end subroutine multi_fam

    subroutine flm_scalar_dense(y, x, t, result, x_test, n_perm, fve_threshold)
        real(dp), intent(in) :: y(:) !! Scalar response vector for functional linear regression.
        real(dp), intent(in) :: x(:,:,:) !! Dense functional predictors shaped subjects-by-grid-by-predictor.
        real(dp), intent(in) :: t(:) !! Common regular support used by all functional predictors.
        type(flm_result), intent(out) :: result !! Intercept, coefficient functions, fitted/predicted responses, R2, and optional permutation p-value.
        real(dp), intent(in), optional :: x_test(:,:,:) !! Optional held-out functional predictors on the same grid and predictor set.
        integer, intent(in), optional :: n_perm !! Optional positive number of response permutations for a global R2 test.
        real(dp), intent(in), optional :: fve_threshold !! FPCA FVE threshold used for every predictor; defaults to 0.95.
        type(fpca_result), allocatable :: fits(:)
        real(dp), allocatable :: beta_score(:,:)
        real(dp), allocatable :: design(:,:)
        real(dp), allocatable :: design_test(:,:)
        real(dp), allocatable :: perm_y(:,:)
        real(dp), allocatable :: pred(:,:)
        real(dp), allocatable :: rhs(:,:)
        real(dp), allocatable :: score_test(:,:)
        real(dp), allocatable :: scores(:,:)
        real(dp) :: fve
        real(dp) :: r2_perm
        real(dp) :: sse
        real(dp) :: sst
        integer, allocatable :: idx(:)
        integer :: b
        integer :: d
        integer :: info
        integer :: j
        integer :: k
        integer :: nb
        integer :: offset
        integer :: total_k

        if (size(y) /= size(x, 1) .or. size(x, 2) /= size(t)) error stop "flm_scalar_dense: incompatible dimensions"
        d = size(x, 3)
        fve = 0.95_dp
        if (present(fve_threshold)) fve = fve_threshold
        allocate(fits(d))
        total_k = 0
        do j = 1, d
            call fpca_dense(x(:, :, j), t, fits(j), fve_threshold=fve, assume_error=.false.)
            total_k = total_k + size(fits(j)%lambda)
        end do
        if (total_k >= size(y)) error stop "flm_scalar_dense: selected score dimension must be smaller than sample size"
        allocate(scores(size(y), total_k), design(size(y), total_k + 1), rhs(size(y), 1))
        offset = 0
        do j = 1, d
            k = size(fits(j)%lambda)
            scores(:, offset + 1:offset + k) = fits(j)%xi_est
            offset = offset + k
        end do
        design(:, 1) = 1.0_dp
        design(:, 2:) = scores
        rhs(:, 1) = y
        call least_squares(design, rhs, beta_score, info, ridge=1.0e-10_dp)
        if (info /= 0) error stop "flm_scalar_dense: least-squares fit failed"
        allocate(result%alpha(1), result%beta(size(t), d), result%y_hat(size(y), 1))
        result%alpha(1) = beta_score(1, 1)
        result%beta = 0.0_dp
        offset = 0
        do j = 1, d
            k = size(fits(j)%lambda)
            result%beta(:, j) = matmul(fits(j)%phi, beta_score(offset + 2:offset + k + 1, 1))
            offset = offset + k
        end do
        result%y_hat = matmul(design, beta_score)
        sst = sum((y - mean_value(y))**2)
        sse = sum((y - result%y_hat(:, 1))**2)
        if (sst > 0.0_dp) result%r2 = max(0.0_dp, min(1.0_dp, 1.0_dp - sse / sst))
        if (present(x_test)) then
            if (size(x_test, 2) /= size(t) .or. size(x_test, 3) /= d) error stop "flm_scalar_dense: x_test shape mismatch"
            allocate(score_test(size(x_test, 1), total_k), design_test(size(x_test, 1), total_k + 1))
            offset = 0
            do j = 1, d
                k = size(fits(j)%lambda)
                call project_scores(x_test(:, :, j), t, fits(j), score_test(:, offset + 1:offset + k))
                offset = offset + k
            end do
            design_test(:, 1) = 1.0_dp
            design_test(:, 2:) = score_test
            result%y_pred = matmul(design_test, beta_score)
        else
            result%y_pred = result%y_hat
        end if
        if (present(n_perm)) then
            nb = n_perm
            if (nb > 0) then
                allocate(idx(size(y)), perm_y(size(y), 1))
                result%p_value = 0.0_dp
                do b = 1, nb
                    call sample_permutation(size(y), idx)
                    perm_y(:, 1) = y(idx)
                    call least_squares(design, perm_y, pred, info, ridge=1.0e-10_dp)
                    if (info /= 0) cycle
                    sse = sum((perm_y(:, 1) - matmul(design, pred(:, 1)))**2)
                    sst = sum((perm_y(:, 1) - mean_value(perm_y(:, 1)))**2)
                    r2_perm = 0.0_dp
                    if (sst > 0.0_dp) r2_perm = 1.0_dp - sse / sst
                    if (r2_perm >= result%r2) result%p_value = result%p_value + 1.0_dp
                end do
                result%p_value = result%p_value / real(nb, dp)
                result%has_p_value = .true.
            end if
        end if
    end subroutine flm_scalar_dense

    subroutine flm_ci_scalar_dense(y, x, t, result, level, n_boot, fve_threshold)
        real(dp), intent(in) :: y(:) !! Scalar response vector for bootstrap functional linear-regression confidence intervals.
        real(dp), intent(in) :: x(:,:,:) !! Dense functional predictors shaped subjects-by-grid-by-predictor.
        real(dp), intent(in) :: t(:) !! Common regular support for columns of x.
        type(flm_ci_result), intent(out) :: result !! Pointwise percentile intervals for intercept and functional coefficients.
        real(dp), intent(in), optional :: level !! Confidence level in [0,1]; defaults to 0.95.
        integer, intent(in), optional :: n_boot !! Positive bootstrap replicate count; defaults to 999.
        real(dp), intent(in), optional :: fve_threshold !! FPCA FVE threshold used by each bootstrap regression; defaults to 0.95.
        type(flm_result) :: fit
        real(dp), allocatable :: alpha_boot(:)
        real(dp), allocatable :: beta_boot(:,:,:)
        real(dp), allocatable :: xb(:,:,:)
        real(dp) :: fve
        real(dp) :: lev
        integer, allocatable :: idx(:)
        integer :: b
        integer :: j
        integer :: k
        integer :: nb

        lev = 0.95_dp
        if (present(level)) lev = level
        if (lev < 0.0_dp .or. lev > 1.0_dp) error stop "flm_ci_scalar_dense: level must lie in [0,1]"
        nb = 999
        if (present(n_boot)) nb = n_boot
        if (nb < 1) error stop "flm_ci_scalar_dense: n_boot must be positive"
        fve = 0.95_dp
        if (present(fve_threshold)) fve = fve_threshold
        allocate(alpha_boot(nb), beta_boot(nb, size(t), size(x, 3)), idx(size(y)))
        allocate(xb(size(x, 1), size(x, 2), size(x, 3)))
        do b = 1, nb
            call sample_indices_with_replacement(size(y), idx)
            xb = x(idx, :, :)
            call flm_scalar_dense(y(idx), xb, t, fit, fve_threshold=fve)
            alpha_boot(b) = fit%alpha(1)
            beta_boot(b, :, :) = fit%beta
        end do
        allocate(result%alpha_lower(1), result%alpha_upper(1))
        allocate(result%beta_lower(size(t), size(x, 3)), result%beta_upper(size(t), size(x, 3)))
        result%alpha_lower(1) = empirical_quantile(alpha_boot, 0.5_dp * (1.0_dp - lev))
        result%alpha_upper(1) = empirical_quantile(alpha_boot, 1.0_dp - 0.5_dp * (1.0_dp - lev))
        do k = 1, size(x, 3)
            do j = 1, size(t)
                result%beta_lower(j, k) = empirical_quantile(beta_boot(:, j, k), 0.5_dp * (1.0_dp - lev))
                result%beta_upper(j, k) = empirical_quantile(beta_boot(:, j, k), 1.0_dp - 0.5_dp * (1.0_dp - lev))
            end do
        end do
        result%level = lev
    end subroutine flm_ci_scalar_dense

    subroutine fc_reg(x, y, t, result, ridge)
        real(dp), intent(in) :: x(:,:,:) !! Dense concurrent functional covariates shaped subjects-by-grid-by-covariate.
        real(dp), intent(in) :: y(:,:) !! Dense functional response shaped subjects-by-grid.
        real(dp), intent(in) :: t(:) !! Common output grid corresponding to the second dimensions of x and y.
        type(fcreg_result), intent(out) :: result !! Pointwise concurrent regression intercepts, slopes, R2 values, and grid.
        real(dp), intent(in), optional :: ridge !! Optional nonnegative pointwise ridge penalty; defaults to 1e-10.
        real(dp), allocatable :: beta(:,:)
        real(dp), allocatable :: design(:,:)
        real(dp), allocatable :: rhs(:,:)
        real(dp) :: r
        real(dp) :: sse
        real(dp) :: sst
        integer :: info
        integer :: j

        if (size(x, 1) /= size(y, 1) .or. size(x, 2) /= size(y, 2) .or. size(y, 2) /= size(t)) then
            error stop "fc_reg: incompatible dense dimensions"
        end if
        r = 1.0e-10_dp
        if (present(ridge)) r = ridge
        allocate(result%beta(size(x, 3), size(t)), result%beta0(size(t)), result%r2(size(t)))
        allocate(design(size(y, 1), size(x, 3) + 1), rhs(size(y, 1), 1))
        design(:, 1) = 1.0_dp
        do j = 1, size(t)
            design(:, 2:) = x(:, j, :)
            rhs(:, 1) = y(:, j)
            call least_squares(design, rhs, beta, info, ridge=r)
            if (info /= 0) error stop "fc_reg: pointwise least-squares fit failed"
            result%beta0(j) = beta(1, 1)
            result%beta(:, j) = beta(2:, 1)
            sse = sum((y(:, j) - matmul(design, beta(:, 1)))**2)
            sst = sum((y(:, j) - mean_value(y(:, j)))**2)
            result%r2(j) = 0.0_dp
            if (sst > 0.0_dp) result%r2(j) = max(0.0_dp, min(1.0_dp, 1.0_dp - sse / sst))
        end do
        result%out_grid = t
    end subroutine fc_reg

    subroutine f_opt_des(y, t, p, ridge, result, response)
        real(dp), intent(in) :: y(:,:) !! Dense regular functional sample used to estimate covariance for optimal design selection.
        real(dp), intent(in) :: t(:) !! Common regular support whose entries are candidate design locations.
        integer, intent(in) :: p !! Positive number of design points selected greedily.
        real(dp), intent(in) :: ridge !! Positive covariance ridge used to regularize selected design submatrices.
        type(fopt_result), intent(out) :: result !! Selected grid indices/locations, coefficient of determination, adjustment, and ridge.
        real(dp), intent(in), optional :: response(:) !! Optional scalar response for prediction design; omission selects for trajectory recovery.
        type(fpca_result) :: fit
        integer, allocatable :: chosen(:)
        logical, allocatable :: used(:)
        real(dp), allocatable :: cross_cov(:)
        real(dp) :: best
        real(dp) :: criterion
        real(dp) :: total_var
        integer :: best_idx
        integer :: c
        integer :: n
        integer :: step

        n = size(y, 1)
        if (size(y, 2) /= size(t) .or. p < 1 .or. p > size(t) .or. ridge <= 0.0_dp) then
            error stop "f_opt_des: invalid inputs"
        end if
        if (present(response)) then
            if (size(response) /= n) error stop "f_opt_des: response length mismatch"
            cross_cov = get_cr_cov_yz(y, response)
        end if
        call fpca_dense(y, t, fit, assume_error=.false.)
        allocate(chosen(p), used(size(t)))
        used = .false.
        do step = 1, p
            best = -huge(1.0_dp)
            best_idx = 0
            do c = 1, size(t)
                if (used(c)) cycle
                chosen(step) = c
                if (present(response)) then
                    criterion = scalar_design_score(fit%fitted_cov, cross_cov, chosen(1:step), ridge)
                else
                    criterion = trajectory_design_score(fit%fitted_cov, chosen(1:step), ridge)
                end if
                if (criterion > best) then
                    best = criterion
                    best_idx = c
                end if
            end do
            chosen(step) = best_idx
            used(best_idx) = .true.
        end do
        result%indices = chosen
        result%opt_des = t(chosen)
        result%ridge = ridge
        if (present(response)) then
            total_var = sum((response - mean_value(response))**2) / real(max(1, n - 1), dp)
            result%r2 = 0.0_dp
            if (total_var > 0.0_dp) result%r2 = scalar_design_score(fit%fitted_cov, cross_cov, chosen, ridge) / total_var
        else
            total_var = sum(fit%lambda)
            result%r2 = 0.0_dp
            if (total_var > 0.0_dp) then
                result%r2 = trajectory_design_score(fit%fitted_cov, chosen, ridge) * &
                    (t(size(t)) - t(1)) / real(max(1, size(t) - 1), dp) / total_var
            end if
        end if
        result%r2 = max(0.0_dp, min(1.0_dp, result%r2))
        if (n > p + 1) then
            result%r2_adj = 1.0_dp - (1.0_dp - result%r2) * real(n - 1, dp) / real(n - p - 1, dp)
        else
            result%r2_adj = result%r2
        end if
    end subroutine f_opt_des

    real(dp) function trajectory_design_score(cov, indices, ridge) result(value)
        real(dp), intent(in) :: cov(:,:) !! Symmetric fitted covariance matrix on the candidate design grid.
        integer, intent(in) :: indices(:) !! One-based selected design indices.
        real(dp), intent(in) :: ridge !! Nonnegative diagonal ridge for the selected covariance matrix.
        real(dp), allocatable :: inv(:,:)
        real(dp), allocatable :: sub(:,:)
        real(dp), allocatable :: cross(:,:)
        integer :: i
        integer :: info
        integer :: j

        allocate(sub(size(indices), size(indices)), cross(size(cov, 1), size(indices)))
        do j = 1, size(indices)
            cross(:, j) = cov(:, indices(j))
            do i = 1, size(indices)
                sub(i, j) = cov(indices(i), indices(j))
            end do
        end do
        call inverse_matrix(sub, inv, info, ridge=ridge)
        if (info /= 0) then
            value = -huge(1.0_dp)
        else
            value = sum(cross * transpose(matmul(inv, transpose(cross))))
        end if
    end function trajectory_design_score

    real(dp) function scalar_design_score(cov, cross_cov, indices, ridge) result(value)
        real(dp), intent(in) :: cov(:,:) !! Symmetric fitted covariance matrix on the candidate design grid.
        real(dp), intent(in) :: cross_cov(:) !! Functional-scalar cross-covariance vector on the same grid.
        integer, intent(in) :: indices(:) !! One-based selected design indices.
        real(dp), intent(in) :: ridge !! Nonnegative diagonal ridge for the selected covariance matrix.
        real(dp), allocatable :: inv(:,:)
        real(dp), allocatable :: sub(:,:)
        real(dp), allocatable :: c(:)
        integer :: i
        integer :: info
        integer :: j

        allocate(sub(size(indices), size(indices)), c(size(indices)))
        do j = 1, size(indices)
            c(j) = cross_cov(indices(j))
            do i = 1, size(indices)
                sub(i, j) = cov(indices(i), indices(j))
            end do
        end do
        call inverse_matrix(sub, inv, info, ridge=ridge)
        if (info /= 0) then
            value = -huge(1.0_dp)
        else
            value = dot_product(c, matmul(inv, c))
        end if
    end function scalar_design_score

    subroutine fpc_quantile(x, t, y, out_q, result, fve_threshold)
        real(dp), intent(in) :: x(:,:) !! Dense functional covariate with subjects in rows.
        real(dp), intent(in) :: t(:) !! Common regular support for columns of x.
        real(dp), intent(in) :: y(:) !! Scalar response paired with rows of x.
        real(dp), intent(in) :: out_q(:) !! Desired conditional quantile probabilities in [0,1].
        type(fpc_quantile_result), intent(out) :: result !! Conditional quantiles, CDF values, logistic coefficient curves, and evaluation grids.
        real(dp), intent(in), optional :: fve_threshold !! FPCA FVE threshold for the functional covariate; defaults to 0.95.
        type(fpca_result) :: fit
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: design(:,:)
        real(dp), allocatable :: raw_cdf(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: ybin(:)
        real(dp), allocatable :: zgrid(:)
        real(dp) :: eta
        real(dp) :: fve
        real(dp) :: lower
        real(dp) :: upper
        real(dp) :: bw
        integer :: i
        integer :: info
        integer :: j
        integer :: k
        integer, parameter :: m = 100
        integer, parameter :: zm = 50

        if (size(y) /= size(x, 1) .or. size(x, 2) /= size(t)) error stop "fpc_quantile: incompatible dimensions"
        if (any(out_q < 0.0_dp) .or. any(out_q > 1.0_dp)) error stop "fpc_quantile: quantiles must lie in [0,1]"
        fve = 0.95_dp
        if (present(fve_threshold)) fve = fve_threshold
        call fpca_dense(x, t, fit, fve_threshold=fve, assume_error=.false.)
        lower = empirical_quantile(y, min(0.05_dp, 5.0_dp / real(max(10, size(y)), dp)))
        upper = empirical_quantile(y, max(0.95_dp, 1.0_dp - 5.0_dp / real(max(10, size(y)), dp)))
        if (upper <= lower) then
            lower = minval(y)
            upper = maxval(y)
        end if
        allocate(zgrid(zm), result%beta(zm, size(fit%xi_est, 2) + 1))
        allocate(design(size(y), size(fit%xi_est, 2) + 1), ybin(size(y)))
        design(:, 1) = 1.0_dp
        design(:, 2:) = fit%xi_est
        do j = 1, zm
            zgrid(j) = lower + real(j - 1, dp) * (upper - lower) / real(zm - 1, dp)
            ybin = merge(1.0_dp, 0.0_dp, y <= zgrid(j))
            call logistic_regression(design, ybin, beta, info, ridge=1.0e-6_dp)
            if (info < 0) error stop "fpc_quantile: logistic fit failed"
            result%beta(j, :) = beta
        end do
        allocate(result%cdf_grid(m), result%pred_cdf(size(y), m), weights(zm), raw_cdf(zm))
        allocate(result%pred_quantile(size(y), size(out_q)), result%quantiles(size(out_q)))
        result%quantiles = out_q
        do j = 1, m
            result%cdf_grid(j) = lower + real(j - 1, dp) * (upper - lower) / real(m - 1, dp)
        end do
        weights = 1.0_dp
        bw = 4.0_dp * (zgrid(2) - zgrid(1))
        do i = 1, size(y)
            do j = 1, zm
                eta = result%beta(j, 1) + dot_product(result%beta(j, 2:), fit%xi_est(i, :))
                if (eta >= 30.0_dp) then
                    raw_cdf(j) = 1.0_dp
                else if (eta <= -30.0_dp) then
                    raw_cdf(j) = 0.0_dp
                else
                    raw_cdf(j) = 1.0_dp / (1.0_dp + exp(-eta))
                end if
            end do
            result%pred_cdf(i, :) = lwls1d(bw, "epan", zgrid, raw_cdf, result%cdf_grid, weights, 1, 0)
            result%pred_cdf(i, :) = max(0.0_dp, min(1.0_dp, result%pred_cdf(i, :)))
            do k = 1, size(out_q)
                result%pred_quantile(i, k) = result%cdf_grid(m)
                do j = 1, m
                    if (result%pred_cdf(i, j) > out_q(k)) then
                        result%pred_quantile(i, k) = result%cdf_grid(j)
                        exit
                    end if
                end do
            end do
        end do
    end subroutine fpc_quantile

    subroutine tvam(obs_t, y, x_obs, grid_t, x_eval, ht, hx, result)
        real(dp), intent(in) :: obs_t(:) !! Pooled longitudinal observation times across all subjects.
        real(dp), intent(in) :: y(:) !! Pooled longitudinal responses corresponding to obs_t rows.
        real(dp), intent(in) :: x_obs(:,:) !! Pooled longitudinal covariates with observations in rows.
        real(dp), intent(in) :: grid_t(:) !! Increasing time grid at which additive surfaces are estimated.
        real(dp), intent(in) :: x_eval(:,:) !! Evaluation points for additive covariate components, one column per covariate.
        real(dp), intent(in) :: ht !! Positive time-window and smoothing bandwidth.
        real(dp), intent(in) :: hx(:) !! Positive additive-component bandwidths, one per covariate.
        type(tvam_result), intent(out) :: result !! Time-varying additive mean and component surfaces on grid_t by x_eval.
        type(sbf_result) :: sbf
        integer, allocatable :: idx(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: center
        integer :: d
        integer :: i
        integer :: j
        integer :: m

        if (size(obs_t) /= size(y) .or. size(x_obs, 1) /= size(y)) error stop "tvam: pooled observation lengths mismatch"
        d = size(x_obs, 2)
        if (size(x_eval, 2) /= d .or. size(hx) /= d .or. ht <= 0.0_dp .or. any(hx <= 0.0_dp)) then
            error stop "tvam: invalid grids or bandwidths"
        end if
        allocate(result%mean_t(size(grid_t)), result%components(size(grid_t), size(x_eval, 1), d))
        do m = 1, size(grid_t)
            idx = pack([(i, i=1, size(y))], abs(obs_t - grid_t(m)) < ht)
            if (size(idx) < max(2, d + 1)) error stop "tvam: time bandwidth leaves too few observations"
            call sb_fitting(y(idx), x_eval, x_obs(idx, :), sbf, h=hx)
            result%mean_t(m) = sbf%mean_y
            result%components(m, :, :) = sbf%fit
        end do
        allocate(weights(size(grid_t)))
        weights = 1.0_dp
        result%mean_t = lwls1d(ht, "epan", grid_t, result%mean_t, grid_t, weights, 0, 0)
        do j = 1, d
            do i = 1, size(x_eval, 1)
                result%components(:, i, j) = lwls1d(ht, "epan", grid_t, result%components(:, i, j), &
                    grid_t, weights, 0, 0)
            end do
        end do
        center = 0.0_dp
        if (size(grid_t) > 0) center = result%mean_t(1)
        result%grid_t = grid_t
        result%x_eval = x_eval
    end subroutine tvam

    subroutine vcam(y, t, x, grid_x, result, bandwidth)
        real(dp), intent(in) :: y(:,:) !! Dense longitudinal responses with subjects in rows and common times in columns.
        real(dp), intent(in) :: t(:) !! Common increasing longitudinal time grid.
        real(dp), intent(in) :: x(:,:) !! Subject-level scalar covariates with subjects in rows.
        real(dp), intent(in) :: grid_x(:,:) !! Evaluation grids for additive covariate transformations, one column per covariate.
        type(vcam_result), intent(out) :: result !! Additive covariate effects, varying coefficients, fitted curves, and evaluation grids.
        real(dp), intent(in), optional :: bandwidth(:) !! Optional positive additive smoothing bandwidths; defaults by sb_fitting.
        type(sbf_result) :: sbf
        real(dp), allocatable :: beta(:,:)
        real(dp), allocatable :: design(:,:)
        real(dp), allocatable :: int_y(:)
        real(dp), allocatable :: phi_subject(:,:)
        real(dp), allocatable :: rhs(:,:)
        real(dp), allocatable :: t_scaled(:)
        integer :: i
        integer :: info
        integer :: j

        if (size(y, 1) /= size(x, 1) .or. size(y, 2) /= size(t) .or. size(grid_x, 2) /= size(x, 2)) then
            error stop "vcam: incompatible dense dimensions"
        end if
        allocate(int_y(size(y, 1)), t_scaled(size(t)))
        t_scaled = (t - t(1)) / (t(size(t)) - t(1))
        do i = 1, size(y, 1)
            int_y(i) = trapz_rcpp(t_scaled, y(i, :))
        end do
        if (present(bandwidth)) then
            call sb_fitting(int_y, grid_x, x, sbf, h=bandwidth)
        else
            call sb_fitting(int_y, grid_x, x, sbf)
        end if
        result%phi_est = sbf%fit
        result%grid_x = grid_x
        allocate(phi_subject(size(x, 1), size(x, 2)), design(size(x, 1), size(x, 2) + 1), rhs(size(x, 1), 1))
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                phi_subject(i, j) = linear_interp(grid_x(:, j), sbf%fit(:, j), x(i, j), .true.)
            end do
        end do
        design(:, 1) = 1.0_dp
        design(:, 2:) = phi_subject
        allocate(result%beta0_est(size(t)), result%beta_est(size(t), size(x, 2)), result%fitted_y(size(y, 1), size(t)))
        do j = 1, size(t)
            rhs(:, 1) = y(:, j)
            call least_squares(design, rhs, beta, info, ridge=1.0e-10_dp)
            if (info /= 0) error stop "vcam: varying-coefficient regression failed"
            result%beta0_est(j) = beta(1, 1)
            result%beta_est(j, :) = beta(2:, 1)
            result%fitted_y(:, j) = matmul(design, beta(:, 1))
        end do
        result%grid_t = t
    end subroutine vcam


    function smooth_unordered(bw, kernel_type, xin, yin, xout, npoly, nder) result(values)
        real(dp), intent(in) :: bw !! Positive local-smoothing bandwidth passed to lwls1d.
        character(len=*), intent(in) :: kernel_type !! Kernel name accepted by lwls1d.
        real(dp), intent(in) :: xin(:) !! Potentially unordered observation coordinates.
        real(dp), intent(in) :: yin(:) !! Observation values paired elementwise with xin.
        real(dp), intent(in) :: xout(:) !! Increasing coordinates at which to evaluate the smoother.
        integer, intent(in) :: npoly !! Local polynomial degree passed to lwls1d.
        integer, intent(in) :: nder !! Derivative order passed to lwls1d.
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: xs(:)
        real(dp), allocatable :: ys(:)
        real(dp), allocatable :: w(:)
        integer :: i
        integer :: j
        integer :: key
        integer, allocatable :: order(:)

        if (size(xin) /= size(yin)) error stop "smooth_unordered: input lengths differ"
        allocate(order(size(xin)), xs(size(xin)), ys(size(yin)), w(size(xin)))
        order = [(i, i=1, size(xin))]
        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (xin(order(j)) <= xin(key)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
        xs = xin(order)
        ys = yin(order)
        w = 1.0_dp
        values = lwls1d(bw, kernel_type, xs, ys, xout, w, npoly, nder)
    end function smooth_unordered

    subroutine project_scores(curves, t, fit, scores)
        real(dp), intent(in) :: curves(:,:) !! New dense curves on the same common support as fit.
        real(dp), intent(in) :: t(:) !! Common support used for numerical score integration.
        type(fpca_result), intent(in) :: fit !! Training FPCA fit defining mean and eigenfunctions.
        real(dp), intent(out) :: scores(:,:) !! Projected scores with one column per retained training eigenfunction.
        integer :: i
        integer :: j

        if (size(curves, 2) /= size(t) .or. size(scores, 1) /= size(curves, 1) .or. &
            size(scores, 2) /= size(fit%phi, 2)) error stop "project_scores: shape mismatch"
        do i = 1, size(curves, 1)
            do j = 1, size(fit%phi, 2)
                scores(i, j) = trapz_rcpp(t, (curves(i, :) - fit%mu) * fit%phi(:, j))
            end do
        end do
    end subroutine project_scores

    subroutine sample_permutation(n, indices)
        integer, intent(in) :: n !! Positive population size whose one-based permutation is requested.
        integer, intent(out) :: indices(:) !! Random permutation of 1:n; output length must equal n.
        real(dp) :: u
        integer :: i
        integer :: j
        integer :: tmp

        if (size(indices) /= n .or. n < 1) error stop "sample_permutation: invalid size"
        indices = [(i, i=1, n)]
        do i = 1, n - 1
            call random_number(u)
            j = i + min(n - i, int(u * real(n - i + 1, dp)))
            tmp = indices(i)
            indices(i) = indices(j)
            indices(j) = tmp
        end do
    end subroutine sample_permutation

end module fdapace_models
