module fdausc_generalized
    use r_kinds, only : dp
    use r_linalg, only : solve_system, inverse_matrix
    use fdausc_extended, only : basis_project_grid, hetero_anova_onefactor
    use fdausc_regression, only : weighted_linear_fit
    use fdausc_statistics, only : fdr_adjusted_pvalue
    implicit none
    private

    public :: basis_scalar_regression, basis_response_regression
    public :: glm_irls, multiclass_glm_classify, additive_backfit
    public :: gls_fit, partially_linear_fit, deviance_smoother_score
    public :: influence_quantile_summary, projected_two_sample_ks
    public :: random_projection_anova
    public :: penalized_glm_irls, multiclass_penalized_glm_classify
    public :: generalized_additive_backfit, multiclass_additive_classify
    public :: forward_glm_select, bootstrap_linear_refits

contains

    subroutine basis_scalar_regression(curves, basis_x, y, coefficients, beta_coefficients, fitted, residuals, hat, &
                                       cross_basis, weights, ridge, center, info)
        real(dp), intent(in) :: curves(:, :) !! Functional predictor curves by rows on a common observation grid.
        real(dp), intent(in) :: basis_x(:, :) !! Predictor basis functions by rows on the same grid as curves.
        real(dp), intent(in) :: y(:) !! Scalar response vector aligned with functional observations.
        real(dp), intent(out) :: coefficients(:) !! Intercept followed by regression coefficients in the transformed basis.
        real(dp), intent(out) :: beta_coefficients(:) !! Functional coefficient-basis coefficients excluding the intercept.
        real(dp), intent(out) :: fitted(:) !! Fitted scalar responses.
        real(dp), intent(out) :: residuals(:) !! Scalar residuals y-fitted.
        real(dp), intent(out) :: hat(:, :) !! Linear-model hat matrix on the observations.
        real(dp), intent(in), optional :: cross_basis(:, :) !! Optional predictor-to-coefficient basis inner-product matrix.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights.
        real(dp), intent(in), optional :: ridge !! Optional nonnegative ridge penalty on transformed basis coefficients.
        logical, intent(in), optional :: center !! Whether to center curves before basis projection; default true.
        integer, intent(out), optional :: info !! Zero on successful projection and regression; nonzero on a linear solve failure.
        real(dp), allocatable :: projected(:, :)
        real(dp), allocatable :: transformed(:, :)
        real(dp), allocatable :: mean_curve(:)
        real(dp), allocatable :: design(:, :)
        integer :: stat

        allocate(projected(size(curves, 1), size(basis_x, 1)))
        allocate(mean_curve(size(curves, 2)))
        call basis_project_grid(curves, basis_x, projected, mean_curve, center, stat)
        if (stat /= 0) then
            coefficients = 0.0_dp
            beta_coefficients = 0.0_dp
            fitted = 0.0_dp
            residuals = y
            hat = 0.0_dp
            if (present(info)) info = stat
            return
        end if
        if (present(cross_basis)) then
            allocate(transformed(size(projected, 1), size(cross_basis, 2)))
            transformed = matmul(projected, cross_basis)
        else
            allocate(transformed(size(projected, 1), size(projected, 2)))
            transformed = projected
        end if
        allocate(design(size(transformed, 1), size(transformed, 2) + 1))
        design(:, 1) = 1.0_dp
        design(:, 2:) = transformed
        call weighted_linear_fit(design, y, coefficients, fitted, residuals, hat, weights, ridge, stat)
        beta_coefficients = coefficients(2:)
        if (present(info)) info = stat
    end subroutine basis_scalar_regression

    subroutine basis_response_regression(x_scores, y_scores, coefficients, fitted_scores, residual_scores, ridge, info)
        real(dp), intent(in) :: x_scores(:, :) !! Predictor basis scores with observations in rows.
        real(dp), intent(in) :: y_scores(:, :) !! Functional-response basis scores with observations in rows.
        real(dp), intent(out) :: coefficients(:, :) !! Intercept and predictor coefficient matrix, one response basis per column.
        real(dp), intent(out) :: fitted_scores(:, :) !! Fitted response-basis scores.
        real(dp), intent(out) :: residual_scores(:, :) !! Response-basis residual scores.
        real(dp), intent(in), optional :: ridge !! Nonnegative ridge penalty applied to non-intercept coefficients.
        integer, intent(out), optional :: info !! Zero on success; nonzero if any response-basis solve fails.
        real(dp), allocatable :: design(:, :)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: solution(:)
        real(dp) :: lambda
        integer :: i
        integer :: j
        integer :: stat
        integer :: final_stat

        lambda = 0.0_dp
        if (present(ridge)) lambda = max(0.0_dp, ridge)
        allocate(design(size(x_scores, 1), size(x_scores, 2) + 1))
        design(:, 1) = 1.0_dp
        design(:, 2:) = x_scores
        gram = matmul(transpose(design), design)
        do i = 2, size(gram, 1)
            gram(i, i) = gram(i, i) + lambda
        end do
        allocate(rhs(size(design, 2)), solution(size(design, 2)))
        final_stat = 0
        do j = 1, size(y_scores, 2)
            rhs = matmul(transpose(design), y_scores(:, j))
            call solve_system(gram, rhs, solution, stat)
            if (stat /= 0) final_stat = stat
            coefficients(:, j) = solution
        end do
        if (final_stat /= 0) then
            fitted_scores = 0.0_dp
            residual_scores = y_scores
        else
            fitted_scores = matmul(design, coefficients)
            residual_scores = y_scores - fitted_scores
        end if
        if (present(info)) info = final_stat
    end subroutine basis_response_regression

    subroutine glm_irls(x, y, family_code, coefficients, fitted, residuals, hat, deviance, iterations, converged, &
                        weights, ridge, max_iterations, tolerance, info)
        real(dp), intent(in) :: x(:, :) !! Numeric design matrix including an intercept column when desired.
        real(dp), intent(in) :: y(:) !! Response values; binomial responses must lie in [0,1] and Poisson responses be nonnegative.
        integer, intent(in) :: family_code !! Family selector: 1 Gaussian identity, 2 binomial logit, or 3 Poisson log.
        real(dp), intent(out) :: coefficients(:) !! Fitted generalized-linear-model coefficient vector.
        real(dp), intent(out) :: fitted(:) !! Fitted response means on the response scale.
        real(dp), intent(out) :: residuals(:) !! Response residuals y-fitted.
        real(dp), intent(out) :: hat(:, :) !! Final IRLS working-weight hat matrix.
        real(dp), intent(out) :: deviance !! Family deviance at the fitted means.
        integer, intent(out) :: iterations !! Number of IRLS iterations performed.
        logical, intent(out) :: converged !! True when the coefficient update satisfies the requested tolerance.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative prior observation weights.
        real(dp), intent(in), optional :: ridge !! Nonnegative ridge penalty on non-intercept coefficients.
        integer, intent(in), optional :: max_iterations !! Positive maximum number of IRLS iterations; default 50.
        real(dp), intent(in), optional :: tolerance !! Positive relative coefficient tolerance; default 1e-8.
        integer, intent(out), optional :: info !! Zero on a successful solve; nonzero when a working normal equation is singular.
        real(dp), allocatable :: eta(:)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: variance(:)
        real(dp), allocatable :: derivative(:)
        real(dp), allocatable :: work_weights(:)
        real(dp), allocatable :: prior_weights(:)
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: wx(:, :)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: beta_new(:)
        real(dp), allocatable :: inv(:, :)
        real(dp), allocatable :: weighted_design(:, :)
        real(dp) :: lambda
        real(dp) :: tol
        real(dp) :: denom
        integer :: i
        integer :: iter
        integer :: maxit
        integer :: stat

        lambda = 0.0_dp
        if (present(ridge)) lambda = max(0.0_dp, ridge)
        tol = 1.0e-8_dp
        if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
        maxit = 50
        if (present(max_iterations)) maxit = max(1, max_iterations)
        allocate(eta(size(y)), mu(size(y)), variance(size(y)), derivative(size(y)))
        allocate(work_weights(size(y)), prior_weights(size(y)), z(size(y)))
        allocate(wx(size(x, 1), size(x, 2)), rhs(size(x, 2)), beta_new(size(x, 2)))
        prior_weights = 1.0_dp
        if (present(weights)) prior_weights = max(weights, 0.0_dp)
        coefficients = 0.0_dp
        if (size(coefficients) >= 1) then
            select case (family_code)
            case (2)
                coefficients(1) = log(max(sum(prior_weights*y), 0.5_dp) &
                    /max(sum(prior_weights*(1.0_dp - y)), 0.5_dp))
            case (3)
                coefficients(1) = log(max(sum(prior_weights*y)/max(sum(prior_weights), tiny(1.0_dp)), 1.0e-8_dp))
            case default
                coefficients(1) = sum(prior_weights*y)/max(sum(prior_weights), tiny(1.0_dp))
            end select
        end if
        converged = .false.
        stat = 0
        do iter = 1, maxit
            eta = matmul(x, coefficients)
            call family_moments(eta, family_code, mu, variance, derivative)
            work_weights = prior_weights*derivative*derivative/max(variance, 1.0e-12_dp)
            z = eta + (y - mu)/max(abs(derivative), 1.0e-12_dp)*sign(1.0_dp, derivative)
            do i = 1, size(y)
                wx(i, :) = work_weights(i)*x(i, :)
            end do
            gram = matmul(transpose(x), wx)
            do i = 2, size(gram, 1)
                gram(i, i) = gram(i, i) + lambda
            end do
            rhs = matmul(transpose(x), work_weights*z)
            call solve_system(gram, rhs, beta_new, stat)
            if (stat /= 0) exit
            denom = max(1.0_dp, maxval(abs(coefficients)))
            if (maxval(abs(beta_new - coefficients))/denom <= tol) then
                coefficients = beta_new
                converged = .true.
                exit
            end if
            coefficients = beta_new
        end do
        iterations = min(iter, maxit)
        eta = matmul(x, coefficients)
        call family_moments(eta, family_code, fitted, variance, derivative)
        residuals = y - fitted
        deviance = family_deviance(y, fitted, prior_weights, family_code)
        if (stat /= 0) then
            hat = 0.0_dp
            if (present(info)) info = stat
            return
        end if
        work_weights = prior_weights*derivative*derivative/max(variance, 1.0e-12_dp)
        allocate(weighted_design(size(x, 1), size(x, 2)))
        do i = 1, size(y)
            weighted_design(i, :) = work_weights(i)*x(i, :)
        end do
        gram = matmul(transpose(x), weighted_design)
        do i = 2, size(gram, 1)
            gram(i, i) = gram(i, i) + lambda
        end do
        call inverse_matrix(gram, inv, stat)
        if (stat == 0) then
            hat = matmul(matmul(x, inv), transpose(weighted_design))
        else
            hat = 0.0_dp
        end if
        if (present(info)) info = stat
    end subroutine glm_irls

    subroutine multiclass_glm_classify(x, labels, classes, probabilities, predicted, ridge, max_iterations, tolerance, info)
        real(dp), intent(in) :: x(:, :) !! Numeric classifier design matrix including any desired intercept column.
        integer, intent(in) :: labels(:) !! Integer training class labels aligned with design rows.
        integer, intent(in) :: classes(:) !! Distinct class labels to fit by one-versus-all binomial GLMs.
        real(dp), intent(out) :: probabilities(:, :) !! One-versus-all fitted probabilities by observation and class.
        integer, intent(out) :: predicted(:) !! Class label having the largest fitted one-versus-all probability.
        real(dp), intent(in), optional :: ridge !! Optional ridge penalty on non-intercept GLM coefficients.
        integer, intent(in), optional :: max_iterations !! Optional maximum IRLS iterations for each class.
        real(dp), intent(in), optional :: tolerance !! Optional IRLS convergence tolerance.
        integer, intent(out), optional :: info !! Zero when every class fit succeeds; otherwise the last nonzero solve status.
        real(dp), allocatable :: binary(:)
        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: residuals(:)
        real(dp), allocatable :: hat(:, :)
        real(dp) :: deviance
        real(dp) :: rowsum
        integer :: best
        integer :: c
        integer :: i
        integer :: iterations
        integer :: stat
        integer :: final_stat
        logical :: converged

        allocate(binary(size(labels)), coefficients(size(x, 2)), residuals(size(labels)))
        allocate(hat(size(labels), size(labels)))
        final_stat = 0
        do c = 1, size(classes)
            binary = merge(1.0_dp, 0.0_dp, labels == classes(c))
            call glm_irls(x, binary, 2, coefficients, probabilities(:, c), residuals, hat, deviance, &
                          iterations, converged, ridge=ridge, max_iterations=max_iterations, tolerance=tolerance, info=stat)
            if (stat /= 0) final_stat = stat
        end do
        do i = 1, size(labels)
            rowsum = sum(probabilities(i, :))
            if (rowsum > tiny(1.0_dp)) probabilities(i, :) = probabilities(i, :)/rowsum
            best = maxloc(probabilities(i, :), dim=1)
            predicted(i) = classes(best)
        end do
        if (present(info)) info = final_stat
    end subroutine multiclass_glm_classify

    pure subroutine additive_backfit(smoothers, y, intercept, components, fitted, residuals, iterations, converged, &
                                     max_iterations, tolerance)
        real(dp), intent(in) :: smoothers(:, :, :) !! One square linear smoother matrix per additive term.
        real(dp), intent(in) :: y(:) !! Scalar response vector aligned with smoother rows and columns.
        real(dp), intent(out) :: intercept !! Estimated additive-model intercept.
        real(dp), intent(out) :: components(:, :) !! Centered fitted contribution of each additive smoother term.
        real(dp), intent(out) :: fitted(:) !! Additive fitted response intercept plus all term contributions.
        real(dp), intent(out) :: residuals(:) !! Response residuals y-fitted.
        integer, intent(out) :: iterations !! Number of backfitting sweeps performed.
        logical, intent(out) :: converged !! True when the largest component update is within tolerance.
        integer, intent(in), optional :: max_iterations !! Positive maximum number of backfitting sweeps; default 100.
        real(dp), intent(in), optional :: tolerance !! Positive maximum absolute component-change tolerance; default 1e-8.
        real(dp), allocatable :: old_components(:, :)
        real(dp), allocatable :: partial(:)
        real(dp) :: tol
        integer :: j
        integer :: maxit

        maxit = 100
        if (present(max_iterations)) maxit = max(1, max_iterations)
        tol = 1.0e-8_dp
        if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
        intercept = sum(y)/real(size(y), dp)
        components = 0.0_dp
        allocate(old_components(size(components, 1), size(components, 2)))
        allocate(partial(size(y)))
        converged = .false.
        do iterations = 1, maxit
            old_components = components
            do j = 1, size(smoothers, 3)
                partial = y - intercept - sum(components, dim=2) + components(:, j)
                components(:, j) = matmul(smoothers(:, :, j), partial)
                components(:, j) = components(:, j) - sum(components(:, j))/real(size(y), dp)
            end do
            if (maxval(abs(components - old_components)) <= tol) then
                converged = .true.
                exit
            end if
        end do
        iterations = min(iterations, maxit)
        fitted = intercept + sum(components, dim=2)
        residuals = y - fitted
    end subroutine additive_backfit

    subroutine gls_fit(x, y, covariance, coefficients, fitted, residuals, hat, ridge, info)
        real(dp), intent(in) :: x(:, :) !! Numeric generalized-least-squares design matrix.
        real(dp), intent(in) :: y(:) !! Scalar response vector aligned with design rows.
        real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-definite observation covariance matrix.
        real(dp), intent(out) :: coefficients(:) !! Generalized least-squares coefficient vector.
        real(dp), intent(out) :: fitted(:) !! Fitted responses X*coefficients.
        real(dp), intent(out) :: residuals(:) !! Response residuals y-fitted.
        real(dp), intent(out) :: hat(:, :) !! GLS linear mapping from observed y to fitted values.
        real(dp), intent(in), optional :: ridge !! Nonnegative ridge penalty on non-intercept coefficients.
        integer, intent(out), optional :: info !! Zero on successful inversions/solve; nonzero otherwise.
        real(dp), allocatable :: precision(:, :)
        real(dp), allocatable :: px(:, :)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: invgram(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp) :: lambda
        integer :: i
        integer :: stat

        lambda = 0.0_dp
        if (present(ridge)) lambda = max(0.0_dp, ridge)
        call inverse_matrix(covariance, precision, stat)
        if (stat /= 0) then
            coefficients = 0.0_dp
            fitted = 0.0_dp
            residuals = y
            hat = 0.0_dp
            if (present(info)) info = stat
            return
        end if
        px = matmul(precision, x)
        gram = matmul(transpose(x), px)
        do i = 2, size(gram, 1)
            gram(i, i) = gram(i, i) + lambda
        end do
        rhs = matmul(transpose(x), matmul(precision, y))
        call solve_system(gram, rhs, coefficients, stat)
        if (stat /= 0) then
            fitted = 0.0_dp
            residuals = y
            hat = 0.0_dp
            if (present(info)) info = stat
            return
        end if
        fitted = matmul(x, coefficients)
        residuals = y - fitted
        call inverse_matrix(gram, invgram, stat)
        if (stat == 0) then
            hat = matmul(matmul(x, invgram), matmul(transpose(x), precision))
        else
            hat = 0.0_dp
        end if
        if (present(info)) info = stat
    end subroutine gls_fit

    subroutine partially_linear_fit(smoother, scalar_design, y, coefficients, nonparametric, fitted, residuals, hat, info)
        real(dp), intent(in) :: smoother(:, :) !! Linear smoother matrix for the functional/nonparametric predictor.
        real(dp), intent(in) :: scalar_design(:, :) !! Parametric scalar-covariate design matrix without an added intercept.
        real(dp), intent(in) :: y(:) !! Scalar response vector aligned with smoother rows.
        real(dp), intent(out) :: coefficients(:) !! Partially-linear coefficients for scalar_design columns.
        real(dp), intent(out) :: nonparametric(:) !! Estimated nonparametric component S(y-X*beta).
        real(dp), intent(out) :: fitted(:) !! Combined parametric and nonparametric fitted response.
        real(dp), intent(out) :: residuals(:) !! Response residuals y-fitted.
        real(dp), intent(out) :: hat(:, :) !! Partially-linear hat matrix from the supplied smoother and scalar design.
        integer, intent(out), optional :: info !! Zero on a successful residualized least-squares solve; nonzero otherwise.
        real(dp), allocatable :: identity(:, :)
        real(dp), allocatable :: residualizer(:, :)
        real(dp), allocatable :: xh(:, :)
        real(dp), allocatable :: yh(:)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: invgram(:, :)
        real(dp), allocatable :: c1(:, :)
        integer :: i
        integer :: n
        integer :: stat

        n = size(y)
        allocate(identity(n, n))
        identity = 0.0_dp
        do i = 1, n
            identity(i, i) = 1.0_dp
        end do
        residualizer = identity - smoother
        yh = matmul(residualizer, y)
        xh = matmul(residualizer, scalar_design)
        gram = matmul(transpose(xh), xh)
        call solve_system(gram, matmul(transpose(xh), yh), coefficients, stat)
        if (stat /= 0) then
            nonparametric = 0.0_dp
            fitted = 0.0_dp
            residuals = y
            hat = smoother
            if (present(info)) info = stat
            return
        end if
        nonparametric = matmul(smoother, y - matmul(scalar_design, coefficients))
        fitted = matmul(scalar_design, coefficients) + nonparametric
        residuals = y - fitted
        call inverse_matrix(gram, invgram, stat)
        if (stat == 0) then
            c1 = matmul(invgram, transpose(xh))
            hat = matmul(matmul(matmul(residualizer, scalar_design), c1), residualizer) + smoother
        else
            hat = smoother
        end if
        if (present(info)) info = stat
    end subroutine partially_linear_fit

    pure subroutine deviance_smoother_score(y, smoother, observed, family_code, offset, offset_df, criterion_code, &
                                            score, weights, trim)
        real(dp), intent(in) :: y(:) !! Working response multiplied by the smoother before inverse-link transformation.
        real(dp), intent(in) :: smoother(:, :) !! Square smoother/hat matrix whose trace supplies effective degrees of freedom.
        real(dp), intent(in) :: observed(:) !! Observed response used to form standardized residuals and deviance score.
        integer, intent(in) :: family_code !! Family selector: 1 Gaussian identity, 2 binomial logit, or 3 Poisson log.
        real(dp), intent(in) :: offset(:) !! Additive linear-predictor offset for each observation.
        real(dp), intent(in) :: offset_df !! Effective degrees of freedom contributed by offset/other terms.
        integer, intent(in) :: criterion_code !! 1 GCV, 2 exponential GCV, 3 Shibata, 4 Rice, or other finite correction.
        real(dp), intent(out) :: score !! Trimmed weighted deviance-style smoothing criterion; infinity if degrees are excessive.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative quadratic-loss weights; defaults to one.
        real(dp), intent(in), optional :: trim !! Fraction of largest standardized residual magnitudes discarded; default zero.
        real(dp), allocatable :: eta(:)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: variance(:)
        real(dp), allocatable :: derivative(:)
        real(dp), allocatable :: abs_residual(:)
        real(dp), allocatable :: sorted(:)
        real(dp), allocatable :: w(:)
        logical, allocatable :: keep(:)
        real(dp) :: cutoff
        real(dp) :: ndf
        real(dp) :: res
        real(dp) :: trim_fraction
        real(dp) :: vv
        real(dp) :: nreal
        integer :: i

        eta = matmul(smoother, y) + offset
        allocate(mu(size(y)), variance(size(y)), derivative(size(y)))
        call family_moments(eta, family_code, mu, variance, derivative)
        allocate(abs_residual(size(y)), sorted(size(y)), w(size(y)), keep(size(y)))
        abs_residual = abs((observed - mu)/sqrt(max(variance, 1.0e-12_dp)))
        w = 1.0_dp
        if (present(weights)) w = max(weights, 0.0_dp)
        keep = .true.
        trim_fraction = 0.0_dp
        if (present(trim)) trim_fraction = max(0.0_dp, min(0.999999_dp, trim))
        if (trim_fraction > 0.0_dp) then
            sorted = abs_residual
            call sort_real_local(sorted)
            cutoff = quantile_type4_local(sorted, 1.0_dp - trim_fraction)
            keep = abs_residual <= cutoff
        end if
        res = 0.0_dp
        ndf = offset_df
        do i = 1, size(y)
            if (keep(i)) then
                res = res + w(i)*(observed(i) - mu(i))**2/max(variance(i), 1.0e-12_dp)
                ndf = ndf + smoother(i, i)
            end if
        end do
        nreal = real(size(y), dp)
        if (ndf > 0.5_dp*nreal) then
            score = huge(1.0_dp)
            return
        end if
        select case (criterion_code)
        case (1)
            vv = (1.0_dp - ndf/nreal)**(-2)
        case (2)
            vv = exp(2.0_dp*ndf/nreal)
        case (3, 4)
            vv = (1.0_dp + ndf/nreal)/(1.0_dp - ndf/nreal)
        case default
            vv = 1.0_dp/(1.0_dp - 2.0_dp*ndf/nreal)
        end select
        score = res*vv/nreal
    end subroutine deviance_smoother_score

    pure subroutine influence_quantile_summary(dcp, dce, dp_stat, bootstrap_dcp, bootstrap_dce, bootstrap_dp, &
                                               quantile_dcp, quantile_dce, quantile_dp)
        real(dp), intent(in) :: dcp(:) !! Original prediction-influence statistics, one per observation.
        real(dp), intent(in) :: dce(:) !! Original estimation-influence statistics, one per observation.
        real(dp), intent(in) :: dp_stat(:) !! Original penalty/influence statistics, one per observation.
        real(dp), intent(in) :: bootstrap_dcp(:) !! Pooled bootstrap prediction-influence statistics.
        real(dp), intent(in) :: bootstrap_dce(:) !! Pooled bootstrap estimation-influence statistics.
        real(dp), intent(in) :: bootstrap_dp(:) !! Pooled bootstrap penalty/influence statistics.
        real(dp), intent(out) :: quantile_dcp(:) !! Empirical bootstrap CDF evaluated at each original DCP statistic.
        real(dp), intent(out) :: quantile_dce(:) !! Empirical bootstrap CDF evaluated at each original DCE statistic.
        real(dp), intent(out) :: quantile_dp(:) !! Empirical bootstrap CDF evaluated at each original DP statistic.
        integer :: i

        do i = 1, size(dcp)
            quantile_dcp(i) = real(count(bootstrap_dcp <= dcp(i)), dp)/real(max(1, size(bootstrap_dcp)), dp)
            quantile_dce(i) = real(count(bootstrap_dce <= dce(i)), dp)/real(max(1, size(bootstrap_dce)), dp)
            quantile_dp(i) = real(count(bootstrap_dp <= dp_stat(i)), dp)/real(max(1, size(bootstrap_dp)), dp)
        end do
    end subroutine influence_quantile_summary

    pure subroutine projected_two_sample_ks(x, y, projections, statistics, pvalues, fdr_pvalue)
        real(dp), intent(in) :: x(:, :) !! First functional sample by rows on a common numeric grid or coefficient basis.
        real(dp), intent(in) :: y(:, :) !! Second functional sample by rows on the same grid or coefficient basis as x.
        real(dp), intent(in) :: projections(:, :) !! Projection directions by rows, each matching the sample column dimension.
        real(dp), intent(out) :: statistics(:) !! Two-sample Kolmogorov-Smirnov statistic for each projection direction.
        real(dp), intent(out) :: pvalues(:) !! Large-sample two-sided KS p-value approximation for each projection.
        real(dp), intent(out) :: fdr_pvalue !! Benjamini-Hochberg adjusted minimum p-value across projections.
        real(dp), allocatable :: xp(:)
        real(dp), allocatable :: yp(:)
        integer :: j

        allocate(xp(size(x, 1)), yp(size(y, 1)))
        do j = 1, size(projections, 1)
            xp = matmul(x, projections(j, :))
            yp = matmul(y, projections(j, :))
            call two_sample_ks(xp, yp, statistics(j), pvalues(j))
        end do
        fdr_pvalue = fdr_adjusted_pvalue(pvalues)
    end subroutine projected_two_sample_ks

    pure subroutine random_projection_anova(curves, projections, groups, statistics, heteroscedastic, min_statistic)
        real(dp), intent(in) :: curves(:, :) !! Functional observations by rows on a common grid or coefficient basis.
        real(dp), intent(in) :: projections(:, :) !! Projection directions by rows matching the curve column dimension.
        integer, intent(in) :: groups(:) !! Positive one-factor group label for each functional observation.
        real(dp), intent(out) :: statistics(:) !! One-way ANOVA statistic for every projection.
        logical, intent(in), optional :: heteroscedastic !! Use the Brunner-Dette-Munk heteroscedastic statistic when true.
        real(dp), intent(out), optional :: min_statistic !! Minimum statistic over the supplied projections.
        real(dp), allocatable :: values(:)
        real(dp) :: df1
        real(dp) :: df2
        logical :: hetero
        integer :: j

        hetero = .true.
        if (present(heteroscedastic)) hetero = heteroscedastic
        allocate(values(size(curves, 1)))
        do j = 1, size(projections, 1)
            values = matmul(curves, projections(j, :))
            if (hetero) then
                call hetero_anova_onefactor(values, groups, statistics(j), df1, df2)
            else
                call ordinary_anova_onefactor(values, groups, statistics(j))
            end if
        end do
        if (present(min_statistic)) min_statistic = minval(statistics)
    end subroutine random_projection_anova

    subroutine penalized_glm_irls(x, y, family_code, penalty, coefficients, fitted, residuals, deviance, iterations, &
                                      converged, weights, max_iterations, tolerance, info)
        real(dp), intent(in) :: x(:, :) !! Expanded numeric design matrix, including smooth-term basis columns and intercept.
        real(dp), intent(in) :: y(:) !! Response values; binomial values must lie in [0,1] and Poisson values be nonnegative.
        integer, intent(in) :: family_code !! Family selector: 1 Gaussian identity, 2 binomial logit, or 3 Poisson log.
        real(dp), intent(in) :: penalty(:, :) !! Symmetric positive-semidefinite coefficient penalty matrix matching design columns.
        real(dp), intent(out) :: coefficients(:) !! Penalized generalized-linear-model coefficients.
        real(dp), intent(out) :: fitted(:) !! Fitted response means on the response scale.
        real(dp), intent(out) :: residuals(:) !! Response residuals y-fitted.
        real(dp), intent(out) :: deviance !! Family deviance at the fitted response means.
        integer, intent(out) :: iterations !! Number of penalized IRLS iterations performed.
        logical, intent(out) :: converged !! True when the relative coefficient update reaches tolerance.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative prior observation weights.
        integer, intent(in), optional :: max_iterations !! Positive maximum penalized IRLS iterations; default 50.
        real(dp), intent(in), optional :: tolerance !! Positive relative coefficient tolerance; default 1e-8.
        integer, intent(out), optional :: info !! Zero on successful normal-equation solves; nonzero on a singular solve.
        real(dp), allocatable :: eta(:)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: variance(:)
        real(dp), allocatable :: derivative(:)
        real(dp), allocatable :: work_weights(:)
        real(dp), allocatable :: prior_weights(:)
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: wx(:, :)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: beta_new(:)
        real(dp) :: denom
        real(dp) :: tol
        integer :: i
        integer :: iter
        integer :: maxit
        integer :: stat

        tol = 1.0e-8_dp
        if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
        maxit = 50
        if (present(max_iterations)) maxit = max(1, max_iterations)
        allocate(eta(size(y)), mu(size(y)), variance(size(y)), derivative(size(y)))
        allocate(work_weights(size(y)), prior_weights(size(y)), z(size(y)))
        allocate(wx(size(x, 1), size(x, 2)), rhs(size(x, 2)), beta_new(size(x, 2)))
        prior_weights = 1.0_dp
        if (present(weights)) prior_weights = max(weights, 0.0_dp)
        coefficients = 0.0_dp
        if (size(coefficients) >= 1) then
            select case (family_code)
            case (2)
                coefficients(1) = log(max(sum(prior_weights*y), 0.5_dp) &
                    /max(sum(prior_weights*(1.0_dp - y)), 0.5_dp))
            case (3)
                coefficients(1) = log(max(sum(prior_weights*y)/max(sum(prior_weights), tiny(1.0_dp)), 1.0e-8_dp))
            case default
                coefficients(1) = sum(prior_weights*y)/max(sum(prior_weights), tiny(1.0_dp))
            end select
        end if
        stat = 0
        converged = .false.
        do iter = 1, maxit
            eta = matmul(x, coefficients)
            call family_moments(eta, family_code, mu, variance, derivative)
            work_weights = prior_weights*derivative*derivative/max(variance, 1.0e-12_dp)
            z = eta + (y - mu)/max(abs(derivative), 1.0e-12_dp)*sign(1.0_dp, derivative)
            do i = 1, size(y)
                wx(i, :) = work_weights(i)*x(i, :)
            end do
            gram = matmul(transpose(x), wx) + penalty
            rhs = matmul(transpose(x), work_weights*z)
            call solve_system(gram, rhs, beta_new, stat)
            if (stat /= 0) exit
            denom = max(1.0_dp, maxval(abs(coefficients)))
            if (maxval(abs(beta_new - coefficients))/denom <= tol) then
                coefficients = beta_new
                converged = .true.
                exit
            end if
            coefficients = beta_new
        end do
        iterations = min(iter, maxit)
        eta = matmul(x, coefficients)
        call family_moments(eta, family_code, fitted, variance, derivative)
        residuals = y - fitted
        deviance = family_deviance(y, fitted, prior_weights, family_code)
        if (present(info)) info = stat
    end subroutine penalized_glm_irls

    subroutine multiclass_penalized_glm_classify(x, labels, classes, penalty, probabilities, predicted, &
                                                  max_iterations, tolerance, info)
        real(dp), intent(in) :: x(:, :) !! Expanded classifier design matrix including smooth-term basis columns.
        integer, intent(in) :: labels(:) !! Integer training class labels aligned with design rows.
        integer, intent(in) :: classes(:) !! Distinct class labels fitted by one-versus-all penalized binomial models.
        real(dp), intent(in) :: penalty(:, :) !! Common coefficient penalty matrix for each one-versus-all model.
        real(dp), intent(out) :: probabilities(:, :) !! Normalized one-versus-all fitted probabilities by observation and class.
        integer, intent(out) :: predicted(:) !! Predicted class label for each observation.
        integer, intent(in), optional :: max_iterations !! Optional maximum penalized IRLS iterations for each class.
        real(dp), intent(in), optional :: tolerance !! Optional penalized IRLS coefficient tolerance.
        integer, intent(out), optional :: info !! Zero when every class fit succeeds; otherwise the last nonzero solve status.
        real(dp), allocatable :: binary(:)
        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: residuals(:)
        real(dp) :: deviance
        real(dp) :: rowsum
        integer :: best
        integer :: c
        integer :: i
        integer :: iterations
        integer :: stat
        integer :: final_stat
        logical :: converged

        allocate(binary(size(labels)), coefficients(size(x, 2)), residuals(size(labels)))
        final_stat = 0
        do c = 1, size(classes)
            binary = merge(1.0_dp, 0.0_dp, labels == classes(c))
            call penalized_glm_irls(x, binary, 2, penalty, coefficients, probabilities(:, c), residuals, deviance, &
                                    iterations, converged, max_iterations=max_iterations, tolerance=tolerance, info=stat)
            if (stat /= 0) final_stat = stat
        end do
        do i = 1, size(labels)
            rowsum = sum(probabilities(i, :))
            if (rowsum > tiny(1.0_dp)) probabilities(i, :) = probabilities(i, :)/rowsum
            best = maxloc(probabilities(i, :), dim=1)
            predicted(i) = classes(best)
        end do
        if (present(info)) info = final_stat
    end subroutine multiclass_penalized_glm_classify

    pure subroutine generalized_additive_backfit(smoothers, y, family_code, intercept, components, fitted, residuals, &
                                                  iterations, converged, max_iterations, tolerance)
        real(dp), intent(in) :: smoothers(:, :, :) !! Fixed linear smoother matrices for additive predictor terms.
        real(dp), intent(in) :: y(:) !! Response vector; binomial responses must lie in [0,1] and Poisson be nonnegative.
        integer, intent(in) :: family_code !! Family selector: 1 Gaussian identity, 2 binomial logit, or 3 Poisson log.
        real(dp), intent(out) :: intercept !! Fitted additive-model intercept on the link scale.
        real(dp), intent(out) :: components(:, :) !! Centered additive component contributions on the link scale.
        real(dp), intent(out) :: fitted(:) !! Fitted response means after applying the inverse link.
        real(dp), intent(out) :: residuals(:) !! Response residuals y-fitted.
        integer, intent(out) :: iterations !! Number of local-scoring/backfitting sweeps performed.
        logical, intent(out) :: converged !! True when the link-scale predictor changes by no more than tolerance.
        integer, intent(in), optional :: max_iterations !! Positive maximum local-scoring sweeps; default 100.
        real(dp), intent(in), optional :: tolerance !! Positive maximum link-scale update tolerance; default 1e-8.
        real(dp), allocatable :: eta(:)
        real(dp), allocatable :: old_eta(:)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: variance(:)
        real(dp), allocatable :: derivative(:)
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: partial(:)
        real(dp) :: tol
        integer :: j
        integer :: maxit

        maxit = 100
        if (present(max_iterations)) maxit = max(1, max_iterations)
        tol = 1.0e-8_dp
        if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
        select case (family_code)
        case (2)
            intercept = log(max(sum(y), 0.5_dp)/max(real(size(y), dp) - sum(y), 0.5_dp))
        case (3)
            intercept = log(max(sum(y)/real(size(y), dp), 1.0e-8_dp))
        case default
            intercept = sum(y)/real(size(y), dp)
        end select
        components = 0.0_dp
        allocate(eta(size(y)), old_eta(size(y)), mu(size(y)), variance(size(y)), derivative(size(y)))
        allocate(z(size(y)), partial(size(y)))
        converged = .false.
        do iterations = 1, maxit
            eta = intercept + sum(components, dim=2)
            old_eta = eta
            call family_moments(eta, family_code, mu, variance, derivative)
            z = eta + (y - mu)/max(abs(derivative), 1.0e-12_dp)*sign(1.0_dp, derivative)
            intercept = sum(z - sum(components, dim=2))/real(size(y), dp)
            do j = 1, size(smoothers, 3)
                partial = z - intercept - sum(components, dim=2) + components(:, j)
                components(:, j) = matmul(smoothers(:, :, j), partial)
                components(:, j) = components(:, j) - sum(components(:, j))/real(size(y), dp)
            end do
            eta = intercept + sum(components, dim=2)
            if (maxval(abs(eta - old_eta)) <= tol) then
                converged = .true.
                exit
            end if
        end do
        iterations = min(iterations, maxit)
        eta = intercept + sum(components, dim=2)
        call family_moments(eta, family_code, fitted, variance, derivative)
        residuals = y - fitted
    end subroutine generalized_additive_backfit

    pure subroutine multiclass_additive_classify(smoothers, labels, classes, probabilities, predicted, &
                                                  max_iterations, tolerance)
        real(dp), intent(in) :: smoothers(:, :, :) !! Fixed linear smoother matrices shared by all one-versus-all class fits.
        integer, intent(in) :: labels(:) !! Integer training class labels aligned with smoother rows.
        integer, intent(in) :: classes(:) !! Distinct class labels fitted by one-versus-all additive logistic models.
        real(dp), intent(out) :: probabilities(:, :) !! Normalized additive logistic fitted probabilities by observation/class.
        integer, intent(out) :: predicted(:) !! Predicted class label for each observation.
        integer, intent(in), optional :: max_iterations !! Optional maximum local-scoring sweeps for each class model.
        real(dp), intent(in), optional :: tolerance !! Optional local-scoring convergence tolerance.
        real(dp), allocatable :: binary(:)
        real(dp), allocatable :: components(:, :)
        real(dp), allocatable :: residuals(:)
        real(dp) :: intercept
        real(dp) :: rowsum
        integer :: best
        integer :: c
        integer :: i
        integer :: iterations
        logical :: converged

        allocate(binary(size(labels)), components(size(labels), size(smoothers, 3)), residuals(size(labels)))
        do c = 1, size(classes)
            binary = merge(1.0_dp, 0.0_dp, labels == classes(c))
            call generalized_additive_backfit(smoothers, binary, 2, intercept, components, probabilities(:, c), residuals, &
                                               iterations, converged, max_iterations, tolerance)
        end do
        do i = 1, size(labels)
            rowsum = sum(probabilities(i, :))
            if (rowsum > tiny(1.0_dp)) probabilities(i, :) = probabilities(i, :)/rowsum
            best = maxloc(probabilities(i, :), dim=1)
            predicted(i) = classes(best)
        end do
    end subroutine multiclass_additive_classify

    subroutine forward_glm_select(x, y, column_groups, family_code, selected_groups, selection_order, selected_count, &
                                  final_deviance, min_improvement, max_iterations, tolerance, info)
        real(dp), intent(in) :: x(:, :) !! Full candidate design matrix; column 1 is always retained as the intercept/base term.
        real(dp), intent(in) :: y(:) !! Response vector aligned with design rows.
        integer, intent(in) :: column_groups(:) !! Group id for each design column; first column should use group 0.
        integer, intent(in) :: family_code !! Family selector passed to the Gaussian/binomial/Poisson IRLS fitter.
        logical, intent(out) :: selected_groups(:) !! Selection flags for group ids 1:size(selected_groups).
        integer, intent(out) :: selection_order(:) !! Group ids in greedy selection order; trailing positions are zero.
        integer, intent(out) :: selected_count !! Number of selected predictor groups.
        real(dp), intent(out) :: final_deviance !! Deviance of the final selected model.
        real(dp), intent(in), optional :: min_improvement !! Required decrease in deviance-plus-2df score; default zero.
        integer, intent(in), optional :: max_iterations !! Optional maximum IRLS iterations in each candidate fit.
        real(dp), intent(in), optional :: tolerance !! Optional IRLS coefficient convergence tolerance.
        integer, intent(out), optional :: info !! Zero when all required candidate fits solve; last nonzero solve code otherwise.
        real(dp), allocatable :: design(:, :)
        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: residuals(:)
        real(dp), allocatable :: hat(:, :)
        real(dp) :: best_score
        real(dp) :: candidate_score
        real(dp) :: current_score
        real(dp) :: deviance
        real(dp) :: candidate_deviance
        real(dp) :: improvement
        integer :: best_group
        integer :: candidate
        integer :: iterations
        integer :: ngroup
        integer :: stat
        integer :: final_stat
        logical :: converged

        selected_groups = .false.
        selection_order = 0
        selected_count = 0
        ngroup = size(selected_groups)
        improvement = 0.0_dp
        if (present(min_improvement)) improvement = max(0.0_dp, min_improvement)
        final_stat = 0
        call build_group_design(x, column_groups, selected_groups, 0, design)
        allocate(coefficients(size(design, 2)), fitted(size(y)), residuals(size(y)), hat(size(y), size(y)))
        call glm_irls(design, y, family_code, coefficients, fitted, residuals, hat, deviance, iterations, converged, &
                      max_iterations=max_iterations, tolerance=tolerance, info=stat)
        if (stat /= 0) final_stat = stat
        current_score = deviance + 2.0_dp*real(size(design, 2), dp)
        final_deviance = deviance
        do while (selected_count < ngroup)
            best_group = 0
            best_score = huge(1.0_dp)
            candidate_deviance = final_deviance
            do candidate = 1, ngroup
                if (selected_groups(candidate)) cycle
                call build_group_design(x, column_groups, selected_groups, candidate, design)
                deallocate(coefficients)
                allocate(coefficients(size(design, 2)))
                call glm_irls(design, y, family_code, coefficients, fitted, residuals, hat, deviance, iterations, converged, &
                              max_iterations=max_iterations, tolerance=tolerance, info=stat)
                if (stat /= 0) then
                    final_stat = stat
                    cycle
                end if
                candidate_score = deviance + 2.0_dp*real(size(design, 2), dp)
                if (candidate_score < best_score) then
                    best_score = candidate_score
                    best_group = candidate
                    candidate_deviance = deviance
                end if
            end do
            if (best_group == 0) exit
            if (current_score - best_score <= improvement) exit
            selected_count = selected_count + 1
            selected_groups(best_group) = .true.
            selection_order(selected_count) = best_group
            current_score = best_score
            final_deviance = candidate_deviance
        end do
        if (present(info)) info = final_stat
    end subroutine forward_glm_select

    subroutine bootstrap_linear_refits(x, fitted, residuals, observation_indices, residual_indices, multipliers, &
                                       coefficient_replicates, prediction_design, prediction_replicates, ridge, info)
        real(dp), intent(in) :: x(:, :) !! Original explicit regression design matrix used for each deterministic bootstrap refit.
        real(dp), intent(in) :: fitted(:) !! Original fitted responses aligned with design rows.
        real(dp), intent(in) :: residuals(:) !! Original residual vector used to construct bootstrap responses.
        integer, intent(in) :: observation_indices(:, :) !! Resampled design-row indices [replicate, observation].
        integer, intent(in) :: residual_indices(:, :) !! Bootstrap residual indices with the same shape as observation_indices.
        real(dp), intent(in) :: multipliers(:, :) !! Residual multipliers/noise factors with the same bootstrap shape.
        real(dp), intent(out) :: coefficient_replicates(:, :) !! Refit coefficient vector for every bootstrap replicate.
        real(dp), intent(in), optional :: prediction_design(:, :) !! Optional design matrix for bootstrap predictions.
        real(dp), intent(out), optional :: prediction_replicates(:, :) !! Predictions by replicate and prediction-design row.
        real(dp), intent(in), optional :: ridge !! Optional nonnegative ridge penalty on non-intercept coefficients.
        integer, intent(out), optional :: info !! Zero if all refit solves succeed; otherwise the last nonzero solve code.
        real(dp), allocatable :: xb(:, :)
        real(dp), allocatable :: yb(:)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: beta(:)
        real(dp) :: lambda
        integer :: b
        integer :: i
        integer :: j
        integer :: stat
        integer :: final_stat

        lambda = 0.0_dp
        if (present(ridge)) lambda = max(0.0_dp, ridge)
        allocate(xb(size(observation_indices, 2), size(x, 2)), yb(size(observation_indices, 2)))
        allocate(rhs(size(x, 2)), beta(size(x, 2)))
        final_stat = 0
        do b = 1, size(observation_indices, 1)
            do i = 1, size(observation_indices, 2)
                xb(i, :) = x(observation_indices(b, i), :)
                yb(i) = fitted(observation_indices(b, i)) &
                    + multipliers(b, i)*residuals(residual_indices(b, i))
            end do
            gram = matmul(transpose(xb), xb)
            do j = 2, size(gram, 1)
                gram(j, j) = gram(j, j) + lambda
            end do
            rhs = matmul(transpose(xb), yb)
            call solve_system(gram, rhs, beta, stat)
            if (stat /= 0) final_stat = stat
            coefficient_replicates(b, :) = beta
            if (present(prediction_design) .and. present(prediction_replicates)) then
                prediction_replicates(b, :) = matmul(prediction_design, beta)
            end if
        end do
        if (present(info)) info = final_stat
    end subroutine bootstrap_linear_refits

    subroutine build_group_design(x, column_groups, selected, candidate, design)
        real(dp), intent(in) :: x(:, :) !! Full candidate design matrix whose selected columns are copied to design.
        integer, intent(in) :: column_groups(:) !! Group id for each full design column; group zero is always retained.
        logical, intent(in) :: selected(:) !! Current selection flags for positive group ids.
        integer, intent(in) :: candidate !! Additional positive group id to include temporarily, or zero for no candidate.
        real(dp), allocatable, intent(out) :: design(:, :) !! Reduced design containing base and selected/candidate group columns.
        integer :: count_columns
        integer :: col
        integer :: out_col
        integer :: group
        logical :: keep

        count_columns = 0
        do col = 1, size(x, 2)
            group = column_groups(col)
            keep = group == 0
            if (group > 0 .and. group <= size(selected)) keep = selected(group) .or. group == candidate
            if (keep) count_columns = count_columns + 1
        end do
        allocate(design(size(x, 1), count_columns))
        out_col = 0
        do col = 1, size(x, 2)
            group = column_groups(col)
            keep = group == 0
            if (group > 0 .and. group <= size(selected)) keep = selected(group) .or. group == candidate
            if (keep) then
                out_col = out_col + 1
                design(:, out_col) = x(:, col)
            end if
        end do
    end subroutine build_group_design

    pure subroutine family_moments(eta, family_code, mu, variance, derivative)
        real(dp), intent(in) :: eta(:) !! Linear predictor values.
        integer, intent(in) :: family_code !! Family selector: 1 Gaussian, 2 binomial-logit, or 3 Poisson-log.
        real(dp), intent(out) :: mu(:) !! Inverse-link mean for each linear predictor.
        real(dp), intent(out) :: variance(:) !! Family variance evaluated at each mean.
        real(dp), intent(out) :: derivative(:) !! Derivative d(mu)/d(eta) used in IRLS.

        select case (family_code)
        case (2)
            mu = 1.0_dp/(1.0_dp + exp(-max(-35.0_dp, min(35.0_dp, eta))))
            mu = max(1.0e-10_dp, min(1.0_dp - 1.0e-10_dp, mu))
            variance = mu*(1.0_dp - mu)
            derivative = variance
        case (3)
            mu = exp(max(-30.0_dp, min(30.0_dp, eta)))
            variance = max(mu, 1.0e-12_dp)
            derivative = mu
        case default
            mu = eta
            variance = 1.0_dp
            derivative = 1.0_dp
        end select
    end subroutine family_moments

    pure real(dp) function family_deviance(y, mu, weights, family_code) result(deviance)
        real(dp), intent(in) :: y(:) !! Observed response values.
        real(dp), intent(in) :: mu(:) !! Fitted response means on the response scale.
        real(dp), intent(in) :: weights(:) !! Nonnegative prior observation weights.
        integer, intent(in) :: family_code !! Family selector: 1 Gaussian, 2 binomial-logit, or 3 Poisson-log.
        real(dp) :: term
        integer :: i

        deviance = 0.0_dp
        select case (family_code)
        case (2)
            do i = 1, size(y)
                term = 0.0_dp
                if (y(i) > 0.0_dp) term = term + y(i)*log(y(i)/max(mu(i), 1.0e-12_dp))
                if (y(i) < 1.0_dp) then
                    term = term + (1.0_dp - y(i))*log((1.0_dp - y(i))/max(1.0_dp - mu(i), 1.0e-12_dp))
                end if
                deviance = deviance + 2.0_dp*weights(i)*term
            end do
        case (3)
            do i = 1, size(y)
                if (y(i) > 0.0_dp) then
                    term = y(i)*log(y(i)/max(mu(i), 1.0e-12_dp)) - (y(i) - mu(i))
                else
                    term = mu(i)
                end if
                deviance = deviance + 2.0_dp*weights(i)*term
            end do
        case default
            deviance = sum(weights*(y - mu)**2)
        end select
    end function family_deviance

    pure subroutine two_sample_ks(x, y, statistic, pvalue)
        real(dp), intent(in) :: x(:) !! First scalar projected sample.
        real(dp), intent(in) :: y(:) !! Second scalar projected sample.
        real(dp), intent(out) :: statistic !! Two-sample Kolmogorov-Smirnov sup-norm empirical CDF difference.
        real(dp), intent(out) :: pvalue !! Two-sided asymptotic KS p-value approximation.
        real(dp), allocatable :: xs(:)
        real(dp), allocatable :: ys(:)
        real(dp) :: cdfx
        real(dp) :: cdfy
        real(dp) :: en
        real(dp) :: lambda
        real(dp) :: term
        real(dp) :: sumq
        real(dp) :: value
        integer :: ix
        integer :: iy
        integer :: k

        xs = x
        ys = y
        call sort_real_local(xs)
        call sort_real_local(ys)
        ix = 0
        iy = 0
        statistic = 0.0_dp
        do while (ix < size(xs) .or. iy < size(ys))
            if (iy >= size(ys)) then
                value = xs(ix + 1)
            else if (ix >= size(xs)) then
                value = ys(iy + 1)
            else
                value = min(xs(ix + 1), ys(iy + 1))
            end if
            do while (ix < size(xs))
                if (xs(ix + 1) > value) exit
                ix = ix + 1
            end do
            do while (iy < size(ys))
                if (ys(iy + 1) > value) exit
                iy = iy + 1
            end do
            cdfx = real(ix, dp)/real(size(xs), dp)
            cdfy = real(iy, dp)/real(size(ys), dp)
            statistic = max(statistic, abs(cdfx - cdfy))
        end do
        en = sqrt(real(size(xs)*size(ys), dp)/real(size(xs) + size(ys), dp))
        lambda = (en + 0.12_dp + 0.11_dp/max(en, tiny(1.0_dp)))*statistic
        sumq = 0.0_dp
        do k = 1, 100
            term = 2.0_dp*(-1.0_dp)**(k - 1)*exp(-2.0_dp*real(k*k, dp)*lambda*lambda)
            sumq = sumq + term
            if (abs(term) < 1.0e-14_dp) exit
        end do
        pvalue = max(0.0_dp, min(1.0_dp, sumq))
    end subroutine two_sample_ks

    pure subroutine ordinary_anova_onefactor(y, groups, statistic)
        real(dp), intent(in) :: y(:) !! Scalar response values for a one-factor ordinary ANOVA.
        integer, intent(in) :: groups(:) !! Positive integer factor level for each response observation.
        real(dp), intent(out) :: statistic !! Classical one-way ANOVA F statistic.
        real(dp), allocatable :: means(:)
        integer, allocatable :: counts(:)
        real(dp) :: grand
        real(dp) :: between
        real(dp) :: within
        integer :: g
        integer :: i
        integer :: j

        g = maxval(groups)
        allocate(means(g), counts(g))
        means = 0.0_dp
        counts = 0
        do i = 1, size(y)
            if (groups(i) >= 1 .and. groups(i) <= g) then
                means(groups(i)) = means(groups(i)) + y(i)
                counts(groups(i)) = counts(groups(i)) + 1
            end if
        end do
        do j = 1, g
            if (counts(j) > 0) means(j) = means(j)/real(counts(j), dp)
        end do
        grand = sum(y)/real(size(y), dp)
        between = 0.0_dp
        within = 0.0_dp
        do j = 1, g
            between = between + real(counts(j), dp)*(means(j) - grand)**2
        end do
        do i = 1, size(y)
            if (groups(i) >= 1 .and. groups(i) <= g) within = within + (y(i) - means(groups(i)))**2
        end do
        if (g <= 1 .or. size(y) <= g .or. within <= tiny(1.0_dp)) then
            statistic = 0.0_dp
        else
            statistic = (between/real(g - 1, dp))/(within/real(size(y) - g, dp))
        end if
    end subroutine ordinary_anova_onefactor

    pure real(dp) function quantile_type4_local(sorted, probability) result(q)
        real(dp), intent(in) :: sorted(:) !! Increasing sample values.
        real(dp), intent(in) :: probability !! Quantile probability in [0,1].
        real(dp) :: h
        real(dp) :: frac
        integer :: j
        integer :: n

        n = size(sorted)
        if (probability <= 0.0_dp) then
            q = sorted(1)
        else if (probability >= 1.0_dp) then
            q = sorted(n)
        else
            h = real(n, dp)*probability
            j = floor(h)
            frac = h - real(j, dp)
            if (j <= 0) then
                q = sorted(1)
            else if (j >= n) then
                q = sorted(n)
            else
                q = (1.0_dp - frac)*sorted(j) + frac*sorted(j + 1)
            end if
        end if
    end function quantile_type4_local

    pure subroutine sort_real_local(x)
        real(dp), intent(inout) :: x(:) !! Real vector sorted into increasing order in place.
        real(dp) :: key
        integer :: i
        integer :: j

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_real_local

end module fdausc_generalized
