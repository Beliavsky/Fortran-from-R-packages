module fdausc_regression
    use r_kinds, only : dp
    use r_linalg, only : solve_system, inverse_matrix
    use fdausc_kernels, only : kernel_value
    implicit none
    private

    public :: weighted_linear_fit, component_regression, nonparametric_regression
    public :: select_bandwidth_cv, regression_diagnostics

contains

    subroutine weighted_linear_fit(x, y, coefficients, fitted, residuals, hat, weights, ridge, info)
        real(dp), intent(in) :: x(:, :) !! Numeric design matrix including any desired intercept column.
        real(dp), intent(in) :: y(:) !! Scalar response vector aligned with design rows.
        real(dp), intent(out) :: coefficients(:) !! Weighted least-squares or ridge coefficient vector.
        real(dp), intent(out) :: fitted(:) !! Fitted response values X*coefficients.
        real(dp), intent(out) :: residuals(:) !! Response residuals y-fitted.
        real(dp), intent(out) :: hat(:, :) !! Linear hat matrix X(X'WX+ridge I)^-1 X'W.
        real(dp), intent(in), optional :: weights(:) !! Nonnegative case weights; all ones if omitted.
        real(dp), intent(in), optional :: ridge !! Nonnegative ridge penalty applied to non-intercept columns, default zero.
        integer, intent(out), optional :: info !! Zero on successful solve/inversion; nonzero value from rfortran-linalg otherwise.
        real(dp), allocatable :: xtwx(:, :)
        real(dp), allocatable :: xtwy(:)
        real(dp), allocatable :: inv(:, :)
        real(dp), allocatable :: wx(:, :)
        real(dp), allocatable :: wy(:)
        real(dp) :: lambda
        integer :: i
        integer :: stat

        lambda = 0.0_dp
        if (present(ridge)) lambda = max(0.0_dp, ridge)
        allocate(wx(size(x, 1), size(x, 2)), wy(size(y)))
        wx = x
        wy = y
        if (present(weights)) then
            do i = 1, size(y)
                wx(i, :) = weights(i)*x(i, :)
                wy(i) = weights(i)*y(i)
            end do
        end if
        xtwx = matmul(transpose(x), wx)
        xtwy = matmul(transpose(x), wy)
        do i = 2, size(xtwx, 1)
            xtwx(i, i) = xtwx(i, i) + lambda
        end do
        call solve_system(xtwx, xtwy, coefficients, stat)
        if (present(info)) info = stat
        if (stat /= 0) then
            fitted = 0.0_dp
            residuals = y
            hat = 0.0_dp
            return
        end if
        fitted = matmul(x, coefficients)
        residuals = y - fitted
        call inverse_matrix(xtwx, inv, stat)
        if (present(info)) info = stat
        if (stat /= 0) then
            hat = 0.0_dp
            return
        end if
        if (present(weights)) then
            hat = matmul(matmul(x, inv), transpose(wx))
        else
            hat = matmul(matmul(x, inv), transpose(x))
        end if
    end subroutine weighted_linear_fit

    subroutine component_regression(scores, y, coefficients, fitted, residuals, hat, weights, ridge, info)
        real(dp), intent(in) :: scores(:, :) !! Selected functional PC or PLS score columns.
        real(dp), intent(in) :: y(:) !! Scalar response vector aligned with score rows.
        real(dp), intent(out) :: coefficients(:) !! Intercept followed by component regression coefficients.
        real(dp), intent(out) :: fitted(:) !! Fitted response values.
        real(dp), intent(out) :: residuals(:) !! Response residuals y-fitted.
        real(dp), intent(out) :: hat(:, :) !! Hat matrix of the component regression.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative case weights.
        real(dp), intent(in), optional :: ridge !! Optional ridge penalty on component coefficients.
        integer, intent(out), optional :: info !! Zero on successful fit; nonzero value from the linear solver otherwise.
        real(dp), allocatable :: x(:, :)
        integer :: stat

        allocate(x(size(scores, 1), size(scores, 2) + 1))
        x(:, 1) = 1.0_dp
        x(:, 2:) = scores
        call weighted_linear_fit(x, y, coefficients, fitted, residuals, hat, weights, ridge, stat)
        if (present(info)) info = stat
    end subroutine component_regression

    pure subroutine nonparametric_regression(dist_train_eval, y, bandwidth, fitted, kernel_code, weights)
        real(dp), intent(in) :: dist_train_eval(:, :) !! Train-to-evaluation distances, with training observations in rows.
        real(dp), intent(in) :: y(:) !! Training responses aligned with distance rows.
        real(dp), intent(in) :: bandwidth !! Positive kernel bandwidth.
        real(dp), intent(out) :: fitted(:) !! Kernel-weighted response prediction for each evaluation column.
        integer, intent(in), optional :: kernel_code !! Symmetric kernel selector; default normal.
        real(dp), intent(in), optional :: weights(:) !! Optional training-observation weights.
        real(dp) :: w
        real(dp) :: total
        real(dp) :: numerator
        integer :: code
        integer :: i
        integer :: j

        code = 1
        if (present(kernel_code)) code = kernel_code
        do j = 1, size(dist_train_eval, 2)
            total = 0.0_dp
            numerator = 0.0_dp
            do i = 1, size(y)
                w = kernel_value(dist_train_eval(i, j)/bandwidth, code)
                if (present(weights)) w = w*weights(i)
                total = total + w
                numerator = numerator + w*y(i)
            end do
            if (total > tiny(1.0_dp)) then
                fitted(j) = numerator/total
            else
                fitted(j) = sum(y)/real(size(y), dp)
            end if
        end do
    end subroutine nonparametric_regression

    pure subroutine select_bandwidth_cv(distances, y, bandwidths, best_bandwidth, scores, kernel_code, weights)
        real(dp), intent(in) :: distances(:, :) !! Square training-distance matrix for leave-one-out bandwidth selection.
        real(dp), intent(in) :: y(:) !! Training responses aligned with distance rows/columns.
        real(dp), intent(in) :: bandwidths(:) !! Positive candidate kernel bandwidths.
        real(dp), intent(out) :: best_bandwidth !! Candidate bandwidth minimizing weighted leave-one-out mean squared error.
        real(dp), intent(out) :: scores(:) !! Leave-one-out mean squared error for each candidate bandwidth.
        integer, intent(in), optional :: kernel_code !! Symmetric kernel selector; default normal.
        real(dp), intent(in), optional :: weights(:) !! Optional case weights in both fitting and validation loss.
        real(dp) :: w
        real(dp) :: total
        real(dp) :: numerator
        real(dp) :: pred
        real(dp) :: loss
        real(dp) :: lossw
        integer :: code
        integer :: h
        integer :: i
        integer :: j
        integer :: best

        code = 1
        if (present(kernel_code)) code = kernel_code
        do h = 1, size(bandwidths)
            loss = 0.0_dp
            lossw = 0.0_dp
            do j = 1, size(y)
                total = 0.0_dp
                numerator = 0.0_dp
                do i = 1, size(y)
                    if (i == j) cycle
                    w = kernel_value(distances(i, j)/bandwidths(h), code)
                    if (present(weights)) w = w*weights(i)
                    total = total + w
                    numerator = numerator + w*y(i)
                end do
                if (total > tiny(1.0_dp)) then
                    pred = numerator/total
                else
                    pred = (sum(y) - y(j))/real(max(1, size(y) - 1), dp)
                end if
                if (present(weights)) then
                    loss = loss + weights(j)*(y(j) - pred)**2
                    lossw = lossw + weights(j)
                else
                    loss = loss + (y(j) - pred)**2
                    lossw = lossw + 1.0_dp
                end if
            end do
            scores(h) = loss/max(lossw, tiny(1.0_dp))
        end do
        best = minloc(scores, dim=1)
        best_bandwidth = bandwidths(best)
    end subroutine select_bandwidth_cv

    pure subroutine regression_diagnostics(y, fitted, hat, residual_variance, r_squared, adjusted_r_squared, df)
        real(dp), intent(in) :: y(:) !! Observed scalar response vector.
        real(dp), intent(in) :: fitted(:) !! Fitted response values aligned with y.
        real(dp), intent(in) :: hat(:, :) !! Hat matrix used to obtain effective degrees of freedom.
        real(dp), intent(out) :: residual_variance !! Residual sum of squares divided by n-trace(hat).
        real(dp), intent(out) :: r_squared !! Coefficient of determination.
        real(dp), intent(out) :: adjusted_r_squared !! Degrees-of-freedom adjusted coefficient of determination.
        real(dp), intent(out) :: df !! Effective degrees of freedom trace(hat).
        real(dp) :: rss
        real(dp) :: tss
        real(dp) :: ym
        integer :: i

        df = 0.0_dp
        do i = 1, min(size(hat, 1), size(hat, 2))
            df = df + hat(i, i)
        end do
        rss = sum((y - fitted)**2)
        ym = sum(y)/real(size(y), dp)
        tss = sum((y - ym)**2)
        residual_variance = rss/max(1.0_dp, real(size(y), dp) - df)
        if (tss <= tiny(1.0_dp)) then
            r_squared = 0.0_dp
            adjusted_r_squared = 0.0_dp
        else
            r_squared = 1.0_dp - rss/tss
            adjusted_r_squared = 1.0_dp - (1.0_dp - r_squared)*real(size(y) - 1, dp) &
                /max(1.0_dp, real(size(y), dp) - df)
        end if
    end subroutine regression_diagnostics

end module fdausc_regression
