module fdapace_numerics
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use fdapace_kinds, only : dp
    use fdapace_math, only : gaussian_solve, mean_value, sort_real
    implicit none
    private

    public :: column_correlation
    public :: empirical_quantile
    public :: inverse_matrix
    public :: least_squares
    public :: logistic_regression
    public :: rank_average
    public :: sample_indices_with_replacement

contains

    function empirical_quantile(x, prob) result(value)
        real(dp), intent(in) :: x(:) !! Finite sample whose empirical linear-interpolation quantile is requested.
        real(dp), intent(in) :: prob !! Quantile probability in the closed interval [0,1].
        real(dp) :: value
        real(dp), allocatable :: work(:)
        real(dp) :: pos
        integer :: hi
        integer :: lo
        integer :: n

        n = size(x)
        if (n < 1) error stop "empirical_quantile: sample is empty"
        if (prob < 0.0_dp .or. prob > 1.0_dp) error stop "empirical_quantile: probability must lie in [0,1]"
        allocate(work(n))
        work = x
        call sort_real(work)
        if (n == 1) then
            value = work(1)
            return
        end if
        pos = 1.0_dp + prob * real(n - 1, dp)
        lo = floor(pos)
        hi = ceiling(pos)
        if (lo == hi) then
            value = work(lo)
        else
            value = work(lo) + (pos - real(lo, dp)) * (work(hi) - work(lo))
        end if
    end function empirical_quantile

    subroutine inverse_matrix(a, ainv, info, ridge)
        real(dp), intent(in) :: a(:,:) !! Square matrix to invert by repeated Gaussian solves.
        real(dp), allocatable, intent(out) :: ainv(:,:) !! Matrix inverse, or NaNs if a solve fails.
        integer, intent(out) :: info !! Zero on success; nonzero for non-square or singular systems.
        real(dp), intent(in), optional :: ridge !! Optional nonnegative diagonal ridge added before inversion.
        real(dp), allocatable :: ar(:,:)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: x(:)
        real(dp) :: r
        integer :: i
        integer :: n
        integer :: stat

        n = size(a, 1)
        allocate(ainv(n, size(a, 2)))
        if (size(a, 2) /= n) then
            ainv = ieee_value(0.0_dp, ieee_quiet_nan)
            info = -1
            return
        end if
        r = 0.0_dp
        if (present(ridge)) r = ridge
        if (r < 0.0_dp) error stop "inverse_matrix: ridge must be nonnegative"
        allocate(ar(n, n), b(n), x(n))
        ar = a
        do i = 1, n
            ar(i, i) = ar(i, i) + r
        end do
        ainv = 0.0_dp
        do i = 1, n
            b = 0.0_dp
            b(i) = 1.0_dp
            call gaussian_solve(ar, b, x, stat)
            if (stat /= 0) then
                ainv = ieee_value(0.0_dp, ieee_quiet_nan)
                info = stat
                return
            end if
            ainv(:, i) = x
        end do
        info = 0
    end subroutine inverse_matrix

    subroutine least_squares(x, y, beta, info, ridge)
        real(dp), intent(in) :: x(:,:) !! Design matrix with observations in rows and predictors in columns.
        real(dp), intent(in) :: y(:,:) !! Response matrix with the same observation count as x.
        real(dp), allocatable, intent(out) :: beta(:,:) !! Least-squares coefficient matrix mapping columns of x to columns of y.
        integer, intent(out) :: info !! Zero on success; nonzero if dimensions are invalid or the normal equations are singular.
        real(dp), intent(in), optional :: ridge !! Optional nonnegative diagonal ridge added to X'X.
        real(dp), allocatable :: gram(:,:)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: sol(:)
        real(dp) :: r
        integer :: i
        integer :: j
        integer :: p
        integer :: stat

        if (size(x, 1) /= size(y, 1) .or. size(x, 1) < 1 .or. size(x, 2) < 1) then
            allocate(beta(max(0, size(x, 2)), max(0, size(y, 2))))
            info = -1
            return
        end if
        p = size(x, 2)
        r = 0.0_dp
        if (present(ridge)) r = ridge
        if (r < 0.0_dp) error stop "least_squares: ridge must be nonnegative"
        allocate(gram(p, p), rhs(p), sol(p), beta(p, size(y, 2)))
        gram = matmul(transpose(x), x)
        do i = 1, p
            gram(i, i) = gram(i, i) + r
        end do
        do j = 1, size(y, 2)
            rhs = matmul(transpose(x), y(:, j))
            call gaussian_solve(gram, rhs, sol, stat)
            if (stat /= 0) then
                beta = ieee_value(0.0_dp, ieee_quiet_nan)
                info = stat
                return
            end if
            beta(:, j) = sol
        end do
        info = 0
    end subroutine least_squares

    pure real(dp) function column_correlation(x, y) result(value)
        real(dp), intent(in) :: x(:) !! First finite observation vector.
        real(dp), intent(in) :: y(:) !! Second finite observation vector of the same length as x.
        real(dp) :: mx
        real(dp) :: my
        real(dp) :: sx
        real(dp) :: sy

        if (size(x) /= size(y) .or. size(x) < 2) then
            value = 0.0_dp
            return
        end if
        mx = sum(x) / real(size(x), dp)
        my = sum(y) / real(size(y), dp)
        sx = sqrt(sum((x - mx)**2))
        sy = sqrt(sum((y - my)**2))
        if (sx <= 0.0_dp .or. sy <= 0.0_dp) then
            value = 0.0_dp
        else
            value = sum((x - mx) * (y - my)) / (sx * sy)
        end if
    end function column_correlation

    pure subroutine rank_average(x, ranks)
        real(dp), intent(in) :: x(:) !! Values to rank in ascending order with average ranks for ties.
        real(dp), intent(out) :: ranks(:) !! One-based average ranks, with the same shape as x.
        integer, allocatable :: idx(:)
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: tmp

        n = size(x)
        if (size(ranks) /= n) error stop "rank_average: output shape mismatch"
        allocate(idx(n))
        idx = [(i, i=1, n)]
        do i = 2, n
            tmp = idx(i)
            j = i - 1
            do while (j >= 1)
                if (x(idx(j)) <= x(tmp)) exit
                idx(j + 1) = idx(j)
                j = j - 1
            end do
            idx(j + 1) = tmp
        end do
        i = 1
        do while (i <= n)
            j = i
            do while (j < n)
                if (x(idx(j + 1)) /= x(idx(i))) exit
                j = j + 1
            end do
            do k = i, j
                ranks(idx(k)) = 0.5_dp * real(i + j, dp)
            end do
            i = j + 1
        end do
    end subroutine rank_average

    subroutine sample_indices_with_replacement(n, indices)
        integer, intent(in) :: n !! Positive population size and number of bootstrap draws.
        integer, intent(out) :: indices(:) !! Random one-based indices sampled uniformly with replacement from 1:n.
        real(dp) :: u
        integer :: i

        if (n < 1 .or. size(indices) < 1) error stop "sample_indices_with_replacement: invalid size"
        do i = 1, size(indices)
            call random_number(u)
            indices(i) = min(n, int(u * real(n, dp)) + 1)
        end do
    end subroutine sample_indices_with_replacement

    subroutine logistic_regression(x, y, beta, info, ridge, max_iter, tolerance)
        real(dp), intent(in) :: x(:,:) !! Logistic-regression design matrix with observations in rows, including any desired intercept column.
        real(dp), intent(in) :: y(:) !! Binary response values encoded as zeros and ones.
        real(dp), allocatable, intent(out) :: beta(:) !! Estimated logistic-regression coefficients.
        integer, intent(out) :: info !! Zero on convergence; positive on iteration limit; negative on invalid input or singular solve.
        real(dp), intent(in), optional :: ridge !! Nonnegative ridge penalty for numerical stability; defaults to 1e-8.
        integer, intent(in), optional :: max_iter !! Maximum Newton/IRLS iterations; defaults to 50.
        real(dp), intent(in), optional :: tolerance !! Maximum coefficient-change convergence tolerance; defaults to 1e-8.
        real(dp), allocatable :: eta(:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: gram(:,:)
        real(dp), allocatable :: grad(:)
        real(dp), allocatable :: step(:)
        real(dp) :: r
        real(dp) :: tol
        real(dp) :: w
        integer :: i
        integer :: iter
        integer :: itmax
        integer :: j
        integer :: k
        integer :: stat

        if (size(x, 1) /= size(y) .or. size(x, 1) < 1 .or. size(x, 2) < 1) then
            allocate(beta(max(0, size(x, 2))))
            info = -1
            return
        end if
        if (any(y < 0.0_dp) .or. any(y > 1.0_dp)) error stop "logistic_regression: y must lie in [0,1]"
        r = 1.0e-8_dp
        if (present(ridge)) r = ridge
        itmax = 50
        if (present(max_iter)) itmax = max_iter
        tol = 1.0e-8_dp
        if (present(tolerance)) tol = tolerance
        allocate(beta(size(x, 2)), eta(size(y)), p(size(y)), gram(size(x, 2), size(x, 2)))
        allocate(grad(size(x, 2)), step(size(x, 2)))
        beta = 0.0_dp
        do iter = 1, itmax
            eta = matmul(x, beta)
            do i = 1, size(y)
                if (eta(i) >= 30.0_dp) then
                    p(i) = 1.0_dp
                else if (eta(i) <= -30.0_dp) then
                    p(i) = 0.0_dp
                else
                    p(i) = 1.0_dp / (1.0_dp + exp(-eta(i)))
                end if
            end do
            gram = 0.0_dp
            grad = matmul(transpose(x), y - p)
            do i = 1, size(y)
                w = max(1.0e-10_dp, p(i) * (1.0_dp - p(i)))
                do j = 1, size(x, 2)
                    do k = 1, size(x, 2)
                        gram(j, k) = gram(j, k) + w * x(i, j) * x(i, k)
                    end do
                end do
            end do
            do j = 1, size(x, 2)
                gram(j, j) = gram(j, j) + r
                grad(j) = grad(j) - r * beta(j)
            end do
            call gaussian_solve(gram, grad, step, stat)
            if (stat /= 0) then
                info = -2
                return
            end if
            beta = beta + step
            if (maxval(abs(step)) <= tol) then
                info = 0
                return
            end if
        end do
        info = 1
    end subroutine logistic_regression

end module fdapace_numerics
