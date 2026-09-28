module car_utils
    use car_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
    implicit none
    private

    public :: chi_square_quantile
    public :: chi_square_survival
    public :: f_quantile
    public :: f_survival
    public :: f_survival_real
    public :: least_squares_residuals
    public :: matrix_determinant
    public :: matrix_inverse
    public :: mean_valid
    public :: median_valid
    public :: normal_quantile
    public :: quantile_type7
    public :: sample_variance
    public :: student_t_quantile
    public :: student_t_two_sided_p


contains

    pure subroutine mean_valid(x, mean_value, n_valid, status)
        real(dp), intent(in) :: x(:) !! Numeric sample; NaN values are omitted.
        real(dp), intent(out) :: mean_value !! Arithmetic mean of the non-NaN values.
        integer, intent(out) :: n_valid !! Number of non-NaN observations used.
        integer, intent(out) :: status !! Zero on success; one when no valid observations exist.
        integer :: i
        real(dp) :: total

        total = 0.0_dp
        n_valid = 0
        do i = 1, size(x)
            if (.not. ieee_is_nan(x(i))) then
                total = total + x(i)
                n_valid = n_valid + 1
            end if
        end do

        if (n_valid == 0) then
            mean_value = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
        else
            mean_value = total / real(n_valid, dp)
            status = 0
        end if
    end subroutine mean_valid

    pure subroutine sample_variance(x, variance, n_valid, status)
        real(dp), intent(in) :: x(:) !! Numeric sample; NaN values are omitted.
        real(dp), intent(out) :: variance !! Unbiased sample variance with denominator n - 1.
        integer, intent(out) :: n_valid !! Number of non-NaN observations used.
        integer, intent(out) :: status !! Zero on success; one when fewer than two valid observations exist.
        real(dp) :: mu
        real(dp) :: ss
        integer :: i
        integer :: mean_status

        call mean_valid(x, mu, n_valid, mean_status)
        if (mean_status /= 0 .or. n_valid < 2) then
            variance = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if

        ss = 0.0_dp
        do i = 1, size(x)
            if (.not. ieee_is_nan(x(i))) ss = ss + (x(i) - mu)**2
        end do
        variance = ss / real(n_valid - 1, dp)
        status = 0
    end subroutine sample_variance

    pure subroutine sort_real_in_place(x)
        real(dp), intent(inout) :: x(:) !! Real vector sorted into ascending order in place.
        integer :: i
        integer :: j
        real(dp) :: key

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
    end subroutine sort_real_in_place

    pure subroutine median_valid(x, median_value, n_valid, status)
        real(dp), intent(in) :: x(:) !! Numeric sample; NaN values are omitted.
        real(dp), intent(out) :: median_value !! Median of the non-NaN observations.
        integer, intent(out) :: n_valid !! Number of non-NaN observations used.
        integer, intent(out) :: status !! Zero on success; one when no valid observations exist.
        real(dp), allocatable :: work(:)
        integer :: i
        integer :: j

        n_valid = count(.not. ieee_is_nan(x))
        if (n_valid == 0) then
            median_value = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if

        allocate(work(n_valid))
        j = 0
        do i = 1, size(x)
            if (.not. ieee_is_nan(x(i))) then
                j = j + 1
                work(j) = x(i)
            end if
        end do
        call sort_real_in_place(work)

        if (mod(n_valid, 2) == 1) then
            median_value = work((n_valid + 1) / 2)
        else
            median_value = 0.5_dp * (work(n_valid / 2) + work(n_valid / 2 + 1))
        end if
        status = 0
    end subroutine median_valid

    pure subroutine quantile_type7(x, probability, quantile, status)
        real(dp), intent(in) :: x(:) !! Numeric sample; NaN values are omitted.
        real(dp), intent(in) :: probability !! Quantile probability in the closed interval [0, 1].
        real(dp), intent(out) :: quantile !! R type-7 sample quantile.
        integer, intent(out) :: status !! Zero on success; nonzero for invalid probability or empty sample.
        real(dp), allocatable :: work(:)
        real(dp) :: h
        real(dp) :: gamma
        integer :: i
        integer :: j
        integer :: n
        integer :: k

        if (probability < 0.0_dp .or. probability > 1.0_dp) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if

        n = count(.not. ieee_is_nan(x))
        if (n == 0) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 2
            return
        end if

        allocate(work(n))
        j = 0
        do i = 1, size(x)
            if (.not. ieee_is_nan(x(i))) then
                j = j + 1
                work(j) = x(i)
            end if
        end do
        call sort_real_in_place(work)

        if (n == 1) then
            quantile = work(1)
            status = 0
            return
        end if

        h = 1.0_dp + real(n - 1, dp) * probability
        k = int(floor(h))
        gamma = h - real(k, dp)
        if (k <= 1) then
            quantile = work(1) + gamma * (work(2) - work(1))
        else if (k >= n) then
            quantile = work(n)
        else
            quantile = work(k) + gamma * (work(k + 1) - work(k))
        end if
        status = 0
    end subroutine quantile_type7

    pure subroutine matrix_determinant(a, determinant, status)
        real(dp), intent(in) :: a(:, :) !! Square real matrix whose determinant is required.
        real(dp), intent(out) :: determinant !! Determinant of a.
        integer, intent(out) :: status !! Zero on success; nonzero for nonsquare or numerically singular input.
        real(dp), allocatable :: work(:, :)
        real(dp) :: pivot
        real(dp) :: scale
        real(dp) :: temp
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: pivot_row
        integer :: sign_det

        n = size(a, 1)
        if (size(a, 2) /= n) then
            determinant = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if
        if (n == 0) then
            determinant = 1.0_dp
            status = 0
            return
        end if

        allocate(work(n, n))
        work = a
        determinant = 1.0_dp
        sign_det = 1
        scale = max(1.0_dp, maxval(abs(work)))

        do k = 1, n
            pivot_row = k
            pivot = abs(work(k, k))
            do i = k + 1, n
                if (abs(work(i, k)) > pivot) then
                    pivot = abs(work(i, k))
                    pivot_row = i
                end if
            end do
            if (pivot <= epsilon(1.0_dp) * scale * real(max(1, n), dp)) then
                determinant = 0.0_dp
                status = 2
                return
            end if
            if (pivot_row /= k) then
                do j = 1, n
                    temp = work(k, j)
                    work(k, j) = work(pivot_row, j)
                    work(pivot_row, j) = temp
                end do
                sign_det = -sign_det
            end if

            determinant = determinant * work(k, k)
            do i = k + 1, n
                temp = work(i, k) / work(k, k)
                work(i, k) = 0.0_dp
                do j = k + 1, n
                    work(i, j) = work(i, j) - temp * work(k, j)
                end do
            end do
        end do
        determinant = real(sign_det, dp) * determinant
        status = 0
    end subroutine matrix_determinant

    pure subroutine matrix_inverse(a, inverse, status)
        real(dp), intent(in) :: a(:, :) !! Square real matrix to invert.
        real(dp), intent(out) :: inverse(:, :) !! Matrix inverse of a; same shape as a.
        integer, intent(out) :: status !! Zero on success; nonzero for shape mismatch or numerical singularity.
        real(dp), allocatable :: aug(:, :)
        real(dp) :: pivot
        real(dp) :: scale
        real(dp) :: temp
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: pivot_row

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(inverse, 1) /= n .or. size(inverse, 2) /= n) then
            status = 1
            if (size(inverse) > 0) inverse = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        allocate(aug(n, 2 * n))
        aug(:, 1:n) = a
        aug(:, n + 1:2 * n) = 0.0_dp
        do i = 1, n
            aug(i, n + i) = 1.0_dp
        end do
        scale = max(1.0_dp, maxval(abs(a)))

        do k = 1, n
            pivot_row = k
            pivot = abs(aug(k, k))
            do i = k + 1, n
                if (abs(aug(i, k)) > pivot) then
                    pivot = abs(aug(i, k))
                    pivot_row = i
                end if
            end do
            if (pivot <= epsilon(1.0_dp) * scale * real(max(1, n), dp)) then
                inverse = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 2
                return
            end if

            if (pivot_row /= k) then
                do j = 1, 2 * n
                    temp = aug(k, j)
                    aug(k, j) = aug(pivot_row, j)
                    aug(pivot_row, j) = temp
                end do
            end if

            temp = aug(k, k)
            aug(k, :) = aug(k, :) / temp
            do i = 1, n
                if (i == k) cycle
                temp = aug(i, k)
                if (abs(temp) > tiny(1.0_dp)) aug(i, :) = aug(i, :) - temp * aug(k, :)
            end do
        end do

        inverse = aug(:, n + 1:2 * n)
        status = 0
    end subroutine matrix_inverse

    pure subroutine solve_linear_system(a, b, x, status)
        real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector with length equal to the order of a.
        real(dp), intent(out) :: x(:) !! Solution vector with the same length as b.
        integer, intent(out) :: status !! Zero on success; nonzero for shape mismatch or singular matrix.
        real(dp), allocatable :: inverse(:, :)
        integer :: n

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
            status = 1
            if (size(x) > 0) x = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        allocate(inverse(n, n))
        call matrix_inverse(a, inverse, status)
        if (status /= 0) then
            x = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        x = matmul(inverse, b)
    end subroutine solve_linear_system

    pure subroutine least_squares_residuals(x, y, beta, residuals, status, weights)
        real(dp), intent(in) :: x(:, :) !! Design matrix with observations in rows and coefficients in columns.
        real(dp), intent(in) :: y(:) !! Response vector with one value per row of x.
        real(dp), intent(out) :: beta(:) !! Least-squares coefficient vector with length size(x, 2).
        real(dp), intent(out) :: residuals(:) !! Residuals on the weighted scale used in the fit.
        integer, intent(out) :: status !! Zero on success; nonzero for invalid shapes, weights, or singular normal equations.
        real(dp), intent(in), optional :: weights(:) !! Optional nonnegative case weights; omitted means unit weights.
        real(dp), allocatable :: xw(:, :)
        real(dp), allocatable :: yw(:)
        real(dp), allocatable :: xtx(:, :)
        real(dp), allocatable :: xty(:)
        real(dp) :: root_w
        integer :: i
        integer :: n
        integer :: p

        n = size(x, 1)
        p = size(x, 2)
        if (size(y) /= n .or. size(beta) /= p .or. size(residuals) /= n) then
            status = 1
            beta = ieee_value(0.0_dp, ieee_quiet_nan)
            residuals = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (present(weights)) then
            if (size(weights) /= n .or. any(weights < 0.0_dp) .or. any(ieee_is_nan(weights))) then
                status = 2
                beta = ieee_value(0.0_dp, ieee_quiet_nan)
                residuals = ieee_value(0.0_dp, ieee_quiet_nan)
                return
            end if
        end if

        allocate(xw(n, p), yw(n), xtx(p, p), xty(p))
        do i = 1, n
            if (present(weights)) then
                root_w = sqrt(weights(i))
            else
                root_w = 1.0_dp
            end if
            xw(i, :) = root_w * x(i, :)
            yw(i) = root_w * y(i)
        end do

        xtx = matmul(transpose(xw), xw)
        xty = matmul(transpose(xw), yw)
        call solve_linear_system(xtx, xty, beta, status)
        if (status /= 0) then
            residuals = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        residuals = yw - matmul(xw, beta)
    end subroutine least_squares_residuals

    pure function beta_continued_fraction(a, b, x) result(cf)
        real(dp), intent(in) :: a !! First beta-shape parameter; must be positive.
        real(dp), intent(in) :: b !! Second beta-shape parameter; must be positive.
        real(dp), intent(in) :: x !! Evaluation point in [0, 1].
        real(dp) :: cf
        integer, parameter :: max_iter = 400
        real(dp), parameter :: fpmin = tiny(1.0_dp) / epsilon(1.0_dp)
        real(dp), parameter :: tol = 8.0_dp * epsilon(1.0_dp)
        real(dp) :: aa
        real(dp) :: c
        real(dp) :: d
        real(dp) :: del
        real(dp) :: h
        real(dp) :: qab
        real(dp) :: qam
        real(dp) :: qap
        integer :: m
        integer :: m2

        qab = a + b
        qap = a + 1.0_dp
        qam = a - 1.0_dp
        c = 1.0_dp
        d = 1.0_dp - qab * x / qap
        if (abs(d) < fpmin) d = fpmin
        d = 1.0_dp / d
        h = d

        do m = 1, max_iter
            m2 = 2 * m
            aa = real(m, dp) * (b - real(m, dp)) * x / &
                ((qam + real(m2, dp)) * (a + real(m2, dp)))
            d = 1.0_dp + aa * d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa / c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            h = h * d * c

            aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
                ((a + real(m2, dp)) * (qap + real(m2, dp)))
            d = 1.0_dp + aa * d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa / c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            del = d * c
            h = h * del
            if (abs(del - 1.0_dp) <= tol) exit
        end do
        cf = h
    end function beta_continued_fraction

    pure function regularized_beta(x, a, b) result(value)
        real(dp), intent(in) :: x !! Evaluation point in [0, 1].
        real(dp), intent(in) :: a !! First positive beta-shape parameter.
        real(dp), intent(in) :: b !! Second positive beta-shape parameter.
        real(dp) :: value
        real(dp) :: bt

        if (x <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (x >= 1.0_dp) then
            value = 1.0_dp
            return
        end if
        if (a <= 0.0_dp .or. b <= 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) + b * log(1.0_dp - x))
        if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
            value = bt * beta_continued_fraction(a, b, x) / a
        else
            value = 1.0_dp - bt * beta_continued_fraction(b, a, 1.0_dp - x) / b
        end if
        value = max(0.0_dp, min(1.0_dp, value))
    end function regularized_beta

    pure function regularized_gamma_q(a, x) result(q)
        real(dp), intent(in) :: a !! Positive gamma-shape parameter.
        real(dp), intent(in) :: x !! Nonnegative evaluation point.
        real(dp) :: q
        integer, parameter :: max_iter = 500
        real(dp), parameter :: tol = 8.0_dp * epsilon(1.0_dp)
        real(dp), parameter :: fpmin = tiny(1.0_dp) / epsilon(1.0_dp)
        real(dp) :: an
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        real(dp) :: del
        real(dp) :: h
        real(dp) :: p
        real(dp) :: sum_series
        real(dp) :: term
        integer :: i

        if (a <= 0.0_dp .or. x < 0.0_dp) then
            q = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (abs(x) <= tiny(1.0_dp)) then
            q = 1.0_dp
            return
        end if

        if (x < a + 1.0_dp) then
            term = 1.0_dp / a
            sum_series = term
            do i = 1, max_iter
                term = term * x / (a + real(i, dp))
                sum_series = sum_series + term
                if (abs(term) <= abs(sum_series) * tol) exit
            end do
            p = sum_series * exp(-x + a * log(x) - log_gamma(a))
            q = 1.0_dp - p
        else
            b = x + 1.0_dp - a
            c = 1.0_dp / fpmin
            d = 1.0_dp / max(fpmin, b)
            h = d
            do i = 1, max_iter
                an = -real(i, dp) * (real(i, dp) - a)
                b = b + 2.0_dp
                d = an * d + b
                if (abs(d) < fpmin) d = fpmin
                c = b + an / c
                if (abs(c) < fpmin) c = fpmin
                d = 1.0_dp / d
                del = d * c
                h = h * del
                if (abs(del - 1.0_dp) <= tol) exit
            end do
            q = exp(-x + a * log(x) - log_gamma(a)) * h
        end if
        q = max(0.0_dp, min(1.0_dp, q))
    end function regularized_gamma_q

    pure function chi_square_survival(x, degrees_freedom) result(probability)
        real(dp), intent(in) :: x !! Chi-square statistic; negative values map to survival probability one.
        integer, intent(in) :: degrees_freedom !! Positive chi-square degrees of freedom.
        real(dp) :: probability

        if (degrees_freedom <= 0) then
            probability = ieee_value(0.0_dp, ieee_quiet_nan)
        else if (x <= 0.0_dp) then
            probability = 1.0_dp
        else
            probability = regularized_gamma_q(0.5_dp * real(degrees_freedom, dp), 0.5_dp * x)
        end if
    end function chi_square_survival

    pure function f_survival(x, numerator_df, denominator_df) result(probability)
        real(dp), intent(in) :: x !! F statistic; negative values map to survival probability one.
        integer, intent(in) :: numerator_df !! Positive integer numerator degrees of freedom.
        integer, intent(in) :: denominator_df !! Positive integer denominator degrees of freedom.
        real(dp) :: probability

        probability = f_survival_real(x, real(numerator_df, dp), real(denominator_df, dp))
    end function f_survival

    pure function f_survival_real(x, numerator_df, denominator_df) result(probability)
        real(dp), intent(in) :: x !! F statistic; negative values map to survival probability one.
        real(dp), intent(in) :: numerator_df !! Positive real-valued numerator degrees of freedom.
        real(dp), intent(in) :: denominator_df !! Positive real-valued denominator degrees of freedom.
        real(dp) :: probability
        real(dp) :: z

        if (numerator_df <= 0.0_dp .or. denominator_df <= 0.0_dp) then
            probability = ieee_value(0.0_dp, ieee_quiet_nan)
        else if (x <= 0.0_dp) then
            probability = 1.0_dp
        else
            z = denominator_df / (denominator_df + numerator_df * x)
            probability = regularized_beta(z, 0.5_dp * denominator_df, 0.5_dp * numerator_df)
        end if
    end function f_survival_real

    pure function student_t_two_sided_p(t_statistic, degrees_freedom) result(probability)
        real(dp), intent(in) :: t_statistic !! Student-t statistic; sign is ignored.
        integer, intent(in) :: degrees_freedom !! Positive Student-t degrees of freedom.
        real(dp) :: probability
        real(dp) :: z

        if (degrees_freedom <= 0) then
            probability = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        z = real(degrees_freedom, dp) / &
            (real(degrees_freedom, dp) + t_statistic * t_statistic)
        probability = regularized_beta(z, 0.5_dp * real(degrees_freedom, dp), 0.5_dp)
    end function student_t_two_sided_p

    pure subroutine chi_square_quantile(probability, degrees_freedom, quantile, status)
        real(dp), intent(in) :: probability !! Lower-tail probability strictly between zero and one.
        integer, intent(in) :: degrees_freedom !! Positive chi-square degrees of freedom.
        real(dp), intent(out) :: quantile !! Chi-square quantile satisfying P(X <= quantile) = probability.
        integer, intent(out) :: status !! Zero on success; nonzero for invalid arguments or bracketing failure.
        real(dp) :: cdf_mid
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        integer :: i

        if (probability <= 0.0_dp .or. probability >= 1.0_dp .or. degrees_freedom <= 0) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if

        lo = 0.0_dp
        hi = max(1.0_dp, real(degrees_freedom, dp))
        do i = 1, 100
            if (1.0_dp - chi_square_survival(hi, degrees_freedom) >= probability) exit
            hi = 2.0_dp * hi
        end do
        if (i > 100) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 2
            return
        end if

        do i = 1, 160
            mid = 0.5_dp * (lo + hi)
            cdf_mid = 1.0_dp - chi_square_survival(mid, degrees_freedom)
            if (cdf_mid < probability) then
                lo = mid
            else
                hi = mid
            end if
            if (abs(hi - lo) <= 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        quantile = 0.5_dp * (lo + hi)
        status = 0
    end subroutine chi_square_quantile


    pure function normal_quantile(probability) result(quantile)
        real(dp), intent(in) :: probability !! Lower-tail standard-normal probability strictly between zero and one.
        real(dp) :: quantile
        real(dp), parameter :: a1 = -3.969683028665376e1_dp
        real(dp), parameter :: a2 = 2.209460984245205e2_dp
        real(dp), parameter :: a3 = -2.759285104469687e2_dp
        real(dp), parameter :: a4 = 1.383577518672690e2_dp
        real(dp), parameter :: a5 = -3.066479806614716e1_dp
        real(dp), parameter :: a6 = 2.506628277459239_dp
        real(dp), parameter :: b1 = -5.447609879822406e1_dp
        real(dp), parameter :: b2 = 1.615858368580409e2_dp
        real(dp), parameter :: b3 = -1.556989798598866e2_dp
        real(dp), parameter :: b4 = 6.680131188771972e1_dp
        real(dp), parameter :: b5 = -1.328068155288572e1_dp
        real(dp), parameter :: c1 = -7.784894002430293e-3_dp
        real(dp), parameter :: c2 = -3.223964580411365e-1_dp
        real(dp), parameter :: c3 = -2.400758277161838_dp
        real(dp), parameter :: c4 = -2.549732539343734_dp
        real(dp), parameter :: c5 = 4.374664141464968_dp
        real(dp), parameter :: c6 = 2.938163982698783_dp
        real(dp), parameter :: d1 = 7.784695709041462e-3_dp
        real(dp), parameter :: d2 = 3.224671290700398e-1_dp
        real(dp), parameter :: d3 = 2.445134137142996_dp
        real(dp), parameter :: d4 = 3.754408661907416_dp
        real(dp), parameter :: p_low = 0.02425_dp
        real(dp), parameter :: p_high = 1.0_dp - p_low
        real(dp) :: q
        real(dp) :: r

        if (probability <= 0.0_dp .or. probability >= 1.0_dp) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
        else if (probability < p_low) then
            q = sqrt(-2.0_dp * log(probability))
            quantile = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
                ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
        else if (probability <= p_high) then
            q = probability - 0.5_dp
            r = q * q
            quantile = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / &
                (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0_dp)
        else
            q = sqrt(-2.0_dp * log(1.0_dp - probability))
            quantile = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
                ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
        end if
    end function normal_quantile

    pure subroutine student_t_quantile(probability, degrees_freedom, quantile, status)
        real(dp), intent(in) :: probability !! Lower-tail probability strictly between zero and one.
        integer, intent(in) :: degrees_freedom !! Positive Student-t degrees of freedom.
        real(dp), intent(out) :: quantile !! Student-t quantile satisfying P(T <= quantile) = probability.
        integer, intent(out) :: status !! Zero on success; nonzero for invalid arguments or bracketing failure.
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        real(dp) :: target_two_sided
        real(dp) :: sign_q
        integer :: i

        if (probability <= 0.0_dp .or. probability >= 1.0_dp .or. degrees_freedom <= 0) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if
        if (abs(probability - 0.5_dp) <= epsilon(1.0_dp)) then
            quantile = 0.0_dp
            status = 0
            return
        end if

        if (probability > 0.5_dp) then
            sign_q = 1.0_dp
            target_two_sided = 2.0_dp * (1.0_dp - probability)
        else
            sign_q = -1.0_dp
            target_two_sided = 2.0_dp * probability
        end if

        lo = 0.0_dp
        hi = 1.0_dp
        do i = 1, 100
            if (student_t_two_sided_p(hi, degrees_freedom) <= target_two_sided) exit
            hi = 2.0_dp * hi
        end do
        if (i > 100) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 2
            return
        end if

        do i = 1, 160
            mid = 0.5_dp * (lo + hi)
            if (student_t_two_sided_p(mid, degrees_freedom) > target_two_sided) then
                lo = mid
            else
                hi = mid
            end if
            if (abs(hi - lo) <= 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        quantile = sign_q * 0.5_dp * (lo + hi)
        status = 0
    end subroutine student_t_quantile

    pure subroutine f_quantile(probability, numerator_df, denominator_df, quantile, status)
        real(dp), intent(in) :: probability !! Lower-tail probability strictly between zero and one.
        integer, intent(in) :: numerator_df !! Positive numerator degrees of freedom.
        integer, intent(in) :: denominator_df !! Positive denominator degrees of freedom.
        real(dp), intent(out) :: quantile !! F quantile satisfying P(F <= quantile) = probability.
        integer, intent(out) :: status !! Zero on success; nonzero for invalid arguments or bracketing failure.
        real(dp) :: cdf_mid
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        integer :: i

        if (probability <= 0.0_dp .or. probability >= 1.0_dp .or. numerator_df <= 0 .or. denominator_df <= 0) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if

        lo = 0.0_dp
        hi = 1.0_dp
        do i = 1, 100
            if (1.0_dp - f_survival(hi, numerator_df, denominator_df) >= probability) exit
            hi = 2.0_dp * hi
        end do
        if (i > 100) then
            quantile = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 2
            return
        end if

        do i = 1, 160
            mid = 0.5_dp * (lo + hi)
            cdf_mid = 1.0_dp - f_survival(mid, numerator_df, denominator_df)
            if (cdf_mid < probability) then
                lo = mid
            else
                hi = mid
            end if
            if (abs(hi - lo) <= 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        quantile = 0.5_dp * (lo + hi)
        status = 0
    end subroutine f_quantile


end module car_utils
