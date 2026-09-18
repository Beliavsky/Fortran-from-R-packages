module otrimle_linalg
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use otrimle_kinds, only : dp
    use otrimle_rng, only : rand_normal
    implicit none
    private

    real(dp), parameter, public :: pi_dp = acos(-1.0_dp)

    public :: symmetric_eigen
    public :: covariance_matrix
    public :: pchisq
    public :: qchisq
    public :: sort_real
    public :: sort_real_index
    public :: count_unique_rows
    public :: mvn_sample
    public :: median_real
    public :: scale_tau2_location_scale

contains

    subroutine symmetric_eigen(a, values, vectors)
        real(dp), intent(in) :: a(:, :) !! Real symmetric matrix to decompose; shape p by p.
        real(dp), intent(out) :: values(:) !! Eigenvalues in ascending order; length p.
        real(dp), intent(out) :: vectors(:, :) !! Orthonormal eigenvectors in columns; shape p by p.
        real(dp), allocatable :: b(:, :)
        real(dp), allocatable :: col(:)
        real(dp) :: app
        real(dp) :: aqq
        real(dp) :: apq
        real(dp) :: c
        real(dp) :: s
        real(dp) :: t
        real(dp) :: tau
        real(dp) :: temp
        real(dp) :: tol
        real(dp) :: tv
        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer :: pidx
        integer :: qidx
        integer :: sweep

        n = size(a, 1)
        if (size(a, 2) /= n) error stop 'symmetric_eigen: matrix must be square'
        if (size(values) /= n) error stop 'symmetric_eigen: bad eigenvalue size'
        if (size(vectors, 1) /= n .or. size(vectors, 2) /= n) error stop 'symmetric_eigen: bad eigenvector shape'
        allocate(b(n, n), col(n))
        b = 0.5_dp * (a + transpose(a))
        vectors = 0.0_dp
        do i = 1, n
            vectors(i, i) = 1.0_dp
        end do
        if (n == 1) then
            values(1) = b(1, 1)
            return
        end if
        tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(b)))
        do sweep = 1, max(40, 60 * n * n)
            apq = 0.0_dp
            pidx = 1
            qidx = 2
            do j = 2, n
                do i = 1, j - 1
                    if (abs(b(i, j)) > abs(apq)) then
                        apq = b(i, j)
                        pidx = i
                        qidx = j
                    end if
                end do
            end do
            if (abs(apq) <= tol) exit
            app = b(pidx, pidx)
            aqq = b(qidx, qidx)
            tau = (aqq - app) / (2.0_dp * apq)
            if (tau >= 0.0_dp) then
                t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
            else
                t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
            end if
            c = 1.0_dp / sqrt(1.0_dp + t * t)
            s = t * c
            do j = 1, n
                if (j /= pidx .and. j /= qidx) then
                    temp = b(j, pidx)
                    b(j, pidx) = c * temp - s * b(j, qidx)
                    b(pidx, j) = b(j, pidx)
                    b(j, qidx) = s * temp + c * b(j, qidx)
                    b(qidx, j) = b(j, qidx)
                end if
            end do
            b(pidx, pidx) = c * c * app - 2.0_dp * c * s * apq + s * s * aqq
            b(qidx, qidx) = s * s * app + 2.0_dp * c * s * apq + c * c * aqq
            b(pidx, qidx) = 0.0_dp
            b(qidx, pidx) = 0.0_dp
            do j = 1, n
                temp = vectors(j, pidx)
                vectors(j, pidx) = c * temp - s * vectors(j, qidx)
                vectors(j, qidx) = s * temp + c * vectors(j, qidx)
            end do
        end do
        do i = 1, n
            values(i) = b(i, i)
        end do
        do i = 1, n - 1
            m = i
            do j = i + 1, n
                if (values(j) < values(m)) m = j
            end do
            if (m /= i) then
                tv = values(i)
                values(i) = values(m)
                values(m) = tv
                col = vectors(:, i)
                vectors(:, i) = vectors(:, m)
                vectors(:, m) = col
            end if
        end do
    end subroutine symmetric_eigen

    subroutine covariance_matrix(x, cov, center, unbiased)
        real(dp), intent(in) :: x(:, :) !! Observations in rows and variables in columns; shape n by p.
        real(dp), intent(out) :: cov(:, :) !! Covariance matrix; shape p by p.
        real(dp), intent(out), optional :: center(:) !! Optional arithmetic mean vector; length p.
        logical, intent(in), optional :: unbiased !! If true divide by n-1, matching R cov(); otherwise divide by n.
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: z(:, :)
        real(dp) :: den
        integer :: n
        integer :: p
        logical :: unb

        n = size(x, 1)
        p = size(x, 2)
        if (size(cov, 1) /= p .or. size(cov, 2) /= p) error stop 'covariance_matrix: shape mismatch'
        allocate(mu(p), z(n, p))
        if (n == 0) then
            cov = 0.0_dp
            if (present(center)) center = 0.0_dp
            return
        end if
        mu = sum(x, dim=1) / real(n, dp)
        z = x - spread(mu, 1, n)
        unb = .false.
        if (present(unbiased)) unb = unbiased
        if (unb .and. n > 1) then
            den = real(n - 1, dp)
        else
            den = real(max(n, 1), dp)
        end if
        cov = matmul(transpose(z), z) / den
        if (present(center)) center = mu
    end subroutine covariance_matrix

    real(dp) function pchisq(x, df) result(p)
        real(dp), intent(in) :: x !! Chi-square variate; negative values have lower-tail probability zero.
        real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.

        if (x <= 0.0_dp) then
            p = 0.0_dp
        else
            p = regularized_gamma_p(0.5_dp * df, 0.5_dp * x)
        end if
    end function pchisq

    real(dp) function qchisq(prob, df) result(x)
        real(dp), intent(in) :: prob !! Lower-tail probability in the closed interval [0,1].
        real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
        real(dp) :: hi
        real(dp) :: lo
        real(dp) :: mid
        integer :: iter

        if (prob <= 0.0_dp) then
            x = 0.0_dp
            return
        end if
        if (prob >= 1.0_dp) then
            x = huge(1.0_dp)
            return
        end if
        lo = 0.0_dp
        hi = max(df, 1.0_dp)
        do while (pchisq(hi, df) < prob)
            hi = 2.0_dp * hi
            if (.not. ieee_is_finite(hi)) exit
        end do
        do iter = 1, 120
            mid = 0.5_dp * (lo + hi)
            if (pchisq(mid, df) < prob) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x = 0.5_dp * (lo + hi)
    end function qchisq

    real(dp) function regularized_gamma_p(a, x) result(pval)
        real(dp), intent(in) :: a !! Positive gamma shape parameter.
        real(dp), intent(in) :: x !! Nonnegative gamma argument.
        real(dp) :: an
        real(dp) :: ap
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        real(dp) :: del
        real(dp) :: h
        real(dp) :: sum_series
        integer :: i

        if (x <= 0.0_dp) then
            pval = 0.0_dp
            return
        end if
        if (x < a + 1.0_dp) then
            ap = a
            sum_series = 1.0_dp / a
            del = sum_series
            do i = 1, 10000
                ap = ap + 1.0_dp
                del = del * x / ap
                sum_series = sum_series + del
                if (abs(del) <= abs(sum_series) * 10.0_dp * epsilon(1.0_dp)) exit
            end do
            pval = sum_series * exp(-x + a * log(x) - log_gamma(a))
        else
            b = x + 1.0_dp - a
            c = 1.0_dp / tiny(1.0_dp)
            d = 1.0_dp / b
            h = d
            do i = 1, 10000
                an = -real(i, dp) * (real(i, dp) - a)
                b = b + 2.0_dp
                d = an * d + b
                if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
                c = b + an / c
                if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
                d = 1.0_dp / d
                del = d * c
                h = h * del
                if (abs(del - 1.0_dp) <= 10.0_dp * epsilon(1.0_dp)) exit
            end do
            pval = 1.0_dp - exp(-x + a * log(x) - log_gamma(a)) * h
        end if
        pval = min(1.0_dp, max(0.0_dp, pval))
    end function regularized_gamma_p

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:) !! Real vector sorted into ascending order in place.
        integer :: i
        integer :: j
        integer :: m
        real(dp) :: t

        do i = 1, size(x) - 1
            m = i
            do j = i + 1, size(x)
                if (x(j) < x(m)) m = j
            end do
            if (m /= i) then
                t = x(i)
                x(i) = x(m)
                x(m) = t
            end if
        end do
    end subroutine sort_real

    subroutine sort_real_index(x, order, ascending)
        real(dp), intent(in) :: x(:) !! Values whose permutation into sorted order is requested.
        integer, intent(out) :: order(:) !! Permutation indices of length size(x).
        logical, intent(in), optional :: ascending !! True for ascending order; false for descending order.
        integer :: i
        integer :: j
        integer :: m
        integer :: t
        logical :: asc

        if (size(order) /= size(x)) error stop 'sort_real_index: size mismatch'
        asc = .true.
        if (present(ascending)) asc = ascending
        do i = 1, size(x)
            order(i) = i
        end do
        do i = 1, size(x) - 1
            m = i
            do j = i + 1, size(x)
                if ((asc .and. x(order(j)) < x(order(m))) .or. &
                    (.not. asc .and. x(order(j)) > x(order(m)))) m = j
            end do
            if (m /= i) then
                t = order(i)
                order(i) = order(m)
                order(m) = t
            end if
        end do
    end subroutine sort_real_index

    integer function count_unique_rows(x, tol) result(nunique)
        real(dp), intent(in) :: x(:, :) !! Matrix whose distinct rows are counted.
        real(dp), intent(in), optional :: tol !! Absolute row-comparison tolerance; default is exact equality.
        real(dp) :: threshold
        integer :: i
        integer :: j
        logical :: seen

        threshold = 0.0_dp
        if (present(tol)) threshold = max(0.0_dp, tol)
        nunique = 0
        do i = 1, size(x, 1)
            seen = .false.
            do j = 1, i - 1
                if (maxval(abs(x(i, :) - x(j, :))) <= threshold) then
                    seen = .true.
                    exit
                end if
            end do
            if (.not. seen) nunique = nunique + 1
        end do
    end function count_unique_rows

    subroutine mvn_sample(mean, cov, x)
        real(dp), intent(in) :: mean(:) !! Mean vector of the multivariate normal distribution.
        real(dp), intent(in) :: cov(:, :) !! Symmetric positive-semidefinite covariance matrix.
        real(dp), intent(out) :: x(:) !! Generated multivariate-normal observation.
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: vec(:, :)
        real(dp), allocatable :: z(:)
        integer :: i
        integer :: p

        p = size(mean)
        allocate(eig(p), vec(p, p), z(p))
        call symmetric_eigen(cov, eig, vec)
        eig = max(eig, 0.0_dp)
        do i = 1, p
            z(i) = rand_normal()
        end do
        x = mean + matmul(vec, sqrt(eig) * z)
    end subroutine mvn_sample

    real(dp) function median_real(x) result(med)
        real(dp), intent(in) :: x(:) !! Real sample whose ordinary median is requested.
        real(dp), allocatable :: y(:)
        integer :: n

        n = size(x)
        if (n == 0) then
            med = 0.0_dp
            return
        end if
        allocate(y(n))
        y = x
        call sort_real(y)
        if (mod(n, 2) == 1) then
            med = y((n + 1) / 2)
        else
            med = 0.5_dp * (y(n / 2) + y(n / 2 + 1))
        end if
    end function median_real


    subroutine scale_tau2_location_scale(x, location, scale, c1, c2, iterations, tolerance)
        ! Numerical translation of robustbase::scaleTau2(..., mu.too=TRUE).
        ! The default constants and consistency correction follow robustbase/R/OGK.R.
        real(dp), intent(in) :: x(:) !! Finite sample to robustly standardize.
        real(dp), intent(out) :: location !! Tau-estimator robust location returned as the ``m`` component.
        real(dp), intent(out) :: scale !! Consistency-corrected tau scale returned as the ``s`` component.
        real(dp), intent(in), optional :: c1 !! Tukey-biweight location tuning constant; default 4.5.
        real(dp), intent(in), optional :: c2 !! Truncated-square scale tuning constant; default 3.0.
        integer, intent(in), optional :: iterations !! Positive number of scale iterations; default 1.
        real(dp), intent(in), optional :: tolerance !! Positive relative scale stopping tolerance; default 1e-7.
        real(dp), allocatable :: absdev(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: standardized(:)
        real(dp) :: b
        real(dp) :: cc1
        real(dp) :: cc2
        real(dp) :: den
        real(dp) :: dnorm_b
        real(dp) :: erho
        real(dp) :: es2
        real(dp) :: mu0
        real(dp) :: pnorm_b
        real(dp) :: s0
        real(dp) :: snew
        real(dp) :: sumw
        real(dp) :: tol
        integer :: it
        integer :: maxit
        real(dp), parameter :: qnorm_075 = 0.6744897501960817_dp

        if (size(x) == 0) then
            location = 0.0_dp
            scale = 0.0_dp
            return
        end if
        if (any(.not. ieee_is_finite(x))) error stop 'scale_tau2_location_scale: x must be finite'
        cc1 = 4.5_dp
        if (present(c1)) cc1 = c1
        cc2 = 3.0_dp
        if (present(c2)) cc2 = c2
        maxit = 1
        if (present(iterations)) maxit = iterations
        tol = 1.0e-7_dp
        if (present(tolerance)) tol = tolerance
        if (maxit < 1) error stop 'scale_tau2_location_scale: iterations must be positive'
        if (tol <= 0.0_dp) error stop 'scale_tau2_location_scale: tolerance must be positive'
        if (cc2 <= 0.0_dp) error stop 'scale_tau2_location_scale: c2 must be positive'

        mu0 = median_real(x)
        allocate(absdev(size(x)), weights(size(x)), standardized(size(x)))
        absdev = abs(x - mu0)
        s0 = median_real(absdev)
        if (s0 <= 0.0_dp) then
            location = mu0
            scale = 0.0_dp
            return
        end if

        b = cc2 * qnorm_075
        pnorm_b = 0.5_dp * erfc(-b / sqrt(2.0_dp))
        dnorm_b = exp(-0.5_dp * b * b) / sqrt(2.0_dp * pi_dp)
        erho = 2.0_dp * ((1.0_dp - b * b) * pnorm_b - b * dnorm_b + b * b) - 1.0_dp
        es2 = erho
        den = real(size(x), dp) * es2

        location = mu0
        scale = s0
        do it = 1, maxit
            if (cc1 > 0.0_dp) then
                weights = max(0.0_dp, 1.0_dp - (absdev / (s0 * cc1)) ** 2) ** 2
                sumw = sum(weights)
                if (sumw > 0.0_dp) then
                    location = sum(x * weights) / sumw
                else
                    location = mu0
                end if
            else
                location = mu0
            end if
            standardized = (x - location) / s0
            standardized = min(standardized * standardized, cc2 * cc2)
            snew = s0 * sqrt(sum(standardized) / den)
            scale = snew
            if (.not. ieee_is_finite(scale)) return
            if (abs(scale - s0) <= tol * scale) exit
            s0 = scale
        end do
    end subroutine scale_tau2_location_scale

end module otrimle_linalg
