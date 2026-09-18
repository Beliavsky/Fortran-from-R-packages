module tclust_linalg
    use tclust_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: symmetric_eigen
    public :: covariance_matrix
    public :: weighted_mean_cov
    public :: log_mvnorm
    public :: mahalanobis_sq
    public :: determinant_sym
    public :: qchisq
    public :: mvn_sample
    public :: sort_real_with_index
    public :: sort_integer_unique

contains

    subroutine symmetric_eigen(a, values, vectors)
        real(dp), intent(in) :: a(:, :) !! Real symmetric matrix to decompose; shape p by p.
        real(dp), intent(out) :: values(:) !! Eigenvalues in ascending order; length p.
        real(dp), intent(out) :: vectors(:, :) !! Corresponding orthonormal eigenvectors in columns; shape p by p.
        real(dp), allocatable :: b(:, :)
        real(dp) :: app
        real(dp) :: aqq
        real(dp) :: apq
        real(dp) :: c
        real(dp) :: s
        real(dp) :: tau
        real(dp) :: t
        real(dp) :: temp
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: p
        integer :: q
        integer :: n
        integer :: sweep
        integer :: max_sweeps

        n = size(a, 1)
        if (size(a, 2) /= n) error stop 'symmetric_eigen: matrix must be square'
        if (size(values) /= n) error stop 'symmetric_eigen: bad values size'
        if (size(vectors, 1) /= n .or. size(vectors, 2) /= n) error stop 'symmetric_eigen: bad vectors size'

        allocate(b(n, n))
        b = 0.5_dp * (a + transpose(a))
        vectors = 0.0_dp
        do i = 1, n
            vectors(i, i) = 1.0_dp
        end do

        tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(b)))
        max_sweeps = max(20, 50 * n * n)
        do sweep = 1, max_sweeps
            apq = 0.0_dp
            p = 1
            q = min(2, n)
            do j = 2, n
                do i = 1, j - 1
                    if (abs(b(i, j)) > abs(apq)) then
                        apq = b(i, j)
                        p = i
                        q = j
                    end if
                end do
            end do
            if (abs(apq) <= tol) exit

            app = b(p, p)
            aqq = b(q, q)
            tau = (aqq - app) / (2.0_dp * apq)
            if (tau >= 0.0_dp) then
                t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
            else
                t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
            end if
            c = 1.0_dp / sqrt(1.0_dp + t * t)
            s = t * c

            do j = 1, n
                if (j /= p .and. j /= q) then
                    temp = b(j, p)
                    b(j, p) = c * temp - s * b(j, q)
                    b(p, j) = b(j, p)
                    b(j, q) = s * temp + c * b(j, q)
                    b(q, j) = b(j, q)
                end if
            end do
            b(p, p) = c * c * app - 2.0_dp * c * s * apq + s * s * aqq
            b(q, q) = s * s * app + 2.0_dp * c * s * apq + c * c * aqq
            b(p, q) = 0.0_dp
            b(q, p) = 0.0_dp

            do j = 1, n
                temp = vectors(j, p)
                vectors(j, p) = c * temp - s * vectors(j, q)
                vectors(j, q) = s * temp + c * vectors(j, q)
            end do
        end do

        do i = 1, n
            values(i) = b(i, i)
        end do
        call sort_eigensystem(values, vectors)
    end subroutine symmetric_eigen

    subroutine sort_eigensystem(values, vectors)
        real(dp), intent(inout) :: values(:) !! Eigenvalues to sort in ascending order in place.
        real(dp), intent(inout) :: vectors(:, :) !! Eigenvector columns permuted consistently with values.
        real(dp) :: tv
        real(dp), allocatable :: col(:)
        integer :: i
        integer :: j
        integer :: m
        integer :: n

        n = size(values)
        allocate(col(size(vectors, 1)))
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
    end subroutine sort_eigensystem

    subroutine covariance_matrix(x, cov, center, unbiased)
        real(dp), intent(in) :: x(:, :) !! Observations in rows and variables in columns; shape n by p.
        real(dp), intent(out) :: cov(:, :) !! Covariance or scatter matrix; shape p by p.
        real(dp), intent(out), optional :: center(:) !! Optional arithmetic mean vector; length p.
        logical, intent(in), optional :: unbiased !! If true, divide by n-1; otherwise divide by n.
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: z(:, :)
        real(dp) :: den
        integer :: n
        integer :: p
        logical :: unb

        n = size(x, 1)
        p = size(x, 2)
        if (size(cov, 1) /= p .or. size(cov, 2) /= p) error stop 'covariance_matrix: bad cov shape'
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

    subroutine weighted_mean_cov(x, w, mean, cov, weight_sum)
        real(dp), intent(in) :: x(:, :) !! Observations in rows and variables in columns; shape n by p.
        real(dp), intent(in) :: w(:) !! Nonnegative observation weights; length n.
        real(dp), intent(out) :: mean(:) !! Weighted mean vector; length p.
        real(dp), intent(out) :: cov(:, :) !! Maximum-likelihood weighted covariance; shape p by p.
        real(dp), intent(out), optional :: weight_sum !! Optional sum of observation weights.
        real(dp), allocatable :: z(:, :)
        real(dp) :: sw
        integer :: i
        integer :: n
        integer :: p

        n = size(x, 1)
        p = size(x, 2)
        if (size(w) /= n) error stop 'weighted_mean_cov: bad weights size'
        sw = sum(w)
        if (sw <= tiny(1.0_dp)) then
            mean = 0.0_dp
            cov = 0.0_dp
            if (present(weight_sum)) weight_sum = sw
            return
        end if
        mean = matmul(transpose(x), w) / sw
        allocate(z(n, p))
        z = x - spread(mean, 1, n)
        cov = 0.0_dp
        do i = 1, n
            cov = cov + w(i) * outer_product(z(i, :), z(i, :))
        end do
        cov = cov / sw
        if (present(weight_sum)) weight_sum = sw
    end subroutine weighted_mean_cov

    pure function outer_product(a, b) result(c)
        real(dp), intent(in) :: a(:) !! Left vector in an outer product.
        real(dp), intent(in) :: b(:) !! Right vector in an outer product.
        real(dp) :: c(size(a), size(b))

        c = spread(a, 2, size(b)) * spread(b, 1, size(a))
    end function outer_product

    subroutine log_mvnorm(x, mean, cov, logdens, zero_tol)
        real(dp), intent(in) :: x(:, :) !! Evaluation points in rows; shape n by p.
        real(dp), intent(in) :: mean(:) !! Gaussian mean vector; length p.
        real(dp), intent(in) :: cov(:, :) !! Symmetric covariance matrix; shape p by p.
        real(dp), intent(out) :: logdens(:) !! Log densities for each observation; length n.
        real(dp), intent(in), optional :: zero_tol !! Eigenvalue floor for nearly singular covariance matrices.
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: vec(:, :)
        real(dp), allocatable :: z(:)
        real(dp) :: floor_ev
        real(dp) :: logdet
        real(dp) :: quad
        integer :: i
        integer :: p

        p = size(mean)
        if (size(x, 2) /= p) error stop 'log_mvnorm: mean dimension mismatch'
        allocate(eig(p), vec(p, p), z(p))
        call symmetric_eigen(cov, eig, vec)
        floor_ev = max(100.0_dp * epsilon(1.0_dp), 1.0e-14_dp)
        if (present(zero_tol)) floor_ev = max(floor_ev, zero_tol)
        eig = max(eig, floor_ev)
        logdet = sum(log(eig))
        do i = 1, size(x, 1)
            z = matmul(transpose(vec), x(i, :) - mean)
            quad = sum((z * z) / eig)
            logdens(i) = -0.5_dp * (real(p, dp) * log(2.0_dp * pi_dp) + logdet + quad)
        end do
    end subroutine log_mvnorm

    real(dp) function mahalanobis_sq(x, mean, cov, zero_tol) result(d2)
        real(dp), intent(in) :: x(:) !! Observation vector whose squared Mahalanobis distance is requested.
        real(dp), intent(in) :: mean(:) !! Center vector with the same length as x.
        real(dp), intent(in) :: cov(:, :) !! Symmetric covariance matrix.
        real(dp), intent(in), optional :: zero_tol !! Eigenvalue floor used when covariance is nearly singular.
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: vec(:, :)
        real(dp), allocatable :: z(:)
        real(dp) :: floor_ev
        integer :: p

        p = size(x)
        allocate(eig(p), vec(p, p), z(p))
        call symmetric_eigen(cov, eig, vec)
        floor_ev = 100.0_dp * epsilon(1.0_dp)
        if (present(zero_tol)) floor_ev = max(floor_ev, zero_tol)
        eig = max(eig, floor_ev)
        z = matmul(transpose(vec), x - mean)
        d2 = sum((z * z) / eig)
    end function mahalanobis_sq

    real(dp) function determinant_sym(a) result(det)
        real(dp), intent(in) :: a(:, :) !! Real symmetric matrix whose determinant is requested.
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: vec(:, :)
        integer :: n

        n = size(a, 1)
        allocate(eig(n), vec(n, n))
        call symmetric_eigen(a, eig, vec)
        det = product(eig)
    end function determinant_sym

    subroutine mvn_sample(mean, cov, x)
        use tclust_rng, only : rand_normal
        real(dp), intent(in) :: mean(:) !! Mean vector of the multivariate normal distribution.
        real(dp), intent(in) :: cov(:, :) !! Symmetric positive-semidefinite covariance matrix.
        real(dp), intent(out) :: x(:) !! One generated multivariate-normal observation.
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

    subroutine sort_real_with_index(x, order, ascending)
        real(dp), intent(in) :: x(:) !! Values to rank without modifying the input vector.
        integer, intent(out) :: order(:) !! Permutation of 1:size(x) giving sorted positions.
        logical, intent(in), optional :: ascending !! True for ascending order; false for descending order.
        integer :: i
        integer :: j
        integer :: m
        integer :: t
        logical :: asc

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
    end subroutine sort_real_with_index

    subroutine sort_integer_unique(x, values, nvalues)
        integer, intent(in) :: x(:) !! Integer labels whose distinct sorted values are requested.
        integer, intent(out) :: values(:) !! Workspace receiving distinct labels in ascending order.
        integer, intent(out) :: nvalues !! Number of distinct labels written to values.
        integer, allocatable :: temp(:)
        integer :: i
        integer :: j
        integer :: m
        integer :: t

        allocate(temp(size(x)))
        temp = x
        do i = 1, size(temp) - 1
            m = i
            do j = i + 1, size(temp)
                if (temp(j) < temp(m)) m = j
            end do
            if (m /= i) then
                t = temp(i)
                temp(i) = temp(m)
                temp(m) = t
            end if
        end do
        nvalues = 0
        do i = 1, size(temp)
            if (i == 1 .or. temp(i) /= temp(i - 1)) then
                nvalues = nvalues + 1
                values(nvalues) = temp(i)
            end if
        end do
    end subroutine sort_integer_unique

    real(dp) function qchisq(prob, df) result(x)
        real(dp), intent(in) :: prob !! Lower-tail probability strictly between zero and one.
        real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
        real(dp) :: lo
        real(dp) :: hi
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
        do while (regularized_gamma_p(0.5_dp * df, 0.5_dp * hi) < prob)
            hi = 2.0_dp * hi
            if (hi > huge(1.0_dp) / 4.0_dp) exit
        end do
        do iter = 1, 120
            mid = 0.5_dp * (lo + hi)
            if (regularized_gamma_p(0.5_dp * df, 0.5_dp * mid) < prob) then
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
        real(dp) :: ap
        real(dp) :: del
        real(dp) :: sum_series
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        real(dp) :: h
        real(dp) :: an
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

end module tclust_linalg
