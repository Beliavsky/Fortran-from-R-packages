module otrimle_density
    use otrimle_kinds, only : dp
    use otrimle_types, only : otrimle_fit, kernel_density_result
    use otrimle_linalg, only : pi_dp, symmetric_eigen, sort_real_index
    implicit none
    private

    public :: kmeanfun
    public :: ksdfun
    public :: kerndensmeasure
    public :: kerndensp
    public :: kerndenscluster

contains

    pure elemental real(dp) function kmeanfun(n) result(out)
        real(dp), intent(in) :: n !! Effective sample size used by the empirical kernel-density calibration formula.

        if (n < 20.0_dp) then
            out = exp(-1.83084_dp) * n ** (-0.55323_dp)
        else
            out = exp(-2.4562342_dp) * n ** (-0.3524332_dp)
        end if
    end function kmeanfun

    pure elemental real(dp) function ksdfun(n) result(out)
        real(dp), intent(in) :: n !! Effective sample size used by the empirical kernel-density scale formula.

        if (n < 12.5_dp) then
            out = exp(-2.1482546_dp) * n ** (-0.4639833_dp)
        else
            out = exp(-2.6575388_dp) * n ** (-0.4547216_dp)
        end if
    end function ksdfun

    subroutine kerndensmeasure(x, result, weights, maxq, kernn)
        real(dp), intent(in) :: x(:) !! One-dimensional standardized scores whose density shape is assessed.
        type(kernel_density_result), intent(out) :: result !! Density grid, paired ordinate envelope, and discrepancy measure.
        real(dp), intent(in), optional :: weights(:) !! Nonnegative case weights; default equal weights.
        real(dp), intent(in), optional :: maxq !! Positive symmetric density-grid endpoint; default qnorm(0.9995).
        integer, intent(in), optional :: kernn !! Even number of density-grid points; default 100.
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: dens(:)
        integer, allocatable :: ord(:)
        real(dp), allocatable :: sorted_dens(:)
        real(dp) :: sw
        real(dp) :: mu
        real(dp) :: h
        real(dp) :: qmax
        real(dp) :: first_err
        real(dp) :: second_err
        integer :: i
        integer :: j
        integer :: m
        integer :: nk

        nk = 100
        if (present(kernn)) nk = kernn
        if (nk < 2 .or. mod(nk, 2) /= 0) error stop 'kerndensmeasure: kernn must be positive and even'
        qmax = 3.2905267314919255_dp
        if (present(maxq)) qmax = maxq
        if (size(x) < 2) error stop 'kerndensmeasure: need at least two observations'
        allocate(w(size(x)), z(size(x)), dens(nk), ord(nk), sorted_dens(nk))
        if (present(weights)) then
            if (size(weights) /= size(x)) error stop 'kerndensmeasure: weights size mismatch'
            w = max(weights, 0.0_dp)
        else
            w = 1.0_dp
        end if
        sw = sum(w)
        if (sw <= 0.0_dp) error stop 'kerndensmeasure: weights must have positive sum'
        w = w / sw
        mu = sum(w * x)
        z = x - mu
        h = bandwidth_nrd0(z)
        call density_r_gaussian(z, w, -qmax, qmax, nk, h, dens)
        call sort_real_index(dens, ord, ascending=.false.)
        do i = 1, nk
            sorted_dens(i) = dens(ord(i))
        end do
        m = nk / 2
        allocate(result%cp(m), result%cpx(nk))
        result%cpx = dens
        do j = 1, m
            result%cp(j) = 0.5_dp * (sorted_dens(2 * j - 1) + sorted_dens(2 * j))
        end do
        first_err = 0.0_dp
        second_err = 0.0_dp
        do j = 1, m
            first_err = first_err + (dens(m - j + 1) - result%cp(j)) ** 2
            second_err = second_err + (dens(m + j) - result%cp(j)) ** 2
        end do
        result%measure = sqrt((first_err + second_err) / real(nk, dp))
    end subroutine kerndensmeasure

    subroutine kerndensp(x, center, cov, measure, component_measures, standardized_measures, weights, maxq, kernn)
        real(dp), intent(in) :: x(:, :) !! Multivariate observations in rows; shape n by p.
        real(dp), intent(in) :: center(:) !! PCA center supplied by the fitted Gaussian component; length p.
        real(dp), intent(in) :: cov(:, :) !! PCA covariance supplied by the fitted Gaussian component; shape p by p.
        real(dp), intent(out) :: measure !! Aggregate positive standardized density-discrepancy measure.
        real(dp), allocatable, intent(out) :: component_measures(:) !! Raw one-dimensional discrepancy for each PC direction.
        real(dp), allocatable, intent(out) :: standardized_measures(:) !! Empirically standardized PC discrepancies.
        real(dp), intent(in), optional :: weights(:) !! Posterior component weights for observations; default equal weights.
        real(dp), intent(in), optional :: maxq !! Density-grid endpoint passed to kerndensmeasure().
        integer, intent(in), optional :: kernn !! Even density-grid size passed to kerndensmeasure().
        real(dp), allocatable :: eig(:)
        real(dp), allocatable :: vec(:, :)
        real(dp), allocatable :: scores(:, :)
        real(dp), allocatable :: w(:)
        real(dp), allocatable :: z(:)
        type(kernel_density_result) :: kd
        real(dp) :: nw
        integer :: i
        integer :: n
        integer :: p

        n = size(x, 1)
        p = size(x, 2)
        if (size(center) /= p) error stop 'kerndensp: center size mismatch'
        if (size(cov, 1) /= p .or. size(cov, 2) /= p) error stop 'kerndensp: covariance shape mismatch'
        allocate(eig(p), vec(p, p), scores(n, p), w(n), z(n))
        allocate(component_measures(p), standardized_measures(p))
        if (present(weights)) then
            if (size(weights) /= n) error stop 'kerndensp: weights size mismatch'
            w = max(weights, 0.0_dp)
        else
            w = 1.0_dp
        end if
        nw = sum(w)
        call symmetric_eigen(cov, eig, vec)
        scores = matmul(x - spread(center, 1, n), vec)
        do i = 1, p
            z = scores(:, i) / sqrt(max(eig(i), tiny(1.0_dp)))
            if (present(maxq) .and. present(kernn)) then
                call kerndensmeasure(z, kd, weights=w, maxq=maxq, kernn=kernn)
            else if (present(maxq)) then
                call kerndensmeasure(z, kd, weights=w, maxq=maxq)
            else if (present(kernn)) then
                call kerndensmeasure(z, kd, weights=w, kernn=kernn)
            else
                call kerndensmeasure(z, kd, weights=w)
            end if
            component_measures(i) = kd%measure
            standardized_measures(i) = (component_measures(i) - kmeanfun(nw)) / ksdfun(nw)
        end do
        measure = sum(max(standardized_measures, 0.0_dp) ** 2) / real(p, dp)
    end subroutine kerndensp

    subroutine kerndenscluster(x, fit, measure, ddpm, maxq, kernn)
        real(dp), intent(in) :: x(:, :) !! Original observations in rows; shape n by p.
        type(otrimle_fit), intent(in) :: fit !! Successful RIMLE/OTRIMLE fit whose Gaussian components are assessed.
        real(dp), intent(out) :: measure !! Mixture-weighted aggregate density diagnostic.
        real(dp), allocatable, intent(out) :: ddpm(:) !! Per-cluster multivariate density-discrepancy measures.
        real(dp), intent(in), optional :: maxq !! Density-grid endpoint used in every principal-component diagnostic.
        integer, intent(in), optional :: kernn !! Even density-grid size used in every principal-component diagnostic.
        real(dp), allocatable :: cm(:)
        real(dp), allocatable :: sm(:)
        integer :: j

        if (.not. allocated(fit%mean) .or. .not. allocated(fit%tau)) error stop 'kerndenscluster: fit is incomplete'
        allocate(ddpm(fit%g))
        do j = 1, fit%g
            if (present(maxq) .and. present(kernn)) then
                call kerndensp(x, fit%mean(:, j), fit%cov(:, :, j), ddpm(j), cm, sm, &
                    weights=fit%tau(:, j + 1), maxq=maxq, kernn=kernn)
            else if (present(maxq)) then
                call kerndensp(x, fit%mean(:, j), fit%cov(:, :, j), ddpm(j), cm, sm, &
                    weights=fit%tau(:, j + 1), maxq=maxq)
            else if (present(kernn)) then
                call kerndensp(x, fit%mean(:, j), fit%cov(:, :, j), ddpm(j), cm, sm, &
                    weights=fit%tau(:, j + 1), kernn=kernn)
            else
                call kerndensp(x, fit%mean(:, j), fit%cov(:, :, j), ddpm(j), cm, sm, weights=fit%tau(:, j + 1))
            end if
        end do
        measure = sqrt(sum(fit%exproportion(2:) * ddpm ** 2))
    end subroutine kerndenscluster

    subroutine density_r_gaussian(x, weights, from, to, n_user, bandwidth, density)
        ! Reproduce the Gaussian branch of R 4.1-era stats::density.default.
        ! R linearly bins observations on an internal grid, applies a zero-padded
        ! circular Gaussian convolution, then linearly interpolates to the user grid.
        real(dp), intent(in) :: x(:) !! Finite observations used by the density estimate.
        real(dp), intent(in) :: weights(:) !! Nonnegative normalized case weights with the same length as x.
        real(dp), intent(in) :: from !! Lower endpoint of the requested density grid.
        real(dp), intent(in) :: to !! Upper endpoint of the requested density grid.
        integer, intent(in) :: n_user !! Number of requested output grid points; at least two.
        real(dp), intent(in) :: bandwidth !! Positive Gaussian-kernel standard deviation.
        real(dp), intent(out) :: density(:) !! Interpolated density ordinates; length n_user.
        real(dp), allocatable :: binned(:)
        real(dp), allocatable :: kernel(:)
        real(dp), allocatable :: conv(:)
        real(dp) :: delta
        real(dp) :: frac
        real(dp) :: kcoord
        real(dp) :: lo
        real(dp) :: pos
        real(dp) :: range
        real(dp) :: up
        real(dp) :: value
        real(dp) :: xdelta
        integer :: i
        integer :: ix
        integer :: j
        integer :: kidx
        integer :: n_grid
        integer :: n_total
        integer :: src

        if (size(weights) /= size(x)) error stop 'density_r_gaussian: weights size mismatch'
        if (size(density) /= n_user) error stop 'density_r_gaussian: density size mismatch'
        if (n_user < 2) error stop 'density_r_gaussian: n_user must be at least two'
        if (bandwidth <= 0.0_dp) error stop 'density_r_gaussian: bandwidth must be positive'
        n_grid = max(n_user, 512)
        if (n_grid > 512) n_grid = next_power_of_two(n_grid)
        n_total = 2 * n_grid
        lo = from - 4.0_dp * bandwidth
        up = to + 4.0_dp * bandwidth
        range = up - lo
        xdelta = range / real(n_grid - 1, dp)
        allocate(binned(n_total), kernel(n_total), conv(n_grid))
        binned = 0.0_dp

        do i = 1, size(x)
            pos = (x(i) - lo) / xdelta
            ix = int(floor(pos))
            frac = pos - real(ix, dp)
            if (ix >= 0 .and. ix <= n_grid - 2) then
                binned(ix + 1) = binned(ix + 1) + (1.0_dp - frac) * weights(i)
                binned(ix + 2) = binned(ix + 2) + frac * weights(i)
            else if (ix == -1) then
                binned(1) = binned(1) + frac * weights(i)
            else if (ix == n_grid - 1) then
                binned(ix + 1) = binned(ix + 1) + (1.0_dp - frac) * weights(i)
            end if
        end do

        delta = 2.0_dp * range / real(n_total - 1, dp)
        do i = 1, n_total
            if (i <= n_grid + 1) then
                kcoord = real(i - 1, dp) * delta
            else
                src = n_total + 2 - i
                kcoord = -real(src - 1, dp) * delta
            end if
            kernel(i) = exp(-0.5_dp * (kcoord / bandwidth) ** 2) / &
                (bandwidth * sqrt(2.0_dp * pi_dp))
        end do

        do i = 1, n_grid
            value = 0.0_dp
            do j = 1, n_total
                kidx = 1 + modulo(j - i, n_total)
                value = value + binned(j) * kernel(kidx)
            end do
            conv(i) = max(0.0_dp, value)
        end do

        do i = 1, n_user
            value = from + (to - from) * real(i - 1, dp) / real(n_user - 1, dp)
            pos = (value - lo) / range * real(n_grid - 1, dp)
            ix = int(floor(pos))
            frac = pos - real(ix, dp)
            if (ix <= 0) then
                density(i) = (1.0_dp - frac) * conv(1) + frac * conv(2)
            else if (ix >= n_grid - 1) then
                density(i) = conv(n_grid)
            else
                density(i) = (1.0_dp - frac) * conv(ix + 1) + frac * conv(ix + 2)
            end if
        end do
    end subroutine density_r_gaussian

    pure integer function next_power_of_two(n) result(out)
        integer, intent(in) :: n !! Positive integer for which the next power of two is required.

        out = 1
        do while (out < n)
            out = 2 * out
        end do
    end function next_power_of_two

    real(dp) function bandwidth_nrd0(x) result(h)
        real(dp), intent(in) :: x(:) !! Sample used by R's bw.nrd0-style normal-reference bandwidth rule.
        real(dp), allocatable :: y(:)
        real(dp) :: mu
        real(dp) :: sd
        real(dp) :: iqr
        real(dp) :: scale
        integer :: n

        n = size(x)
        if (n <= 1) error stop 'bandwidth_nrd0: need at least two observations'
        allocate(y(n))
        y = x
        call sort_values(y)
        mu = sum(x) / real(n, dp)
        sd = sqrt(sum((x - mu) ** 2) / real(n - 1, dp))
        iqr = quantile_type7_sorted(y, 0.75_dp) - quantile_type7_sorted(y, 0.25_dp)
        scale = min(sd, iqr / 1.34_dp)
        if (scale <= 0.0_dp) scale = sd
        if (scale <= 0.0_dp) scale = abs(x(1))
        if (scale <= 0.0_dp) scale = 1.0_dp
        h = 0.9_dp * scale * real(n, dp) ** (-0.2_dp)
    end function bandwidth_nrd0

    real(dp) function quantile_type7_sorted(x, prob) result(q)
        real(dp), intent(in) :: x(:) !! Ascending sorted sample.
        real(dp), intent(in) :: prob !! Probability in [0,1] for R type-7 interpolation.
        real(dp) :: h
        real(dp) :: frac
        integer :: j
        integer :: n

        n = size(x)
        if (n == 1) then
            q = x(1)
            return
        end if
        h = 1.0_dp + real(n - 1, dp) * min(1.0_dp, max(0.0_dp, prob))
        j = int(floor(h))
        frac = h - real(j, dp)
        if (j >= n) then
            q = x(n)
        else
            q = (1.0_dp - frac) * x(j) + frac * x(j + 1)
        end if
    end function quantile_type7_sorted

    subroutine sort_values(x)
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
    end subroutine sort_values

end module otrimle_density
