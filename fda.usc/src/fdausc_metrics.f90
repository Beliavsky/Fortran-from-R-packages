module fdausc_metrics
    use r_kinds, only : dp
    use r_distributions, only : r_pt
    use r_linalg, only : inverse_matrix
    use fdausc_integration, only : integrate_curve
    implicit none
    private

    public :: lp_distance_matrix, ordinary_distance_matrix, hausdorff_distance_matrix
    public :: kl_distance_matrix, dtw_distance_matrix, wdtw_distance_matrix, twed_distance_matrix
    public :: distance_correlation, bias_corrected_distance_correlation, distance_correlation_test
    public :: combine_distance_matrices

contains

    pure subroutine lp_distance_matrix(x, curves1, curves2, distances, p, weights, method, dscale)
        real(dp), intent(in) :: x(:) !! Increasing functional argument grid shared by both curve matrices.
        real(dp), intent(in) :: curves1(:, :) !! First functional sample, one curve per row.
        real(dp), intent(in) :: curves2(:, :) !! Second functional sample, one curve per row.
        real(dp), intent(out) :: distances(:, :) !! Pairwise Lp distances with shape [size(curves1,1),size(curves2,1)].
        real(dp), intent(in), optional :: p !! Lp exponent; zero selects the supremum norm, default two.
        real(dp), intent(in), optional :: weights(:) !! Optional pointwise integration weights on x.
        integer, intent(in), optional :: method !! Integration selector: 1 trapezoid, 2 composite Simpson, 3 extended Simpson.
        real(dp), intent(in), optional :: dscale !! Positive divisor applied to every computed distance, default one.
        integer :: i
        integer :: j
        integer :: meth
        real(dp) :: pp
        real(dp) :: scale
        real(dp), allocatable :: f(:)

        pp = 2.0_dp
        if (present(p)) pp = p
        meth = 2
        if (present(method)) meth = method
        scale = 1.0_dp
        if (present(dscale)) scale = dscale
        if (abs(scale) <= tiny(1.0_dp)) scale = 1.0_dp
        allocate(f(size(x)))
        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                if (abs(pp) <= tiny(1.0_dp)) then
                    distances(i, j) = maxval(abs(curves1(i, :) - curves2(j, :)))/scale
                else
                    f = abs(curves1(i, :) - curves2(j, :))**pp
                    if (present(weights)) f = f*weights
                    distances(i, j) = max(0.0_dp, integrate_curve(x, f, meth))**(1.0_dp/pp)/scale
                end if
            end do
        end do
    end subroutine lp_distance_matrix

    subroutine ordinary_distance_matrix(x, y, distances, method_code, p, dscale, info)
        real(dp), intent(in) :: x(:, :) !! First multivariate sample, observations in rows and variables in columns.
        real(dp), intent(in) :: y(:, :) !! Second multivariate sample with the same number of variables as x.
        real(dp), intent(out) :: distances(:, :) !! Pairwise row distances with shape [size(x,1),size(y,1)].
        integer, intent(in), optional :: method_code !! Distance-method code 1 through 7; default Euclidean.
        real(dp), intent(in), optional :: p !! Minkowski exponent for method 6, default two.
        real(dp), intent(in), optional :: dscale !! Nonzero divisor applied to distances, default one.
        integer, intent(out), optional :: info !! Zero on success; positive if covariance inversion for Mahalanobis fails.
        integer :: i
        integer :: j
        integer :: k
        integer :: code
        integer :: stat
        integer :: nall
        real(dp) :: pp
        real(dp) :: scale
        real(dp) :: denom
        real(dp) :: val
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: z(:, :)
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: invcov(:, :)
        real(dp), allocatable :: meanv(:)

        code = 1
        if (present(method_code)) code = method_code
        pp = 2.0_dp
        if (present(p)) pp = p
        scale = 1.0_dp
        if (present(dscale)) scale = dscale
        if (abs(scale) <= tiny(1.0_dp)) scale = 1.0_dp
        if (present(info)) info = 0
        allocate(d(size(x, 2)))
        if (code == 7) then
            nall = size(x, 1) + size(y, 1)
            allocate(z(nall, size(x, 2)), meanv(size(x, 2)))
            z(:size(x, 1), :) = x
            z(size(x, 1) + 1:, :) = y
            meanv = sum(z, dim=1)/real(nall, dp)
            allocate(cov(size(x, 2), size(x, 2)))
            cov = 0.0_dp
            do i = 1, nall
                d = z(i, :) - meanv
                do j = 1, size(d)
                    do k = 1, size(d)
                        cov(j, k) = cov(j, k) + d(j)*d(k)
                    end do
                end do
            end do
            cov = cov/real(max(1, nall - 1), dp)
            call inverse_matrix(cov, invcov, stat)
            if (stat /= 0) then
                distances = huge(1.0_dp)
                if (present(info)) info = stat
                return
            end if
        end if

        do i = 1, size(x, 1)
            do j = 1, size(y, 1)
                d = abs(x(i, :) - y(j, :))
                select case (code)
                case (2)
                    val = maxval(d)
                case (3)
                    val = sum(d)
                case (4)
                    val = 0.0_dp
                    do k = 1, size(d)
                        denom = abs(x(i, k)) + abs(y(j, k))
                        if (denom > tiny(1.0_dp)) val = val + d(k)/denom
                    end do
                case (5)
                    denom = 0.0_dp
                    val = 0.0_dp
                    do k = 1, size(d)
                        if (abs(x(i, k)) > tiny(1.0_dp) .or. abs(y(j, k)) > tiny(1.0_dp)) then
                            denom = denom + 1.0_dp
                            if ((abs(x(i, k)) > tiny(1.0_dp)) .neqv. &
                                (abs(y(j, k)) > tiny(1.0_dp))) val = val + 1.0_dp
                        end if
                    end do
                    if (denom > 0.0_dp) val = val/denom
                case (6)
                    val = sum(d**pp)**(1.0_dp/pp)
                case (7)
                    d = x(i, :) - y(j, :)
                    val = sqrt(max(0.0_dp, dot_product(d, matmul(invcov, d))))
                case default
                    val = sqrt(sum(d*d))
                end select
                distances(i, j) = val/scale
            end do
        end do
    end subroutine ordinary_distance_matrix

    pure subroutine hausdorff_distance_matrix(t, curves1, curves2, distances)
        real(dp), intent(in) :: t(:) !! Common curve argument grid used as the first graph coordinate.
        real(dp), intent(in) :: curves1(:, :) !! First functional sample, curves in rows.
        real(dp), intent(in) :: curves2(:, :) !! Second functional sample, curves in rows.
        real(dp), intent(out) :: distances(:, :) !! Pairwise Hausdorff distances between discretized curve graphs.
        integer :: i
        integer :: j

        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                distances(i, j) = curve_hausdorff(t, curves1(i, :), curves2(j, :))
            end do
        end do
    end subroutine hausdorff_distance_matrix

    pure subroutine kl_distance_matrix(x, curves1, curves2, distances, symmetric, eps, method)
        real(dp), intent(in) :: x(:) !! Increasing grid on which density-like curves are observed.
        real(dp), intent(in) :: curves1(:, :) !! First nonnegative functional sample, curves in rows.
        real(dp), intent(in) :: curves2(:, :) !! Second nonnegative functional sample, curves in rows.
        real(dp), intent(out) :: distances(:, :) !! Pairwise Kullback-Leibler divergences or symmetrized divergences.
        logical, intent(in), optional :: symmetric !! If true, average KL(f,g) and KL(g,f); default true.
        real(dp), intent(in), optional :: eps !! Positive lower clipping level before logarithms; default sqrt(epsilon).
        integer, intent(in), optional :: method !! Integration selector passed to integrate_curve.
        logical :: sym
        integer :: meth
        integer :: i
        integer :: j
        real(dp) :: ee
        real(dp) :: z1
        real(dp) :: z2
        real(dp) :: kl12
        real(dp) :: kl21
        real(dp), allocatable :: f(:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: q(:)

        sym = .true.
        if (present(symmetric)) sym = symmetric
        ee = sqrt(epsilon(1.0_dp))
        if (present(eps)) ee = max(eps, tiny(1.0_dp))
        meth = 2
        if (present(method)) meth = method
        allocate(f(size(x)), g(size(x)), q(size(x)))
        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                f = max(curves1(i, :), ee)
                g = max(curves2(j, :), ee)
                z1 = integrate_curve(x, f, meth)
                z2 = integrate_curve(x, g, meth)
                if (z1 > 0.0_dp) f = f/z1
                if (z2 > 0.0_dp) g = g/z2
                q = f*log(f/g)
                kl12 = integrate_curve(x, q, meth)
                if (sym) then
                    q = g*log(g/f)
                    kl21 = integrate_curve(x, q, meth)
                    distances(i, j) = 0.5_dp*(kl12 + kl21)
                else
                    distances(i, j) = kl12
                end if
            end do
        end do
    end subroutine kl_distance_matrix

    pure subroutine dtw_distance_matrix(curves1, curves2, distances, p, window)
        real(dp), intent(in) :: curves1(:, :) !! First time-series sample, one series per row.
        real(dp), intent(in) :: curves2(:, :) !! Second time-series sample, one series per row.
        real(dp), intent(out) :: distances(:, :) !! Pairwise dynamic-time-warping terminal costs.
        real(dp), intent(in), optional :: p !! Power applied to local absolute differences, default two.
        integer, intent(in), optional :: window !! Sakoe-Chiba half-window in index units; default max series length.
        integer :: i
        integer :: j
        integer :: w
        real(dp) :: pp

        pp = 2.0_dp
        if (present(p)) pp = p
        w = max(size(curves1, 2), size(curves2, 2))
        if (present(window)) w = window
        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                distances(i, j) = dtw_pair(curves1(i, :), curves2(j, :), pp, w)
            end do
        end do
    end subroutine dtw_distance_matrix

    pure subroutine wdtw_distance_matrix(curves1, curves2, distances, p, window, wmax, g)
        real(dp), intent(in) :: curves1(:, :) !! First time-series sample, one series per row.
        real(dp), intent(in) :: curves2(:, :) !! Second time-series sample, one series per row.
        real(dp), intent(out) :: distances(:, :) !! Pairwise weighted dynamic-time-warping terminal costs.
        real(dp), intent(in), optional :: p !! Power applied to weighted local differences, default two.
        integer, intent(in), optional :: window !! Sakoe-Chiba half-window in index units; default max series length.
        real(dp), intent(in), optional :: wmax !! Maximum logistic weight, default one.
        real(dp), intent(in), optional :: g !! Logistic steepness parameter, default 0.05.
        integer :: i
        integer :: j
        integer :: w
        real(dp) :: pp
        real(dp) :: wm
        real(dp) :: gg

        pp = 2.0_dp
        if (present(p)) pp = p
        w = max(size(curves1, 2), size(curves2, 2))
        if (present(window)) w = window
        wm = 1.0_dp
        if (present(wmax)) wm = wmax
        gg = 0.05_dp
        if (present(g)) gg = g
        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                distances(i, j) = wdtw_pair(curves1(i, :), curves2(j, :), pp, w, wm, gg)
            end do
        end do
    end subroutine wdtw_distance_matrix

    pure subroutine twed_distance_matrix(curves1, curves2, distances, p, lambda, nu)
        real(dp), intent(in) :: curves1(:, :) !! First equally spaced time-series sample, one series per row.
        real(dp), intent(in) :: curves2(:, :) !! Second equally spaced time-series sample, one series per row.
        real(dp), intent(out) :: distances(:, :) !! Pairwise time-warp-edit distances.
        real(dp), intent(in), optional :: p !! Power used in value differences, default two.
        real(dp), intent(in), optional :: lambda !! Edit penalty, default one.
        real(dp), intent(in), optional :: nu !! Stiffness penalty per index displacement, default 0.05.
        integer :: i
        integer :: j
        real(dp) :: pp
        real(dp) :: lam
        real(dp) :: stiffness

        pp = 2.0_dp
        if (present(p)) pp = p
        lam = 1.0_dp
        if (present(lambda)) lam = lambda
        stiffness = 0.05_dp
        if (present(nu)) stiffness = nu
        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                distances(i, j) = twed_pair(curves1(i, :), curves2(j, :), pp, lam, stiffness)
            end do
        end do
    end subroutine twed_distance_matrix

    pure real(dp) function distance_correlation(d1, d2) result(r)
        real(dp), intent(in) :: d1(:, :) !! First square distance matrix.
        real(dp), intent(in) :: d2(:, :) !! Second square distance matrix of the same order as d1.
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:, :)
        real(dp) :: xy
        real(dp) :: xx
        real(dp) :: yy

        call double_center(d1, a)
        call double_center(d2, b)
        xy = sqrt(max(0.0_dp, sum(a*b)))/real(size(d1, 1), dp)
        xx = sqrt(max(0.0_dp, sum(a*a)))/real(size(d1, 1), dp)
        yy = sqrt(max(0.0_dp, sum(b*b)))/real(size(d1, 1), dp)
        if (xx*yy <= tiny(1.0_dp)) then
            r = 0.0_dp
        else
            r = xy/sqrt(xx*yy)
        end if
    end function distance_correlation

    pure real(dp) function bias_corrected_distance_correlation(d1, d2) result(r)
        real(dp), intent(in) :: d1(:, :) !! First square distance matrix with sample size at least four.
        real(dp), intent(in) :: d2(:, :) !! Second square distance matrix of the same order as d1.
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: b(:, :)
        integer :: n
        integer :: i
        real(dp) :: xy
        real(dp) :: xx
        real(dp) :: yy
        real(dp) :: diagxy
        real(dp) :: diagxx
        real(dp) :: diagyy

        n = size(d1, 1)
        if (n <= 3) then
            r = 0.0_dp
            return
        end if
        call astar2(d1, a)
        call astar2(d2, b)
        diagxy = 0.0_dp
        diagxx = 0.0_dp
        diagyy = 0.0_dp
        do i = 1, n
            diagxy = diagxy + a(i, i)*b(i, i)
            diagxx = diagxx + a(i, i)*a(i, i)
            diagyy = diagyy + b(i, i)*b(i, i)
        end do
        xy = sum(a*b) - real(n, dp)/real(n - 2, dp)*diagxy
        xx = sum(a*a) - real(n, dp)/real(n - 2, dp)*diagxx
        yy = sum(b*b) - real(n, dp)/real(n - 2, dp)*diagyy
        if (xx <= 0.0_dp .or. yy <= 0.0_dp) then
            r = 0.0_dp
        else
            r = xy/sqrt(xx*yy)
        end if
    end function bias_corrected_distance_correlation

    pure subroutine distance_correlation_test(d1, d2, r, t_stat, df, p_value)
        real(dp), intent(in) :: d1(:, :) !! First square distance matrix with at least four observations.
        real(dp), intent(in) :: d2(:, :) !! Second square distance matrix of the same order as d1.
        real(dp), intent(out) :: r !! Bias-corrected distance correlation estimate.
        real(dp), intent(out) :: t_stat !! Szekely-Rizzo transformed distance-correlation statistic.
        real(dp), intent(out) :: df !! Approximate Student-t degrees of freedom n(n-3)/2-1.
        real(dp), intent(out) :: p_value !! Upper-tail approximate Student-t p-value.
        real(dp) :: m
        real(dp) :: rc

        r = bias_corrected_distance_correlation(d1, d2)
        m = real(size(d1, 1)*(size(d1, 1) - 3), dp)/2.0_dp
        df = m - 1.0_dp
        rc = max(-1.0_dp + 10.0_dp*epsilon(1.0_dp), min(1.0_dp - 10.0_dp*epsilon(1.0_dp), r))
        t_stat = sqrt(max(0.0_dp, m - 1.0_dp))*rc/sqrt(max(tiny(1.0_dp), 1.0_dp - rc*rc))
        p_value = 1.0_dp - r_pt(t_stat, df)
    end subroutine distance_correlation_test

    pure real(dp) function curve_hausdorff(t, a, b) result(h)
        real(dp), intent(in) :: t(:) !! Common curve argument coordinates.
        real(dp), intent(in) :: a(:) !! First curve ordinates on t.
        real(dp), intent(in) :: b(:) !! Second curve ordinates on t.
        real(dp) :: hab
        real(dp) :: hba
        real(dp) :: d
        real(dp) :: best
        integer :: i
        integer :: j

        hab = 0.0_dp
        do i = 1, size(t)
            best = huge(1.0_dp)
            do j = 1, size(t)
                d = sqrt((t(i) - t(j))**2 + (a(i) - b(j))**2)
                best = min(best, d)
            end do
            hab = max(hab, best)
        end do
        hba = 0.0_dp
        do i = 1, size(t)
            best = huge(1.0_dp)
            do j = 1, size(t)
                d = sqrt((t(i) - t(j))**2 + (b(i) - a(j))**2)
                best = min(best, d)
            end do
            hba = max(hba, best)
        end do
        h = max(hab, hba)
    end function curve_hausdorff

    pure real(dp) function dtw_pair(a, b, p, window) result(distance)
        real(dp), intent(in) :: a(:) !! First time series.
        real(dp), intent(in) :: b(:) !! Second time series.
        real(dp), intent(in) :: p !! Power applied to local absolute differences.
        integer, intent(in) :: window !! Sakoe-Chiba half-window in index units.
        real(dp), allocatable :: d(:, :)
        integer :: i
        integer :: j
        integer :: j1
        integer :: j2

        allocate(d(0:size(a), 0:size(b)))
        d = huge(1.0_dp)
        d(0, 0) = 0.0_dp
        do i = 1, size(a)
            j1 = max(1, i - window)
            j2 = min(size(b), i + window)
            do j = j1, j2
                d(i, j) = abs(a(i) - b(j))**p + min(d(i - 1, j), d(i, j - 1), d(i - 1, j - 1))
            end do
        end do
        distance = d(size(a), size(b))
    end function dtw_pair

    pure real(dp) function wdtw_pair(a, b, p, window, wmax, g) result(distance)
        real(dp), intent(in) :: a(:) !! First time series.
        real(dp), intent(in) :: b(:) !! Second time series.
        real(dp), intent(in) :: p !! Power applied to weighted local absolute differences.
        integer, intent(in) :: window !! Sakoe-Chiba half-window in index units.
        real(dp), intent(in) :: wmax !! Maximum logistic warping weight.
        real(dp), intent(in) :: g !! Logistic steepness parameter.
        real(dp), allocatable :: d(:, :)
        integer :: i
        integer :: j
        integer :: j1
        integer :: j2
        integer :: nref
        real(dp) :: weight

        nref = max(size(a), size(b))
        allocate(d(0:size(a), 0:size(b)))
        d = huge(1.0_dp)
        d(0, 0) = 0.0_dp
        do i = 1, size(a)
            j1 = max(1, i - window)
            j2 = min(size(b), i + window)
            do j = j1, j2
                weight = wmax/(1.0_dp + exp(-g*(real(abs(i - j), dp) - 0.5_dp*real(nref, dp))))
                d(i, j) = abs(weight*(a(i) - b(j)))**p + min(d(i - 1, j), d(i, j - 1), d(i - 1, j - 1))
            end do
        end do
        distance = d(size(a), size(b))
    end function wdtw_pair

    pure real(dp) function twed_pair(a, b, p, lambda, nu) result(distance)
        real(dp), intent(in) :: a(:) !! First equally spaced time series.
        real(dp), intent(in) :: b(:) !! Second equally spaced time series.
        real(dp), intent(in) :: p !! Power applied to value differences.
        real(dp), intent(in) :: lambda !! Edit penalty added to insertion and deletion steps.
        real(dp), intent(in) :: nu !! Stiffness penalty multiplying index displacement.
        real(dp), allocatable :: d(:, :)
        real(dp) :: ai0
        real(dp) :: bj0
        real(dp) :: del_a
        real(dp) :: del_b
        real(dp) :: match
        integer :: i
        integer :: j

        allocate(d(0:size(a), 0:size(b)))
        d = huge(1.0_dp)
        d(0, 0) = 0.0_dp
        do i = 1, size(a)
            ai0 = a(max(1, i - 1))
            if (i == 1) ai0 = 0.0_dp
            d(i, 0) = d(i - 1, 0) + abs(a(i) - ai0)**p + nu + lambda
        end do
        do j = 1, size(b)
            bj0 = b(max(1, j - 1))
            if (j == 1) bj0 = 0.0_dp
            d(0, j) = d(0, j - 1) + abs(b(j) - bj0)**p + nu + lambda
        end do
        do i = 1, size(a)
            ai0 = a(max(1, i - 1))
            if (i == 1) ai0 = 0.0_dp
            do j = 1, size(b)
                bj0 = b(max(1, j - 1))
                if (j == 1) bj0 = 0.0_dp
                del_a = d(i - 1, j) + abs(a(i) - ai0)**p + nu + lambda
                del_b = d(i, j - 1) + abs(b(j) - bj0)**p + nu + lambda
                match = d(i - 1, j - 1) + abs(a(i) - b(j))**p + abs(ai0 - bj0)**p &
                    + 2.0_dp*nu*abs(real(i - j, dp))
                d(i, j) = min(del_a, del_b, match)
            end do
        end do
        distance = d(size(a), size(b))
    end function twed_pair

    pure subroutine double_center(d, a)
        real(dp), intent(in) :: d(:, :) !! Square distance matrix to double-center.
        real(dp), allocatable, intent(out) :: a(:, :) !! Double-centered distance matrix.
        real(dp), allocatable :: rowm(:)
        real(dp), allocatable :: colm(:)
        real(dp) :: grand
        integer :: i
        integer :: j
        integer :: n

        n = size(d, 1)
        allocate(a(n, n), rowm(n), colm(n))
        rowm = sum(d, dim=2)/real(n, dp)
        colm = sum(d, dim=1)/real(n, dp)
        grand = sum(d)/real(n*n, dp)
        do i = 1, n
            do j = 1, n
                a(i, j) = d(i, j) - rowm(i) - colm(j) + grand
            end do
        end do
    end subroutine double_center

    pure subroutine astar2(d, a)
        real(dp), intent(in) :: d(:, :) !! Square distance matrix for upstream fda.usc bias correction.
        real(dp), allocatable, intent(out) :: a(:, :) !! Bias-corrected centered matrix used by bcdcor.dist.
        real(dp), allocatable :: m(:)
        real(dp) :: grand
        integer :: i
        integer :: j
        integer :: n

        n = size(d, 1)
        allocate(a(n, n), m(n))
        m = sum(d, dim=2)/real(n, dp)
        grand = sum(d)/real(n*n, dp)
        do i = 1, n
            do j = 1, n
                a(i, j) = d(i, j) - m(i) - m(j) + grand - d(i, j)/real(n, dp)
            end do
        end do
        do i = 1, n
            a(i, i) = m(i) - grand
        end do
        a = real(n, dp)/real(n - 1, dp)*a
    end subroutine astar2

    pure subroutine combine_distance_matrices(distance_components, component_weights, method, combined)
        real(dp), intent(in) :: distance_components(:, :, :) !! Component distance matrices with shape (n1,n2,ncomponents).
        real(dp), intent(in) :: component_weights(:) !! Nonnegative component weights aligned with the third array dimension.
        integer, intent(in) :: method !! Combination code: 1 Euclidean, 2 maximum, 3 Manhattan, or 4 Minkowski with p=2.
        real(dp), intent(out) :: combined(:, :) !! Combined distance matrix with shape (n1,n2).
        integer :: k

        select case (method)
        case (1, 4)
            combined = 0.0_dp
            do k = 1, size(distance_components, 3)
                combined = combined + component_weights(k)*distance_components(:, :, k)**2
            end do
            combined = sqrt(max(combined, 0.0_dp))
        case (2)
            combined = 0.0_dp
            do k = 1, size(distance_components, 3)
                combined = max(combined, component_weights(k)*distance_components(:, :, k))
            end do
        case (3)
            combined = 0.0_dp
            do k = 1, size(distance_components, 3)
                combined = combined + component_weights(k)*abs(distance_components(:, :, k))
            end do
        case default
            combined = 0.0_dp
        end select
    end subroutine combine_distance_matrices

end module fdausc_metrics
