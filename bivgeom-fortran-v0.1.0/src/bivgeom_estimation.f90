! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
module bivgeom_estimation
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use bivgeom_kinds, only : dp
    use bivgeom_math, only : solve_3x3, upper_string
    use bivgeom_distribution, only : feasible_roy, log_dbivgeom_roy, empirical_survival_roy
    implicit none
    private

    type, public :: bivgeom_fit
        real(dp) :: theta(3) = 0.0_dp
        real(dp) :: nll = huge(1.0_dp)
        integer :: iterations = 0
        logical :: converged = .false.
    end type bivgeom_fit

    public :: negative_loglik_roy
    public :: fit_bivgeom_ml
    public :: estbivgeom_roy
    public :: estimate_ls_roy
    public :: estimate_mmp_roy
    public :: estimate_mm1_roy
    public :: estimate_mm2_roy
    public :: estimate_mm3_roy
    public :: estimate_mm4_roy

contains

    real(dp) function negative_loglik_roy(theta, x, y) result(nll)
        real(dp), intent(in) :: theta(3)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: lp
        integer :: i

        if (size(x) /= size(y) .or. size(x) == 0) then
            nll = huge(1.0_dp)
            return
        end if
        if (.not. feasible_roy(theta(1), theta(2), theta(3))) then
            nll = huge(1.0_dp)
            return
        end if

        nll = 0.0_dp
        do i = 1, size(x)
            lp = log_dbivgeom_roy(x(i), y(i), theta(1), theta(2), theta(3))
            if (lp <= -huge(1.0_dp) / 2.0_dp) then
                nll = huge(1.0_dp)
                return
            end if
            nll = nll - lp
        end do
    end function negative_loglik_roy

    function fit_bivgeom_ml(x, y, start, max_iter, tol) result(fit)
        integer, intent(in) :: x(:), y(:)
        real(dp), intent(in), optional :: start(3)
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        type(bivgeom_fit) :: fit
        real(dp) :: simplex(3, 4), f(4), centroid(3), xr(3), xe(3), xc(3)
        real(dp) :: fr, fe, fc, alpha, gamma, rho, sigma, eps, tolerance
        real(dp) :: spread_x, spread_f
        integer :: j, iter, imax

        if (size(x) /= size(y) .or. size(x) == 0) return
        if (any(x < 0) .or. any(y < 0)) return

        eps = 1.0e-8_dp
        tolerance = 1.0e-9_dp
        if (present(tol)) tolerance = max(tol, 100.0_dp * epsilon(1.0_dp))
        imax = 4000
        if (present(max_iter)) imax = max(1, max_iter)

        simplex(:, 1) = [0.5_dp, 0.5_dp, 1.0_dp]
        if (present(start)) simplex(:, 1) = project_theta(start, eps)
        do j = 1, 3
            simplex(:, j + 1) = simplex(:, 1)
            if (j < 3) then
                simplex(j, j + 1) = min(1.0_dp - eps, simplex(j, j + 1) + 0.08_dp)
            else
                simplex(j, j + 1) = max(eps, simplex(j, j + 1) - 0.08_dp)
            end if
            if (.not. feasible_roy(simplex(1, j + 1), simplex(2, j + 1), simplex(3, j + 1))) then
                simplex(:, j + 1) = feasible_nudge(simplex(:, j + 1), eps)
            end if
        end do

        do j = 1, 4
            f(j) = penalized_nll(simplex(:, j), x, y)
        end do

        alpha = 1.0_dp
        gamma = 2.0_dp
        rho = 0.5_dp
        sigma = 0.5_dp

        do iter = 1, imax
            call sort_simplex(simplex, f)
            spread_f = maxval(abs(f - f(1)))
            spread_x = 0.0_dp
            do j = 2, 4
                spread_x = max(spread_x, maxval(abs(simplex(:, j) - simplex(:, 1))))
            end do
            if (spread_f <= tolerance * (1.0_dp + abs(f(1))) .and. spread_x <= sqrt(tolerance)) then
                fit%converged = .true.
                exit
            end if

            centroid = 0.0_dp
            do j = 1, 3
                centroid = centroid + simplex(:, j)
            end do
            centroid = centroid / 3.0_dp

            xr = centroid + alpha * (centroid - simplex(:, 4))
            xr = project_theta(xr, eps)
            fr = penalized_nll(xr, x, y)

            if (fr < f(1)) then
                xe = centroid + gamma * (xr - centroid)
                xe = project_theta(xe, eps)
                fe = penalized_nll(xe, x, y)
                if (fe < fr) then
                    simplex(:, 4) = xe
                    f(4) = fe
                else
                    simplex(:, 4) = xr
                    f(4) = fr
                end if
            else if (fr < f(3)) then
                simplex(:, 4) = xr
                f(4) = fr
            else
                if (fr < f(4)) then
                    xc = centroid + rho * (xr - centroid)
                else
                    xc = centroid + rho * (simplex(:, 4) - centroid)
                end if
                xc = project_theta(xc, eps)
                fc = penalized_nll(xc, x, y)
                if (fc < min(fr, f(4))) then
                    simplex(:, 4) = xc
                    f(4) = fc
                else
                    do j = 2, 4
                        simplex(:, j) = simplex(:, 1) + sigma * (simplex(:, j) - simplex(:, 1))
                        simplex(:, j) = project_theta(simplex(:, j), eps)
                        f(j) = penalized_nll(simplex(:, j), x, y)
                    end do
                end if
            end if
        end do

        call sort_simplex(simplex, f)
        fit%theta = simplex(:, 1)
        fit%nll = f(1)
        fit%iterations = min(iter, imax)
        if (.not. fit%converged) then
            fit%converged = isfinite_nll(f(1)) .and. &
                maxval(abs(simplex(:, 4) - simplex(:, 1))) < 1.0e-5_dp
        end if
    end function fit_bivgeom_ml

    function estbivgeom_roy(x, y, method) result(theta)
        integer, intent(in) :: x(:), y(:)
        character(len=*), intent(in) :: method
        real(dp) :: theta(3)
        character(len=:), allocatable :: m
        type(bivgeom_fit) :: fit

        m = trim(upper_string(trim(method)))
        select case (m)
        case ('ML')
            fit = fit_bivgeom_ml(x, y)
            theta = fit%theta
        case ('LS')
            theta = estimate_ls_roy(x, y)
        case ('MMP')
            theta = estimate_mmp_roy(x, y)
        case ('MM1', 'M1')
            theta = estimate_mm1_roy(x, y)
        case ('MM2', 'M2')
            theta = estimate_mm2_roy(x, y)
        case ('MM3', 'M3')
            theta = estimate_mm3_roy(x, y)
        case ('MM4', 'M4')
            theta = estimate_mm4_roy(x, y)
        case default
            theta = nan3()
        end select
    end function estbivgeom_roy

    function estimate_ls_roy(x, y) result(theta)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: theta(3)
        real(dp) :: xtx(3, 3), xtz(3), beta(3), row(3), s
        logical :: ok
        integer :: i, j, k

        if (.not. valid_sample(x, y)) then
            theta = nan3()
            return
        end if

        xtx = 0.0_dp
        xtz = 0.0_dp
        do i = 1, size(x)
            s = empirical_survival_roy(x(i), y(i), x, y)
            if (s <= 0.0_dp) then
                theta = nan3()
                return
            end if
            row = [real(x(i), dp), real(y(i), dp), real(x(i) * y(i), dp)]
            do j = 1, 3
                xtz(j) = xtz(j) + row(j) * log(s)
                do k = 1, 3
                    xtx(j, k) = xtx(j, k) + row(j) * row(k)
                end do
            end do
        end do

        call solve_3x3(xtx, xtz, beta, ok)
        if (ok) then
            theta = exp(beta)
        else
            theta = nan3()
        end if
    end function estimate_ls_roy

    function estimate_mmp_roy(x, y) result(theta)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: theta(3)
        real(dp) :: mx, my, p11
        integer :: i, count

        if (.not. valid_sample(x, y)) then
            theta = nan3()
            return
        end if
        mx = mean_int(x)
        my = mean_int(y)
        theta(1) = mx / (1.0_dp + mx)
        theta(2) = my / (1.0_dp + my)
        count = 0
        do i = 1, size(x)
            if (x(i) >= 1 .and. y(i) >= 1) count = count + 1
        end do
        p11 = real(count, dp) / real(size(x), dp)
        if (theta(1) <= 0.0_dp .or. theta(2) <= 0.0_dp) then
            theta(3) = ieee_value(1.0_dp, ieee_quiet_nan)
        else
            theta(3) = p11 / (theta(1) * theta(2))
        end if
    end function estimate_mmp_roy

    function estimate_mm1_roy(x, y) result(theta)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: theta(3)
        real(dp) :: ex

        if (.not. valid_sample(x, y)) then
            theta = nan3()
            return
        end if
        theta = marginal_thetas(x, y)
        ex = conditional_positive_mean(x, x, y)
        if (.not. (ex > 0.0_dp) .or. theta(1) <= 0.0_dp) then
            theta(3) = ieee_value(1.0_dp, ieee_quiet_nan)
        else
            theta(3) = (ex - 1.0_dp) / (theta(1) * ex)
        end if
    end function estimate_mm1_roy

    function estimate_mm2_roy(x, y) result(theta)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: theta(3)
        real(dp) :: ey

        if (.not. valid_sample(x, y)) then
            theta = nan3()
            return
        end if
        theta = marginal_thetas(x, y)
        ey = conditional_positive_mean(y, x, y)
        if (.not. (ey > 0.0_dp) .or. theta(2) <= 0.0_dp) then
            theta(3) = ieee_value(1.0_dp, ieee_quiet_nan)
        else
            theta(3) = (ey - 1.0_dp) / (theta(2) * ey)
        end if
    end function estimate_mm2_roy

    function estimate_mm3_roy(x, y) result(theta)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: theta(3)
        real(dp) :: ex, ey, den

        if (.not. valid_sample(x, y)) then
            theta = nan3()
            return
        end if
        theta = marginal_thetas(x, y)
        ex = conditional_positive_mean(x, x, y)
        ey = conditional_positive_mean(y, x, y)
        den = theta(1) * ex + theta(2) * ey
        if (.not. (den > 0.0_dp)) then
            theta(3) = ieee_value(1.0_dp, ieee_quiet_nan)
        else
            theta(3) = (ex + ey - 2.0_dp) / den
        end if
    end function estimate_mm3_roy

    function estimate_mm4_roy(x, y) result(theta)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: theta(3)
        real(dp) :: ex, ey, s, disc, den

        if (.not. valid_sample(x, y)) then
            theta = nan3()
            return
        end if
        theta = marginal_thetas(x, y)
        ex = conditional_positive_mean(x, x, y)
        ey = conditional_positive_mean(y, x, y)
        s = ex + ey
        den = 2.0_dp * s * theta(1) * theta(2)
        disc = s * (s - 2.0_dp) * (theta(1) - theta(2))**2 + &
            (theta(1) + theta(2))**2
        if (.not. (den > 0.0_dp) .or. disc < 0.0_dp) then
            theta(3) = ieee_value(1.0_dp, ieee_quiet_nan)
        else
            theta(3) = ((theta(1) + theta(2)) * (s - 1.0_dp) - sqrt(disc)) / den
        end if
    end function estimate_mm4_roy

    pure function marginal_thetas(x, y) result(theta)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: theta(3)
        real(dp) :: mx, my

        mx = mean_int(x)
        my = mean_int(y)
        theta(1) = mx / (1.0_dp + mx)
        theta(2) = my / (1.0_dp + my)
        theta(3) = ieee_value(1.0_dp, ieee_quiet_nan)
    end function marginal_thetas

    pure real(dp) function conditional_positive_mean(values, x, y) result(m)
        integer, intent(in) :: values(:), x(:), y(:)
        integer :: i, count
        real(dp) :: total

        total = 0.0_dp
        count = 0
        do i = 1, size(values)
            if (x(i) > 0 .and. y(i) > 0) then
                total = total + real(values(i), dp)
                count = count + 1
            end if
        end do
        if (count == 0) then
            m = ieee_value(1.0_dp, ieee_quiet_nan)
        else
            m = total / real(count, dp)
        end if
    end function conditional_positive_mean

    pure real(dp) function mean_int(x) result(m)
        integer, intent(in) :: x(:)
        integer :: i
        real(dp) :: total

        total = 0.0_dp
        do i = 1, size(x)
            total = total + real(x(i), dp)
        end do
        m = total / real(size(x), dp)
    end function mean_int

    pure logical function valid_sample(x, y) result(ok)
        integer, intent(in) :: x(:), y(:)

        ok = size(x) == size(y)
        if (.not. ok) return
        ok = size(x) > 0
        if (.not. ok) return
        ok = .not. any(x < 0) .and. .not. any(y < 0)
    end function valid_sample

    pure function nan3() result(x)
        real(dp) :: x(3)
        real(dp) :: nan

        nan = ieee_value(1.0_dp, ieee_quiet_nan)
        x = nan
    end function nan3

    real(dp) function penalized_nll(theta, x, y) result(v)
        real(dp), intent(in) :: theta(3)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: lower, violation

        if (theta(1) <= 0.0_dp .or. theta(1) >= 1.0_dp .or. &
            theta(2) <= 0.0_dp .or. theta(2) >= 1.0_dp .or. &
            theta(3) <= 0.0_dp .or. theta(3) > 1.0_dp) then
            v = 1.0e100_dp
            return
        end if
        lower = max(0.0_dp, (theta(1) + theta(2) - 1.0_dp) / (theta(1) * theta(2)))
        violation = lower - theta(3)
        if (violation > 0.0_dp) then
            v = 1.0e50_dp * (1.0_dp + violation * violation)
        else
            v = negative_loglik_roy(theta, x, y)
        end if
    end function penalized_nll

    pure function project_theta(theta, eps) result(out)
        real(dp), intent(in) :: theta(3), eps
        real(dp) :: out(3), lower

        out(1) = min(1.0_dp - eps, max(eps, theta(1)))
        out(2) = min(1.0_dp - eps, max(eps, theta(2)))
        out(3) = min(1.0_dp, max(eps, theta(3)))
        lower = max(0.0_dp, (out(1) + out(2) - 1.0_dp) / (out(1) * out(2)))
        if (out(3) < lower) out(3) = min(1.0_dp, lower + eps)
    end function project_theta

    pure function feasible_nudge(theta, eps) result(out)
        real(dp), intent(in) :: theta(3), eps
        real(dp) :: out(3)

        out = project_theta(theta, eps)
    end function feasible_nudge

    subroutine sort_simplex(simplex, f)
        real(dp), intent(inout) :: simplex(3, 4), f(4)
        real(dp) :: xtmp(3), ftmp
        integer :: i, j

        do i = 2, 4
            xtmp = simplex(:, i)
            ftmp = f(i)
            j = i - 1
            do while (j >= 1)
                if (f(j) <= ftmp) exit
                simplex(:, j + 1) = simplex(:, j)
                f(j + 1) = f(j)
                j = j - 1
            end do
            simplex(:, j + 1) = xtmp
            f(j + 1) = ftmp
        end do
    end subroutine sort_simplex

    pure logical function isfinite_nll(x) result(ok)
        real(dp), intent(in) :: x

        ok = x < huge(1.0_dp) / 10.0_dp
    end function isfinite_nll

end module bivgeom_estimation
