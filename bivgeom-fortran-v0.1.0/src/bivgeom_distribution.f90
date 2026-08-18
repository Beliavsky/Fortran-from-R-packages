! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
module bivgeom_distribution
    use bivgeom_kinds, only : dp
    use bivgeom_math, only : qgeom_theta
    implicit none
    private

    public :: feasible_roy
    public :: dbivgeom_roy
    public :: log_dbivgeom_roy
    public :: fbivgeom_roy
    public :: sbivgeom_roy
    public :: fyxbivgeom_roy
    public :: eyxbivgeom_roy
    public :: lambda1_roy
    public :: lambda2_roy
    public :: corbivgeom_roy
    public :: relbivgeom_roy
    public :: rbivgeom_roy
    public :: empirical_survival_roy

    interface dbivgeom_roy
        module procedure dbivgeom_roy_i
        module procedure dbivgeom_roy_r
    end interface dbivgeom_roy

    interface log_dbivgeom_roy
        module procedure log_dbivgeom_roy_i
        module procedure log_dbivgeom_roy_r
    end interface log_dbivgeom_roy

contains

    pure logical function feasible_roy(theta1, theta2, theta3) result(ok)
        real(dp), intent(in) :: theta1, theta2, theta3
        real(dp) :: lower

        ok = .false.
        if (theta1 <= 0.0_dp .or. theta1 >= 1.0_dp) return
        if (theta2 <= 0.0_dp .or. theta2 >= 1.0_dp) return
        if (theta3 <= 0.0_dp .or. theta3 > 1.0_dp) return
        lower = (theta1 + theta2 - 1.0_dp) / (theta1 * theta2)
        if (theta3 < max(0.0_dp, lower)) return
        ok = .true.
    end function feasible_roy

    pure real(dp) function dbivgeom_roy_i(x, y, theta1, theta2, theta3) result(p)
        integer, intent(in) :: x, y
        real(dp), intent(in) :: theta1, theta2, theta3
        real(dp) :: lp

        lp = log_dbivgeom_roy_i(x, y, theta1, theta2, theta3)
        if (lp <= -huge(1.0_dp) / 2.0_dp) then
            p = 0.0_dp
        else
            p = exp(lp)
        end if
    end function dbivgeom_roy_i

    pure real(dp) function dbivgeom_roy_r(x, y, theta1, theta2, theta3) result(p)
        real(dp), intent(in) :: x, y
        real(dp), intent(in) :: theta1, theta2, theta3
        integer :: ix, iy

        if (x < 0.0_dp .or. y < 0.0_dp) then
            p = 0.0_dp
            return
        end if
        ix = nint(x)
        iy = nint(y)
        if (abs(x - real(ix, dp)) > 32.0_dp * epsilon(1.0_dp) .or. &
            abs(y - real(iy, dp)) > 32.0_dp * epsilon(1.0_dp)) then
            p = 0.0_dp
            return
        end if
        p = dbivgeom_roy_i(ix, iy, theta1, theta2, theta3)
    end function dbivgeom_roy_r

    pure real(dp) function log_dbivgeom_roy_i(x, y, theta1, theta2, theta3) result(lp)
        integer, intent(in) :: x, y
        real(dp), intent(in) :: theta1, theta2, theta3
        real(dp) :: bracket

        if (x < 0 .or. y < 0 .or. .not. feasible_roy(theta1, theta2, theta3)) then
            lp = -huge(1.0_dp)
            return
        end if

        bracket = 1.0_dp - theta1 * theta3**y - theta2 * theta3**x + &
            theta1 * theta2 * theta3**(x + y + 1)
        if (bracket <= 0.0_dp) then
            lp = -huge(1.0_dp)
            return
        end if

        lp = real(x, dp) * log(theta1) + real(y, dp) * log(theta2) + &
            real(x * y, dp) * log(theta3) + log(bracket)
    end function log_dbivgeom_roy_i

    pure real(dp) function log_dbivgeom_roy_r(x, y, theta1, theta2, theta3) result(lp)
        real(dp), intent(in) :: x, y
        real(dp), intent(in) :: theta1, theta2, theta3
        integer :: ix, iy

        if (x < 0.0_dp .or. y < 0.0_dp) then
            lp = -huge(1.0_dp)
            return
        end if
        ix = nint(x)
        iy = nint(y)
        if (abs(x - real(ix, dp)) > 32.0_dp * epsilon(1.0_dp) .or. &
            abs(y - real(iy, dp)) > 32.0_dp * epsilon(1.0_dp)) then
            lp = -huge(1.0_dp)
            return
        end if
        lp = log_dbivgeom_roy_i(ix, iy, theta1, theta2, theta3)
    end function log_dbivgeom_roy_r

    pure real(dp) function fbivgeom_roy(x, y, theta1, theta2, theta3) result(cdf)
        real(dp), intent(in) :: x, y
        real(dp), intent(in) :: theta1, theta2, theta3
        integer :: ix, iy

        if (.not. feasible_roy(theta1, theta2, theta3)) then
            cdf = -1.0_dp
            return
        end if
        if (x < 0.0_dp .or. y < 0.0_dp) then
            cdf = 0.0_dp
            return
        end if
        ix = floor(x)
        iy = floor(y)
        cdf = 1.0_dp - theta1**(ix + 1) - theta2**(iy + 1) + &
            theta1**(ix + 1) * theta2**(iy + 1) * theta3**((ix + 1) * (iy + 1))
        cdf = min(1.0_dp, max(0.0_dp, cdf))
    end function fbivgeom_roy

    pure real(dp) function sbivgeom_roy(x, y, theta1, theta2, theta3) result(s)
        real(dp), intent(in) :: x, y
        real(dp), intent(in) :: theta1, theta2, theta3
        integer :: a, b

        if (.not. feasible_roy(theta1, theta2, theta3)) then
            s = -1.0_dp
            return
        end if
        a = ceiling(max(x, 0.0_dp))
        b = ceiling(max(y, 0.0_dp))
        s = theta1**a * theta2**b * theta3**(a * b)
    end function sbivgeom_roy

    pure real(dp) function fyxbivgeom_roy(y, theta1, theta2, theta3, x) result(cdf)
        real(dp), intent(in) :: y, theta1, theta2, theta3
        integer, intent(in) :: x
        real(dp) :: a, b

        if (.not. feasible_roy(theta1, theta2, theta3) .or. x < 0) then
            cdf = -1.0_dp
            return
        end if
        if (y < 0.0_dp) then
            cdf = 0.0_dp
            return
        end if
        a = theta2 * theta3**x
        b = theta2 * theta3**(x + 1)
        cdf = (1.0_dp - a**(y + 1.0_dp) - theta1 * (1.0_dp - b**(y + 1.0_dp))) / &
            (1.0_dp - theta1)
        cdf = min(1.0_dp, max(0.0_dp, cdf))
    end function fyxbivgeom_roy

    pure real(dp) function eyxbivgeom_roy(theta1, theta2, theta3, x) result(mean_y)
        real(dp), intent(in) :: theta1, theta2, theta3
        integer, intent(in) :: x
        real(dp) :: a, b

        if (.not. feasible_roy(theta1, theta2, theta3) .or. x < 0) then
            mean_y = -1.0_dp
            return
        end if
        a = theta2 * theta3**x
        b = theta2 * theta3**(x + 1)
        mean_y = a / (1.0_dp - a) * &
            (1.0_dp - theta1 * theta3 * (1.0_dp - a) / (1.0_dp - b)) / &
            (1.0_dp - theta1)
    end function eyxbivgeom_roy

    pure real(dp) function lambda1_roy(y, theta1, theta3) result(lambda)
        integer, intent(in) :: y
        real(dp), intent(in) :: theta1, theta3

        lambda = 1.0_dp - theta1 * theta3**y
    end function lambda1_roy

    pure real(dp) function lambda2_roy(x, theta2, theta3) result(lambda)
        integer, intent(in) :: x
        real(dp), intent(in) :: theta2, theta3

        lambda = 1.0_dp - theta2 * theta3**x
    end function lambda2_roy

    real(dp) function corbivgeom_roy(theta1, theta2, theta3, alpha) result(corr)
        real(dp), intent(in) :: theta1, theta2, theta3
        real(dp), intent(in), optional :: alpha
        real(dp) :: a, exy, mx, my, vx, vy
        integer :: k, kmax

        if (.not. feasible_roy(theta1, theta2, theta3)) then
            corr = -huge(1.0_dp)
            return
        end if
        if (abs(theta3 - 1.0_dp) <= 16.0_dp * epsilon(1.0_dp)) then
            corr = 0.0_dp
            return
        end if

        a = 1.0e-5_dp
        if (present(alpha)) a = alpha
        kmax = max(2 * qgeom_theta(1.0_dp - a, theta1), &
            2 * qgeom_theta(1.0_dp - a, theta2))
        exy = 0.0_dp
        do k = 0, kmax
            exy = exy + theta2 * (theta1 * theta3)**(k + 1) / &
                (1.0_dp - theta2 * theta3**(k + 1))
        end do

        mx = theta1 / (1.0_dp - theta1)
        my = theta2 / (1.0_dp - theta2)
        vx = theta1 / (1.0_dp - theta1)**2
        vy = theta2 / (1.0_dp - theta2)**2
        corr = (exy - mx * my) / sqrt(vx * vy)
    end function corbivgeom_roy

    real(dp) function relbivgeom_roy(theta1, theta2, theta3, alpha) result(rel)
        real(dp), intent(in) :: theta1, theta2, theta3
        real(dp), intent(in), optional :: alpha
        real(dp) :: a, total
        integer :: k, kmax

        if (.not. feasible_roy(theta1, theta2, theta3)) then
            rel = -huge(1.0_dp)
            return
        end if

        a = 1.0e-5_dp
        if (present(alpha)) a = alpha
        kmax = max(2 * qgeom_theta(1.0_dp - a, theta1), &
            2 * qgeom_theta(1.0_dp - a, theta2))
        total = 0.0_dp
        do k = 0, kmax
            total = total + theta2**k * (theta1 * theta3**k)**(k + 1) * &
                (1.0_dp - theta3**(k + 1) * theta2)
        end do
        rel = 1.0_dp - total
    end function relbivgeom_roy

    subroutine rbivgeom_roy(n, theta1, theta2, theta3, sample, ok)
        integer, intent(in) :: n
        real(dp), intent(in) :: theta1, theta2, theta3
        integer, allocatable, intent(out) :: sample(:, :)
        logical, intent(out), optional :: ok
        real(dp) :: u, v
        integer :: i, x, y
        logical :: valid

        valid = n >= 0 .and. feasible_roy(theta1, theta2, theta3)
        if (present(ok)) ok = valid
        if (.not. valid) then
            allocate(sample(0, 2))
            return
        end if

        allocate(sample(n, 2))
        do i = 1, n
            call random_number(u)
            call random_number(v)
            x = qgeom_theta(min(u, 1.0_dp - epsilon(1.0_dp)), theta1)
            y = conditional_quantile(v, theta1, theta2, theta3, x)
            sample(i, 1) = x
            sample(i, 2) = y
        end do
    end subroutine rbivgeom_roy

    integer function conditional_quantile(prob, theta1, theta2, theta3, x) result(q)
        real(dp), intent(in) :: prob, theta1, theta2, theta3
        integer, intent(in) :: x
        integer :: lo, hi, mid

        if (prob <= 0.0_dp) then
            q = 0
            return
        end if

        lo = -1
        hi = max(1, qgeom_theta(min(prob, 1.0_dp - 1.0e-12_dp), theta2) + 2)
        do while (fyxbivgeom_roy(real(hi, dp), theta1, theta2, theta3, x) < prob)
            if (hi > 1000000000) exit
            hi = 2 * hi + 1
        end do

        do while (hi - lo > 1)
            mid = lo + (hi - lo) / 2
            if (fyxbivgeom_roy(real(mid, dp), theta1, theta2, theta3, x) >= prob) then
                hi = mid
            else
                lo = mid
            end if
        end do
        q = hi
    end function conditional_quantile

    pure real(dp) function empirical_survival_roy(xq, yq, x, y) result(s)
        integer, intent(in) :: xq, yq
        integer, intent(in) :: x(:), y(:)
        integer :: i, count

        if (size(x) /= size(y) .or. size(x) == 0) then
            s = -1.0_dp
            return
        end if
        count = 0
        do i = 1, size(x)
            if (x(i) >= xq .and. y(i) >= yq) count = count + 1
        end do
        s = real(count, dp) / real(size(x), dp)
    end function empirical_survival_roy

end module bivgeom_distribution
