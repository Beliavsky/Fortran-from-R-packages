module bivgeo_moments
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use bivgeo_kinds, only : dp
    use bivgeo_types, only : bivgeo_params, valid_bivgeo_params
    implicit none
    private

    public :: cfbivgeo
    public :: covbivgeo
    public :: corbivgeo
    public :: mean_bivgeo
    public :: variance_bivgeo
    public :: mombivgeo

contains

    pure function cfbivgeo(theta) result(value)
        type(bivgeo_params), intent(in) :: theta
        real(dp) :: value
        real(dp) :: p1, p2

        if (.not. valid_bivgeo_params(theta)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        p1 = 1.0_dp - theta%theta1 * theta%theta2 * theta%theta3**2
        p2 = (1.0_dp - theta%theta1 * theta%theta3)
        p2 = p2 * (1.0_dp - theta%theta2 * theta%theta3)
        p2 = p2 * (1.0_dp - theta%theta1 * theta%theta2 * theta%theta3)
        value = p1 / p2
    end function cfbivgeo

    pure function covbivgeo(theta) result(value)
        type(bivgeo_params), intent(in) :: theta
        real(dp) :: value
        real(dp) :: p1, p2

        if (.not. valid_bivgeo_params(theta)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        p1 = (1.0_dp - theta%theta3) * theta%theta1 * theta%theta2 * theta%theta3
        p2 = (1.0_dp - theta%theta1 * theta%theta3)
        p2 = p2 * (1.0_dp - theta%theta2 * theta%theta3)
        p2 = p2 * (1.0_dp - theta%theta1 * theta%theta2 * theta%theta3)
        value = p1 / p2
    end function covbivgeo

    pure function corbivgeo(theta) result(value)
        type(bivgeo_params), intent(in) :: theta
        real(dp) :: value
        real(dp) :: denom

        if (.not. valid_bivgeo_params(theta)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        denom = 1.0_dp - theta%theta1 * theta%theta2 * theta%theta3
        value = (1.0_dp - theta%theta3) * sqrt(theta%theta1 * theta%theta2) / denom
    end function corbivgeo

    pure function mean_bivgeo(theta) result(mu)
        type(bivgeo_params), intent(in) :: theta
        real(dp) :: mu(2)

        if (.not. valid_bivgeo_params(theta)) then
            mu = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        mu(1) = 1.0_dp / (1.0_dp - theta%theta1 * theta%theta3)
        mu(2) = 1.0_dp / (1.0_dp - theta%theta2 * theta%theta3)
    end function mean_bivgeo

    pure function variance_bivgeo(theta) result(var)
        type(bivgeo_params), intent(in) :: theta
        real(dp) :: var(2)
        real(dp) :: q1, q2

        if (.not. valid_bivgeo_params(theta)) then
            var = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        q1 = theta%theta1 * theta%theta3
        q2 = theta%theta2 * theta%theta3
        var(1) = q1 / (1.0_dp - q1)**2
        var(2) = q2 / (1.0_dp - q2)**2
    end function variance_bivgeo

    subroutine mombivgeo(x, y, estimate, ok)
        integer, intent(in) :: x(:), y(:)
        type(bivgeo_params), intent(out) :: estimate
        logical, intent(out), optional :: ok
        real(dp) :: xbar, ybar, zbar
        integer :: i, n
        logical :: good

        n = size(x)
        good = n > 0 .and. size(y) == n
        if (.not. good) then
            estimate%theta1 = ieee_value(0.0_dp, ieee_quiet_nan)
            estimate%theta2 = ieee_value(0.0_dp, ieee_quiet_nan)
            estimate%theta3 = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(ok)) ok = .false.
            return
        end if

        xbar = sum(real(x, dp)) / real(n, dp)
        ybar = sum(real(y, dp)) / real(n, dp)
        zbar = 0.0_dp
        do i = 1, n
            zbar = zbar + real(min(x(i), y(i)), dp)
        end do
        zbar = zbar / real(n, dp)

        good = xbar > 1.0_dp .and. ybar > 1.0_dp .and. zbar > 1.0_dp
        if (good) then
            estimate%theta1 = ybar * (1.0_dp - zbar) / (zbar * (1.0_dp - ybar))
            estimate%theta2 = xbar * (zbar - 1.0_dp) / (zbar * (xbar - 1.0_dp))
            estimate%theta3 = zbar * (xbar - 1.0_dp) * (ybar - 1.0_dp)
            estimate%theta3 = estimate%theta3 / (xbar * ybar * (zbar - 1.0_dp))
            good = valid_bivgeo_params(estimate)
        end if

        if (.not. good) then
            estimate%theta1 = ieee_value(0.0_dp, ieee_quiet_nan)
            estimate%theta2 = ieee_value(0.0_dp, ieee_quiet_nan)
            estimate%theta3 = ieee_value(0.0_dp, ieee_quiet_nan)
        end if
        if (present(ok)) ok = good
    end subroutine mombivgeo

end module bivgeo_moments
