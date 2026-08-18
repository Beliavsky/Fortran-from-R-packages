module bivgeo_distribution
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use bivgeo_kinds, only : dp
    use bivgeo_types, only : bivgeo_params, valid_bivgeo_params
    implicit none
    private

    public :: dbivgeo1
    public :: dbivgeo2
    public :: pbivgeo
    public :: sbivgeo

contains

    pure elemental function dbivgeo1(x, y, theta, log_p) result(value)
        integer, intent(in) :: x, y
        type(bivgeo_params), intent(in) :: theta
        logical, intent(in), optional :: log_p
        real(dp) :: value
        real(dp) :: p1, p2, p3, p4, pmf
        integer :: z1, z2, z3, z4
        logical :: give_log

        give_log = .false.
        if (present(log_p)) give_log = log_p

        if (.not. valid_bivgeo_params(theta)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        if (x < 1 .or. y < 1) then
            if (give_log) then
                value = -huge(1.0_dp)
            else
                value = 0.0_dp
            end if
            return
        end if

        z1 = max(x - 1, y - 1)
        z2 = max(x, y - 1)
        z3 = max(x - 1, y)
        z4 = max(x, y)

        p1 = theta%theta1**(x - 1) * theta%theta2**(y - 1) * theta%theta3**z1
        p2 = theta%theta1**x * theta%theta2**(y - 1) * theta%theta3**z2
        p3 = theta%theta1**(x - 1) * theta%theta2**y * theta%theta3**z3
        p4 = theta%theta1**x * theta%theta2**y * theta%theta3**z4
        pmf = max(0.0_dp, p1 - p2 - p3 + p4)

        if (give_log) then
            if (pmf > 0.0_dp) then
                value = log(pmf)
            else
                value = -huge(1.0_dp)
            end if
        else
            value = pmf
        end if
    end function dbivgeo1

    pure elemental function dbivgeo2(x, y, theta, log_p) result(value)
        integer, intent(in) :: x, y
        type(bivgeo_params), intent(in) :: theta
        logical, intent(in), optional :: log_p
        real(dp) :: value
        real(dp) :: gamma1, gamma2, gamma3, pmf, logpmf
        logical :: give_log

        give_log = .false.
        if (present(log_p)) give_log = log_p

        if (.not. valid_bivgeo_params(theta)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        if (x < 1 .or. y < 1) then
            if (give_log) then
                value = -huge(1.0_dp)
            else
                value = 0.0_dp
            end if
            return
        end if

        gamma1 = theta%theta1 * theta%theta3
        gamma2 = theta%theta2 * theta%theta3
        gamma3 = theta%theta1 * theta%theta2 * theta%theta3

        if (x < y) then
            logpmf = real(x - 1, dp) * log(theta%theta1)
            logpmf = logpmf + real(y - 1, dp) * log(gamma2)
            logpmf = logpmf + log(1.0_dp - gamma2) + log(1.0_dp - theta%theta1)
        else if (x > y) then
            logpmf = real(y - 1, dp) * log(theta%theta2)
            logpmf = logpmf + real(x - 1, dp) * log(gamma1)
            logpmf = logpmf + log(1.0_dp - gamma1) + log(1.0_dp - theta%theta2)
        else
            pmf = 1.0_dp - gamma1 - gamma2 + gamma3
            if (pmf <= 0.0_dp) then
                if (give_log) then
                    value = -huge(1.0_dp)
                else
                    value = 0.0_dp
                end if
                return
            end if
            logpmf = real(x - 1, dp) * log(gamma3) + log(pmf)
        end if

        if (give_log) then
            value = logpmf
        else
            value = exp(logpmf)
        end if
    end function dbivgeo2

    pure elemental function sbivgeo(x, y, theta) result(value)
        integer, intent(in) :: x, y
        type(bivgeo_params), intent(in) :: theta
        real(dp) :: value
        integer :: xx, yy

        if (.not. valid_bivgeo_params(theta)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        xx = max(x, 0)
        yy = max(y, 0)
        value = theta%theta1**xx * theta%theta2**yy * theta%theta3**max(xx, yy)
    end function sbivgeo

    pure elemental function pbivgeo(x, y, theta, lower_tail) result(value)
        integer, intent(in) :: x, y
        type(bivgeo_params), intent(in) :: theta
        logical, intent(in), optional :: lower_tail
        real(dp) :: value
        real(dp) :: sf, sfx, sfy
        logical :: lower

        lower = .true.
        if (present(lower_tail)) lower = lower_tail

        if (.not. valid_bivgeo_params(theta)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        if (.not. lower) then
            value = sbivgeo(x, y, theta)
            return
        end if

        if (x < 1 .or. y < 1) then
            value = 0.0_dp
            return
        end if

        sf = sbivgeo(x, y, theta)
        sfx = (theta%theta1 * theta%theta3)**x
        sfy = (theta%theta2 * theta%theta3)**y
        value = 1.0_dp - sfx - sfy + sf
        value = min(1.0_dp, max(0.0_dp, value))
    end function pbivgeo

end module bivgeo_distribution
