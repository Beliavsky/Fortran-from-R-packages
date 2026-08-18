! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
module bivgeom_compat
    use bivgeom_kinds, only : dp
    use bivgeom_distribution, only : dbivgeomroy => dbivgeom_roy, &
        fbivgeomroy => fbivgeom_roy, sbivgeomroy => sbivgeom_roy, &
        fyxbivgeomroy => fyxbivgeom_roy, eyxbivgeomroy => eyxbivgeom_roy, &
        corbivgeomroy => corbivgeom_roy, relbivgeomroy => relbivgeom_roy, &
        rbivgeomroy => rbivgeom_roy, empirical_survival_roy
    use bivgeom_estimation, only : estbivgeomroy => estbivgeom_roy, &
        loglikgeomroy => negative_loglik_roy
    implicit none
    private

    public :: dbivgeomroy
    public :: fbivgeomroy
    public :: sbivgeomroy
    public :: fyxbivgeomroy
    public :: eyxbivgeomroy
    public :: corbivgeomroy
    public :: relbivgeomroy
    public :: rbivgeomroy
    public :: estbivgeomroy
    public :: loglikgeomroy
    public :: minuslogroy
    public :: lambda1roy
    public :: lambda2roy
    public :: s_n

contains

    real(dp) function minuslogroy(x, y, theta1, theta2, theta3) result(nll)
        integer, intent(in) :: x(:), y(:)
        real(dp), intent(in), optional :: theta1, theta2, theta3
        real(dp) :: theta(3)

        theta = [0.5_dp, 0.5_dp, 1.0_dp]
        if (present(theta1)) theta(1) = theta1
        if (present(theta2)) theta(2) = theta2
        if (present(theta3)) theta(3) = theta3
        nll = loglikgeomroy(theta, x, y)
    end function minuslogroy

    pure real(dp) function lambda1roy(x, y, theta1, theta2, theta3) result(lambda)
        integer, intent(in) :: x, y
        real(dp), intent(in) :: theta1, theta2, theta3

        lambda = 1.0_dp - theta1 * theta3**y + &
            0.0_dp * (real(x, dp) + theta2)
    end function lambda1roy

    pure real(dp) function lambda2roy(x, y, theta1, theta2, theta3) result(lambda)
        integer, intent(in) :: x, y
        real(dp), intent(in) :: theta1, theta2, theta3

        lambda = 1.0_dp - theta2 * theta3**x + &
            0.0_dp * (real(y, dp) + theta1)
    end function lambda2roy

    pure real(dp) function s_n(xq, yq, x, y) result(s)
        integer, intent(in) :: xq, yq
        integer, intent(in) :: x(:), y(:)

        s = empirical_survival_roy(xq, yq, x, y)
    end function s_n

end module bivgeom_compat
