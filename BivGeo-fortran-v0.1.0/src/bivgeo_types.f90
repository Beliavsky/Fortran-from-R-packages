module bivgeo_types
    use bivgeo_kinds, only : dp
    implicit none
    private

    type, public :: bivgeo_params
        real(dp) :: theta1 = 0.5_dp
        real(dp) :: theta2 = 0.5_dp
        real(dp) :: theta3 = 1.0_dp
    end type bivgeo_params

    public :: make_bivgeo_params
    public :: valid_bivgeo_params

contains

    pure function make_bivgeo_params(theta1, theta2, theta3) result(theta)
        real(dp), intent(in) :: theta1, theta2, theta3
        type(bivgeo_params) :: theta

        theta%theta1 = theta1
        theta%theta2 = theta2
        theta%theta3 = theta3
    end function make_bivgeo_params

    pure logical function valid_bivgeo_params(theta) result(ok)
        type(bivgeo_params), intent(in) :: theta

        ok = theta%theta1 > 0.0_dp .and. theta%theta1 < 1.0_dp
        ok = ok .and. theta%theta2 > 0.0_dp .and. theta%theta2 < 1.0_dp
        ok = ok .and. theta%theta3 > 0.0_dp .and. theta%theta3 <= 1.0_dp
    end function valid_bivgeo_params

end module bivgeo_types
