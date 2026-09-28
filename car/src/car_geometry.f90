module car_geometry
    use car_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    implicit none
    private

    public :: ellipse_coordinates

contains

    pure subroutine ellipse_coordinates(center, shape, radius, segments, coordinates, status, log_x, log_y)
        real(dp), intent(in) :: center(2) !! Ellipse center in x-y coordinates, or log coordinates when the corresponding log flag is true.
        real(dp), intent(in) :: shape(2, 2) !! Symmetric positive-definite 2 by 2 shape matrix.
        real(dp), intent(in) :: radius !! Nonnegative ellipse radius multiplier.
        integer, intent(in) :: segments !! Number of polygon segments; output contains segments+1 points and must be at least three.
        real(dp), allocatable, intent(out) :: coordinates(:, :) !! Ellipse coordinates with shape (segments+1, 2).
        integer, intent(out) :: status !! Zero on success; nonzero for invalid radius, segment count, symmetry, or positive definiteness.
        logical, intent(in), optional :: log_x !! If true exponentiate the generated x coordinates, matching log="x".
        logical, intent(in), optional :: log_y !! If true exponentiate the generated y coordinates, matching log="y".
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: angle
        real(dp) :: det_shape
        real(dp) :: max_shape
        real(dp) :: u11
        real(dp) :: u12
        real(dp) :: u22
        logical :: exponentiate_x
        logical :: exponentiate_y
        integer :: i

        if (segments < 3 .or. radius < 0.0_dp) then
            allocate(coordinates(0, 2))
            status = 1
            return
        end if

        max_shape = maxval(abs(shape))
        if (max_shape <= 0.0_dp) then
            allocate(coordinates(0, 2))
            status = 2
            return
        end if
        if (maxval(abs(shape - transpose(shape))) / max_shape > 1.0e-10_dp) then
            allocate(coordinates(0, 2))
            status = 2
            return
        end if

        det_shape = shape(1, 1) * shape(2, 2) - shape(1, 2) * shape(2, 1)
        if (shape(1, 1) <= 0.0_dp .or. det_shape <= 0.0_dp) then
            allocate(coordinates(0, 2))
            status = 3
            return
        end if

        u11 = sqrt(shape(1, 1))
        u12 = shape(1, 2) / u11
        if (shape(2, 2) - u12 * u12 <= 0.0_dp) then
            allocate(coordinates(0, 2))
            status = 3
            return
        end if
        u22 = sqrt(shape(2, 2) - u12 * u12)

        exponentiate_x = .false.
        exponentiate_y = .false.
        if (present(log_x)) exponentiate_x = log_x
        if (present(log_y)) exponentiate_y = log_y

        allocate(coordinates(segments + 1, 2))
        do i = 0, segments
            angle = 2.0_dp * pi * real(i, dp) / real(segments, dp)
            coordinates(i + 1, 1) = center(1) + radius * (cos(angle) * u11)
            coordinates(i + 1, 2) = center(2) + radius * (cos(angle) * u12 + sin(angle) * u22)
        end do
        if (exponentiate_x) coordinates(:, 1) = exp(coordinates(:, 1))
        if (exponentiate_y) coordinates(:, 2) = exp(coordinates(:, 2))
        status = 0
    end subroutine ellipse_coordinates

end module car_geometry
