module bivgeo_random
    use bivgeo_kinds, only : dp
    use bivgeo_types, only : bivgeo_params, valid_bivgeo_params
    implicit none
    private

    public :: bivgeo_seed
    public :: rbivgeo1
    public :: rbivgeo2

contains

    subroutine bivgeo_seed(seed)
        integer, intent(in) :: seed
        integer :: i, n
        integer, allocatable :: put(:)

        call random_seed(size=n)
        allocate(put(n))
        do i = 1, n
            put(i) = modulo(abs(seed) + 104729 * i, 2147483646)
            if (put(i) == 0) put(i) = i + 1
        end do
        call random_seed(put=put)
    end subroutine bivgeo_seed

    integer function geometric_one(q) result(k)
        real(dp), intent(in) :: q
        real(dp) :: u

        if (q <= 0.0_dp) then
            k = 1
            return
        end if
        if (q >= 1.0_dp) then
            k = huge(1)
            return
        end if

        call random_number(u)
        u = min(u, 1.0_dp - epsilon(1.0_dp))
        k = 1 + int(log(1.0_dp - u) / log(q))
        k = max(1, k)
    end function geometric_one

    subroutine rbivgeo2(n, theta, sample, ok)
        integer, intent(in) :: n
        type(bivgeo_params), intent(in) :: theta
        integer, intent(out) :: sample(:, :)
        logical, intent(out), optional :: ok
        integer :: i, r1, r2, r3
        logical :: good

        good = valid_bivgeo_params(theta)
        good = good .and. n >= 0
        good = good .and. size(sample, 1) >= n .and. size(sample, 2) >= 2
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if

        do i = 1, n
            r1 = geometric_one(theta%theta1)
            r2 = geometric_one(theta%theta2)
            if (theta%theta3 < 1.0_dp) then
                r3 = geometric_one(theta%theta3)
                sample(i, 1) = min(r1, r3)
                sample(i, 2) = min(r2, r3)
            else
                sample(i, 1) = r1
                sample(i, 2) = r2
            end if
        end do
        if (present(ok)) ok = .true.
    end subroutine rbivgeo2

    subroutine rbivgeo1(n, theta, sample, ok)
        integer, intent(in) :: n
        type(bivgeo_params), intent(in) :: theta
        integer, intent(out) :: sample(:, :)
        logical, intent(out), optional :: ok
        real(dp) :: u, cdf_before, diag_mass, cdf_diag
        real(dp) :: gamma1, gamma2, gamma3
        integer :: i, x, y
        logical :: good

        good = valid_bivgeo_params(theta)
        good = good .and. n >= 0
        good = good .and. size(sample, 1) >= n .and. size(sample, 2) >= 2
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if

        gamma1 = theta%theta1 * theta%theta3
        gamma2 = theta%theta2 * theta%theta3
        gamma3 = theta%theta1 * theta%theta2 * theta%theta3

        do i = 1, n
            x = geometric_one(gamma1)
            call random_number(u)

            if (x > 1) then
                cdf_before = 1.0_dp - theta%theta2**(x - 1)
            else
                cdf_before = 0.0_dp
            end if

            diag_mass = theta%theta2**(x - 1)
            diag_mass = diag_mass * (1.0_dp - gamma1 - gamma2 + gamma3)
            diag_mass = diag_mass / (1.0_dp - gamma1)
            cdf_diag = min(1.0_dp, cdf_before + diag_mass)

            if (u < cdf_before) then
                y = 1 + int(log(1.0_dp - u) / log(theta%theta2))
                y = min(y, x - 1)
                y = max(1, y)
            else if (u < cdf_diag .or. cdf_diag >= 1.0_dp) then
                y = x
            else
                y = x + geometric_one(gamma2)
            end if

            sample(i, 1) = x
            sample(i, 2) = y
        end do
        if (present(ok)) ok = .true.
    end subroutine rbivgeo1

end module bivgeo_random
