! SPDX-License-Identifier: GPL-2.0-only
module lmoments_core
    use lmoments_utils, only: dp, sort_real, choose_real
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private

    public :: lmoments_sample
    public :: lmoments_matrix
    public :: lcoefs_sample
    public :: lmom_cov
    public :: t1_lmoments
    public :: shifted_legendre
    public :: hosking_lmoments

contains

    subroutine lmoments_sample(data, lmom, info)
        real(dp), intent(in) :: data(:)
        real(dp), intent(out) :: lmom(:)
        integer, intent(out), optional :: info
        real(dp), allocatable :: x(:), beta(:), w(:)
        integer :: n, rmax, r, k
        real(dp) :: coef

        call set_info(info, 0)
        n = size(data)
        rmax = size(lmom)
        if (rmax < 1 .or. rmax > n) then
            lmom = ieee_value(0.0_dp, ieee_quiet_nan)
            call set_info(info, 1)
            return
        end if

        allocate(x(n), beta(rmax), w(n))
        x = data
        call sort_real(x)
        beta = 0.0_dp
        lmom = 0.0_dp

        beta(1) = sum(x) / real(n, dp)
        if (rmax >= 2) then
            do k = 1, n
                w(k) = real(k - 1, dp) / real(n - 1, dp)
            end do
            beta(2) = sum(w * x) / real(n, dp)
        end if

        do r = 2, rmax - 1
            do k = 1, n
                w(k) = w(k) * real(k - r, dp) / real(n - r, dp)
            end do
            beta(r + 1) = sum(w * x) / real(n, dp)
        end do

        lmom(1) = beta(1)
        do r = 1, rmax - 1
            do k = 0, r
                coef = (-1.0_dp) ** (r - k)
                coef = coef * gamma(real(r + k + 1, dp))
                coef = coef / gamma(real(k + 1, dp)) ** 2
                coef = coef / gamma(real(r - k + 1, dp))
                lmom(r + 1) = lmom(r + 1) + coef * beta(k + 1)
            end do
        end do
    end subroutine lmoments_sample


    subroutine lmoments_matrix(data, lmom, info)
        real(dp), intent(in) :: data(:, :)
        real(dp), intent(out) :: lmom(:, :)
        integer, intent(out), optional :: info
        integer :: j, stat

        call set_info(info, 0)
        if (size(lmom, 1) /= size(data, 2)) then
            lmom = ieee_value(0.0_dp, ieee_quiet_nan)
            call set_info(info, 1)
            return
        end if
        do j = 1, size(data, 2)
            call lmoments_sample(data(:, j), lmom(j, :), stat)
            if (stat /= 0) then
                call set_info(info, stat)
                return
            end if
        end do
    end subroutine lmoments_matrix


    subroutine lcoefs_sample(data, coefs, info)
        real(dp), intent(in) :: data(:)
        real(dp), intent(out) :: coefs(:)
        integer, intent(out), optional :: info
        integer :: r, stat

        call lmoments_sample(data, coefs, stat)
        call set_info(info, stat)
        if (stat /= 0) return
        if (size(coefs) >= 3) then
            if (abs(coefs(2)) <= tiny(1.0_dp)) then
                coefs(3:) = ieee_value(0.0_dp, ieee_quiet_nan)
                call set_info(info, 2)
                return
            end if
            do r = 3, size(coefs)
                coefs(r) = coefs(r) / coefs(2)
            end do
        end if
    end subroutine lcoefs_sample


    subroutine shifted_legendre(rmax, c, info)
        integer, intent(in) :: rmax
        real(dp), intent(out) :: c(:, :)
        integer, intent(out), optional :: info
        real(dp), allocatable :: work(:, :), same(:), up(:)
        integer :: k, kn

        call set_info(info, 0)
        if (rmax <= 0 .or. size(c, 1) /= rmax .or. size(c, 2) /= rmax) then
            c = ieee_value(0.0_dp, ieee_quiet_nan)
            call set_info(info, 1)
            return
        end if

        allocate(work(rmax, rmax))
        work = 0.0_dp
        work(1, 1) = 1.0_dp
        if (rmax > 1) then
            work(1, 2) = -1.0_dp
            work(2, 2) = 2.0_dp
            do k = 3, rmax
                kn = k - 2
                allocate(same(k - 1), up(k - 1))
                same = (-(2.0_dp * real(kn, dp) + 1.0_dp) * work(1:k - 1, k - 1) &
                    - real(kn, dp) * work(1:k - 1, k - 2)) / real(kn + 1, dp)
                up = 2.0_dp * (2.0_dp * real(kn, dp) + 1.0_dp) &
                    / real(kn + 1, dp) * work(1:k - 1, k - 1)
                work(1:k - 1, k) = same
                work(2:k, k) = work(2:k, k) + up
                deallocate(same, up)
            end do
        end if
        c = transpose(work)
    end subroutine shifted_legendre


    subroutine lmom_cov(data, cov, info)
        real(dp), intent(in) :: data(:)
        real(dp), intent(out) :: cov(:, :)
        integer, intent(out), optional :: info
        real(dp), allocatable :: x(:), c(:, :), bcoef(:, :), b1(:, :)
        real(dp), allocatable :: b20(:, :), b2(:, :, :), beta(:), theta(:, :)
        integer :: n, rmax, r, k, l, i, j, stat
        real(dp) :: joint, term1, term2

        call set_info(info, 0)
        n = size(data)
        if (size(cov, 1) /= size(cov, 2)) then
            cov = ieee_value(0.0_dp, ieee_quiet_nan)
            call set_info(info, 1)
            return
        end if
        rmax = size(cov, 1)
        if (rmax <= 1 .or. rmax > n / 2) then
            cov = ieee_value(0.0_dp, ieee_quiet_nan)
            call set_info(info, 2)
            return
        end if

        allocate(x(n), c(rmax, rmax), bcoef(n, rmax), b1(n, rmax))
        allocate(b20(n, rmax), b2(n, rmax, rmax), beta(rmax), theta(rmax, rmax))
        x = data
        call sort_real(x)
        call shifted_legendre(rmax, c, stat)
        if (stat /= 0) then
            cov = ieee_value(0.0_dp, ieee_quiet_nan)
            call set_info(info, stat)
            return
        end if

        bcoef = 0.0_dp
        b1 = 0.0_dp
        b20 = 0.0_dp
        b2 = 0.0_dp
        beta = 0.0_dp
        theta = 0.0_dp

        do i = 1, n
            bcoef(i, 1) = real(i - 1, dp) / real(n - 1, dp)
        end do
        beta(1) = sum(x) / real(n, dp)
        beta(2) = sum(bcoef(:, 1) * x) / real(n, dp)

        do i = 1, n
            b20(i, 1) = -1.0_dp / real(n - 2, dp) &
                + real(i - 1, dp) * (1.0_dp + 1.0_dp / real(n - 2, dp)) / real(n - 1, dp)
            b1(i, 1) = real(i - 1, dp) / real(n - 2, dp)
        end do

        do r = 2, rmax - 1
            do i = 1, n
                bcoef(i, r) = bcoef(i, r - 1) * (-(real(r - 1, dp)) / real(n - r, dp) &
                    + real(i - 1, dp) * (1.0_dp + real(r - 1, dp) / real(n - r, dp)) &
                    / real(n - 1, dp))
            end do
            beta(r + 1) = sum(bcoef(:, r) * x) / real(n, dp)
        end do

        do k = 1, rmax - 1
            if (k > 1) then
                do i = 1, n
                    b20(i, k) = b20(i, k - 1) * (-real(k, dp) / real(n - 1 - k, dp) &
                        + real(i - 1, dp) * (1.0_dp + real(k, dp) / real(n - 1 - k, dp)) &
                        / real(n - 1, dp))
                    b1(i, k) = b1(i, k - 1) * (real(-k + 1, dp) / real(n - 1 - k, dp) &
                        + real(i - 1, dp) * real(n - 1, dp) / real(n - 1 - k, dp) &
                        / real(n - 1, dp))
                end do
            end if
            do i = 1, n
                b2(i, k, 1) = real(-k - 1, dp) / real(n - k - 2, dp) &
                    + real(i - 1, dp) * (1.0_dp - real(-k - 1, dp) / real(n - k - 2, dp)) &
                    / real(n - 1, dp)
            end do
            do l = 2, rmax - 1
                do i = 1, n
                    b2(i, k, l) = b2(i, k, l - 1) * (real(-k - l, dp) / real(n - k - l - 1, dp) &
                        + real(i - 1, dp) * (1.0_dp - real(-k - l, dp) / real(n - k - l - 1, dp)) &
                        / real(n - 1, dp))
                end do
            end do
        end do

        do k = 0, rmax - 1
            do l = 0, rmax - 1
                joint = 0.0_dp
                do i = 1, n - 1
                    do j = i + 1, n
                        if (k > 0 .and. l > 0) then
                            term1 = b1(i, k) * b2(j, k, l)
                            term2 = b1(i, l) * b2(j, l, k)
                        else if (k == 0 .and. l > 0) then
                            term1 = b20(j, l)
                            term2 = b1(i, l)
                        else if (k > 0 .and. l == 0) then
                            term1 = b1(i, k)
                            term2 = b20(j, k)
                        else
                            term1 = 1.0_dp
                            term2 = 1.0_dp
                        end if
                        joint = joint + (term1 + term2) * x(i) * x(j)
                    end do
                end do
                joint = joint / real(n * (n - 1), dp)
                theta(k + 1, l + 1) = beta(k + 1) * beta(l + 1) - joint
            end do
        end do

        cov = matmul(c, matmul(theta, transpose(c)))
    end subroutine lmom_cov


    subroutine t1_lmoments(data, tlmom, info)
        real(dp), intent(in) :: data(:)
        real(dp), intent(out) :: tlmom(:)
        integer, intent(out), optional :: info
        real(dp), allocatable :: x(:)
        integer :: n, rmax, r, k, i, m, j
        real(dp) :: est, weight, denom, signcoef

        call set_info(info, 0)
        n = size(data)
        rmax = size(tlmom)
        if (rmax < 1 .or. rmax > 4 .or. n < rmax + 2) then
            tlmom = ieee_value(0.0_dp, ieee_quiet_nan)
            call set_info(info, 1)
            return
        end if

        allocate(x(n))
        x = data
        call sort_real(x)
        tlmom = 0.0_dp

        do r = 1, rmax
            m = r + 2
            denom = choose_real(n, m)
            do k = 0, r - 1
                j = r + 1 - k
                est = 0.0_dp
                do i = 1, n
                    weight = choose_real(i - 1, j - 1) * choose_real(n - i, m - j) / denom
                    est = est + weight * x(i)
                end do
                signcoef = (-1.0_dp) ** k * choose_real(r - 1, k)
                tlmom(r) = tlmom(r) + signcoef * est / real(r, dp)
            end do
        end do
    end subroutine t1_lmoments


    ! Modernized reference path based on Hosking's SAMLMR unbiased PWM logic.
    ! Returns raw L-moments lambda_r, matching the R package, not T3/T4 ratios.
    subroutine hosking_lmoments(data, lmom, info)
        real(dp), intent(in) :: data(:)
        real(dp), intent(out) :: lmom(:)
        integer, intent(out), optional :: info
        real(dp), allocatable :: x(:), beta(:)
        integer :: n, rmax, i, j, k
        real(dp) :: z, term, y, p0, p, temp, ak, ai

        call set_info(info, 0)
        n = size(data)
        rmax = size(lmom)
        if (rmax < 1 .or. rmax > n) then
            lmom = ieee_value(0.0_dp, ieee_quiet_nan)
            call set_info(info, 1)
            return
        end if

        allocate(x(n), beta(rmax))
        x = data
        call sort_real(x)
        beta = 0.0_dp

        do i = 1, n
            z = real(i, dp)
            term = x(i)
            beta(1) = beta(1) + term
            do j = 2, rmax
                z = z - 1.0_dp
                term = term * z
                beta(j) = beta(j) + term
            end do
        end do
        z = real(n, dp)
        beta(1) = beta(1) / z
        y = real(n, dp)
        do j = 2, rmax
            y = y - 1.0_dp
            z = z * y
            beta(j) = beta(j) / z
        end do

        k = rmax
        p0 = 1.0_dp
        if (mod(rmax, 2) == 1) p0 = -1.0_dp
        do j = 2, rmax
            ak = real(k, dp)
            p0 = -p0
            p = p0
            temp = p * beta(1)
            do i = 1, k - 1
                ai = real(i, dp)
                p = -p * (ak + ai - 1.0_dp) * (ak - ai) / (ai * ai)
                temp = temp + p * beta(i + 1)
            end do
            beta(k) = temp
            k = k - 1
        end do
        lmom = beta
    end subroutine hosking_lmoments


    subroutine set_info(info, value)
        integer, intent(out), optional :: info
        integer, intent(in) :: value
        if (present(info)) info = value
    end subroutine set_info

end module lmoments_core
