! SPDX-License-Identifier: GPL-2.0-only
module lmoments_quantile_mixtures
    use lmoments_utils, only: dp, pi, standard_normal_quantile
    use lmoments_core, only: lmoments_sample, lmom_cov, t1_lmoments
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private

    public :: lmom2normpoly4, lmom2normpoly6
    public :: data2normpoly4, data2normpoly6
    public :: qnormpoly, pnormpoly, dnormpoly, rnormpoly
    public :: normpoly_inv, normpoly_cdf, normpoly_pdf, normpoly_rnd
    public :: covnormpoly4
    public :: t1lmom2cauchypoly4, data2cauchypoly4
    public :: qcauchypoly, pcauchypoly, dcauchypoly, rcauchypoly
    public :: cauchypoly_inv, cauchypoly_cdf, cauchypoly_pdf, cauchypoly_rnd

contains

    pure function lmom2normpoly4(lmom) result(param)
        real(dp), intent(in) :: lmom(:)
        real(dp) :: param(4)

        if (size(lmom) < 4) then
            param = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        param(1) = lmom(1) - 3.0_dp * lmom(2) + 5.0_dp * lmom(3) &
            + 24.46947735493778_dp * lmom(4)
        param(2) = 6.0_dp * lmom(2) - 30.0_dp * lmom(3) &
            - 48.93895470987556_dp * lmom(4)
        param(3) = 30.0_dp * lmom(3)
        param(4) = 14.457006455887834_dp * lmom(4)
    end function lmom2normpoly4


    pure function lmom2normpoly6(lmom) result(param)
        real(dp), intent(in) :: lmom(:)
        real(dp) :: param(6)

        if (size(lmom) < 6) then
            param = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        param(1) = lmom(1) - 3.0_dp * lmom(2) + 5.0_dp * lmom(3) - 7.0_dp * lmom(4) &
            + 9.0_dp * lmom(5) + 88.36715690214288_dp * lmom(6)
        param(2) = 6.0_dp * lmom(2) - 30.0_dp * lmom(3) + 84.0_dp * lmom(4) &
            - 180.0_dp * lmom(5) - 373.2962367553012_dp * lmom(6)
        param(3) = 30.0_dp * lmom(3) - 210.0_dp * lmom(4) + 810.0_dp * lmom(5) &
            + 589.6857688530462_dp * lmom(6)
        param(4) = 140.0_dp * lmom(4) - 1260.0_dp * lmom(5) &
            - 393.1238459020308_dp * lmom(6)
        param(5) = 630.0_dp * lmom(5)
        param(6) = 40.595671272636515_dp * lmom(6)
    end function lmom2normpoly6


    function data2normpoly4(data) result(param)
        real(dp), intent(in) :: data(:)
        real(dp) :: param(4), lmom(4)
        integer :: info

        call lmoments_sample(data, lmom, info)
        if (info /= 0) then
            param = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            param = lmom2normpoly4(lmom)
        end if
    end function data2normpoly4


    function data2normpoly6(data) result(param)
        real(dp), intent(in) :: data(:)
        real(dp) :: param(6), lmom(6)
        integer :: info

        call lmoments_sample(data, lmom, info)
        if (info /= 0) then
            param = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            param = lmom2normpoly6(lmom)
        end if
    end function data2normpoly6


    pure real(dp) function qnormpoly(cp, param) result(x)
        real(dp), intent(in) :: cp
        real(dp), intent(in) :: param(:)
        integer :: k, poly_order

        if (size(param) < 2 .or. cp < 0.0_dp .or. cp > 1.0_dp) then
            x = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        poly_order = size(param) - 1
        x = 0.0_dp
        do k = 1, poly_order
            x = x + param(k) * cp ** (k - 1)
        end do
        x = x + param(poly_order + 1) * standard_normal_quantile(cp)
    end function qnormpoly


    pure real(dp) function pnormpoly(x, param) result(cp)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: param(:)
        integer :: j
        real(dp) :: diff

        if (size(param) < 2) then
            cp = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        cp = 0.5_dp
        do j = -2, -50, -1
            diff = x - qnormpoly(cp, param)
            if (diff > 0.0_dp) then
                cp = cp + 2.0_dp ** j
            else if (diff < 0.0_dp) then
                cp = cp - 2.0_dp ** j
            end if
        end do
    end function pnormpoly


    pure real(dp) function dnormpoly(x, param) result(density)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: param(:)
        integer :: k, poly_order
        real(dp) :: cp, denom, z

        if (size(param) < 2) then
            density = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        poly_order = size(param) - 1
        cp = pnormpoly(x, param)
        denom = 0.0_dp
        do k = 2, poly_order
            denom = denom + param(k) * real(k - 1, dp) * cp ** (k - 2)
        end do
        z = standard_normal_quantile(cp)
        denom = denom + param(poly_order + 1) * sqrt(2.0_dp * pi) * exp(0.5_dp * z * z)
        density = 1.0_dp / denom
        if (abs(param(poly_order + 1)) <= tiny(1.0_dp)) then
            if (x > sum(param) .or. x < param(1)) density = 0.0_dp
        end if
    end function dnormpoly


    subroutine rnormpoly(param, x)
        real(dp), intent(in) :: param(:)
        real(dp), intent(out) :: x(:)
        real(dp), allocatable :: u(:)
        integer :: i

        allocate(u(size(x)))
        call random_number(u)
        do i = 1, size(x)
            x(i) = qnormpoly(u(i), param)
        end do
    end subroutine rnormpoly


    pure real(dp) function normpoly_inv(cp, param) result(x)
        real(dp), intent(in) :: cp
        real(dp), intent(in) :: param(:)
        x = qnormpoly(cp, param)
    end function normpoly_inv


    pure real(dp) function normpoly_cdf(x, param) result(cp)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: param(:)
        cp = pnormpoly(x, param)
    end function normpoly_cdf


    pure real(dp) function normpoly_pdf(x, param) result(density)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: param(:)
        density = dnormpoly(x, param)
    end function normpoly_pdf


    subroutine normpoly_rnd(param, x)
        real(dp), intent(in) :: param(:)
        real(dp), intent(out) :: x(:)
        call rnormpoly(param, x)
    end subroutine normpoly_rnd


    function covnormpoly4(data) result(covparam)
        real(dp), intent(in) :: data(:)
        real(dp) :: covparam(4, 4), covlmom(4, 4), a(4, 4)
        real(dp) :: eta2, eta4
        integer :: info

        call lmom_cov(data, covlmom, info)
        if (info /= 0) then
            covparam = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        eta2 = 1.0_dp / sqrt(pi)
        eta4 = (30.0_dp * atan(sqrt(2.0_dp)) / pi - 9.0_dp) * eta2
        a = 0.0_dp
        a(1, :) = [1.0_dp, -3.0_dp, 5.0_dp, 3.0_dp * eta2 / eta4]
        a(2, :) = [0.0_dp, 6.0_dp, -30.0_dp, -6.0_dp * eta2 / eta4]
        a(3, :) = [0.0_dp, 0.0_dp, 30.0_dp, 0.0_dp]
        a(4, :) = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp / eta4]
        covparam = matmul(a, matmul(covlmom, transpose(a)))
    end function covnormpoly4


    pure function t1lmom2cauchypoly4(t1lmom) result(param)
        real(dp), intent(in) :: t1lmom(:)
        real(dp) :: param(4)
        real(dp), parameter :: v1 = 0.6978_dp
        real(dp), parameter :: v2 = 0.3428_dp * v1

        if (size(t1lmom) < 4) then
            param = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        param(1) = (25.0_dp * t1lmom(4) * v1 + 5.0_dp * t1lmom(1) * v2 &
            - 25.0_dp * t1lmom(2) * v2 + 63.0_dp * t1lmom(3) * v2) / (5.0_dp * v2)
        param(2) = 10.0_dp * t1lmom(2) - 63.0_dp * t1lmom(3) &
            - 10.0_dp * t1lmom(4) * v1 / v2
        param(3) = 63.0_dp * t1lmom(3)
        param(4) = t1lmom(4) / v2
    end function t1lmom2cauchypoly4


    function data2cauchypoly4(data) result(param)
        real(dp), intent(in) :: data(:)
        real(dp) :: param(4), tlmom(4)
        integer :: info

        call t1_lmoments(data, tlmom, info)
        if (info /= 0) then
            param = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            param = t1lmom2cauchypoly4(tlmom)
        end if
    end function data2cauchypoly4


    pure real(dp) function qcauchypoly(cp, param) result(x)
        real(dp), intent(in) :: cp
        real(dp), intent(in) :: param(:)
        integer :: k, poly_order

        if (size(param) < 2 .or. cp < 0.0_dp .or. cp > 1.0_dp) then
            x = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        poly_order = size(param) - 1
        x = 0.0_dp
        do k = 1, poly_order
            x = x + param(k) * cp ** (k - 1)
        end do
        x = x + param(poly_order + 1) * tan(pi * (cp - 0.5_dp))
    end function qcauchypoly


    pure real(dp) function pcauchypoly(x, param) result(cp)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: param(:)
        integer :: j
        real(dp) :: diff

        if (size(param) < 2) then
            cp = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        cp = 0.5_dp
        do j = -2, -50, -1
            diff = x - qcauchypoly(cp, param)
            if (diff > 0.0_dp) then
                cp = cp + 2.0_dp ** j
            else if (diff < 0.0_dp) then
                cp = cp - 2.0_dp ** j
            end if
        end do
    end function pcauchypoly


    pure real(dp) function dcauchypoly(x, param) result(density)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: param(:)
        integer :: k, poly_order
        real(dp) :: cp, denom

        if (size(param) < 2) then
            density = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        poly_order = size(param) - 1
        cp = pcauchypoly(x, param)
        denom = 0.0_dp
        do k = 2, poly_order
            denom = denom + param(k) * real(k - 1, dp) * cp ** (k - 2)
        end do
        denom = denom + param(poly_order + 1) * pi &
            / cos(pi * (cp - 0.5_dp)) ** 2
        density = 1.0_dp / denom
        if (abs(param(poly_order + 1)) <= tiny(1.0_dp)) then
            if (x > sum(param) .or. x < param(1)) density = 0.0_dp
        end if
    end function dcauchypoly


    subroutine rcauchypoly(param, x)
        real(dp), intent(in) :: param(:)
        real(dp), intent(out) :: x(:)
        real(dp), allocatable :: u(:)
        integer :: i

        allocate(u(size(x)))
        call random_number(u)
        do i = 1, size(x)
            x(i) = qcauchypoly(u(i), param)
        end do
    end subroutine rcauchypoly


    pure real(dp) function cauchypoly_inv(cp, param) result(x)
        real(dp), intent(in) :: cp
        real(dp), intent(in) :: param(:)
        x = qcauchypoly(cp, param)
    end function cauchypoly_inv


    pure real(dp) function cauchypoly_cdf(x, param) result(cp)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: param(:)
        cp = pcauchypoly(x, param)
    end function cauchypoly_cdf


    pure real(dp) function cauchypoly_pdf(x, param) result(density)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: param(:)
        density = dcauchypoly(x, param)
    end function cauchypoly_pdf


    subroutine cauchypoly_rnd(param, x)
        real(dp), intent(in) :: param(:)
        real(dp), intent(out) :: x(:)
        call rcauchypoly(param, x)
    end subroutine cauchypoly_rnd

end module lmoments_quantile_mixtures
