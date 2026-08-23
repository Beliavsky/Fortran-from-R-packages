! SPDX-License-Identifier: GPL-3.0-only
module ecpdist_moments
    use ecpdist_kinds, only: dp
    use ecpdist_math, only: nan_dp
    use ecpdist_distribution, only: qecp, pecp, ecp_valid_parameters
    implicit none
    private

    type, public :: ecp_integral_result
        real(dp) :: estimate = 0.0_dp
        real(dp) :: abs_error = 0.0_dp
        integer :: status = 0
    end type ecp_integral_result

    public :: ecp_kmoment, ecp_kmoment_cond, ecp_mrl, ecp_shape

    real(dp), parameter :: xgk(8) = [ &
        0.991455371120812639206854697526329_dp, &
        0.949107912342758524526189684047851_dp, &
        0.864864423359769072789712788640926_dp, &
        0.741531185599394439863864773280788_dp, &
        0.586087235467691130294144838258730_dp, &
        0.405845151377397166906606412076961_dp, &
        0.207784955007898467600689403773245_dp, &
        0.0_dp ]
    real(dp), parameter :: wgk(8) = [ &
        0.022935322010529224963732008058970_dp, &
        0.063092092629978553290700663189204_dp, &
        0.104790010322250183839876322541518_dp, &
        0.140653259715525918745189590510238_dp, &
        0.169004726639267902826583426598550_dp, &
        0.190350578064785409913256402421014_dp, &
        0.204432940075298892414161999234649_dp, &
        0.209482141084727828012999174891714_dp ]
    real(dp), parameter :: wg(4) = [ &
        0.129484966168869693270611432679082_dp, &
        0.279705391489276667901467771423780_dp, &
        0.381830050505118944950369775488975_dp, &
        0.417959183673469387755102040816327_dp ]

contains

    function ecp_kmoment(k, lambda, gamma, phi, rel_tol, abs_tol) result(out)
        integer, intent(in) :: k
        real(dp), intent(in) :: lambda, gamma, phi
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(ecp_integral_result) :: out
        real(dp) :: rt, at

        rt = 1.0e-8_dp
        at = 1.0e-10_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        if (k < 1 .or. .not. ecp_valid_parameters(lambda, gamma, phi)) then
            out%estimate = nan_dp()
            out%abs_error = nan_dp()
            out%status = 1
            return
        end if
        call adaptive_quantile_integral(0.0_dp, 1.0_dp, k, lambda, gamma, phi, &
            rt, at, 0, out%estimate, out%abs_error, out%status)
        if (out%abs_error <= max(at, rt*abs(out%estimate))) out%status = 0
    end function ecp_kmoment

    function ecp_kmoment_cond(x, k, lambda, gamma, phi, rel_tol, abs_tol) result(out)
        real(dp), intent(in) :: x
        integer, intent(in) :: k
        real(dp), intent(in) :: lambda, gamma, phi
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(ecp_integral_result) :: out
        real(dp) :: rt, at, fx, sx, val, err
        integer :: istat

        rt = 1.0e-8_dp
        at = 1.0e-10_dp
        if (present(rel_tol)) rt = rel_tol
        if (present(abs_tol)) at = abs_tol
        if (x < 0.0_dp .or. k < 1 .or. .not. ecp_valid_parameters(lambda, gamma, phi)) then
            out%estimate = nan_dp()
            out%abs_error = nan_dp()
            out%status = 1
            return
        end if

        fx = pecp(x, lambda, gamma, phi)
        sx = pecp(x, lambda, gamma, phi, lower_tail=.false.)
        if (sx <= tiny(1.0_dp)) then
            out%estimate = nan_dp()
            out%abs_error = nan_dp()
            out%status = 2
            return
        end if
        call adaptive_quantile_integral(fx, 1.0_dp, k, lambda, gamma, phi, &
            rt, at*sx, 0, val, err, istat)
        out%estimate = val/sx
        out%abs_error = err/sx
        out%status = istat
        if (out%abs_error <= max(at, rt*abs(out%estimate))) out%status = 0
    end function ecp_kmoment_cond

    function ecp_mrl(x, lambda, gamma, phi, rel_tol, abs_tol) result(out)
        real(dp), intent(in) :: x, lambda, gamma, phi
        real(dp), intent(in), optional :: rel_tol, abs_tol
        type(ecp_integral_result) :: out

        out = ecp_kmoment_cond(x, 1, lambda, gamma, phi, rel_tol, abs_tol)
        if (out%status /= 1) out%estimate = out%estimate - x
    end function ecp_mrl

    real(dp) function ecp_shape(lambda, gamma, phi, measure) result(ans)
        real(dp), intent(in) :: lambda, gamma, phi
        character(len=*), intent(in) :: measure
        real(dp) :: q(7), p(7)
        integer :: i

        if (.not. ecp_valid_parameters(lambda, gamma, phi)) then
            ans = nan_dp()
            return
        end if
        do i = 1, 7
            p(i) = real(i, dp)/8.0_dp
            q(i) = qecp(p(i), lambda, gamma, phi)
        end do
        select case (trim(adjustl(measure)))
        case ('bowley', 'Bowley', 'BOWLEY')
            ans = (q(2) - 2.0_dp*q(4) + q(6))/(q(6) - q(2))
        case ('moors', 'Moors', 'MOORS')
            ans = (q(7) - q(5) - q(3) + q(1))/(q(6) - q(2))
        case default
            ans = nan_dp()
        end select
    end function ecp_shape

    recursive subroutine adaptive_quantile_integral(a, b, k, lambda, gamma, phi, &
            rel_tol, abs_tol, depth, value, error, status)
        real(dp), intent(in) :: a, b, lambda, gamma, phi, rel_tol, abs_tol
        integer, intent(in) :: k, depth
        real(dp), intent(out) :: value, error
        integer, intent(out) :: status
        real(dp) :: whole, err0, left, right, err_l, err_r, mid
        integer :: stat_l, stat_r

        call gk15_quantile(a, b, k, lambda, gamma, phi, whole, err0)
        if (err0 <= max(abs_tol, rel_tol*abs(whole)) .or. depth >= 20) then
            value = whole
            error = err0
            if (depth >= 20 .and. err0 > max(abs_tol, rel_tol*abs(whole))) then
                status = 2
            else
                status = 0
            end if
            return
        end if

        mid = 0.5_dp*(a + b)
        call adaptive_quantile_integral(a, mid, k, lambda, gamma, phi, rel_tol, &
            0.5_dp*abs_tol, depth + 1, left, err_l, stat_l)
        call adaptive_quantile_integral(mid, b, k, lambda, gamma, phi, rel_tol, &
            0.5_dp*abs_tol, depth + 1, right, err_r, stat_r)
        value = left + right
        error = err_l + err_r
        status = max(stat_l, stat_r)
    end subroutine adaptive_quantile_integral

    subroutine gk15_quantile(a, b, k, lambda, gamma, phi, result, error)
        real(dp), intent(in) :: a, b, lambda, gamma, phi
        integer, intent(in) :: k
        real(dp), intent(out) :: result, error
        real(dp) :: center, half, fc, f1, f2, resg, resk, absc
        integer :: j

        center = 0.5_dp*(a + b)
        half = 0.5_dp*(b - a)
        fc = qecp(center, lambda, gamma, phi)**k
        resg = wg(4)*fc
        resk = wgk(8)*fc

        do j = 1, 7
            absc = half*xgk(j)
            f1 = qecp(center - absc, lambda, gamma, phi)**k
            f2 = qecp(center + absc, lambda, gamma, phi)**k
            resk = resk + wgk(j)*(f1 + f2)
            select case (j)
            case (2)
                resg = resg + wg(1)*(f1 + f2)
            case (4)
                resg = resg + wg(2)*(f1 + f2)
            case (6)
                resg = resg + wg(3)*(f1 + f2)
            end select
        end do
        result = resk*half
        error = abs((resk - resg)*half)
    end subroutine gk15_quantile

end module ecpdist_moments
