module actuar_special_v02
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use actuar_kinds, only: dp
    use actuar_special, only: reg_beta
    implicit none
    private
    public :: betaint_raw, choose_int
contains
    pure real(dp) function choose_int(n, k) result(v)
        integer, intent(in) :: n, k
        integer :: i, kk
        if (k < 0 .or. k > n) then
            v = 0.0_dp
            return
        end if
        kk = min(k, n-k)
        v = 1.0_dp
        do i = 1, kk
            v = v * real(n-kk+i, dp) / real(i, dp)
        end do
    end function choose_int

    pure real(dp) function betaint_raw(x, a, b) result(value)
        real(dp), intent(in) :: x, a, b
        real(dp) :: r, ap, bp, lx, lx1m, x1, c, tmp, sumv, ratio, ix
        integer :: i, ir
        if (x <= 0.0_dp .or. x >= 1.0_dp .or. a <= 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (b > 0.0_dp) then
            value = gamma(a) * gamma(b) * reg_beta(x, a, b)
            return
        end if
        r = floor(-b)
        ir = int(r)
        if (abs(b - real(nint(b),dp)) < 64.0_dp*epsilon(1.0_dp) .or. a-r-1.0_dp <= 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        ap = a
        bp = b
        lx = log(x)
        lx1m = log(1.0_dp-x)
        x1 = exp(lx1m-lx)
        ap = ap - 1.0_dp
        c = exp(ap*lx + bp*lx1m) / bp
        sumv = c
        ratio = 1.0_dp / bp
        bp = bp + 1.0_dp
        do i = 0, ir-1
            tmp = ap / bp
            c = tmp * (c*x1)
            sumv = sumv + c
            ratio = ratio * tmp
            ap = ap - 1.0_dp
            bp = bp + 1.0_dp
        end do
        ix = reg_beta(x, ap, bp)
        value = -gamma(a+b)*sumv
        if (ix > 0.0_dp) then
            value = value + (ratio*ap) * exp(log_gamma(ap)+log_gamma(bp)+log(ix))
        end if
    end function betaint_raw
end module actuar_special_v02
