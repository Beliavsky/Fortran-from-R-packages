! SPDX-License-Identifier: GPL-2.0-only
!
! IIR filter design translated from CRAN signal / Octave-Forge algorithms.
module signal_iir_design
    use signal_kinds, only : dp, signal_pi
    use signal_types, only : arma_filter, zpg_filter, filter_order
    use signal_filters, only : arma_from_zpg, make_filter_order, make_zpg
    implicit none
    private

    public :: bilinear_transform
    public :: splane_frequency_transform
    public :: butter_filter
    public :: butter_order
    public :: cheby1_filter
    public :: cheby2_filter
    public :: cheby1_order
    public :: ellip_filter
    public :: ellip_order

contains

    pure function bilinear_transform(s_filter, sample_period) result(z_filter)
        type(zpg_filter), intent(in) :: s_filter !! Continuous-time zero-pole-gain filter in the s plane.
        real(dp), intent(in) :: sample_period !! Sampling period T used by the bilinear transform.
        type(zpg_filter) :: z_filter
        complex(dp), allocatable :: num_terms(:)
        complex(dp), allocatable :: den_terms(:)
        complex(dp), allocatable :: zeros(:)
        integer :: p
        integer :: z

        p = size(s_filter%pole)
        z = size(s_filter%zero)
        if (p == 0 .or. z > p) then
            allocate(z_filter%zero(0), z_filter%pole(0))
            z_filter%gain = 0.0_dp
            return
        end if
        allocate(num_terms(z), den_terms(p))
        if (z > 0) num_terms = (2.0_dp - s_filter%zero * sample_period) / sample_period
        den_terms = (2.0_dp - s_filter%pole * sample_period) / sample_period
        z_filter%gain = real(s_filter%gain * complex_product(num_terms) / complex_product(den_terms), dp)
        z_filter%pole = (2.0_dp + s_filter%pole * sample_period) / &
            (2.0_dp - s_filter%pole * sample_period)
        allocate(zeros(p))
        zeros = cmplx(-1.0_dp, 0.0_dp, dp)
        if (z > 0) zeros(1:z) = (2.0_dp + s_filter%zero * sample_period) / &
            (2.0_dp - s_filter%zero * sample_period)
        z_filter%zero = zeros
    end function bilinear_transform

    pure function splane_frequency_transform(s_filter, w, stop) result(out_filter)
        type(zpg_filter), intent(in) :: s_filter !! Prototype low-pass zero-pole-gain filter in the s plane.
        real(dp), intent(in) :: w(:) !! One cutoff frequency or two band-edge frequencies in the s plane.
        logical, intent(in), optional :: stop !! True requests high-pass or band-stop transformation.
        type(zpg_filter) :: out_filter
        complex(dp), allocatable :: pnew(:)
        complex(dp), allocatable :: znew(:)
        complex(dp), allocatable :: b(:)
        complex(dp) :: root_term
        complex(dp) :: extension(2)
        real(dp) :: gain
        real(dp) :: fl
        real(dp) :: fh
        real(dp) :: fc
        logical :: is_stop
        integer :: p
        integer :: z
        integer :: i
        integer :: k

        p = size(s_filter%pole)
        z = size(s_filter%zero)
        is_stop = .false.
        if (present(stop)) is_stop = stop
        if (p == 0 .or. z > p .or. (size(w) /= 1 .and. size(w) /= 2)) then
            allocate(out_filter%zero(0), out_filter%pole(0))
            out_filter%gain = 0.0_dp
            return
        end if
        gain = s_filter%gain
        if (size(w) == 2) then
            fl = w(1)
            fh = w(2)
            if (is_stop) then
                gain = gain * real(complex_product(-s_filter%zero) / complex_product(-s_filter%pole), dp)
                allocate(b(p), pnew(2 * p))
                b = ((fh - fl) / 2.0_dp) / s_filter%pole
                do i = 1, p
                    root_term = sqrt(b(i) * b(i) - fh * fl)
                    pnew(i) = b(i) + root_term
                    pnew(p + i) = b(i) - root_term
                end do
                extension(1) = sqrt(cmplx(-fh * fl, 0.0_dp, dp))
                extension(2) = -extension(1)
                if (z == 0) then
                    allocate(znew(2 * p))
                    do i = 1, p
                        znew(2 * i - 1:2 * i) = extension
                    end do
                else
                    allocate(b(z), znew(2 * z + 2 * (p - z)))
                    b = ((fh - fl) / 2.0_dp) / s_filter%zero
                    do i = 1, z
                        root_term = sqrt(b(i) * b(i) - fh * fl)
                        znew(i) = b(i) + root_term
                        znew(z + i) = b(i) - root_term
                    end do
                    k = 2 * z
                    do i = 1, p - z
                        znew(k + 2 * i - 1:k + 2 * i) = extension
                    end do
                end if
            else
                gain = gain * (1.0_dp / (fh - fl)) ** (z - p)
                allocate(b(p), pnew(2 * p))
                b = s_filter%pole * (fh - fl) / 2.0_dp
                do i = 1, p
                    root_term = sqrt(b(i) * b(i) - fh * fl)
                    pnew(i) = b(i) + root_term
                    pnew(p + i) = b(i) - root_term
                end do
                if (z == 0) then
                    allocate(znew(p))
                    znew = cmplx(0.0_dp, 0.0_dp, dp)
                else
                    allocate(b(z), znew(2 * z + p - z))
                    b = s_filter%zero * (fh - fl) / 2.0_dp
                    do i = 1, z
                        root_term = sqrt(b(i) * b(i) - fh * fl)
                        znew(i) = b(i) + root_term
                        znew(z + i) = b(i) - root_term
                    end do
                    if (p > z) znew(2 * z + 1:) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if
        else
            fc = w(1)
            allocate(pnew(p))
            if (is_stop) then
                gain = gain * real(complex_product(-s_filter%zero) / complex_product(-s_filter%pole), dp)
                pnew = fc / s_filter%pole
                allocate(znew(p))
                znew = cmplx(0.0_dp, 0.0_dp, dp)
                if (z > 0) znew(1:z) = fc / s_filter%zero
            else
                gain = gain * (1.0_dp / fc) ** (z - p)
                pnew = fc * s_filter%pole
                allocate(znew(z))
                if (z > 0) znew = fc * s_filter%zero
            end if
        end if
        out_filter%pole = pnew
        out_filter%zero = znew
        out_filter%gain = gain
    end function splane_frequency_transform

    pure function butter_filter(n, w, filter_type, digital) result(filt)
        integer, intent(in) :: n !! Positive Butterworth prototype order.
        real(dp), intent(in) :: w(:) !! One or two critical frequencies; digital values are normalised to Nyquist.
        character(len=*), intent(in), optional :: filter_type !! low, high, stop, or pass; defaults to low.
        logical, intent(in), optional :: digital !! True for a z-plane design; defaults to true.
        type(arma_filter) :: filt
        type(zpg_filter) :: proto
        type(zpg_filter) :: transformed
        real(dp), allocatable :: work_w(:)
        character(len=8) :: kind
        logical :: is_digital
        logical :: is_stop
        integer :: i

        kind = 'low'
        if (present(filter_type)) kind = adjustl(filter_type)
        is_digital = .true.
        if (present(digital)) is_digital = digital
        is_stop = kind == 'stop' .or. kind == 'high'
        allocate(proto%pole(max(0, n)), proto%zero(0))
        do i = 1, n
            proto%pole(i) = exp(cmplx(0.0_dp, signal_pi * real(2 * i + n - 1, dp) / real(2 * n, dp), dp))
        end do
        if (modulo(n, 2) == 1) proto%pole((n + 1) / 2) = cmplx(-1.0_dp, 0.0_dp, dp)
        proto%gain = 1.0_dp
        work_w = w
        if (is_digital) work_w = tan(signal_pi * work_w / 2.0_dp)
        transformed = splane_frequency_transform(proto, work_w, is_stop)
        if (is_digital) transformed = bilinear_transform(transformed, 2.0_dp)
        filt = arma_from_zpg(transformed)
    end function butter_filter

    pure function cheby1_filter(n, rp, w, filter_type, digital) result(filt)
        integer, intent(in) :: n !! Positive Chebyshev type-I prototype order.
        real(dp), intent(in) :: rp !! Pass-band ripple in decibels.
        real(dp), intent(in) :: w(:) !! One or two critical frequencies; digital values are normalised to Nyquist.
        character(len=*), intent(in), optional :: filter_type !! low, high, stop, or pass; defaults to low.
        logical, intent(in), optional :: digital !! True for a z-plane design; defaults to true.
        type(arma_filter) :: filt
        type(zpg_filter) :: proto
        type(zpg_filter) :: transformed
        real(dp), allocatable :: work_w(:)
        real(dp) :: epsilon_ripple
        real(dp) :: v0
        real(dp) :: angle
        character(len=8) :: kind
        logical :: is_digital
        logical :: is_stop
        integer :: i

        kind = 'low'
        if (present(filter_type)) kind = adjustl(filter_type)
        is_digital = .true.
        if (present(digital)) is_digital = digital
        is_stop = kind == 'stop' .or. kind == 'high'
        epsilon_ripple = sqrt(10.0_dp ** (rp / 10.0_dp) - 1.0_dp)
        v0 = asinh(1.0_dp / epsilon_ripple) / real(n, dp)
        allocate(proto%pole(max(0, n)), proto%zero(0))
        do i = 1, n
            angle = signal_pi * real(-(n - 1) + 2 * (i - 1), dp) / real(2 * n, dp)
            proto%pole(i) = cmplx(-sinh(v0) * cos(angle), cosh(v0) * sin(angle), dp)
        end do
        proto%gain = real(complex_product(-proto%pole), dp)
        if (modulo(n, 2) == 0) proto%gain = proto%gain / 10.0_dp ** (rp / 20.0_dp)
        work_w = w
        if (is_digital) work_w = tan(signal_pi * work_w / 2.0_dp)
        transformed = splane_frequency_transform(proto, work_w, is_stop)
        if (is_digital) transformed = bilinear_transform(transformed, 2.0_dp)
        filt = arma_from_zpg(transformed)
    end function cheby1_filter

    pure function cheby2_filter(n, rs, w, filter_type, digital) result(filt)
        integer, intent(in) :: n !! Positive Chebyshev type-II prototype order.
        real(dp), intent(in) :: rs !! Stop-band attenuation in decibels, named Rp in upstream signal::cheby2.
        real(dp), intent(in) :: w(:) !! One or two critical frequencies; digital values are normalised to Nyquist.
        character(len=*), intent(in), optional :: filter_type !! low, high, stop, or pass; defaults to low.
        logical, intent(in), optional :: digital !! True for a z-plane design; defaults to true.
        type(arma_filter) :: filt
        type(zpg_filter) :: proto
        type(zpg_filter) :: transformed
        real(dp), allocatable :: work_w(:)
        real(dp) :: lambda
        real(dp) :: phi
        real(dp) :: theta
        real(dp) :: alpha
        real(dp) :: beta
        character(len=8) :: kind
        logical :: is_digital
        logical :: is_stop
        integer :: i
        integer :: iz

        kind = 'low'
        if (present(filter_type)) kind = adjustl(filter_type)
        is_digital = .true.
        if (present(digital)) is_digital = digital
        is_stop = kind == 'stop' .or. kind == 'high'
        lambda = 10.0_dp ** (rs / 20.0_dp)
        phi = log(lambda + sqrt(lambda * lambda - 1.0_dp)) / real(n, dp)
        allocate(proto%pole(max(0, n)), proto%zero(n - modulo(n, 2)))
        iz = 0
        do i = 1, n
            theta = signal_pi * (real(i, dp) - 0.5_dp) / real(n, dp)
            alpha = -sinh(phi) * sin(theta)
            beta = cosh(phi) * cos(theta)
            proto%pole(i) = cmplx(alpha, -beta, dp) / (alpha * alpha + beta * beta)
            if (.not. (modulo(n, 2) == 1 .and. i == (n + 1) / 2)) then
                iz = iz + 1
                proto%zero(iz) = cmplx(0.0_dp, 1.0_dp / cos(theta), dp)
            end if
        end do
        proto%gain = abs(real(complex_product(proto%pole) / complex_product(proto%zero), dp))
        work_w = w
        if (is_digital) work_w = tan(signal_pi * work_w / 2.0_dp)
        transformed = splane_frequency_transform(proto, work_w, is_stop)
        if (is_digital) transformed = bilinear_transform(transformed, 2.0_dp)
        filt = arma_from_zpg(transformed)
    end function cheby2_filter

    pure function butter_order(wp, ws, rp, rs) result(spec)
        real(dp), intent(in) :: wp(:) !! Pass-band edge or pair of edges normalised to Nyquist.
        real(dp), intent(in) :: ws(:) !! Stop-band edge or pair of edges normalised to Nyquist.
        real(dp), intent(in) :: rp !! Maximum pass-band loss in decibels.
        real(dp), intent(in) :: rs !! Minimum stop-band attenuation in decibels.
        type(filter_order) :: spec
        real(dp), allocatable :: wpw(:)
        real(dp), allocatable :: wsw(:)
        real(dp), allocatable :: wc(:)
        real(dp) :: qs
        real(dp) :: qp
        real(dp) :: candidate
        integer :: n
        integer :: i
        character(len=8) :: kind

        wpw = wp
        wsw = ws
        kind = 'low'
        do i = 1, min(size(wpw), size(wsw))
            if (wpw(i) > wsw(i)) then
                wpw(i) = 1.0_dp - wpw(i)
                wsw(i) = 1.0_dp - wsw(i)
            end if
        end do
        if (size(wp) == 2) then
            if (any(wp > ws)) then
                kind = 'stop'
            else
                kind = 'pass'
            end if
        else if (wp(1) > ws(1)) then
            kind = 'high'
        end if
        wpw = tan(signal_pi * wpw / 2.0_dp)
        wsw = tan(signal_pi * wsw / 2.0_dp)
        qs = log(10.0_dp ** (rs / 10.0_dp) - 1.0_dp)
        qp = log(10.0_dp ** (rp / 10.0_dp) - 1.0_dp)
        candidate = 0.0_dp
        do i = 1, size(wpw)
            candidate = max(candidate, 0.5_dp * (qs - qp) / log(wsw(i) / wpw(i)))
        end do
        n = ceiling(candidate / real(size(wpw), dp))
        allocate(wc(size(wpw)))
        wc = exp(log(wpw) - qp / (2.0_dp * real(n, dp)))
        wc = 2.0_dp * atan(wc) / signal_pi
        do i = 1, min(size(wc), size(wp))
            if (wp(i) > ws(i)) wc(i) = 1.0_dp - wc(i)
        end do
        spec = make_filter_order(n, wc, kind, rp=rp, rs=rs)
    end function butter_order

    pure function cheby1_order(wp, ws, rp, rs) result(spec)
        real(dp), intent(in) :: wp(:) !! Pass-band edge or pair of edges normalised to Nyquist.
        real(dp), intent(in) :: ws(:) !! Stop-band edge or pair of edges normalised to Nyquist.
        real(dp), intent(in) :: rp !! Pass-band ripple in decibels.
        real(dp), intent(in) :: rs !! Stop-band attenuation in decibels.
        type(filter_order) :: spec
        real(dp), allocatable :: wpw(:)
        real(dp), allocatable :: wsw(:)
        real(dp), allocatable :: wa(:)
        real(dp) :: ratio
        integer :: n
        character(len=8) :: kind

        wpw = tan(signal_pi * wp / 2.0_dp)
        wsw = tan(signal_pi * ws / 2.0_dp)
        if (size(wp) == 1) then
            if (wp(1) < ws(1)) then
                kind = 'low'
                allocate(wa(1))
                wa = wsw / wpw
            else
                kind = 'high'
                allocate(wa(1))
                wa = wpw / wsw
            end if
        else if (wp(1) > ws(1)) then
            kind = 'pass'
            allocate(wa(size(ws)))
            wa = (wsw * wsw - wpw(1) * wpw(2)) / (wsw * (wpw(1) - wpw(2)))
        else
            kind = 'stop'
            n = 0
            spec = make_filter_order(n, wp, kind, rp=rp, rs=rs)
            return
        end if
        ratio = minval(abs(wa))
        n = ceiling(acosh(sqrt((10.0_dp ** (abs(rs) / 10.0_dp) - 1.0_dp) / &
            (10.0_dp ** (abs(rp) / 10.0_dp) - 1.0_dp))) / acosh(ratio))
        spec = make_filter_order(n, wp, kind, rp=rp, rs=rs)
    end function cheby1_order

    pure function ellip_order(wp, ws, rp, rs) result(spec)
        real(dp), intent(in) :: wp(:) !! Pass-band edge or pair of edges normalised to Nyquist.
        real(dp), intent(in) :: ws(:) !! Stop-band edge or pair of edges normalised to Nyquist.
        real(dp), intent(in) :: rp !! Pass-band ripple in decibels.
        real(dp), intent(in) :: rs !! Stop-band attenuation in decibels.
        type(filter_order) :: spec
        real(dp), allocatable :: wpw(:)
        real(dp), allocatable :: wsw(:)
        real(dp) :: wp0
        real(dp) :: ws0
        real(dp) :: w02
        real(dp) :: w3
        real(dp) :: w4
        real(dp) :: k
        real(dp) :: k1
        real(dp) :: q0
        real(dp) :: q
        real(dp) :: d
        integer :: n
        character(len=8) :: kind

        wpw = tan(signal_pi * wp / 2.0_dp)
        wsw = tan(signal_pi * ws / 2.0_dp)
        if (size(wp) == 2 .and. size(ws) == 2) then
            kind = 'pass'
            wp0 = 1.0_dp
            w02 = wpw(1) * wpw(2)
            w3 = w02 / wsw(2)
            w4 = w02 / wsw(1)
            if (w3 > wsw(1)) then
                ws0 = (wsw(2) - w3) / (wpw(2) - wpw(1))
            else if (w4 < wsw(2)) then
                ws0 = (w4 - wsw(1)) / (wpw(2) - wpw(1))
            else
                ws0 = (wsw(2) - wsw(1)) / (wpw(2) - wpw(1))
            end if
        else if (size(wp) == 2 .and. size(ws) == 1) then
            kind = 'pass'
            wp0 = 1.0_dp
            w02 = wpw(1) * wpw(2)
            if (wsw(1) > wpw(2)) then
                w3 = w02 / wsw(1)
                ws0 = (wsw(1) - w3) / (wpw(2) - wpw(1))
            else
                w4 = w02 / wsw(1)
                ws0 = (w4 - wsw(1)) / (wpw(2) - wpw(1))
            end if
        else
            if (wpw(1) < wsw(1)) then
                kind = 'low'
            else
                kind = 'high'
            end if
            wp0 = wpw(1)
            ws0 = wsw(1)
        end if
        k = wp0 / ws0
        k1 = sqrt(max(0.0_dp, 1.0_dp - k * k))
        q0 = 0.5_dp * (1.0_dp - sqrt(k1)) / (1.0_dp + sqrt(k1))
        q = q0 + 2.0_dp * q0 ** 5 + 15.0_dp * q0 ** 9 + 150.0_dp * q0 ** 13
        d = (10.0_dp ** (0.1_dp * rs) - 1.0_dp) / (10.0_dp ** (0.1_dp * rp) - 1.0_dp)
        n = ceiling(log10(16.0_dp * d) / log10(1.0_dp / q))
        spec = make_filter_order(n, wp, kind, rp=rp, rs=rs)
    end function ellip_order

    pure function ellip_filter(n, rp, rs, w, filter_type, digital) result(filt)
        integer, intent(in) :: n !! Positive elliptic prototype order.
        real(dp), intent(in) :: rp !! Pass-band ripple in decibels.
        real(dp), intent(in) :: rs !! Stop-band attenuation in decibels.
        real(dp), intent(in) :: w(:) !! One or two critical frequencies; digital values are normalised to Nyquist.
        character(len=*), intent(in), optional :: filter_type !! low, high, stop, or pass; defaults to low.
        logical, intent(in), optional :: digital !! True for a z-plane design; defaults to true.
        type(arma_filter) :: filt
        type(zpg_filter) :: proto
        type(zpg_filter) :: transformed
        real(dp), allocatable :: work_w(:)
        character(len=8) :: kind
        logical :: is_digital
        logical :: is_stop

        kind = 'low'
        if (present(filter_type)) kind = adjustl(filter_type)
        is_digital = .true.
        if (present(digital)) is_digital = digital
        is_stop = kind == 'stop' .or. kind == 'high'
        proto = ncauer(rp, rs, n)
        work_w = w
        if (is_digital) work_w = tan(signal_pi * work_w / 2.0_dp)
        transformed = splane_frequency_transform(proto, work_w, is_stop)
        if (is_digital) transformed = bilinear_transform(transformed, 2.0_dp)
        filt = arma_from_zpg(transformed)
    end function ellip_filter

    pure function ncauer(rp, rs, n) result(proto)
        real(dp), intent(in) :: rp !! Pass-band ripple in decibels.
        real(dp), intent(in) :: rs !! Stop-band attenuation in decibels.
        integer, intent(in) :: n !! Positive elliptic prototype order.
        type(zpg_filter) :: proto
        real(dp), allocatable :: wi(:)
        real(dp), allocatable :: vi(:)
        real(dp), allocatable :: a0i(:)
        real(dp), allocatable :: b0i(:)
        real(dp) :: ws
        real(dp) :: k
        real(dp) :: k1
        real(dp) :: q0
        real(dp) :: q
        real(dp) :: lambda_l
        real(dp) :: sig01
        real(dp) :: sig02
        real(dp) :: sig0
        real(dp) :: omega
        real(dp) :: mu
        real(dp) :: sum1
        real(dp) :: sum2
        real(dp) :: t0
        complex(dp), allocatable :: zeros(:)
        complex(dp), allocatable :: poles(:)
        integer :: r
        integer :: i
        integer :: m

        ws = elliptic_stop_edge(n, rp, rs)
        k = 1.0_dp / ws
        k1 = sqrt(max(0.0_dp, 1.0_dp - k * k))
        q0 = 0.5_dp * (1.0_dp - sqrt(k1)) / (1.0_dp + sqrt(k1))
        q = q0 + 2.0_dp * q0 ** 5 + 15.0_dp * q0 ** 9 + 150.0_dp * q0 ** 13
        lambda_l = log((10.0_dp ** (0.05_dp * rp) + 1.0_dp) / &
            (10.0_dp ** (0.05_dp * rp) - 1.0_dp)) / real(2 * n, dp)
        sig01 = 0.0_dp
        sig02 = 0.0_dp
        do m = 0, 30
            sig01 = sig01 + (-1.0_dp) ** m * q ** (m * (m + 1)) * sinh(real(2 * m + 1, dp) * lambda_l)
        end do
        do m = 1, 30
            sig02 = sig02 + (-1.0_dp) ** m * q ** (m * m) * cosh(real(2 * m, dp) * lambda_l)
        end do
        sig0 = abs(2.0_dp * q ** 0.25_dp * sig01 / (1.0_dp + 2.0_dp * sig02))
        omega = sqrt((1.0_dp + k * sig0 * sig0) * (1.0_dp + sig0 * sig0 / k))
        r = (n - modulo(n, 2)) / 2
        allocate(wi(r), vi(r), a0i(r), b0i(r))
        do i = 1, r
            mu = real(i, dp) - real(1 - modulo(n, 2), dp) / 2.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            do m = 0, 30
                sum1 = sum1 + 2.0_dp * q ** 0.25_dp * (-1.0_dp) ** m * q ** (m * (m + 1)) * &
                    sin(real(2 * m + 1, dp) * signal_pi * mu / real(n, dp))
            end do
            do m = 1, 30
                sum2 = sum2 + 2.0_dp * (-1.0_dp) ** m * q ** (m * m) * &
                    cos(real(2 * m, dp) * signal_pi * mu / real(n, dp))
            end do
            wi(i) = sum1 / (1.0_dp + sum2)
        end do
        vi = sqrt((1.0_dp - k * wi * wi) * (1.0_dp - wi * wi / k))
        a0i = 1.0_dp / (wi * wi)
        b0i = ((sig0 * vi) ** 2 + (omega * wi) ** 2) / (1.0_dp + sig0 * sig0 * wi * wi) ** 2
        if (modulo(n, 2) == 1) then
            t0 = sig0 * product(b0i / a0i) * sqrt(ws)
        else
            t0 = 10.0_dp ** (-0.05_dp * rp) * product(b0i / a0i)
        end if
        allocate(zeros(2 * r), poles(2 * r + modulo(n, 2)))
        do i = 1, r
            zeros(i) = cmplx(0.0_dp, sqrt(ws) / wi(i), dp)
            zeros(r + i) = -zeros(i)
            poles(i) = sqrt(ws) * cmplx(-sig0 * vi(i), wi(i) * omega, dp) / &
                (1.0_dp + sig0 * sig0 * wi(i) * wi(i))
            poles(r + i) = conjg(poles(i))
        end do
        if (modulo(n, 2) == 1) poles(2 * r + 1) = cmplx(-sqrt(ws) * sig0, 0.0_dp, dp)
        proto = make_zpg(zeros, poles, t0)
    end function ncauer

    pure function elliptic_stop_edge(n, rp, rs) result(ws)
        integer, intent(in) :: n !! Elliptic filter order.
        real(dp), intent(in) :: rp !! Pass-band ripple in decibels.
        real(dp), intent(in) :: rs !! Stop-band attenuation in decibels.
        real(dp) :: ws
        real(dp) :: kl0
        real(dp) :: x_target
        real(dp) :: lo
        real(dp) :: hi
        real(dp) :: c
        real(dp) :: d
        real(dp) :: fc
        real(dp) :: fd
        real(dp), parameter :: phi = 0.6180339887498948482_dp
        integer :: iter

        kl0 = (10.0_dp ** (0.1_dp * rp) - 1.0_dp) / (10.0_dp ** (0.1_dp * rs) - 1.0_dp)
        x_target = real(n, dp) * elliptic_k(kl0) / elliptic_k(1.0_dp - kl0)
        lo = epsilon(1.0_dp)
        hi = 1.0_dp - epsilon(1.0_dp)
        c = hi - phi * (hi - lo)
        d = lo + phi * (hi - lo)
        fc = abs(elliptic_k(c) / elliptic_k(1.0_dp - c) - x_target)
        fd = abs(elliptic_k(d) / elliptic_k(1.0_dp - d) - x_target)
        do iter = 1, 160
            if (fc < fd) then
                hi = d
                d = c
                fd = fc
                c = hi - phi * (hi - lo)
                fc = abs(elliptic_k(c) / elliptic_k(1.0_dp - c) - x_target)
            else
                lo = c
                c = d
                fc = fd
                d = lo + phi * (hi - lo)
                fd = abs(elliptic_k(d) / elliptic_k(1.0_dp - d) - x_target)
            end if
        end do
        ws = sqrt(1.0_dp / ((lo + hi) / 2.0_dp))
    end function elliptic_stop_edge

    pure elemental function elliptic_k(m) result(kval)
        real(dp), intent(in) :: m !! Elliptic parameter m, normally in [0,1).
        real(dp) :: kval
        real(dp) :: a
        real(dp) :: b
        real(dp) :: next_a
        real(dp) :: transform
        real(dp) :: mm
        integer :: iter

        if (m >= 1.0_dp) then
            kval = huge(1.0_dp)
            return
        end if
        transform = 1.0_dp
        mm = m
        if (m < 0.0_dp) then
            transform = 1.0_dp / sqrt(1.0_dp - m)
            mm = -m / (1.0_dp - m)
        end if
        a = 1.0_dp
        b = sqrt(max(0.0_dp, 1.0_dp - mm))
        do iter = 1, 64
            next_a = 0.5_dp * (a + b)
            b = sqrt(a * b)
            a = next_a
            if (abs(a - b) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(a))) exit
        end do
        kval = transform * signal_pi / (2.0_dp * a)
    end function elliptic_k

    pure function complex_product(x) result(value)
        complex(dp), intent(in) :: x(:) !! Complex vector whose elements are multiplied; an empty vector has product one.
        complex(dp) :: value
        integer :: i

        value = cmplx(1.0_dp, 0.0_dp, dp)
        do i = 1, size(x)
            value = value * x(i)
        end do
    end function complex_product

end module signal_iir_design
