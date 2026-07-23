! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

module longmemo
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use longmemo_kinds, only : dp, pi, two_pi
    use longmemo_fft, only : fft_forward, fft_inverse_raw
    use longmemo_linalg, only : inverse_spd
    use longmemo_optimize, only : minimize_scalar, nelder_mead
    use longmemo_stats, only : normal_sf, student_t_cdf, sample_sd, &
                               random_normal_vector, random_exponential_vector
    implicit none
    private

    type, public :: spectrum_result
        real(dp), allocatable :: frequency(:)
        real(dp), allocatable :: spectrum(:)
        real(dp) :: theta1 = 1.0_dp
        real(dp) :: hurst = 0.5_dp
        integer :: p = 0
        integer :: q = 0
        real(dp), allocatable :: eta(:)
    end type spectrum_result

    type, public :: qeta_result
        integer :: n = 0
        real(dp) :: hurst = 0.5_dp
        real(dp), allocatable :: eta(:)
        real(dp) :: a = 0.0_dp
        real(dp) :: b = 0.0_dp
        real(dp) :: tn = 0.0_dp
        real(dp) :: z = 0.0_dp
        real(dp) :: p_value = 0.0_dp
        real(dp) :: theta1 = 0.0_dp
        real(dp), allocatable :: spectrum(:)
    end type qeta_result

    type, public :: whittle_result
        character(len=8) :: model = "fGn"
        integer :: n = 0
        integer :: p = 0
        integer :: q = 0
        logical :: converged = .false.
        real(dp) :: theta1 = 0.0_dp
        real(dp), allocatable :: eta(:)
        real(dp), allocatable :: std_error(:)
        real(dp), allocatable :: z_value(:)
        real(dp), allocatable :: p_value(:)
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: periodogram(:)
        real(dp), allocatable :: spectrum(:)
        type(qeta_result) :: goodness_of_fit
    end type whittle_result

    type, public :: fexp_result
        integer :: n = 0
        integer :: order_poly = 0
        integer :: max_order_poly = 0
        logical :: early_stop = .false.
        real(dp) :: max_p_value = 0.0_dp
        real(dp) :: hurst = 0.5_dp
        real(dp), allocatable :: coefficients(:, :)
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: periodogram(:)
        real(dp), allocatable :: spectrum(:)
    end type fexp_result

    type :: whittle_context
        character(len=8) :: model = "fGn"
        integer :: n = 0
        integer :: p = 0
        integer :: q = 0
        real(dp), allocatable :: periodogram(:)
    end type whittle_context

    public :: b_spec_fgn, spec_fgn, spec_arima
    public :: ck_fgn0, ck_arma0
    public :: sim_gauss, sim_fgn0, sim_arma0, sim_fgn_fft
    public :: periodogram, fourier_frequencies
    public :: qeta, ceta_fgn, ceta_arima
    public :: whittle_estimate, fexp_estimate

contains

    subroutine b_spec_fgn(lambda, hurst, b, k_approx, adjust, nsum)
        real(dp), intent(in) :: lambda(:), hurst
        real(dp), allocatable, intent(out) :: b(:)
        integer, intent(in), optional :: k_approx, nsum
        logical, intent(in), optional :: adjust
        real(dp) :: d, d1, a1, b1, a3, b3, aj, bj, anext, bnext
        real(dp) :: correction
        integer :: i, j, k, ns
        logical :: do_adjust

        k = 3
        if (present(k_approx)) k = k_approx
        ns = 200
        if (present(nsum)) ns = nsum
        do_adjust = k == 3
        if (present(adjust)) do_adjust = adjust
        if (hurst <= 0.0_dp) error stop "b_spec_fgn: hurst must be positive"

        allocate(b(size(lambda)))
        d = -2.0_dp*hurst - 1.0_dp
        d1 = -2.0_dp*hurst

        if (k <= 0) then
            if (ns < 5) error stop "b_spec_fgn: nsum must be at least 5"
            do i = 1, size(lambda)
                b(i) = 0.0_dp
                do j = 1, ns
                    b(i) = b(i) + abs(two_pi*real(j, dp) + lambda(i))**d + &
                                    abs(two_pi*real(j, dp) - lambda(i))**d
                end do
            end do
            return
        end if

        do i = 1, size(lambda)
            a1 = two_pi + lambda(i)
            b1 = two_pi - lambda(i)
            if (k == 3) then
                a3 = 3.0_dp*two_pi + lambda(i)
                b3 = 3.0_dp*two_pi - lambda(i)
                b(i) = a1**d + b1**d + (a1 + two_pi)**d + (b1 + two_pi)**d + &
                       a3**d + b3**d + &
                       (a3**d1 + b3**d1 + (a3 + two_pi)**d1 + (b3 + two_pi)**d1)/ &
                       (8.0_dp*pi*hurst)
                if (do_adjust) then
                    correction = 2.0_dp**(-7.65_dp*hurst - 7.4_dp)
                    b(i) = (1.0002_dp - 0.000134_dp*lambda(i))*(b(i) - correction)
                end if
            else
                aj = a1
                bj = b1
                b(i) = aj**d + bj**d
                anext = aj + two_pi
                bnext = bj + two_pi
                do j = 1, k - 1
                    aj = anext
                    bj = bnext
                    b(i) = b(i) + aj**d + bj**d
                    anext = aj + two_pi
                    bnext = bj + two_pi
                end do
                b(i) = b(i) + (aj**d1 + bj**d1 + anext**d1 + bnext**d1)/ &
                                (8.0_dp*pi*hurst)
            end if
        end do
    end subroutine b_spec_fgn


    subroutine spec_fgn(eta, m, result, k_approx, adjust, nsum)
        real(dp), intent(in) :: eta(:)
        integer, intent(in) :: m
        type(spectrum_result), intent(out) :: result
        integer, intent(in), optional :: k_approx, nsum
        logical, intent(in), optional :: adjust

        if (size(eta) < 1) error stop "spec_fgn: eta must contain H"
        call spec_fgn_values(eta(1), m, result%frequency, result%spectrum, result%theta1, &
                             k_approx, adjust, nsum)
        result%hurst = eta(1)
        result%p = 0
        result%q = 0
        allocate(result%eta(size(eta)))
        result%eta = eta
    end subroutine spec_fgn


    subroutine spec_fgn_values(hurst, m, frequency, spectrum, theta1, k_approx, adjust, nsum)
        real(dp), intent(in) :: hurst
        integer, intent(in) :: m
        real(dp), allocatable, intent(out) :: frequency(:), spectrum(:)
        real(dp), intent(out) :: theta1
        integer, intent(in), optional :: k_approx, nsum
        logical, intent(in), optional :: adjust
        real(dp), allocatable :: b(:)
        real(dp) :: constant
        integer :: m2, j

        if (m < 2) error stop "spec_fgn: m must be at least 2"
        m2 = (m - 1)/2
        allocate(frequency(m2), spectrum(m2))
        do j = 1, m2
            frequency(j) = two_pi*real(j, dp)/real(m, dp)
        end do
        if (m2 == 0) then
            theta1 = 1.0_dp
            return
        end if

        call b_spec_fgn(frequency, hurst, b, k_approx, adjust, nsum)
        constant = sin(pi*hurst)/pi*gamma(2.0_dp*hurst + 1.0_dp)
        spectrum = constant*(1.0_dp - cos(frequency))* &
                   (frequency**(-2.0_dp*hurst - 1.0_dp) + b)
        theta1 = exp(2.0_dp/real(m, dp)*sum(log(spectrum)))
        spectrum = spectrum/theta1
    end subroutine spec_fgn_values


    subroutine spec_arima(eta, p, q, m, result)
        real(dp), intent(in) :: eta(:)
        integer, intent(in) :: p, q, m
        type(spectrum_result), intent(out) :: result

        call spec_arima_values(eta, p, q, m, result%frequency, result%spectrum)
        result%theta1 = 1.0_dp/(2.0_dp*pi)
        result%hurst = eta(1)
        result%p = p
        result%q = q
        allocate(result%eta(size(eta)))
        result%eta = eta
    end subroutine spec_arima


    subroutine spec_arima_values(eta, p, q, m, frequency, spectrum)
        real(dp), intent(in) :: eta(:)
        integer, intent(in) :: p, q, m
        real(dp), allocatable, intent(out) :: frequency(:), spectrum(:)
        real(dp) :: ar_re, ar_im, ma_re, ma_im, far, fma, angle
        integer :: m2, i, j

        if (p < 0 .or. q < 0) error stop "spec_arima: p and q must be nonnegative"
        if (size(eta) /= 1 + p + q) error stop "spec_arima: eta size must equal 1+p+q"
        if (m < 2) error stop "spec_arima: m must be at least 2"

        m2 = (m - 1)/2
        allocate(frequency(m2), spectrum(m2))
        do i = 1, m2
            frequency(i) = two_pi*real(i, dp)/real(m, dp)
            ar_re = 0.0_dp
            ar_im = 0.0_dp
            do j = 1, p
                angle = frequency(i)*real(j, dp)
                ar_re = ar_re + eta(1 + j)*cos(angle)
                ar_im = ar_im + eta(1 + j)*sin(angle)
            end do
            far = (1.0_dp - ar_re)**2 + ar_im**2

            ma_re = 0.0_dp
            ma_im = 0.0_dp
            do j = 1, q
                angle = frequency(i)*real(j, dp)
                ma_re = ma_re + eta(1 + p + j)*cos(angle)
                ma_im = ma_im + eta(1 + p + j)*sin(angle)
            end do
            fma = (1.0_dp + ma_re)**2 + ma_im**2
            spectrum(i) = fma/far*sqrt(2.0_dp - 2.0_dp*cos(frequency(i)))** &
                          (1.0_dp - 2.0_dp*eta(1))
        end do
    end subroutine spec_arima_values


    subroutine ck_fgn0(n, hurst, autocov)
        integer, intent(in) :: n
        real(dp), intent(in) :: hurst
        real(dp), allocatable, intent(out) :: autocov(:)
        real(dp) :: k
        integer :: i

        if (n < 1) error stop "ck_fgn0: n must be positive"
        allocate(autocov(n))
        do i = 1, n
            k = real(i - 1, dp)
            autocov(i) = 0.5_dp*(abs(k - 1.0_dp)**(2.0_dp*hurst) - &
                                  2.0_dp*k**(2.0_dp*hurst) + &
                                  (k + 1.0_dp)**(2.0_dp*hurst))
        end do
    end subroutine ck_fgn0


    subroutine ck_arma0(n, hurst, autocov)
        integer, intent(in) :: n
        real(dp), intent(in) :: hurst
        real(dp), allocatable, intent(out) :: autocov(:)
        real(dp) :: d, g1d, gd, k
        integer :: i, k_exact

        if (n < 1) error stop "ck_arma0: n must be positive"
        d = hurst - 0.5_dp
        if (d <= 0.0_dp .or. d >= 0.5_dp) then
            error stop "ck_arma0: H must be strictly between 0.5 and 1"
        end if

        allocate(autocov(n))
        g1d = gamma(1.0_dp - d)
        gd = pi/(sin(pi*d)*g1d)
        autocov(1) = gamma(1.0_dp - 2.0_dp*d)/(g1d*g1d)
        k_exact = min(50, n - 1)
        do i = 1, k_exact
            k = real(i, dp)
            autocov(i + 1) = autocov(1)*gamma(k + d)*g1d/ &
                             (gamma(k - d + 1.0_dp)*gd)
        end do
        do i = 51, n - 1
            k = real(i, dp)
            autocov(i + 1) = autocov(1)*g1d/gd*k**(2.0_dp*hurst - 2.0_dp)
        end do
    end subroutine ck_arma0


    subroutine sim_gauss(autocov, x)
        real(dp), intent(in) :: autocov(:)
        real(dp), allocatable, intent(out) :: x(:)
        real(dp), allocatable :: covariance_circle(:), gk(:), zr0(:), zi0(:), zr(:), zi(:)
        complex(dp), allocatable :: cseq(:), transformed(:), zc(:), ztime(:)
        integer :: n, m, i

        n = size(autocov)
        if (n < 2) error stop "sim_gauss: need at least two autocovariances"
        m = 2*(n - 1)
        allocate(covariance_circle(m), cseq(m), gk(m))
        covariance_circle(1:n) = autocov
        if (n >= 3) then
            do i = 1, n - 2
                covariance_circle(n + i) = autocov(n - i)
            end do
        end if
        cseq = cmplx(covariance_circle, 0.0_dp, dp)
        call fft_inverse_raw(cseq, transformed)
        gk = real(transformed, dp)
        where (gk > -100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(gk))))
            gk = max(gk, 0.0_dp)
        end where
        if (any(gk < 0.0_dp)) error stop "sim_gauss: autocovariance embedding is not nonnegative"

        allocate(zr0(n), zi0(max(0, n - 2)), zr(m), zi(m), zc(m))
        call random_normal_vector(zr0)
        if (n > 2) call random_normal_vector(zi0)
        zr0(1) = sqrt(2.0_dp)*zr0(1)
        zr0(n) = sqrt(2.0_dp)*zr0(n)
        zr(1:n) = zr0
        do i = 1, n - 2
            zr(n + i) = zr0(n - i)
        end do
        zi = 0.0_dp
        if (n > 2) then
            zi(2:n - 1) = zi0
            do i = 1, n - 2
                zi(n + i) = -zi0(n - 1 - i)
            end do
        end if
        zc = cmplx(zr, zi, dp)*sqrt(gk)
        call fft_inverse_raw(zc, ztime)
        allocate(x(n))
        x = real(ztime(1:n), dp)/(2.0_dp*sqrt(real(n - 1, dp)))
    end subroutine sim_gauss


    subroutine sim_fgn0(n, hurst, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: hurst
        real(dp), allocatable, intent(out) :: x(:)
        real(dp), allocatable :: autocov(:)

        call ck_fgn0(n, hurst, autocov)
        call sim_gauss(autocov, x)
    end subroutine sim_fgn0


    subroutine sim_arma0(n, hurst, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: hurst
        real(dp), allocatable, intent(out) :: x(:)
        real(dp), allocatable :: autocov(:)

        call ck_arma0(n, hurst, autocov)
        call sim_gauss(autocov, x)
    end subroutine sim_arma0


    subroutine sim_fgn_fft(n, hurst, x, k_approx, adjust, nsum)
        integer, intent(in) :: n
        real(dp), intent(in) :: hurst
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(in), optional :: k_approx, nsum
        logical, intent(in), optional :: adjust
        real(dp), allocatable :: lambda(:), b(:), f(:), expo(:), phase(:)
        complex(dp), allocatable :: z(:), full_z(:), time_z(:)
        integer :: nh, i

        if (n < 2 .or. modulo(n, 2) /= 0) error stop "sim_fgn_fft: n must be positive and even"
        nh = n/2
        allocate(lambda(nh), f(nh), expo(nh), phase(nh), z(nh), full_z(n))
        do i = 1, nh
            lambda(i) = pi*real(i, dp)/real(nh, dp)
        end do
        call b_spec_fgn(lambda, hurst, b, k_approx, adjust, nsum)
        f = 2.0_dp*sin(pi*hurst)*gamma(2.0_dp*hurst + 1.0_dp)* &
            (1.0_dp - cos(lambda))*(lambda**(-2.0_dp*hurst - 1.0_dp) + b)
        call random_exponential_vector(expo)
        call random_number(phase)
        f = f*expo
        z = sqrt(f)*cmplx(cos(two_pi*phase), sin(two_pi*phase), dp)
        z(nh) = cmplx(abs(z(nh)), 0.0_dp, dp)
        full_z(1) = (0.0_dp, 0.0_dp)
        full_z(2:nh + 1) = z
        do i = 2, nh
            full_z(nh + i) = conjg(z(nh - i + 2))
        end do
        call fft_inverse_raw(full_z, time_z)
        allocate(x(n))
        x = real(time_z, dp)
    end subroutine sim_fgn_fft


    subroutine periodogram(z, values)
        real(dp), intent(in) :: z(:)
        real(dp), allocatable, intent(out) :: values(:)
        complex(dp), allocatable :: input(:), output(:)
        integer :: n, nout

        n = size(z)
        if (n < 1) error stop "periodogram: input cannot be empty"
        nout = n/2 + 1
        allocate(input(n))
        input = cmplx(z, 0.0_dp, dp)
        call fft_forward(input, output)
        allocate(values(nout))
        values = abs(output(1:nout))**2/(two_pi*real(n, dp))
    end subroutine periodogram


    subroutine fourier_frequencies(n, frequency, full)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: frequency(:)
        logical, intent(in), optional :: full
        integer :: n2, j
        logical :: use_full

        if (n < 1) error stop "fourier_frequencies: n must be positive"
        use_full = .false.
        if (present(full)) use_full = full
        if (use_full) then
            n2 = n/2
        else
            n2 = (n - 1)/2
        end if
        allocate(frequency(n2))
        do j = 1, n2
            frequency(j) = two_pi*real(j, dp)/real(n, dp)
        end do
    end subroutine fourier_frequencies


    subroutine qeta(eta, model, n, yper, result, p, q, give_b_only, b_only)
        real(dp), intent(in) :: eta(:), yper(:)
        character(len=*), intent(in) :: model
        integer, intent(in) :: n
        type(qeta_result), intent(out), optional :: result
        integer, intent(in), optional :: p, q
        logical, intent(in), optional :: give_b_only
        real(dp), intent(out), optional :: b_only
        real(dp), allocatable :: frequency(:), spectrum(:), yf(:)
        real(dp) :: theta_dummy, a_value, b_value
        integer :: pp, qq
        logical :: only_b

        if (size(yper) /= (n - 1)/2) error stop "qeta: yper has the wrong length"
        only_b = .false.
        if (present(give_b_only)) only_b = give_b_only
        pp = 0
        qq = 0
        if (present(p)) pp = p
        if (present(q)) qq = q

        select case (trim(adjustl(model)))
        case ("fGn", "fgn", "FGN")
            call spec_fgn_values(eta(1), n, frequency, spectrum, theta_dummy)
        case ("fARIMA", "farima", "FARIMA")
            call spec_arima_values(eta, pp, qq, n, frequency, spectrum)
        case default
            error stop "qeta: model must be fGn or fARIMA"
        end select

        if (any(spectrum <= 0.0_dp) .or. any(.not. ieee_is_finite(spectrum))) then
            b_value = huge(1.0_dp)
            if (present(b_only)) b_only = b_value
            return
        end if
        allocate(yf(size(yper)))
        yf = yper/spectrum
        b_value = 2.0_dp*(two_pi/real(n, dp))*sum(yf)
        if (present(b_only)) b_only = b_value
        if (only_b .or. .not. present(result)) return

        a_value = 2.0_dp*(two_pi/real(n, dp))*sum(yf*yf)
        result%n = n
        result%hurst = eta(1)
        allocate(result%eta(size(eta)), result%spectrum(size(spectrum)))
        result%eta = eta
        result%spectrum = spectrum
        result%a = a_value
        result%b = b_value
        result%tn = a_value/(b_value*b_value)
        result%z = sqrt(real(n, dp))*(pi*result%tn - 1.0_dp)/sqrt(2.0_dp)
        result%p_value = normal_sf(result%z)
        result%theta1 = b_value/two_pi
    end subroutine qeta


    subroutine ceta_fgn(eta, covariance, m, delta)
        real(dp), intent(in) :: eta(:)
        real(dp), allocatable, intent(out) :: covariance(:, :)
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: delta
        real(dp), allocatable :: f0(:), fj(:), freq(:), lf(:, :), gram(:, :), inv(:, :), etaj(:)
        real(dp) :: theta, theta_j, h
        integer :: mm, m2, j, info

        if (size(eta) < 1) error stop "ceta_fgn: eta cannot be empty"
        mm = 10000
        if (present(m)) mm = m
        h = 1.0e-9_dp
        if (present(delta)) h = delta
        m2 = (mm - 1)/2
        allocate(lf(m2, size(eta)), etaj(size(eta)))
        call spec_fgn_values(eta(1), mm, freq, f0, theta)
        do j = 1, size(eta)
            etaj = eta
            etaj(j) = etaj(j) + h
            call spec_fgn_values(etaj(1), mm, freq, fj, theta_j)
            lf(:, j) = log(fj/f0)/h
        end do
        allocate(gram(size(eta), size(eta)))
        gram = matmul(transpose(lf), lf)
        call inverse_spd(gram, inv, info)
        if (info /= 0) error stop "ceta_fgn: singular derivative cross-product"
        allocate(covariance(size(eta), size(eta)))
        covariance = real(mm, dp)*inv
    end subroutine ceta_fgn


    subroutine ceta_arima(eta, p, q, covariance, m, delta)
        real(dp), intent(in) :: eta(:)
        integer, intent(in) :: p, q
        real(dp), allocatable, intent(out) :: covariance(:, :)
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: delta
        real(dp), allocatable :: f0(:), fj(:), freq(:), lf(:, :), gram(:, :), inv(:, :), etaj(:)
        real(dp) :: h
        integer :: mm, m2, j, info

        if (size(eta) /= 1 + p + q) error stop "ceta_arima: eta size must equal 1+p+q"
        mm = 10000
        if (present(m)) mm = m
        h = 1.0e-9_dp
        if (present(delta)) h = delta
        m2 = (mm - 1)/2
        allocate(lf(m2, size(eta)), etaj(size(eta)))
        call spec_arima_values(eta, p, q, mm, freq, f0)
        do j = 1, size(eta)
            etaj = eta
            etaj(j) = etaj(j) + h
            call spec_arima_values(etaj, p, q, mm, freq, fj)
            lf(:, j) = log(fj/f0)/h
        end do
        allocate(gram(size(eta), size(eta)))
        gram = matmul(transpose(lf), lf)
        call inverse_spd(gram, inv, info)
        if (info /= 0) error stop "ceta_arima: singular derivative cross-product"
        allocate(covariance(size(eta), size(eta)))
        covariance = real(mm, dp)*inv
    end subroutine ceta_arima


    subroutine whittle_estimate(x, model, result, p, q, start_eta, scale_data, covariance_m, delta)
        real(dp), intent(in) :: x(:)
        character(len=*), intent(in) :: model
        type(whittle_result), intent(out) :: result
        integer, intent(in), optional :: p, q, covariance_m
        real(dp), intent(in), optional :: start_eta(:), delta
        logical, intent(in), optional :: scale_data
        real(dp), allocatable :: work(:), full_per(:), eta0(:), eta_hat(:), covariance(:, :)
        real(dp) :: sd_x, fmin, hdelta
        integer :: n, pp, qq, m_cov, i
        logical :: do_scale, converged
        type(qeta_result) :: qres
        type(whittle_context) :: context

        n = size(x)
        if (n < 4) error stop "whittle_estimate: x must contain at least four observations"
        pp = 0
        qq = 0
        if (present(p)) pp = p
        if (present(q)) qq = q
        do_scale = .false.
        if (present(scale_data)) do_scale = scale_data
        m_cov = 10000
        if (present(covariance_m)) m_cov = covariance_m
        hdelta = 1.0e-9_dp
        if (present(delta)) hdelta = delta

        allocate(work(n))
        work = x
        if (do_scale) then
            sd_x = sample_sd(work)
            if (sd_x <= 0.0_dp) error stop "whittle_estimate: cannot scale a constant series"
            work = work/sd_x
        end if
        call periodogram(work, full_per)
        allocate(result%periodogram((n - 1)/2))
        result%periodogram = full_per(2:1 + (n - 1)/2)
        context%n = n
        context%p = pp
        context%q = qq
        allocate(context%periodogram(size(result%periodogram)))
        context%periodogram = result%periodogram

        select case (trim(adjustl(model)))
        case ("fGn", "fgn", "FGN")
            pp = 0
            qq = 0
            allocate(eta0(1))
            eta0 = 0.5_dp
            if (present(start_eta)) then
                if (size(start_eta) /= 1) error stop "whittle_estimate: fGn start_eta must have length one"
                eta0 = start_eta
            end if
            eta0(1) = max(0.2_dp, min(eta0(1), 0.9_dp))
            allocate(eta_hat(1))
            context%model = "fGn"
            context%p = 0
            context%q = 0
            call minimize_scalar(whittle_objective_scalar, context, 0.1_dp, 0.99_dp, eta_hat(1), fmin)
            converged = .true.
            result%model = "fGn"
        case ("fARIMA", "farima", "FARIMA")
            if (pp < 0 .or. qq < 0) error stop "whittle_estimate: p and q must be nonnegative"
            allocate(eta0(1 + pp + qq))
            eta0 = 0.0_dp
            eta0(1) = 0.5_dp
            if (present(start_eta)) then
                if (size(start_eta) /= size(eta0)) then
                    error stop "whittle_estimate: start_eta must have length 1+p+q"
                end if
                eta0 = start_eta
            end if
            eta0(1) = max(0.2_dp, min(eta0(1), 0.9_dp))
            context%model = "fARIMA"
            context%p = pp
            context%q = qq
            call nelder_mead(whittle_objective_vector, context, eta0, eta_hat, fmin, converged)
            result%model = "fARIMA"
        case default
            error stop "whittle_estimate: model must be fGn or fARIMA"
        end select

        call qeta(eta_hat, trim(result%model), n, result%periodogram, qres, p=pp, q=qq)
        if (trim(result%model) == "fGn") then
            call ceta_fgn(eta_hat, covariance, m=m_cov, delta=hdelta)
        else
            call ceta_arima(eta_hat, pp, qq, covariance, m=m_cov, delta=hdelta)
        end if
        covariance = covariance/real(n, dp)

        result%n = n
        result%p = pp
        result%q = qq
        result%converged = converged
        result%theta1 = qres%theta1
        allocate(result%eta(size(eta_hat)), result%covariance(size(eta_hat), size(eta_hat)))
        allocate(result%std_error(size(eta_hat)), result%z_value(size(eta_hat)), result%p_value(size(eta_hat)))
        allocate(result%spectrum(size(qres%spectrum)))
        result%eta = eta_hat
        result%covariance = covariance
        result%std_error = sqrt(max(0.0_dp, [(covariance(i, i), i=1,size(eta_hat))]))
        result%z_value = result%eta/result%std_error
        result%p_value = 2.0_dp*normal_sf(abs(result%z_value))
        result%spectrum = result%theta1*qres%spectrum
        result%goodness_of_fit = qres

    end subroutine whittle_estimate


    function whittle_objective_scalar(hurst, generic_context) result(value)
        real(dp), intent(in) :: hurst
        class(*), intent(in) :: generic_context
        real(dp) :: value
        real(dp) :: eta_one(1)

        select type (context => generic_context)
        type is (whittle_context)
            eta_one(1) = hurst
            call qeta(eta_one, "fGn", context%n, context%periodogram, p=0, q=0, &
                      give_b_only=.true., b_only=value)
        class default
            value = huge(1.0_dp)
        end select
        if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
    end function whittle_objective_scalar


    function whittle_objective_vector(eta_trial, generic_context) result(value)
        real(dp), intent(in) :: eta_trial(:)
        class(*), intent(in) :: generic_context
        real(dp) :: value

        select type (context => generic_context)
        type is (whittle_context)
            call qeta(eta_trial, "fARIMA", context%n, context%periodogram, &
                      p=context%p, q=context%q, give_b_only=.true., b_only=value)
        class default
            value = huge(1.0_dp)
        end select
        if (.not. ieee_is_finite(value)) value = huge(1.0_dp)
    end function whittle_objective_vector


    subroutine fexp_estimate(x, order_poly, pvalmax, result, max_iter, tolerance)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: order_poly
        real(dp), intent(in) :: pvalmax
        type(fexp_result), intent(out) :: result
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tolerance
        real(dp), allocatable :: full_per(:), yper(:), freq(:), xlong(:), design(:, :)
        real(dp), allocatable :: beta(:), covariance(:, :), std_error(:), t_value(:), p_value(:), fitted(:)
        real(dp), allocatable :: last_beta(:), last_covariance(:, :), last_std_error(:)
        real(dp), allocatable :: last_t_value(:), last_p_value(:), last_fitted(:)
        real(dp) :: max_p
        integer :: n, nfreq, j, selected, info, niter
        real(dp) :: tol
        logical :: early

        if (size(x) < 4) error stop "fexp_estimate: x must contain at least four observations"
        if (order_poly < 0) error stop "fexp_estimate: order_poly must be nonnegative"
        if (pvalmax < 0.0_dp .or. pvalmax > 1.0_dp) error stop "fexp_estimate: pvalmax must be in [0,1]"
        niter = 100
        if (present(max_iter)) niter = max_iter
        tol = 1.0e-10_dp
        if (present(tolerance)) tol = tolerance

        n = size(x)
        nfreq = (n - 1)/2
        call periodogram(x, full_per)
        allocate(yper(nfreq), freq(nfreq), xlong(nfreq))
        yper = full_per(2:1 + nfreq)
        call fourier_frequencies(n, freq)
        xlong = log(sqrt((1.0_dp - cos(freq))**2 + sin(freq)**2))

        early = .false.
        max_p = 0.0_dp
        selected = 0
        do j = 0, order_poly
            call make_fexp_design(freq, xlong, j, design)
            call fit_gamma_log(yper, design, beta, covariance, std_error, t_value, p_value, fitted, &
                               info, niter, tol)
            if (info /= 0) error stop "fexp_estimate: Gamma GLM failed"
            if (j >= 1) then
                max_p = maxval(p_value(2:))
                if (max_p > pvalmax) then
                    early = .true.
                    selected = j - 1
                    exit
                end if
            end if
            selected = j
            call move_alloc(beta, last_beta)
            call move_alloc(covariance, last_covariance)
            call move_alloc(std_error, last_std_error)
            call move_alloc(t_value, last_t_value)
            call move_alloc(p_value, last_p_value)
            call move_alloc(fitted, last_fitted)
        end do

        result%n = n
        result%order_poly = selected
        result%max_order_poly = order_poly
        result%early_stop = early
        result%max_p_value = max_p
        result%hurst = (1.0_dp - last_beta(2))/2.0_dp
        allocate(result%coefficients(size(last_beta), 4))
        result%coefficients(:, 1) = last_beta
        result%coefficients(:, 2) = last_std_error
        result%coefficients(:, 3) = last_t_value
        result%coefficients(:, 4) = last_p_value
        allocate(result%covariance(size(last_covariance, 1), size(last_covariance, 2)))
        result%covariance = last_covariance
        allocate(result%periodogram(size(yper)), result%spectrum(size(last_fitted)))
        result%periodogram = yper
        result%spectrum = last_fitted
    end subroutine fexp_estimate


    subroutine make_fexp_design(freq, xlong, degree, design)
        real(dp), intent(in) :: freq(:), xlong(:)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: design(:, :)
        real(dp), allocatable :: basis(:, :)
        integer :: n

        n = size(freq)
        allocate(design(n, 2 + degree))
        design(:, 1) = 1.0_dp
        design(:, 2) = xlong
        if (degree == 1) then
            design(:, 3) = freq
        else if (degree >= 2) then
            call orthogonal_polynomial_basis(freq, degree, basis)
            design(:, 3:) = basis
        end if
    end subroutine make_fexp_design


    subroutine orthogonal_polynomial_basis(x, degree, basis)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: basis(:, :)
        real(dp), allocatable :: q0(:), v(:)
        real(dp) :: norm_v
        integer :: k, j

        allocate(basis(size(x), degree), q0(size(x)), v(size(x)))
        q0 = 1.0_dp/sqrt(real(size(x), dp))
        do k = 1, degree
            v = x**k
            v = v - dot_product(q0, v)*q0
            do j = 1, k - 1
                v = v - dot_product(basis(:, j), v)*basis(:, j)
            end do
            norm_v = sqrt(dot_product(v, v))
            if (norm_v <= 100.0_dp*epsilon(1.0_dp)) then
                error stop "orthogonal_polynomial_basis: rank-deficient polynomial basis"
            end if
            basis(:, k) = v/norm_v
        end do
    end subroutine orthogonal_polynomial_basis


    subroutine fit_gamma_log(y, design, beta, covariance, std_error, t_value, p_value, fitted, &
                             info, max_iter, tolerance)
        real(dp), intent(in) :: y(:), design(:, :)
        real(dp), allocatable, intent(out) :: beta(:), covariance(:, :), std_error(:)
        real(dp), allocatable, intent(out) :: t_value(:), p_value(:), fitted(:)
        integer, intent(out) :: info
        integer, intent(in) :: max_iter
        real(dp), intent(in) :: tolerance
        real(dp), allocatable :: xtx(:, :), xtx_inv(:, :), z(:), eta(:), mu(:), beta_new(:)
        real(dp) :: dispersion, df
        integer :: n, k, iter, j

        n = size(y)
        k = size(design, 2)
        info = 0
        if (size(design, 1) /= n .or. n <= k .or. any(y <= 0.0_dp)) then
            info = -1
            allocate(beta(0), covariance(0, 0), std_error(0), t_value(0), p_value(0), fitted(0))
            return
        end if

        allocate(xtx(k, k))
        xtx = matmul(transpose(design), design)
        call inverse_spd(xtx, xtx_inv, info)
        if (info /= 0) then
            allocate(beta(0), covariance(0, 0), std_error(0), t_value(0), p_value(0), fitted(0))
            return
        end if

        allocate(beta(k), beta_new(k), z(n), eta(n), mu(n), fitted(n))
        beta = 0.0_dp
        beta(1) = log(sum(y)/real(n, dp))
        do iter = 1, max_iter
            eta = matmul(design, beta)
            eta = max(-700.0_dp, min(700.0_dp, eta))
            mu = exp(eta)
            z = eta + y/mu - 1.0_dp
            beta_new = matmul(xtx_inv, matmul(transpose(design), z))
            if (maxval(abs(beta_new - beta)) <= tolerance*(1.0_dp + maxval(abs(beta)))) then
                beta = beta_new
                exit
            end if
            beta = beta_new
        end do
        if (iter > max_iter) info = 1

        eta = matmul(design, beta)
        eta = max(-700.0_dp, min(700.0_dp, eta))
        fitted = exp(eta)
        df = real(n - k, dp)
        dispersion = sum(((y - fitted)/fitted)**2)/df
        allocate(covariance(k, k), std_error(k), t_value(k), p_value(k))
        covariance = dispersion*xtx_inv
        do j = 1, k
            std_error(j) = sqrt(max(0.0_dp, covariance(j, j)))
        end do
        t_value = beta/std_error
        do j = 1, k
            p_value(j) = 2.0_dp*(1.0_dp - student_t_cdf(abs(t_value(j)), df))
        end do
    end subroutine fit_gamma_log

end module longmemo
