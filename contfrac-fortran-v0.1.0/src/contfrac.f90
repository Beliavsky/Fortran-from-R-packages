! SPDX-License-Identifier: GPL-2.0-only
module contfrac
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use contfrac_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: tiny_cf = 1.0e-30_dp
    real(dp), parameter :: lentz_eps = 2.22044604925031e-16_dp

    type, public :: contfrac_info
        logical :: converged = .false.
        integer :: iterations = 0
        real(dp) :: residual = huge(1.0_dp)
        real(dp) :: tolerance = epsilon(1.0_dp)
    end type contfrac_info

    public :: dp
    public :: cf, gcf
    public :: convergents, gconvergents
    public :: as_cf

    interface cf
        module procedure cf_real
        module procedure cf_complex
    end interface cf

    interface gcf
        module procedure gcf_real
        module procedure gcf_complex
    end interface gcf

    interface convergents
        module procedure convergents_real
        module procedure convergents_complex
    end interface convergents

    interface gconvergents
        module procedure gconvergents_real
        module procedure gconvergents_complex
    end interface gconvergents

contains

    function gcf_real(a, b, b0, finite, tol, info) result(value)
        real(dp), intent(in) :: a(:), b(:)
        real(dp), intent(in), optional :: b0, tol
        logical, intent(in), optional :: finite
        type(contfrac_info), intent(out), optional :: info
        real(dp) :: value

        real(dp) :: c, d, delta, f, b0_use, tol_use, residual
        logical :: finite_use, internal_converged, accepted
        integer :: j, n, iterations

        n = size(a)
        if (n < 1) error stop "gcf: a and b must be nonempty"
        if (size(b) /= n) error stop "gcf: a and b must have the same size"

        b0_use = 0.0_dp
        if (present(b0)) b0_use = b0
        finite_use = .false.
        if (present(finite)) finite_use = finite
        tol_use = epsilon(1.0_dp)
        if (present(tol)) then
            if (tol > 0.0_dp) tol_use = tol
        end if

        f = tiny_cf
        c = f
        d = 0.0_dp
        delta = 0.0_dp
        internal_converged = .false.
        iterations = 0

        do j = 1, n
            d = b(j) + a(j) * d
            if (.not. (abs(d) > 0.0_dp)) d = tiny_cf

            c = b(j) + a(j) / c
            if (.not. (abs(c) > 0.0_dp)) c = tiny_cf

            d = 1.0_dp / d
            delta = c * d
            f = f * delta
            iterations = j

            if (abs(delta - 1.0_dp) <= lentz_eps) then
                internal_converged = .true.
                exit
            end if
        end do

        residual = abs(delta - 1.0_dp)
        accepted = residual <= tol_use .or. finite_use
        if (accepted) then
            value = b0_use + f
        else
            value = ieee_value(value, ieee_quiet_nan)
        end if

        if (present(info)) then
            info%converged = internal_converged .or. residual <= tol_use
            info%iterations = iterations
            info%residual = residual
            info%tolerance = tol_use
        end if
    end function gcf_real

    function gcf_complex(a, b, b0, finite, tol, info) result(value)
        complex(dp), intent(in) :: a(:), b(:)
        complex(dp), intent(in), optional :: b0
        logical, intent(in), optional :: finite
        real(dp), intent(in), optional :: tol
        type(contfrac_info), intent(out), optional :: info
        complex(dp) :: value

        complex(dp) :: c, d, delta, f, b0_use
        real(dp) :: tol_use, residual
        logical :: finite_use, internal_converged, accepted
        integer :: j, n, iterations

        n = size(a)
        if (n < 1) error stop "gcf: a and b must be nonempty"
        if (size(b) /= n) error stop "gcf: a and b must have the same size"

        b0_use = cmplx(0.0_dp, 0.0_dp, kind=dp)
        if (present(b0)) b0_use = b0
        finite_use = .false.
        if (present(finite)) finite_use = finite
        tol_use = epsilon(1.0_dp)
        if (present(tol)) then
            if (tol > 0.0_dp) tol_use = tol
        end if

        f = cmplx(tiny_cf, 0.0_dp, kind=dp)
        c = f
        d = cmplx(0.0_dp, 0.0_dp, kind=dp)
        delta = cmplx(0.0_dp, 0.0_dp, kind=dp)
        internal_converged = .false.
        iterations = 0

        do j = 1, n
            d = b(j) + a(j) * d
            if (.not. (abs(d) > 0.0_dp)) d = cmplx(tiny_cf, 0.0_dp, kind=dp)

            c = b(j) + a(j) / c
            if (.not. (abs(c) > 0.0_dp)) c = cmplx(tiny_cf, 0.0_dp, kind=dp)

            d = 1.0_dp / d
            delta = c * d
            f = f * delta
            iterations = j

            if (abs(real(delta, dp) - 1.0_dp) <= lentz_eps .and. &
                abs(aimag(delta)) <= lentz_eps) then
                internal_converged = .true.
                exit
            end if
        end do

        residual = abs(delta - cmplx(1.0_dp, 0.0_dp, kind=dp))
        accepted = residual <= tol_use .or. finite_use
        if (accepted) then
            value = b0_use + f
        else
            value = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), &
                          ieee_value(0.0_dp, ieee_quiet_nan), kind=dp)
        end if

        if (present(info)) then
            info%converged = internal_converged .or. residual <= tol_use
            info%iterations = iterations
            info%residual = residual
            info%tolerance = tol_use
        end if
    end function gcf_complex

    function cf_real(a, finite, tol, info) result(value)
        real(dp), intent(in) :: a(:)
        logical, intent(in), optional :: finite
        real(dp), intent(in), optional :: tol
        type(contfrac_info), intent(out), optional :: info
        real(dp) :: value

        real(dp), allocatable :: num(:), den(:), work(:)
        real(dp) :: tol_use
        logical :: finite_use
        integer :: first_bad, j, n
        type(contfrac_info) :: local_info

        if (size(a) < 2) error stop "cf: at least two coefficients are required"

        first_bad = 0
        do j = 1, size(a)
            if (.not. ieee_is_finite(a(j))) then
                first_bad = j
                exit
            end if
        end do

        finite_use = .false.
        if (present(finite)) finite_use = finite
        if (first_bad > 0) then
            if (first_bad <= 1) error stop "cf: first coefficient may not be infinite"
            n = first_bad - 1
            finite_use = .true.
        else
            n = size(a)
        end if
        if (n < 2) error stop "cf: finite continued fraction must contain at least two coefficients"

        allocate(work(n))
        work = a(1:n)
        allocate(num(n - 1), den(n - 1))
        num = 1.0_dp
        den = work(2:n)

        tol_use = epsilon(1.0_dp)
        if (present(tol)) then
            if (tol > 0.0_dp) tol_use = tol
        end if

        value = gcf_real(num, den, b0=work(1), finite=finite_use, tol=tol_use, info=local_info)
        if (present(info)) info = local_info
    end function cf_real

    function cf_complex(a, finite, tol, info) result(value)
        complex(dp), intent(in) :: a(:)
        logical, intent(in), optional :: finite
        real(dp), intent(in), optional :: tol
        type(contfrac_info), intent(out), optional :: info
        complex(dp) :: value

        complex(dp), allocatable :: num(:), den(:), work(:)
        real(dp) :: tol_use
        logical :: finite_use
        integer :: first_bad, j, n
        type(contfrac_info) :: local_info

        if (size(a) < 2) error stop "cf: at least two coefficients are required"

        first_bad = 0
        do j = 1, size(a)
            if (.not. ieee_is_finite(real(a(j), dp)) .or. .not. ieee_is_finite(aimag(a(j)))) then
                first_bad = j
                exit
            end if
        end do

        finite_use = .false.
        if (present(finite)) finite_use = finite
        if (first_bad > 0) then
            if (first_bad <= 1) error stop "cf: first coefficient may not be infinite"
            n = first_bad - 1
            finite_use = .true.
        else
            n = size(a)
        end if
        if (n < 2) error stop "cf: finite continued fraction must contain at least two coefficients"

        allocate(work(n))
        work = a(1:n)
        allocate(num(n - 1), den(n - 1))
        num = cmplx(1.0_dp, 0.0_dp, kind=dp)
        den = work(2:n)

        tol_use = epsilon(1.0_dp)
        if (present(tol)) then
            if (tol > 0.0_dp) tol_use = tol
        end if

        value = gcf_complex(num, den, b0=work(1), finite=finite_use, tol=tol_use, info=local_info)
        if (present(info)) info = local_info
    end function cf_complex

    subroutine gconvergents_real(a, b, b0, numerators, denominators)
        real(dp), intent(in) :: a(:), b(:), b0
        real(dp), allocatable, intent(out) :: numerators(:), denominators(:)

        integer :: j, n

        n = size(a)
        if (n < 1) error stop "gconvergents: a and b must be nonempty"
        if (size(b) /= n) error stop "gconvergents: a and b must have the same size"

        allocate(numerators(n + 1), denominators(n + 1))
        numerators(1) = b0
        denominators(1) = 1.0_dp
        numerators(2) = b(1) * numerators(1) + a(1)
        denominators(2) = b(1) * denominators(1)

        do j = 3, n + 1
            numerators(j) = b(j - 1) * numerators(j - 1) + a(j - 1) * numerators(j - 2)
            denominators(j) = b(j - 1) * denominators(j - 1) + a(j - 1) * denominators(j - 2)
        end do
    end subroutine gconvergents_real

    subroutine gconvergents_complex(a, b, b0, numerators, denominators)
        complex(dp), intent(in) :: a(:), b(:), b0
        complex(dp), allocatable, intent(out) :: numerators(:), denominators(:)

        integer :: j, n

        n = size(a)
        if (n < 1) error stop "gconvergents: a and b must be nonempty"
        if (size(b) /= n) error stop "gconvergents: a and b must have the same size"

        allocate(numerators(n + 1), denominators(n + 1))
        numerators(1) = b0
        denominators(1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        numerators(2) = b(1) * numerators(1) + a(1)
        denominators(2) = b(1) * denominators(1)

        do j = 3, n + 1
            numerators(j) = b(j - 1) * numerators(j - 1) + a(j - 1) * numerators(j - 2)
            denominators(j) = b(j - 1) * denominators(j - 1) + a(j - 1) * denominators(j - 2)
        end do
    end subroutine gconvergents_complex

    subroutine convergents_real(a, numerators, denominators)
        real(dp), intent(in) :: a(:)
        real(dp), allocatable, intent(out) :: numerators(:), denominators(:)

        real(dp), allocatable :: num(:)

        if (size(a) < 2) error stop "convergents: at least two coefficients are required"
        allocate(num(size(a) - 1))
        num = 1.0_dp
        call gconvergents_real(num, a(2:), a(1), numerators, denominators)
    end subroutine convergents_real

    subroutine convergents_complex(a, numerators, denominators)
        complex(dp), intent(in) :: a(:)
        complex(dp), allocatable, intent(out) :: numerators(:), denominators(:)

        complex(dp), allocatable :: num(:)

        if (size(a) < 2) error stop "convergents: at least two coefficients are required"
        allocate(num(size(a) - 1))
        num = cmplx(1.0_dp, 0.0_dp, kind=dp)
        call gconvergents_complex(num, a(2:), a(1), numerators, denominators)
    end subroutine convergents_complex

    function as_cf(x, n, n_used) result(terms)
        real(dp), intent(in) :: x
        integer, intent(in), optional :: n
        integer, intent(out), optional :: n_used
        real(dp), allocatable :: terms(:)

        real(dp) :: frac, work, whole, nan_value
        integer :: i, n_use, used

        if (.not. ieee_is_finite(x)) error stop "as_cf: x must be finite"
        n_use = 10
        if (present(n)) n_use = n
        if (n_use < 1) error stop "as_cf: n must be positive"

        allocate(terms(n_use))
        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        work = x
        used = 0

        do i = 1, n_use
            whole = work - modulo(work, 1.0_dp)
            terms(i) = whole
            used = i
            frac = work - whole
            if (.not. (abs(frac) > 0.0_dp)) then
                if (i < n_use) terms(i + 1:n_use) = nan_value
                exit
            end if
            work = 1.0_dp / frac
        end do

        if (present(n_used)) n_used = used
    end function as_cf

end module contfrac
