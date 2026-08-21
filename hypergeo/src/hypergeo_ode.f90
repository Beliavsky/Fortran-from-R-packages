! SPDX-License-Identifier: GPL-2.0-only
module hypergeo_ode
    use hypergeo_kinds, only : dp, pi, ci
    use hypergeo_special, only : nan_complex, finite_complex
    use hypergeo_gauss_core, only : hypergeo_core
    use desolve, only : ode, ode_result
    implicit none
    private

    abstract interface
        function complex_path(t) result(z)
            import dp
            real(dp), intent(in) :: t
            complex(dp) :: z
        end function complex_path
    end interface

    complex(dp), save :: ctx_a, ctx_b, ctx_c, ctx_start, ctx_target
    procedure(complex_path), pointer, save :: ctx_u => null(), ctx_udash => null()

    public :: hypergeo_ode_continue, hypergeo_press, f15_5_1
    public :: semicircle, semidash, straight, straightdash
    public :: to_real, to_complex, complex_path

contains

    function to_real(z) result(r)
        complex(dp), intent(in) :: z(:)
        real(dp) :: r(2 * size(z))
        integer :: i
        do i = 1, size(z)
            r(2 * i - 1) = real(z(i), dp)
            r(2 * i) = aimag(z(i))
        end do
    end function to_real

    function to_complex(r) result(z)
        real(dp), intent(in) :: r(:)
        complex(dp) :: z(size(r) / 2)
        integer :: i
        if (mod(size(r), 2) /= 0) error stop 'to_complex: real vector length must be even'
        do i = 1, size(z)
            z(i) = cmplx(r(2 * i - 1), r(2 * i), kind=dp)
        end do
    end function to_complex

    function hypergeo_ode_continue(a, b, c, z, method, rtol, atol) result(value)
        complex(dp), intent(in) :: a, b, c, z
        character(len=*), intent(in), optional :: method
        real(dp), intent(in), optional :: rtol, atol
        complex(dp) :: value, startz, f0, fp0
        real(dp) :: times(2), y0(4), rt, at
        type(ode_result) :: sol
        character(len=16) :: meth

        if (real(z, dp) <= 0.0_dp) then
            startz = cmplx(-0.5_dp, 0.0_dp, kind=dp)
        else if (real(z, dp) <= 0.5_dp) then
            startz = cmplx(0.5_dp, 0.0_dp, kind=dp)
        else if (aimag(z) >= 0.0_dp) then
            startz = cmplx(0.0_dp, 0.5_dp, kind=dp)
        else
            startz = cmplx(0.0_dp, -0.5_dp, kind=dp)
        end if
        f0 = hypergeo_core(a, b, c, startz)
        fp0 = (a * b / c) * hypergeo_core(a + 1.0_dp, b + 1.0_dp, c + 1.0_dp, startz)
        if (.not. finite_complex(f0) .or. .not. finite_complex(fp0)) then
            value = nan_complex()
            return
        end if
        ctx_a = a; ctx_b = b; ctx_c = c
        ctx_start = startz; ctx_target = z
        nullify(ctx_u, ctx_udash)
        y0 = to_real([f0, fp0])
        times = [0.0_dp, 1.0_dp]
        meth = 'lsoda'
        if (present(method)) meth = method
        rt = 1.0e-10_dp; at = 1.0e-12_dp
        if (present(rtol)) rt = rtol
        if (present(atol)) at = atol
        sol = ode(hypergeo_rhs, y0, times, method=trim(meth), rtol=rt, atol=at)
        if (sol%status < 0) then
            value = nan_complex()
        else
            value = cmplx(sol%y(1, size(times)), sol%y(2, size(times)), kind=dp)
        end if
    end function hypergeo_ode_continue

    function hypergeo_press(a, b, c, z, method, rtol, atol) result(value)
        complex(dp), intent(in) :: a, b, c, z
        character(len=*), intent(in), optional :: method
        real(dp), intent(in), optional :: rtol, atol
        complex(dp) :: value
        value = hypergeo_ode_continue(a, b, c, z, method=method, rtol=rtol, atol=atol)
    end function hypergeo_press

    function f15_5_1(a, b, c, startz, u, udash, method, rtol, atol) result(value)
        complex(dp), intent(in) :: a, b, c, startz
        procedure(complex_path) :: u, udash
        character(len=*), intent(in), optional :: method
        real(dp), intent(in), optional :: rtol, atol
        complex(dp) :: value, f0, fp0
        real(dp) :: times(2), y0(4), rt, at
        type(ode_result) :: sol
        character(len=16) :: meth

        f0 = hypergeo_core(a, b, c, startz)
        fp0 = (a * b / c) * hypergeo_core(a + 1.0_dp, b + 1.0_dp, c + 1.0_dp, startz)
        ctx_a = a; ctx_b = b; ctx_c = c
        ctx_start = startz; ctx_target = u(1.0_dp)
        ctx_u => u; ctx_udash => udash
        y0 = to_real([f0, fp0])
        times = [0.0_dp, 1.0_dp]
        meth = 'lsoda'; if (present(method)) meth = method
        rt = 1.0e-10_dp; at = 1.0e-12_dp
        if (present(rtol)) rt = rtol
        if (present(atol)) at = atol
        sol = ode(hypergeo_rhs, y0, times, method=trim(meth), rtol=rt, atol=at)
        nullify(ctx_u, ctx_udash)
        if (sol%status < 0) then
            value = nan_complex()
        else
            value = cmplx(sol%y(1, 2), sol%y(2, 2), kind=dp)
        end if
    end function f15_5_1

    subroutine hypergeo_rhs(t, y, dydt)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        complex(dp) :: state(2), deriv(2), z, dz
        state = to_complex(y)
        if (associated(ctx_u)) then
            z = ctx_u(t)
            dz = ctx_udash(t)
        else
            z = ctx_start + (ctx_target - ctx_start) * t
            dz = ctx_target - ctx_start
        end if
        deriv(1) = dz * state(2)
        deriv(2) = dz * (ctx_a * ctx_b * state(1) &
            - (ctx_c - (ctx_a + ctx_b + 1.0_dp) * z) * state(2)) / (z * (1.0_dp - z))
        dydt = to_real(deriv)
    end subroutine hypergeo_rhs

    pure function semicircle(t, z0, z1, clockwise) result(z)
        real(dp), intent(in) :: t
        complex(dp), intent(in) :: z0, z1
        logical, intent(in), optional :: clockwise
        complex(dp) :: z, center
        real(dp) :: m
        m = -1.0_dp
        if (present(clockwise)) then
            if (.not. clockwise) m = 1.0_dp
        end if
        center = (z0 + z1) / 2.0_dp
        z = center + (z0 - center) * exp(ci * t * pi * m)
    end function semicircle

    pure function semidash(t, z0, z1, clockwise) result(z)
        real(dp), intent(in) :: t
        complex(dp), intent(in) :: z0, z1
        logical, intent(in), optional :: clockwise
        complex(dp) :: z, center
        real(dp) :: m
        m = -1.0_dp
        if (present(clockwise)) then
            if (.not. clockwise) m = 1.0_dp
        end if
        center = (z0 + z1) / 2.0_dp
        z = (z0 - center) * (ci * pi * m) * exp(ci * t * pi * m)
    end function semidash

    pure function straight(t, z0, z1) result(z)
        real(dp), intent(in) :: t
        complex(dp), intent(in) :: z0, z1
        complex(dp) :: z
        z = z0 + t * (z1 - z0)
    end function straight

    pure function straightdash(t, z0, z1) result(z)
        real(dp), intent(in) :: t
        complex(dp), intent(in) :: z0, z1
        complex(dp) :: z
        if (t < -huge(1.0_dp)) continue
        z = z1 - z0
    end function straightdash

end module hypergeo_ode
