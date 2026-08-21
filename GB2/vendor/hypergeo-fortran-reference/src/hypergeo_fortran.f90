! SPDX-License-Identifier: GPL-2.0-only
module hypergeo_fortran
    use hypergeo_kinds, only : dp, pi, ci
    use hypergeo_types, only : hypergeo_info
    use hypergeo_special, only : complex_gamma, complex_log_gamma, complex_factorial, complex_digamma, &
        is_near_integer, is_nonpos_integer, is_zero_parameter, finite_complex, nan_complex
    use hypergeo_generalized, only : genhypergeo, genhypergeo_series, genhypergeo_contfrac_single, &
        genhypergeo_contfrac, genhypergeo_shanks, shanks_transform
    use hypergeo_gauss_core
    use hypergeo_buhring_mod, only : lpham, buhring_eqn11, buhring_eqn12, buhring_eqn5_factors, &
        buhring_eqn5_series, hypergeo_buhring
    use hypergeo_ode, only : hypergeo_ode_continue, hypergeo_press, f15_5_1, &
        semicircle, semidash, straight, straightdash, to_real, to_complex, complex_path
    use elliptic_numeric, only : residue
    implicit none
    private

    complex(dp), save :: residue_a = (0.0_dp, 0.0_dp)
    complex(dp), save :: residue_b = (0.0_dp, 0.0_dp)
    complex(dp), save :: residue_c = (1.0_dp, 0.0_dp)
    real(dp), save :: residue_tol = 0.0_dp
    integer, save :: residue_maxiter = 2000

    public :: dp, pi, ci, hypergeo_info
    public :: hypergeo, hypergeo_scalar, hypergeo_array
    public :: complex_gamma, complex_log_gamma, complex_factorial, complex_digamma, lanczos
    public :: is_near_integer, is_nonpos_integer, is_zero_parameter
    public :: genhypergeo, genhypergeo_series, genhypergeo_contfrac_single, genhypergeo_contfrac
    public :: genhypergeo_shanks, hypergeo_shanks, shanks
    public :: hypergeo_contfrac, hypergeo_residue_general
    public :: hypergeo_residue_close_to_crit_single, hypergeo_residue_close_to_crit_multiple
    public :: lpham, buhring_eqn11, buhring_eqn12, buhring_eqn5_factors, buhring_eqn5_series, hypergeo_buhring
    public :: hypergeo_ode_continue, hypergeo_press, f15_5_1
    public :: semicircle, semidash, straight, straightdash, to_real, to_complex, complex_path

    ! Re-export the computational A&S/Wolfram methods from hypergeo_gauss_core.
    public :: hypergeo_core, hypergeo_powerseries, hypergeo_general, hypergeo_taylor
    public :: hypergeo_a_nonpos_int, hypergeo_aorb_nonpos_int, hypergeo_gosper
    public :: f15_1_1, f15_3_1, f15_3_3, f15_3_4, f15_3_5, f15_3_6, f15_3_7, f15_3_8, f15_3_9
    public :: i15_3_6, i15_3_7, i15_3_8, i15_3_9, j15_3_6, j15_3_7, j15_3_8, j15_3_9
    public :: thingfun, crit_points, hypergeo_cover1, hypergeo_cover2, hypergeo_cover3
    public :: f15_3_10, f15_3_10_a, f15_3_10_b
    public :: f15_3_11, f15_3_11_bit1, f15_3_11_bit2_a, f15_3_11_bit2_b
    public :: f15_3_12, f15_3_12_bit1, f15_3_12_bit2_a, f15_3_12_bit2_b
    public :: f15_3_13, f15_3_13_a, f15_3_13_b
    public :: f15_3_14, f15_3_14_bit1_a, f15_3_14_bit1_b, f15_3_14_bit2
    public :: f15_3_10_11_12, f15_3_13_14
    public :: w07_23_06_0029_01, w07_23_06_0031_01, w07_23_06_0026_01
    public :: w07_23_06_0031_01_bit1, w07_23_06_0031_01_bit2
    public :: w07_23_06_0026_01_bit1, w07_23_06_0026_01_bit2
    public :: w07_23_06_0026_01_bit3_a, w07_23_06_0026_01_bit3_b, w07_23_06_0026_01_bit3_c

    interface hypergeo
        module procedure hypergeo_scalar
        module procedure hypergeo_array
    end interface hypergeo

contains

    function hypergeo_scalar(a, b, c, z, tol, maxiter, info) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        type(hypergeo_info), intent(out), optional :: info
        complex(dp) :: value
        real(dp) :: t
        integer :: nmax

        t = 0.0_dp; if (present(tol)) t = tol
        nmax = 2000; if (present(maxiter)) nmax = maxiter

        if (minval(abs(z - crit_points())) < 0.1_dp) then
            value = hypergeo_gosper(a, b, c, z, tol=t, maxiter=nmax)
            if (finite_complex(value)) then
                call set_info(info, .true., 0, 'gosper', 1)
                return
            end if
        else
            value = hypergeo_core(a, b, c, z, tol=t, maxiter=nmax)
            if (finite_complex(value)) then
                call set_info(info, .true., 0, 'analytic', 1)
                return
            end if
        end if

        value = hypergeo_contfrac(a, b, c, z, tol=t, maxiter=nmax)
        if (finite_complex(value)) then
            call set_info(info, .true., 0, 'continued-fraction', 2)
            return
        end if

        if (real(b, dp) > 0.0_dp .and. real(c - b, dp) > 0.0_dp) then
            value = f15_3_1(a, b, c, z)
            if (finite_complex(value)) then
                call set_info(info, .true., 0, 'euler-integral', 3)
                return
            end if
        end if

        value = hypergeo_ode_continue(a, b, c, z)
        if (finite_complex(value)) then
            call set_info(info, .true., 0, 'ode', 4)
            return
        end if
        value = nan_complex()
        call set_info(info, .false., nmax, 'failed', -1)
    end function hypergeo_scalar

    function hypergeo_array(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c
        complex(dp), intent(in) :: z(:)
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value(size(z))
        integer :: i
        do i = 1, size(z)
            value(i) = hypergeo_scalar(a, b, c, z(i), tol=tol, maxiter=maxiter)
        end do
    end function hypergeo_array

    function hypergeo_contfrac(a, b, c, z, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        value = genhypergeo_contfrac_single([a, b], [c], z, tol=tol, maxiter=maxiter)
    end function hypergeo_contfrac

    function shanks(last, this, next) result(value)
        complex(dp), intent(in) :: last, this, next
        complex(dp) :: value
        value = shanks_transform(last, this, next)
    end function shanks

    function hypergeo_shanks(a, b, c, z, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        integer, intent(in), optional :: maxiter
        complex(dp) :: value
        value = genhypergeo_shanks([a, b], [c], z, maxiter=maxiter)
    end function hypergeo_shanks

    function lanczos(z, log_value) result(value)
        complex(dp), intent(in) :: z
        logical, intent(in), optional :: log_value
        complex(dp) :: value
        logical :: lg
        lg = .false.; if (present(log_value)) lg = log_value
        if (lg) then
            value = complex_log_gamma(z)
        else
            value = complex_gamma(z)
        end if
    end function lanczos

    function hypergeo_residue_general(a, b, c, z, r, center, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        real(dp), intent(in), optional :: r, tol
        complex(dp), intent(in), optional :: center
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, o
        real(dp) :: rr

        rr = 0.15_dp
        if (present(r)) rr = r
        o = z
        if (present(center)) o = center
        residue_a = a
        residue_b = b
        residue_c = c
        residue_tol = 0.0_dp
        if (present(tol)) residue_tol = tol
        residue_maxiter = 2000
        if (present(maxiter)) residue_maxiter = maxiter
        value = residue(residue_integrand, z, rr, center=o)
    end function hypergeo_residue_general

    function residue_integrand(w) result(value)
        complex(dp), intent(in) :: w
        complex(dp) :: value
        value = hypergeo_scalar(residue_a, residue_b, residue_c, w, &
            tol=residue_tol, maxiter=residue_maxiter)
    end function residue_integrand

    function hypergeo_residue_close_to_crit_single(a, b, c, z, strategy, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z
        character(len=*), intent(in), optional :: strategy
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value, cp(2), critical, center
        character(len=1) :: s
        cp = crit_points()
        if (abs(z - cp(1)) <= 0.1_dp) then
            critical = cp(1)
        else if (abs(z - cp(2)) <= 0.1_dp) then
            critical = cp(2)
        else
            value = nan_complex(); return
        end if
        s = 'A'; if (present(strategy)) s = strategy(1:1)
        if (s == 'A' .or. s == 'a') then
            center = critical
        else
            center = z
        end if
        value = hypergeo_residue_general(a, b, c, z, r=0.15_dp, center=center, tol=tol, maxiter=maxiter)
    end function hypergeo_residue_close_to_crit_single

    function hypergeo_residue_close_to_crit_multiple(a, b, c, z, strategy, tol, maxiter) result(value)
        complex(dp), intent(in) :: a, b, c, z(:)
        character(len=*), intent(in), optional :: strategy
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        complex(dp) :: value(size(z))
        integer :: i
        do i = 1, size(z)
            value(i) = hypergeo_residue_close_to_crit_single(a, b, c, z(i), strategy, tol, maxiter)
        end do
    end function hypergeo_residue_close_to_crit_multiple

    subroutine set_info(info, ok, iterations, method, code)
        type(hypergeo_info), intent(out), optional :: info
        logical, intent(in) :: ok
        integer, intent(in) :: iterations, code
        character(len=*), intent(in) :: method
        if (.not. present(info)) return
        info%converged = ok
        info%iterations = iterations
        info%method = method
        info%method_code = code
        info%residual = 0.0_dp
    end subroutine set_info

end module hypergeo_fortran
