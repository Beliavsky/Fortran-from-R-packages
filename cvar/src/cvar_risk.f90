! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of cvar 0.6 by Georgi N. Boshnakov.
module cvar_risk
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use cvar_kinds, only : dp
    use cvar_status, only : cvar_ok, cvar_invalid_probability, cvar_invalid_scale, &
                            cvar_invalid_sample, cvar_bracket_failure, cvar_nonconvergence
    implicit none
    private

    abstract interface
        function cvar_scalar_function(x) result(value)
            import dp
            real(dp), intent(in) :: x
            real(dp) :: value
        end function cvar_scalar_function
    end interface

    interface var_qf
        module procedure var_qf_scalar
        module procedure var_qf_vector
    end interface var_qf

    interface var_cdf
        module procedure var_cdf_scalar
        module procedure var_cdf_vector
    end interface var_cdf

    interface var_sample
        module procedure var_sample_scalar
        module procedure var_sample_vector
        module procedure var_sample_matrix_scalar
        module procedure var_sample_matrix_vector
    end interface var_sample

    interface es_qf
        module procedure es_qf_scalar
        module procedure es_qf_vector
    end interface es_qf

    interface es_cdf
        module procedure es_cdf_scalar
        module procedure es_cdf_vector
    end interface es_cdf

    interface es_pdf
        module procedure es_pdf_scalar
        module procedure es_pdf_vector
    end interface es_pdf

    interface es_sample
        module procedure es_sample_scalar
        module procedure es_sample_vector
        module procedure es_sample_matrix_scalar
        module procedure es_sample_matrix_vector
    end interface es_sample

    public :: cvar_scalar_function
    public :: var_qf, var_cdf, var_sample
    public :: es_qf, es_cdf, es_pdf, es_sample
    public :: empirical_quantile, quantile_from_cdf

    real(dp), parameter :: default_tol = sqrt(epsilon(1.0_dp))
    integer, parameter :: max_adapt_depth = 20

contains

    function var_qf_scalar(qf, p_loss, intercept, slope, transf, status) result(value)
        procedure(cvar_scalar_function) :: qf
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp) :: value
        real(dp) :: a, b, q
        logical :: do_transf
        integer :: istat

        call risk_arguments(p_loss, intercept, slope, transf, a, b, do_transf, istat)
        if (istat /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = istat
            return
        end if

        q = a + b * qf(p_loss)
        if (do_transf) then
            value = 1.0_dp - exp(q)
        else
            value = -q
        end if
        if (.not. ieee_is_finite(value)) istat = cvar_nonconvergence
        if (present(status)) status = istat
    end function var_qf_scalar

    function var_qf_vector(qf, p_loss, intercept, slope, transf, status) result(values)
        procedure(cvar_scalar_function) :: qf
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: i, istat, one_status

        allocate(values(size(p_loss)))
        istat = cvar_ok
        do i = 1, size(p_loss)
            values(i) = var_qf_scalar(qf, p_loss(i), intercept, slope, transf, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function var_qf_vector

    function var_cdf_scalar(cdf, p_loss, intercept, slope, transf, tol, status) result(value)
        procedure(cvar_scalar_function) :: cdf
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope, tol
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp) :: value
        real(dp) :: a, b, q, use_tol
        logical :: do_transf
        integer :: istat

        call risk_arguments(p_loss, intercept, slope, transf, a, b, do_transf, istat)
        if (istat /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = istat
            return
        end if
        use_tol = default_tol
        if (present(tol)) use_tol = max(tol, epsilon(1.0_dp))
        q = quantile_from_cdf(cdf, p_loss, use_tol, istat)
        if (istat == cvar_ok) then
            q = a + b * q
            if (do_transf) then
                value = 1.0_dp - exp(q)
            else
                value = -q
            end if
        else
            value = ieee_value(1.0_dp, ieee_quiet_nan)
        end if
        if (present(status)) status = istat
    end function var_cdf_scalar

    function var_cdf_vector(cdf, p_loss, intercept, slope, transf, tol, status) result(values)
        procedure(cvar_scalar_function) :: cdf
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope, tol
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: i, istat, one_status

        allocate(values(size(p_loss)))
        istat = cvar_ok
        do i = 1, size(p_loss)
            values(i) = var_cdf_scalar(cdf, p_loss(i), intercept, slope, transf, tol, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function var_cdf_vector

    function var_sample_scalar(sample, p_loss, intercept, slope, transf, status) result(value)
        real(dp), intent(in) :: sample(:)
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp) :: value
        real(dp), allocatable :: clean(:)
        real(dp) :: a, b, q
        logical :: do_transf
        integer :: istat

        call risk_arguments(p_loss, intercept, slope, transf, a, b, do_transf, istat)
        if (istat == cvar_ok) call finite_sample(sample, clean, istat)
        if (istat /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = istat
            return
        end if

        q = a + b * empirical_quantile(clean, p_loss)
        if (do_transf) then
            value = 1.0_dp - exp(q)
        else
            value = -q
        end if
        if (present(status)) status = cvar_ok
    end function var_sample_scalar

    function var_sample_vector(sample, p_loss, intercept, slope, transf, status) result(values)
        real(dp), intent(in) :: sample(:)
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: i, istat, one_status

        allocate(values(size(p_loss)))
        istat = cvar_ok
        do i = 1, size(p_loss)
            values(i) = var_sample_scalar(sample, p_loss(i), intercept, slope, transf, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function var_sample_vector

    function var_sample_matrix_scalar(sample, p_loss, intercept, slope, transf, status) result(values)
        real(dp), intent(in) :: sample(:, :)
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: j, istat, one_status

        allocate(values(size(sample, 2)))
        istat = cvar_ok
        do j = 1, size(sample, 2)
            values(j) = var_sample_scalar(sample(:, j), p_loss, intercept, slope, &
                                          transf, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function var_sample_matrix_scalar

    function var_sample_matrix_vector(sample, p_loss, intercept, slope, transf, status) result(values)
        real(dp), intent(in) :: sample(:, :)
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:, :)
        integer :: j, istat, one_status

        allocate(values(size(p_loss), size(sample, 2)))
        istat = cvar_ok
        do j = 1, size(sample, 2)
            values(:, j) = var_sample_vector(sample(:, j), p_loss, intercept, slope, &
                                             transf, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function var_sample_matrix_vector

    function es_qf_scalar(qf, p_loss, intercept, slope, transf, tol, status) result(value)
        procedure(cvar_scalar_function) :: qf
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope, tol
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp) :: value
        real(dp) :: a, b, use_tol, error
        logical :: do_transf
        integer :: istat

        call risk_arguments(p_loss, intercept, slope, transf, a, b, do_transf, istat)
        if (istat /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = istat
            return
        end if
        use_tol = 1.0e-9_dp
        if (present(tol)) use_tol = max(tol, 100.0_dp * epsilon(1.0_dp))
        call adaptive_qf(qf, 0.0_dp, 1.0_dp, p_loss, a, b, do_transf, &
                         use_tol, 0, value, error, istat)
        if (present(status)) status = istat
    end function es_qf_scalar

    function es_qf_vector(qf, p_loss, intercept, slope, transf, tol, status) result(values)
        procedure(cvar_scalar_function) :: qf
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope, tol
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: i, istat, one_status

        allocate(values(size(p_loss)))
        istat = cvar_ok
        do i = 1, size(p_loss)
            values(i) = es_qf_scalar(qf, p_loss(i), intercept, slope, transf, tol, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function es_qf_vector

    function es_cdf_scalar(cdf, p_loss, intercept, slope, transf, tol, status) result(value)
        procedure(cvar_scalar_function) :: cdf
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope, tol
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp) :: value
        real(dp) :: a, b, use_tol, error
        logical :: do_transf
        integer :: istat

        call risk_arguments(p_loss, intercept, slope, transf, a, b, do_transf, istat)
        if (istat /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = istat
            return
        end if
        use_tol = 1.0e-8_dp
        if (present(tol)) use_tol = max(tol, 100.0_dp * epsilon(1.0_dp))
        call adaptive_cdf(cdf, 0.0_dp, 1.0_dp, p_loss, a, b, do_transf, &
                          use_tol, 0, value, error, istat)
        if (present(status)) status = istat
    end function es_cdf_scalar

    function es_cdf_vector(cdf, p_loss, intercept, slope, transf, tol, status) result(values)
        procedure(cvar_scalar_function) :: cdf
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope, tol
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: i, istat, one_status

        allocate(values(size(p_loss)))
        istat = cvar_ok
        do i = 1, size(p_loss)
            values(i) = es_cdf_scalar(cdf, p_loss(i), intercept, slope, transf, tol, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function es_cdf_vector

    function es_pdf_scalar(pdf, qf, p_loss, intercept, slope, transf, tol, status) result(value)
        procedure(cvar_scalar_function) :: pdf, qf
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope, tol
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp) :: value
        real(dp) :: a, b, use_tol, error, upper
        logical :: do_transf
        integer :: istat

        call risk_arguments(p_loss, intercept, slope, transf, a, b, do_transf, istat)
        if (istat /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = istat
            return
        end if
        use_tol = 1.0e-9_dp
        if (present(tol)) use_tol = max(tol, 100.0_dp * epsilon(1.0_dp))
        upper = qf(p_loss)
        call adaptive_pdf(pdf, 0.0_dp, 1.0_dp, upper, p_loss, a, b, do_transf, &
                           use_tol, 0, value, error, istat)
        if (present(status)) status = istat
    end function es_pdf_scalar

    function es_pdf_vector(pdf, qf, p_loss, intercept, slope, transf, tol, status) result(values)
        procedure(cvar_scalar_function) :: pdf, qf
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope, tol
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: i, istat, one_status

        allocate(values(size(p_loss)))
        istat = cvar_ok
        do i = 1, size(p_loss)
            values(i) = es_pdf_scalar(pdf, qf, p_loss(i), intercept, slope, transf, tol, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function es_pdf_vector

    function es_sample_scalar(sample, p_loss, intercept, slope, transf, status) result(value)
        real(dp), intent(in) :: sample(:)
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp) :: value
        real(dp), allocatable :: clean(:)
        real(dp) :: a, b, cutoff
        logical :: do_transf
        integer :: istat, i, count_tail

        call risk_arguments(p_loss, intercept, slope, transf, a, b, do_transf, istat)
        if (istat == cvar_ok) call finite_sample(sample, clean, istat)
        if (istat /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = istat
            return
        end if

        cutoff = empirical_quantile(clean, p_loss)
        value = 0.0_dp
        count_tail = 0
        do i = 1, size(clean)
            if (clean(i) <= cutoff) then
                count_tail = count_tail + 1
                if (do_transf) then
                    value = value - (exp(a + b * clean(i)) - 1.0_dp)
                else
                    value = value - (a + b * clean(i))
                end if
            end if
        end do
        if (count_tail > 0) then
            value = value / real(count_tail, dp)
            istat = cvar_ok
        else
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            istat = cvar_invalid_sample
        end if
        if (present(status)) status = istat
    end function es_sample_scalar

    function es_sample_vector(sample, p_loss, intercept, slope, transf, status) result(values)
        real(dp), intent(in) :: sample(:)
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: i, istat, one_status

        allocate(values(size(p_loss)))
        istat = cvar_ok
        do i = 1, size(p_loss)
            values(i) = es_sample_scalar(sample, p_loss(i), intercept, slope, transf, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function es_sample_vector

    function es_sample_matrix_scalar(sample, p_loss, intercept, slope, transf, status) result(values)
        real(dp), intent(in) :: sample(:, :)
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:)
        integer :: j, istat, one_status

        allocate(values(size(sample, 2)))
        istat = cvar_ok
        do j = 1, size(sample, 2)
            values(j) = es_sample_scalar(sample(:, j), p_loss, intercept, slope, &
                                         transf, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function es_sample_matrix_scalar

    function es_sample_matrix_vector(sample, p_loss, intercept, slope, transf, status) result(values)
        real(dp), intent(in) :: sample(:, :)
        real(dp), intent(in) :: p_loss(:)
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        integer, intent(out), optional :: status
        real(dp), allocatable :: values(:, :)
        integer :: j, istat, one_status

        allocate(values(size(p_loss), size(sample, 2)))
        istat = cvar_ok
        do j = 1, size(sample, 2)
            values(:, j) = es_sample_vector(sample(:, j), p_loss, intercept, slope, &
                                            transf, one_status)
            if (one_status /= cvar_ok .and. istat == cvar_ok) istat = one_status
        end do
        if (present(status)) status = istat
    end function es_sample_matrix_vector

    function quantile_from_cdf(cdf, p, tol, status) result(value)
        procedure(cvar_scalar_function) :: cdf
        real(dp), intent(in) :: p
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: status
        real(dp) :: value
        real(dp) :: lo, hi, mid, flo, fhi, fmid, use_tol
        integer :: iter, istat

        istat = cvar_ok
        if (p <= 0.0_dp .or. p >= 1.0_dp) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = cvar_invalid_probability
            return
        end if
        use_tol = default_tol
        if (present(tol)) use_tol = max(tol, epsilon(1.0_dp))

        lo = -1.0_dp
        hi = 1.0_dp
        flo = cdf(lo)
        fhi = cdf(hi)
        do iter = 1, 200
            if (.not. ieee_is_finite(flo) .or. .not. ieee_is_finite(fhi)) exit
            if (flo <= p .and. fhi >= p) exit
            if (flo > p) then
                hi = lo
                fhi = flo
                lo = 2.0_dp * lo
                flo = cdf(lo)
            else
                lo = hi
                flo = fhi
                hi = 2.0_dp * hi
                fhi = cdf(hi)
            end if
        end do
        if (.not. (flo <= p .and. fhi >= p)) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            if (present(status)) status = cvar_bracket_failure
            return
        end if

        do iter = 1, 300
            mid = 0.5_dp * (lo + hi)
            fmid = cdf(mid)
            if (.not. ieee_is_finite(fmid)) then
                istat = cvar_nonconvergence
                exit
            end if
            if (fmid < p) then
                lo = mid
            else
                hi = mid
            end if
            if (abs(hi - lo) <= use_tol * (1.0_dp + abs(mid))) exit
        end do
        if (iter > 300) istat = cvar_nonconvergence
        value = 0.5_dp * (lo + hi)
        if (present(status)) status = istat
    end function quantile_from_cdf

    function empirical_quantile(sample, p) result(value)
        real(dp), intent(in) :: sample(:)
        real(dp), intent(in) :: p
        real(dp) :: value
        real(dp), allocatable :: sorted(:)
        real(dp) :: h, gamma
        integer :: j, n

        n = size(sample)
        if (n == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        sorted = sample
        call quicksort(sorted, 1, n)
        if (n == 1 .or. p <= 0.0_dp) then
            value = sorted(1)
        else if (p >= 1.0_dp) then
            value = sorted(n)
        else
            h = 1.0_dp + real(n - 1, dp) * p
            j = int(floor(h))
            gamma = h - real(j, dp)
            if (j >= n) then
                value = sorted(n)
            else
                value = (1.0_dp - gamma) * sorted(j) + gamma * sorted(j + 1)
            end if
        end if
    end function empirical_quantile

    subroutine risk_arguments(p_loss, intercept, slope, transf, a, b, do_transf, status)
        real(dp), intent(in) :: p_loss
        real(dp), intent(in), optional :: intercept, slope
        logical, intent(in), optional :: transf
        real(dp), intent(out) :: a, b
        logical, intent(out) :: do_transf
        integer, intent(out) :: status

        a = 0.0_dp
        b = 1.0_dp
        do_transf = .false.
        if (present(intercept)) a = intercept
        if (present(slope)) b = slope
        if (present(transf)) do_transf = transf

        if (p_loss <= 0.0_dp .or. p_loss >= 1.0_dp) then
            status = cvar_invalid_probability
        else if (b <= 0.0_dp .or. .not. ieee_is_finite(b)) then
            status = cvar_invalid_scale
        else
            status = cvar_ok
        end if
    end subroutine risk_arguments

    subroutine finite_sample(sample, clean, status)
        real(dp), intent(in) :: sample(:)
        real(dp), allocatable, intent(out) :: clean(:)
        integer, intent(out) :: status
        integer :: i, n

        n = count(ieee_is_finite(sample))
        if (n < 1) then
            allocate(clean(0))
            status = cvar_invalid_sample
            return
        end if
        allocate(clean(n))
        n = 0
        do i = 1, size(sample)
            if (ieee_is_finite(sample(i))) then
                n = n + 1
                clean(n) = sample(i)
            end if
        end do
        status = cvar_ok
    end subroutine finite_sample

    recursive subroutine quicksort(x, left, right)
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: left, right
        integer :: i, j
        real(dp) :: pivot, temp

        if (left >= right) return
        i = left
        j = right
        pivot = x((left + right) / 2)
        do
            do while (x(i) < pivot)
                i = i + 1
            end do
            do while (x(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                temp = x(i)
                x(i) = x(j)
                x(j) = temp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do
        if (left < j) call quicksort(x, left, j)
        if (i < right) call quicksort(x, i, right)
    end subroutine quicksort

    recursive subroutine adaptive_qf(qf, a0, b0, p_loss, intercept, slope, transf, &
                                      tol, depth, value, error, status)
        procedure(cvar_scalar_function) :: qf
        real(dp), intent(in) :: a0, b0, p_loss, intercept, slope, tol
        logical, intent(in) :: transf
        integer, intent(in) :: depth
        real(dp), intent(out) :: value, error
        integer, intent(out) :: status
        real(dp) :: local_value, local_error, left_value, right_value, left_error, right_error, mid
        integer :: left_status, right_status

        call qk15_qf(qf, a0, b0, p_loss, intercept, slope, transf, local_value, local_error, status)
        if (status /= cvar_ok) then
            value = local_value
            error = local_error
            return
        end if
        if (local_error <= tol * (1.0_dp + abs(local_value)) .or. depth >= max_adapt_depth) then
            value = local_value
            error = local_error
            if (depth >= max_adapt_depth .and. local_error > 10.0_dp * tol * (1.0_dp + abs(local_value))) &
                status = cvar_nonconvergence
            return
        end if
        mid = 0.5_dp * (a0 + b0)
        call adaptive_qf(qf, a0, mid, p_loss, intercept, slope, transf, 0.5_dp * tol, &
                         depth + 1, left_value, left_error, left_status)
        call adaptive_qf(qf, mid, b0, p_loss, intercept, slope, transf, 0.5_dp * tol, &
                         depth + 1, right_value, right_error, right_status)
        value = left_value + right_value
        error = left_error + right_error
        status = left_status
        if (status == cvar_ok) status = right_status
    end subroutine adaptive_qf

    recursive subroutine adaptive_cdf(cdf, a0, b0, p_loss, intercept, slope, transf, &
                                       tol, depth, value, error, status)
        procedure(cvar_scalar_function) :: cdf
        real(dp), intent(in) :: a0, b0, p_loss, intercept, slope, tol
        logical, intent(in) :: transf
        integer, intent(in) :: depth
        real(dp), intent(out) :: value, error
        integer, intent(out) :: status
        real(dp) :: local_value, local_error, left_value, right_value, left_error, right_error, mid
        integer :: left_status, right_status

        call qk15_cdf(cdf, a0, b0, p_loss, intercept, slope, transf, local_value, local_error, status)
        if (status /= cvar_ok) then
            value = local_value
            error = local_error
            return
        end if
        if (local_error <= tol * (1.0_dp + abs(local_value)) .or. depth >= max_adapt_depth) then
            value = local_value
            error = local_error
            if (depth >= max_adapt_depth .and. local_error > 10.0_dp * tol * (1.0_dp + abs(local_value))) &
                status = cvar_nonconvergence
            return
        end if
        mid = 0.5_dp * (a0 + b0)
        call adaptive_cdf(cdf, a0, mid, p_loss, intercept, slope, transf, 0.5_dp * tol, &
                          depth + 1, left_value, left_error, left_status)
        call adaptive_cdf(cdf, mid, b0, p_loss, intercept, slope, transf, 0.5_dp * tol, &
                          depth + 1, right_value, right_error, right_status)
        value = left_value + right_value
        error = left_error + right_error
        status = left_status
        if (status == cvar_ok) status = right_status
    end subroutine adaptive_cdf

    recursive subroutine adaptive_pdf(pdf, a0, b0, upper, p_loss, intercept, slope, transf, &
                                       tol, depth, value, error, status)
        procedure(cvar_scalar_function) :: pdf
        real(dp), intent(in) :: a0, b0, upper, p_loss, intercept, slope, tol
        logical, intent(in) :: transf
        integer, intent(in) :: depth
        real(dp), intent(out) :: value, error
        integer, intent(out) :: status
        real(dp) :: local_value, local_error, left_value, right_value, left_error, right_error, mid
        integer :: left_status, right_status

        call qk15_pdf(pdf, a0, b0, upper, p_loss, intercept, slope, transf, local_value, local_error, status)
        if (status /= cvar_ok) then
            value = local_value
            error = local_error
            return
        end if
        if (local_error <= tol * (1.0_dp + abs(local_value)) .or. depth >= max_adapt_depth) then
            value = local_value
            error = local_error
            if (depth >= max_adapt_depth .and. local_error > 10.0_dp * tol * (1.0_dp + abs(local_value))) &
                status = cvar_nonconvergence
            return
        end if
        mid = 0.5_dp * (a0 + b0)
        call adaptive_pdf(pdf, a0, mid, upper, p_loss, intercept, slope, transf, 0.5_dp * tol, &
                          depth + 1, left_value, left_error, left_status)
        call adaptive_pdf(pdf, mid, b0, upper, p_loss, intercept, slope, transf, 0.5_dp * tol, &
                          depth + 1, right_value, right_error, right_status)
        value = left_value + right_value
        error = left_error + right_error
        status = left_status
        if (status == cvar_ok) status = right_status
    end subroutine adaptive_pdf

    subroutine qk15_qf(qf, a0, b0, p_loss, intercept, slope, transf, value, error, status)
        procedure(cvar_scalar_function) :: qf
        real(dp), intent(in) :: a0, b0, p_loss, intercept, slope
        logical, intent(in) :: transf
        real(dp), intent(out) :: value, error
        integer, intent(out) :: status
        call qk15_core(1, qf, a0, b0, p_loss, intercept, slope, transf, 0.0_dp, value, error, status)
    end subroutine qk15_qf

    subroutine qk15_cdf(cdf, a0, b0, p_loss, intercept, slope, transf, value, error, status)
        procedure(cvar_scalar_function) :: cdf
        real(dp), intent(in) :: a0, b0, p_loss, intercept, slope
        logical, intent(in) :: transf
        real(dp), intent(out) :: value, error
        integer, intent(out) :: status
        call qk15_core(2, cdf, a0, b0, p_loss, intercept, slope, transf, 0.0_dp, value, error, status)
    end subroutine qk15_cdf

    subroutine qk15_pdf(pdf, a0, b0, upper, p_loss, intercept, slope, transf, value, error, status)
        procedure(cvar_scalar_function) :: pdf
        real(dp), intent(in) :: a0, b0, upper, p_loss, intercept, slope
        logical, intent(in) :: transf
        real(dp), intent(out) :: value, error
        integer, intent(out) :: status
        call qk15_core(3, pdf, a0, b0, p_loss, intercept, slope, transf, upper, value, error, status)
    end subroutine qk15_pdf

    subroutine qk15_core(mode, fun, a0, b0, p_loss, intercept, slope, transf, upper, &
                          value, error, status)
        integer, intent(in) :: mode
        procedure(cvar_scalar_function) :: fun
        real(dp), intent(in) :: a0, b0, p_loss, intercept, slope, upper
        logical, intent(in) :: transf
        real(dp), intent(out) :: value, error
        integer, intent(out) :: status
        real(dp), parameter :: xgk(8) = [ &
            0.9914553711208126_dp, 0.9491079123427585_dp, &
            0.8648644233597691_dp, 0.7415311855993944_dp, &
            0.5860872354676911_dp, 0.4058451513773972_dp, &
            0.2077849550078985_dp, 0.0_dp ]
        real(dp), parameter :: wgk(8) = [ &
            0.02293532201052922_dp, 0.06309209262997855_dp, &
            0.1047900103222502_dp, 0.1406532597155259_dp, &
            0.1690047266392679_dp, 0.1903505780647854_dp, &
            0.2044329400752989_dp, 0.2094821410847278_dp ]
        real(dp), parameter :: wg(4) = [ &
            0.1294849661688697_dp, 0.2797053914892767_dp, &
            0.3818300505051189_dp, 0.4179591836734694_dp ]
        real(dp) :: center, half_length, fc, f1, f2, resg, resk, absc
        integer :: j, one_status

        center = 0.5_dp * (a0 + b0)
        half_length = 0.5_dp * (b0 - a0)
        fc = integrand_value(mode, fun, center, p_loss, intercept, slope, transf, upper, one_status)
        if (one_status /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            error = huge(1.0_dp)
            status = one_status
            return
        end if
        resg = wg(4) * fc
        resk = wgk(8) * fc

        do j = 1, 7
            absc = half_length * xgk(j)
            f1 = integrand_value(mode, fun, center - absc, p_loss, intercept, slope, transf, upper, one_status)
            if (one_status /= cvar_ok) exit
            f2 = integrand_value(mode, fun, center + absc, p_loss, intercept, slope, transf, upper, one_status)
            if (one_status /= cvar_ok) exit
            resk = resk + wgk(j) * (f1 + f2)
            select case (j)
            case (2)
                resg = resg + wg(1) * (f1 + f2)
            case (4)
                resg = resg + wg(2) * (f1 + f2)
            case (6)
                resg = resg + wg(3) * (f1 + f2)
            end select
        end do
        if (one_status /= cvar_ok) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            error = huge(1.0_dp)
            status = one_status
            return
        end if
        value = resk * half_length
        error = abs((resk - resg) * half_length)
        status = cvar_ok
    end subroutine qk15_core

    function integrand_value(mode, fun, u, p_loss, intercept, slope, transf, upper, status) result(value)
        integer, intent(in) :: mode
        procedure(cvar_scalar_function) :: fun
        real(dp), intent(in) :: u, p_loss, intercept, slope, upper
        logical, intent(in) :: transf
        integer, intent(out) :: status
        real(dp) :: value
        real(dp) :: x, density, q
        integer :: qstatus

        status = cvar_ok
        select case (mode)
        case (1)
            q = fun(p_loss * u * u)
            x = intercept + slope * q
            if (transf) then
                value = -2.0_dp * u * (exp(x) - 1.0_dp)
            else
                value = -2.0_dp * u * x
            end if
        case (2)
            q = quantile_from_cdf(fun, p_loss * u * u, 1.0e-11_dp, qstatus)
            if (qstatus /= cvar_ok) then
                value = ieee_value(1.0_dp, ieee_quiet_nan)
                status = qstatus
                return
            end if
            x = intercept + slope * q
            if (transf) then
                value = -2.0_dp * u * (exp(x) - 1.0_dp)
            else
                value = -2.0_dp * u * x
            end if
        case (3)
            x = upper + 1.0_dp - 1.0_dp / u
            density = fun(x)
            if (transf) then
                value = -(exp(intercept + slope * x) - 1.0_dp) * density / (p_loss * u * u)
            else
                value = -(intercept + slope * x) * density / (p_loss * u * u)
            end if
        case default
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            status = cvar_nonconvergence
            return
        end select
        if (.not. ieee_is_finite(value)) status = cvar_nonconvergence
    end function integrand_value

end module cvar_risk
