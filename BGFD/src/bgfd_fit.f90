! SPDX-License-Identifier: GPL-2.0-or-later
module bgfd_fit
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use adequacy_kinds, only: dp
    use adequacy_gof, only: goodness_result, goodness_fit
    use bgfd_core, only: bgfd_npar, bgfd_pdf, bgfd_cdf
    use bgfd_core, only: id_e => bgfd_e, id_ee => bgfd_ee, id_w => bgfd_w, id_ew => bgfd_ew
    use bgfd_core, only: id_f => bgfd_f, id_l => bgfd_l, id_burr => bgfd_b, id_bx => bgfd_bx
    implicit none
    private

    type, public :: bgfd_fit_result
        real(dp), allocatable :: params(:)
        real(dp), allocatable :: se(:)
        real(dp) :: log_likelihood = 0.0_dp
        real(dp) :: neg2loglik = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: aicc = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: hqic = 0.0_dp
        real(dp) :: cramer_von_mises = 0.0_dp
        real(dp) :: anderson_darling = 0.0_dp
        real(dp) :: ks_statistic = 0.0_dp
        real(dp) :: ks_pvalue = 1.0_dp
        integer :: convergence = 0
        logical :: pdf_integral_ok = .false.
        logical :: cdf_endpoints_ok = .false.
    end type bgfd_fit_result

    integer, save :: active_family = 0
    logical, save :: active_complementary = .false.

    public :: fit_bgfd
    public :: m_bell_e, m_bell_ee, m_bell_w, m_bell_ew
    public :: m_bell_f, m_bell_l, m_bell_b, m_bell_bx
    public :: m_cbell_e, m_cbell_ee, m_cbell_w, m_cbell_ew
    public :: m_cbell_f, m_cbell_l, m_cbell_b, m_cbell_bx

contains

    subroutine fit_bgfd(data, family, complementary, starts, result, method, max_iter)
        real(dp), intent(in) :: data(:), starts(:)
        integer, intent(in) :: family
        logical, intent(in) :: complementary
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        type(goodness_result) :: gof, gof_refined
        real(dp), allocatable :: log_starts(:), lower(:), upper(:)
        real(dp) :: domain(2)
        character(len=20) :: meth
        integer :: npar, maxit

        npar = bgfd_npar(family)
        if (npar <= 0) error stop 'fit_bgfd: unknown family'
        if (size(starts) /= npar) error stop 'fit_bgfd: wrong number of starting parameters'
        if (any(starts <= 0.0_dp)) error stop 'fit_bgfd: starting parameters must be positive'
        if (size(data) < 2 .or. any(data < 0.0_dp)) error stop 'fit_bgfd: data must be nonnegative'

        meth = 'BFGS'
        if (present(method)) meth = adjustl(method)
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        log_starts = log(starts)
        domain = [0.0_dp, ieee_value(1.0_dp, ieee_positive_inf)]
        active_family = family
        active_complementary = complementary

        if (trim(meth) == 'P' .or. trim(meth) == 'p' .or. trim(meth) == 'PSO' .or. trim(meth) == 'pso') then
            allocate(lower(npar), upper(npar))
            lower = log(1.0e-6_dp)
            upper = log(1.0e6_dp)
            call goodness_fit(logpar_pdf, logpar_cdf, log_starts, data, gof, method=meth, &
                domain=domain, lower=lower, upper=upper, max_iter=maxit)
        else
            call goodness_fit(logpar_pdf, logpar_cdf, log_starts, data, gof, method=meth, &
                domain=domain, max_iter=maxit)
        end if

        if ((trim(meth) == 'B' .or. trim(meth) == 'b' .or. trim(meth) == 'BFGS' .or. &
            trim(meth) == 'bfgs') .and. gof%convergence /= 0) then
            call goodness_fit(logpar_pdf, logpar_cdf, gof%mle, data, gof_refined, method='N', &
                domain=domain, max_iter=maxit)
            if (gof_refined%value <= gof%value) gof = gof_refined
        end if

        call copy_fit_result(gof, result)
        active_family = 0
        active_complementary = .false.
    end subroutine fit_bgfd

    function logpar_pdf(logpar, x) result(value)
        real(dp), intent(in) :: logpar(:), x
        real(dp) :: value
        real(dp) :: par(size(logpar))
        par = exp(logpar)
        value = bgfd_pdf(active_family, active_complementary, x, par)
    end function logpar_pdf

    function logpar_cdf(logpar, x) result(value)
        real(dp), intent(in) :: logpar(:), x
        real(dp) :: value
        real(dp) :: par(size(logpar))
        par = exp(logpar)
        value = bgfd_cdf(active_family, active_complementary, x, par)
    end function logpar_cdf

    subroutine copy_fit_result(gof, result)
        type(goodness_result), intent(in) :: gof
        type(bgfd_fit_result), intent(out) :: result
        integer :: n

        n = size(gof%mle)
        allocate(result%params(n), result%se(n))
        result%params = exp(gof%mle)
        if (size(gof%se) == n) then
            result%se = result%params * gof%se
        else
            result%se = 0.0_dp
        end if
        result%log_likelihood = gof%log_likelihood
        result%neg2loglik = -2.0_dp*gof%log_likelihood
        result%aic = gof%aic
        result%aicc = gof%aicc
        result%bic = gof%bic
        result%hqic = gof%hqic
        result%cramer_von_mises = gof%w
        result%anderson_darling = gof%a
        result%ks_statistic = gof%ks
        result%ks_pvalue = gof%ks_pvalue
        result%convergence = gof%convergence
        result%pdf_integral_ok = gof%pdf_integral_ok
        result%cdf_endpoints_ok = gof%cdf_endpoints_ok
    end subroutine copy_fit_result

    subroutine m_bell_e(x, alpha, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_e, .false., [alpha, lambda], result, meth, maxit)
    end subroutine m_bell_e

    subroutine m_bell_ee(x, alpha, beta, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha, beta, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_ee, .false., [alpha, beta, lambda], result, meth, maxit)
    end subroutine m_bell_ee

    subroutine m_bell_w(x, alpha, beta, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha, beta, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_w, .false., [alpha, beta, lambda], result, meth, maxit)
    end subroutine m_bell_w

    subroutine m_bell_ew(x, alpha, beta, theta, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha, beta, theta, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_ew, .false., [alpha, beta, theta, lambda], result, meth, maxit)
    end subroutine m_bell_ew

    subroutine m_bell_f(x, a, b, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: a, b, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_f, .false., [a, b, lambda], result, meth, maxit)
    end subroutine m_bell_f

    subroutine m_bell_l(x, b, q, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: b, q, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_l, .false., [b, q, lambda], result, meth, maxit)
    end subroutine m_bell_l

    subroutine m_bell_b(x, a, b, k, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: a, b, k, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_burr, .false., [a, b, k, lambda], result, meth, maxit)
    end subroutine m_bell_b

    subroutine m_bell_bx(x, a, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: a, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_bx, .false., [a, lambda], result, meth, maxit)
    end subroutine m_bell_bx

    subroutine m_cbell_e(x, alpha, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_e, .true., [alpha, lambda], result, meth, maxit)
    end subroutine m_cbell_e

    subroutine m_cbell_ee(x, alpha, beta, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha, beta, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_ee, .true., [alpha, beta, lambda], result, meth, maxit)
    end subroutine m_cbell_ee

    subroutine m_cbell_w(x, alpha, beta, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha, beta, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_w, .true., [alpha, beta, lambda], result, meth, maxit)
    end subroutine m_cbell_w

    subroutine m_cbell_ew(x, alpha, beta, theta, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: alpha, beta, theta, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_ew, .true., [alpha, beta, theta, lambda], result, meth, maxit)
    end subroutine m_cbell_ew

    subroutine m_cbell_f(x, a, b, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: a, b, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_f, .true., [a, b, lambda], result, meth, maxit)
    end subroutine m_cbell_f

    subroutine m_cbell_l(x, b, q, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: b, q, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_l, .true., [b, q, lambda], result, meth, maxit)
    end subroutine m_cbell_l

    subroutine m_cbell_b(x, a, b, k, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: a, b, k, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_burr, .true., [a, b, k, lambda], result, meth, maxit)
    end subroutine m_cbell_b

    subroutine m_cbell_bx(x, a, lambda, result, method, max_iter)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: a, lambda
        type(bgfd_fit_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_iter
        character(len=20) :: meth
        integer :: maxit
        meth = 'BFGS'
        if (present(method)) meth = method
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        call fit_bgfd(x, id_bx, .true., [a, lambda], result, meth, maxit)
    end subroutine m_cbell_bx

end module bgfd_fit
