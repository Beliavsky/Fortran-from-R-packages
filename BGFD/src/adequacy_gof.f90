! SPDX-License-Identifier: GPL-2.0-or-later
module adequacy_gof
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use adequacy_kinds, only: dp
    use adequacy_interfaces, only: density_fn, cdf_fn
    use adequacy_math, only: mean_real, sample_variance, sort_real, normal_cdf, normal_quantile
    use adequacy_math, only: invert_matrix, kolmogorov_pvalue
    use adequacy_optim, only: optimize_result, pso_optimize, nelder_mead_optimize
    use adequacy_optim, only: bfgs_optimize, cg_optimize, sann_optimize
    implicit none
    private

    type, public :: goodness_result
        real(dp) :: w = 0.0_dp
        real(dp) :: a = 0.0_dp
        real(dp) :: ks = 0.0_dp
        real(dp) :: ks_pvalue = 1.0_dp
        real(dp), allocatable :: mle(:)
        real(dp), allocatable :: se(:)
        real(dp) :: log_likelihood = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: aicc = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: hqic = 0.0_dp
        real(dp) :: value = 0.0_dp
        real(dp) :: pdf_integral = 0.0_dp
        real(dp) :: cdf_lower = 0.0_dp
        real(dp) :: cdf_upper = 0.0_dp
        integer :: convergence = 0
        logical :: pdf_integral_ok = .true.
        logical :: cdf_endpoints_ok = .true.
    end type goodness_result

    procedure(density_fn), pointer, save :: active_pdf => null()

    public :: goodness_fit, goodness_from_mle

contains

    subroutine goodness_fit(pdf, cdf, starts, data, result, method, domain, lower, upper, &
                            swarm_size, max_iter)
        procedure(density_fn) :: pdf
        procedure(cdf_fn) :: cdf
        real(dp), intent(in) :: starts(:), data(:)
        type(goodness_result), intent(out) :: result
        character(len=*), intent(in), optional :: method
        real(dp), intent(in), optional :: domain(2), lower(:), upper(:)
        integer, intent(in), optional :: swarm_size, max_iter
        type(optimize_result) :: opt
        character(len=20) :: meth
        real(dp) :: dom(2)
        integer :: maxit, ss

        meth = 'PSO'
        if (present(method)) meth = adjustl(method)
        dom = [0.0_dp, ieee_value(1.0_dp, ieee_positive_inf)]
        if (present(domain)) dom = domain
        maxit = 5000
        if (present(max_iter)) maxit = max_iter
        ss = 350
        if (present(swarm_size)) ss = swarm_size

        active_pdf => pdf
        select case (trim(meth))
        case ('PSO', 'P', 'pso', 'p')
            if (.not. present(lower) .or. .not. present(upper)) then
                error stop 'goodness_fit: PSO requires lower and upper bounds'
            end if
            call pso_optimize(negloglik_adapter, data, lower, upper, opt, &
                              swarm_size=ss, max_iter=maxit)
        case ('Nelder-Mead', 'N', 'nelder-mead', 'n')
            call nelder_mead_optimize(negloglik_adapter, data, starts, opt, max_iter=maxit)
        case ('BFGS', 'B', 'bfgs', 'b')
            call bfgs_optimize(negloglik_adapter, data, starts, opt, max_iter=maxit)
        case ('CG', 'C', 'cg', 'c')
            call cg_optimize(negloglik_adapter, data, starts, opt, max_iter=maxit)
        case ('SANN', 'S', 'sann', 's')
            call sann_optimize(negloglik_adapter, data, starts, opt, max_iter=maxit)
        case default
            error stop 'goodness_fit: unknown optimization method'
        end select
        nullify(active_pdf)

        call goodness_from_mle(pdf, cdf, opt%par, data, result, dom)
        result%value = opt%value
        result%convergence = opt%convergence
        call standard_errors(opt%hessian, result%se)
    end subroutine goodness_fit

    subroutine goodness_from_mle(pdf, cdf, mle, data, result, domain)
        procedure(density_fn) :: pdf
        procedure(cdf_fn) :: cdf
        real(dp), intent(in) :: mle(:), data(:)
        type(goodness_result), intent(out) :: result
        real(dp), intent(in), optional :: domain(2)
        real(dp), allocatable :: sx(:), v(:), y(:), u(:)
        real(dp) :: dom(2), vary, epsp, dplus, dminus
        integer :: i, n, p

        n = size(data)
        p = size(mle)
        if (n < 2) error stop 'goodness_from_mle: at least two observations are required'
        dom = [0.0_dp, ieee_value(1.0_dp, ieee_positive_inf)]
        if (present(domain)) dom = domain
        allocate(result%mle(p), result%se(0))
        result%mle = mle

        result%cdf_lower = cdf(mle, dom(1))
        if (ieee_is_finite(dom(2))) then
            result%cdf_upper = cdf(mle, dom(2))
        else
            result%cdf_upper = 1.0_dp
        end if
        result%cdf_endpoints_ok = abs(result%cdf_lower) <= 1.0e-6_dp .and. &
                                  abs(result%cdf_upper - 1.0_dp) <= 1.0e-6_dp
        result%pdf_integral = integrate_density(pdf, mle, dom(1), dom(2))
        result%pdf_integral_ok = ieee_is_finite(result%pdf_integral) .and. &
                                 abs(result%pdf_integral - 1.0_dp) <= 0.1_dp

        sx = data
        call sort_real(sx)
        allocate(v(n), y(n), u(n))
        epsp = 10.0_dp * epsilon(1.0_dp)
        do i = 1, n
            v(i) = min(1.0_dp-epsp, max(epsp, cdf(mle, sx(i))))
            y(i) = normal_quantile(v(i))
        end do
        vary = sample_variance(y)
        if (vary <= tiny(1.0_dp)) then
            u = 0.5_dp
        else
            u = normal_cdf((y - mean_real(y)) / sqrt(vary))
        end if
        u = min(1.0_dp-epsp, max(epsp, u))

        result%w = 1.0_dp / (12.0_dp*real(n,dp))
        result%a = -real(n, dp)
        do i = 1, n
            result%w = result%w + (u(i) - real(2*i-1,dp)/(2.0_dp*real(n,dp)))**2
            result%a = result%a - (real(2*i-1,dp)*log(u(i)) + &
                       real(2*n+1-2*i,dp)*log(1.0_dp-u(i))) / real(n,dp)
        end do
        result%w = result%w * (1.0_dp + 0.5_dp/real(n,dp))
        result%a = result%a * (1.0_dp + 0.75_dp/real(n,dp) + 2.25_dp/real(n*n,dp))

        dplus = 0.0_dp
        dminus = 0.0_dp
        do i = 1, n
            dplus = max(dplus, real(i,dp)/real(n,dp) - v(i))
            dminus = max(dminus, v(i) - real(i-1,dp)/real(n,dp))
        end do
        result%ks = max(dplus, dminus)
        result%ks_pvalue = kolmogorov_pvalue(result%ks, n)

        result%log_likelihood = loglikelihood(pdf, mle, data)
        result%aic = -2.0_dp*result%log_likelihood + 2.0_dp*real(p,dp)
        if (n > p + 1) then
            result%aicc = result%aic + 2.0_dp*real(p*(p+1),dp)/real(n-p-1,dp)
        else
            result%aicc = ieee_value(1.0_dp, ieee_quiet_nan)
        end if
        result%bic = -2.0_dp*result%log_likelihood + real(p,dp)*log(real(n,dp))
        result%hqic = -2.0_dp*result%log_likelihood + 2.0_dp*log(log(real(n,dp)))*real(p,dp)
        result%value = -result%log_likelihood
    end subroutine goodness_from_mle

    function negloglik_adapter(par, data) result(value)
        real(dp), intent(in) :: par(:), data(:)
        real(dp) :: value
        integer :: i
        real(dp) :: d

        if (.not. associated(active_pdf)) error stop 'negloglik_adapter: no active density'
        value = 0.0_dp
        do i = 1, size(data)
            d = active_pdf(par, data(i))
            if (.not. ieee_is_finite(d) .or. d <= 0.0_dp) then
                value = huge(1.0_dp) / 1000.0_dp
                return
            end if
            value = value - log(d)
        end do
    end function negloglik_adapter

    function loglikelihood(pdf, par, data) result(value)
        procedure(density_fn) :: pdf
        real(dp), intent(in) :: par(:), data(:)
        real(dp) :: value, d
        integer :: i

        value = 0.0_dp
        do i = 1, size(data)
            d = pdf(par, data(i))
            if (.not. ieee_is_finite(d) .or. d <= 0.0_dp) then
                value = -huge(1.0_dp) / 1000.0_dp
                return
            end if
            value = value + log(d)
        end do
    end function loglikelihood

    subroutine standard_errors(hessian, se)
        real(dp), intent(in) :: hessian(:, :)
        real(dp), allocatable, intent(out) :: se(:)
        real(dp), allocatable :: inv(:, :)
        logical :: ok
        integer :: i, n

        n = size(hessian, 1)
        allocate(se(n), inv(n,n))
        call invert_matrix(hessian, inv, ok)
        if (.not. ok) then
            se = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        do i = 1, n
            if (inv(i,i) >= 0.0_dp) then
                se(i) = sqrt(inv(i,i))
            else
                se(i) = ieee_value(1.0_dp, ieee_quiet_nan)
            end if
        end do
    end subroutine standard_errors

    function integrate_density(pdf, par, lo, hi) result(value)
        procedure(density_fn) :: pdf
        real(dp), intent(in) :: par(:), lo, hi
        real(dp) :: value
        real(dp), parameter :: delta = 1.0e-8_dp
        integer :: mode

        if (ieee_is_finite(lo) .and. ieee_is_finite(hi)) then
            mode = 0
            value = adaptive_simpson(pdf, par, lo, hi, mode, lo, hi, 1.0e-8_dp, 16)
        else if (ieee_is_finite(lo)) then
            mode = 1
            value = adaptive_simpson(pdf, par, delta, 1.0_dp-delta, mode, lo, hi, 1.0e-8_dp, 18)
        else if (ieee_is_finite(hi)) then
            mode = 2
            value = adaptive_simpson(pdf, par, delta, 1.0_dp-delta, mode, lo, hi, 1.0e-8_dp, 18)
        else
            mode = 3
            value = adaptive_simpson(pdf, par, delta, 1.0_dp-delta, mode, lo, hi, 1.0e-8_dp, 18)
        end if
    end function integrate_density

    recursive function adaptive_simpson(pdf, par, a, b, mode, lo, hi, tol, depth) result(ans)
        procedure(density_fn) :: pdf
        real(dp), intent(in) :: par(:), a, b, lo, hi, tol
        integer, intent(in) :: mode, depth
        real(dp) :: ans, fa, fm, fb, whole

        fa = transformed_density(pdf, par, a, mode, lo, hi)
        fm = transformed_density(pdf, par, 0.5_dp*(a+b), mode, lo, hi)
        fb = transformed_density(pdf, par, b, mode, lo, hi)
        whole = (b-a)*(fa + 4.0_dp*fm + fb)/6.0_dp
        ans = adaptive_simpson_rec(pdf, par, a, b, fa, fm, fb, whole, mode, lo, hi, tol, depth)
    end function adaptive_simpson

    recursive function adaptive_simpson_rec(pdf, par, a, b, fa, fm, fb, whole, &
                                            mode, lo, hi, tol, depth) result(ans)
        procedure(density_fn) :: pdf
        real(dp), intent(in) :: par(:), a, b, fa, fm, fb, whole, lo, hi, tol
        integer, intent(in) :: mode, depth
        real(dp) :: ans, m, lm, rm, flm, frm, left, right

        m = 0.5_dp*(a+b)
        lm = 0.5_dp*(a+m)
        rm = 0.5_dp*(m+b)
        flm = transformed_density(pdf, par, lm, mode, lo, hi)
        frm = transformed_density(pdf, par, rm, mode, lo, hi)
        left = (m-a)*(fa + 4.0_dp*flm + fm)/6.0_dp
        right = (b-m)*(fm + 4.0_dp*frm + fb)/6.0_dp
        if (depth <= 0 .or. abs(left+right-whole) <= 15.0_dp*tol) then
            ans = left + right + (left + right - whole)/15.0_dp
        else
            ans = adaptive_simpson_rec(pdf, par, a, m, fa, flm, fm, left, &
                                       mode, lo, hi, 0.5_dp*tol, depth-1) + &
                  adaptive_simpson_rec(pdf, par, m, b, fm, frm, fb, right, &
                                       mode, lo, hi, 0.5_dp*tol, depth-1)
        end if
    end function adaptive_simpson_rec

    function transformed_density(pdf, par, t, mode, lo, hi) result(y)
        procedure(density_fn) :: pdf
        real(dp), intent(in) :: par(:), t, lo, hi
        integer, intent(in) :: mode
        real(dp) :: y, x, jac, angle

        select case (mode)
        case (0)
            x = t
            jac = 1.0_dp
        case (1)
            x = lo + t/(1.0_dp-t)
            jac = 1.0_dp/(1.0_dp-t)**2
        case (2)
            x = hi - t/(1.0_dp-t)
            jac = 1.0_dp/(1.0_dp-t)**2
        case default
            angle = acos(-1.0_dp)*(t-0.5_dp)
            x = tan(angle)
            jac = acos(-1.0_dp)/(cos(angle)**2)
        end select
        y = pdf(par, x)*jac
        if (.not. ieee_is_finite(y)) y = 0.0_dp
    end function transformed_density

end module adequacy_gof
