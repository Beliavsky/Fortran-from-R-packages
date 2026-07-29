! SPDX-License-Identifier: GPL-2.0-or-later
module evir_fitting
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use evir_kinds, only : dp
    use evir_types, only : gev_fit_result, gpd_fit_result, pot_fit_result, profile_result, decluster_result, &
        evir_ok, evir_invalid_input, evir_no_exceedances, evir_optimization_failed, &
        evir_singular_hessian, evir_domain_error
    use evir_math, only : mean_value, sample_variance, normal_quantile, &
        chi_square_quantile_df1, sort_ascending, safe_nan
    use evir_data, only : find_threshold, block_maxima, decluster
    use evir_distributions, only : qgev
    use evir_optimize, only : nelder_mead, numerical_hessian, objective_function
    implicit none
    private

    type :: vector_context
        real(dp), allocatable :: data(:)
    end type vector_context

    type :: pot_context
        real(dp), allocatable :: data(:)
        real(dp) :: threshold = 0.0_dp
        real(dp) :: span = 1.0_dp
    end type pot_context

    type :: gev_profile_context
        real(dp), allocatable :: data(:)
        real(dp) :: pp = 0.05_dp
        real(dp) :: return_level = 0.0_dp
    end type gev_profile_context

    type :: gpd_profile_context
        real(dp), allocatable :: excess(:)
        real(dp) :: a = 0.01_dp
        real(dp) :: threshold = 0.0_dp
        real(dp) :: target = 0.0_dp
        logical :: shortfall = .false.
    end type gpd_profile_context

    public :: fit_gev, fit_gumbel, fit_gpd, fit_pot
    public :: gev_return_level, gev_return_level_profile
    public :: gpd_quantile_estimate, gpd_quantile_wald
    public :: gpd_quantile_profile, gpd_shortfall_estimate, gpd_shortfall_profile
    public :: risk_measures
    public :: gev_negloglik_value, gpd_negloglik_value, pot_negloglik_value

contains

    pure real(dp) function gev_negloglik_value(theta, data) result(f)
        real(dp), intent(in) :: theta(:), data(:)
        real(dp) :: xi, sigma, mu, z, t
        integer :: i
        if (size(theta) /= 3 .or. size(data) == 0) then
            f = huge(1.0_dp)/100.0_dp
            return
        end if
        xi = theta(1)
        sigma = theta(2)
        mu = theta(3)
        if (sigma <= 0.0_dp) then
            f = huge(1.0_dp)/100.0_dp
            return
        end if
        f = real(size(data), dp)*log(sigma)
        if (abs(xi) <= 1.0e-7_dp) then
            do i = 1, size(data)
                z = (data(i)-mu)/sigma
                f = f + z + exp(-z)
            end do
        else
            do i = 1, size(data)
                t = 1.0_dp + xi*(data(i)-mu)/sigma
                if (t <= 0.0_dp) then
                    f = huge(1.0_dp)/100.0_dp
                    return
                end if
                f = f + (1.0_dp+1.0_dp/xi)*log(t) + t**(-1.0_dp/xi)
            end do
        end if
    end function gev_negloglik_value

    pure real(dp) function gumbel_negloglik_value(theta, data) result(f)
        real(dp), intent(in) :: theta(:), data(:)
        real(dp) :: sigma, mu, z
        integer :: i
        if (size(theta) /= 2 .or. size(data) == 0) then
            f = huge(1.0_dp)/100.0_dp
            return
        end if
        sigma = theta(1)
        mu = theta(2)
        if (sigma <= 0.0_dp) then
            f = huge(1.0_dp)/100.0_dp
            return
        end if
        f = real(size(data), dp)*log(sigma)
        do i = 1, size(data)
            z = (data(i)-mu)/sigma
            f = f + z + exp(-z)
        end do
    end function gumbel_negloglik_value

    pure real(dp) function gpd_negloglik_value(theta, excess) result(f)
        real(dp), intent(in) :: theta(:), excess(:)
        real(dp) :: xi, beta, t
        integer :: i
        if (size(theta) /= 2 .or. size(excess) == 0) then
            f = huge(1.0_dp)/100.0_dp
            return
        end if
        xi = theta(1)
        beta = theta(2)
        if (beta <= 0.0_dp) then
            f = huge(1.0_dp)/100.0_dp
            return
        end if
        f = real(size(excess), dp)*log(beta)
        if (abs(xi) <= 1.0e-7_dp) then
            f = f + sum(excess)/beta
        else
            do i = 1, size(excess)
                t = 1.0_dp+xi*excess(i)/beta
                if (t <= 0.0_dp) then
                    f = huge(1.0_dp)/100.0_dp
                    return
                end if
                f = f + (1.0_dp+xi)*log(t)/xi
            end do
        end if
    end function gpd_negloglik_value

    pure real(dp) function pot_negloglik_value(theta, exceedances, threshold, span) result(f)
        real(dp), intent(in) :: theta(:), exceedances(:), threshold, span
        real(dp) :: xi, sigma, mu, t, z
        integer :: i
        if (size(theta) /= 3 .or. size(exceedances) == 0 .or. span <= 0.0_dp) then
            f = huge(1.0_dp)/100.0_dp
            return
        end if
        xi = theta(1)
        sigma = theta(2)
        mu = theta(3)
        if (sigma <= 0.0_dp) then
            f = huge(1.0_dp)/100.0_dp
            return
        end if
        if (abs(xi) <= 1.0e-7_dp) then
            f = span*exp(-(threshold-mu)/sigma) + real(size(exceedances), dp)*log(sigma)
            f = f + sum((exceedances-mu)/sigma)
        else
            t = 1.0_dp+xi*(threshold-mu)/sigma
            if (t <= 0.0_dp) then
                f = huge(1.0_dp)/100.0_dp
                return
            end if
            f = span*t**(-1.0_dp/xi) + real(size(exceedances), dp)*log(sigma)
            do i = 1, size(exceedances)
                z = 1.0_dp+xi*(exceedances(i)-mu)/sigma
                if (z <= 0.0_dp) then
                    f = huge(1.0_dp)/100.0_dp
                    return
                end if
                f = f + (1.0_dp+1.0_dp/xi)*log(z)
            end do
        end if
    end function pot_negloglik_value

    function fit_gev(data, block_size) result(out)
        real(dp), intent(in) :: data(:)
        integer, intent(in), optional :: block_size
        type(gev_fit_result) :: out
        type(vector_context) :: ctx
        real(dp), allocatable :: work(:)
        real(dp) :: x0(3), xbest(3), fbest, sigma0, mu0
        real(dp) :: hess(3,3), cov(3,3)
        logical :: conv
        integer :: hstatus

        if (size(data) < 3) then
            out%status = evir_invalid_input
            return
        end if
        out%n_all = size(data)
        if (present(block_size)) then
            if (block_size <= 0) then
                out%status = evir_invalid_input
                return
            end if
            work = block_maxima(data, block_size)
            out%block_size = block_size
        else
            allocate(work(size(data)))
            work = data
            out%n_all = 0
        end if
        out%data = work
        out%n = size(work)
        sigma0 = sqrt(max(6.0_dp*sample_variance(work), epsilon(1.0_dp)))/acos(-1.0_dp)
        mu0 = mean_value(work)-0.57722_dp*sigma0
        x0 = [0.1_dp, sigma0, mu0]
        ctx%data = work
        call nelder_mead(gev_objective, ctx, x0, xbest, fbest, conv, initial_step=0.08_dp)
        out%xi = xbest(1)
        out%sigma = xbest(2)
        out%mu = xbest(3)
        out%nllh = fbest
        out%converged = conv
        if (.not. conv) out%status = evir_optimization_failed
        call numerical_hessian(gev_objective, ctx, xbest, hess, cov, hstatus)
        if (hstatus == 0) then
            out%varcov = cov
            out%se = sqrt(max(diagonal3(cov), 0.0_dp))
        else
            out%status = evir_singular_hessian
            out%varcov = safe_nan()
            out%se = safe_nan()
        end if
    end function fit_gev

    function fit_gumbel(data, block_size) result(out)
        real(dp), intent(in) :: data(:)
        integer, intent(in), optional :: block_size
        type(gev_fit_result) :: out
        type(vector_context) :: ctx
        real(dp), allocatable :: work(:)
        real(dp) :: x0(2), xbest(2), fbest, sigma0, mu0
        real(dp) :: hess(2,2), cov(2,2)
        logical :: conv
        integer :: hstatus

        if (size(data) < 3) then
            out%status = evir_invalid_input
            return
        end if
        out%gumbel = .true.
        out%n_all = size(data)
        if (present(block_size)) then
            if (block_size <= 0) then
                out%status = evir_invalid_input
                return
            end if
            work = block_maxima(data, block_size)
            out%block_size = block_size
        else
            allocate(work(size(data)))
            work = data
            out%n_all = 0
        end if
        out%data = work
        out%n = size(work)
        sigma0 = sqrt(max(6.0_dp*sample_variance(work), epsilon(1.0_dp)))/acos(-1.0_dp)
        mu0 = mean_value(work)-0.57722_dp*sigma0
        x0 = [sigma0, mu0]
        ctx%data = work
        call nelder_mead(gumbel_objective, ctx, x0, xbest, fbest, conv, initial_step=0.08_dp)
        out%xi = 0.0_dp
        out%sigma = xbest(1)
        out%mu = xbest(2)
        out%nllh = fbest
        out%converged = conv
        if (.not. conv) out%status = evir_optimization_failed
        call numerical_hessian(gumbel_objective, ctx, xbest, hess, cov, hstatus)
        if (hstatus == 0) then
            out%varcov = 0.0_dp
            out%varcov(2:3,2:3) = cov
            out%se = 0.0_dp
            out%se(2:3) = sqrt(max([cov(1,1), cov(2,2)], 0.0_dp))
        else
            out%status = evir_singular_hessian
            out%varcov = safe_nan()
            out%se = safe_nan()
        end if
    end function fit_gumbel

    function fit_gpd(data, threshold, nextremes, method, information) result(out)
        real(dp), intent(in) :: data(:)
        real(dp), intent(in), optional :: threshold
        integer, intent(in), optional :: nextremes
        character(len=*), intent(in), optional :: method, information
        type(gpd_fit_result) :: out
        type(vector_context) :: ctx
        real(dp), allocatable :: excess(:), sorted(:)
        real(dp) :: u, xbar, s2, xi0, beta0, x0(2), xbest(2), fbest
        real(dp) :: hess(2,2), cov(2,2), a0, a1, gamma, delta, denom
        real(dp) :: one, two, cv
        integer :: i, nu, hstatus, istat
        logical :: conv
        character(len=8) :: meth, info

        meth = 'ml'
        if (present(method)) meth = lowercase(method)
        info = 'observed'
        if (present(information)) info = lowercase(information)
        if (size(data) < 3 .or. (present(threshold) .eqv. present(nextremes))) then
            out%status = evir_invalid_input
            return
        end if
        if (present(threshold)) then
            u = threshold
        else
            u = find_threshold(data, nextremes, istat)
            if (istat /= evir_ok) then
                out%status = istat
                return
            end if
        end if
        nu = count(data > u)
        if (nu < 2) then
            out%status = evir_no_exceedances
            return
        end if
        allocate(excess(nu), out%exceedances(nu))
        excess = pack(data-u, data > u)
        out%exceedances = excess+u
        out%n = size(data)
        out%n_exceed = nu
        out%threshold = u
        out%p_less_threshold = 1.0_dp-real(nu,dp)/real(size(data),dp)
        out%method = meth
        out%information = info
        xbar = mean_value(excess)

        if (trim(meth) == 'pwm') then
            sorted = excess
            call sort_ascending(sorted)
            gamma = -0.35_dp
            delta = 0.0_dp
            a0 = xbar
            a1 = 0.0_dp
            do i = 1, nu
                a1 = a1 + sorted(i)*(1.0_dp-(real(i,dp)+gamma)/(real(nu,dp)+delta))
            end do
            a1 = a1/real(nu,dp)
            out%xi = 2.0_dp-a0/(a0-2.0_dp*a1)
            out%beta = (2.0_dp*a0*a1)/(a0-2.0_dp*a1)
            denom = real(nu,dp)*(1.0_dp-2.0_dp*out%xi)*(3.0_dp-2.0_dp*out%xi)
            if (out%xi > 0.5_dp .or. abs(denom) <= epsilon(1.0_dp)) then
                out%varcov = safe_nan()
                out%se = safe_nan()
            else
                one = (1.0_dp-out%xi)*(1.0_dp-out%xi+2.0_dp*out%xi**2)*(2.0_dp-out%xi)**2
                two = (7.0_dp-18.0_dp*out%xi+11.0_dp*out%xi**2-2.0_dp*out%xi**3)*out%beta**2
                cv = out%beta*(2.0_dp-out%xi)*(2.0_dp-6.0_dp*out%xi+7.0_dp*out%xi**2-2.0_dp*out%xi**3)
                out%varcov = reshape([one,cv,cv,two],[2,2])/denom
                out%se = sqrt(max([out%varcov(1,1),out%varcov(2,2)],0.0_dp))
            end if
            out%information = 'expected'
            out%converged = .true.
            out%nllh = safe_nan()
            return
        else if (trim(meth) /= 'ml') then
            out%status = evir_invalid_input
            return
        end if

        s2 = sample_variance(excess)
        if (s2 > epsilon(1.0_dp)) then
            xi0 = -0.5_dp*((xbar*xbar)/s2-1.0_dp)
            beta0 = 0.5_dp*xbar*((xbar*xbar)/s2+1.0_dp)
        else
            xi0 = 0.1_dp
            beta0 = max(xbar, epsilon(1.0_dp))
        end if
        if (beta0 <= 0.0_dp) beta0 = max(xbar, 1.0e-3_dp)
        x0 = [xi0,beta0]
        ctx%data = excess
        call nelder_mead(gpd_objective, ctx, x0, xbest, fbest, conv, initial_step=0.08_dp)
        out%xi = xbest(1)
        out%beta = xbest(2)
        out%nllh = fbest
        out%converged = conv
        if (.not. conv) out%status = evir_optimization_failed
        if (trim(info) == 'expected') then
            out%varcov(1,1) = (1.0_dp+out%xi)**2/real(nu,dp)
            out%varcov(2,2) = 2.0_dp*(1.0_dp+out%xi)*out%beta**2/real(nu,dp)
            out%varcov(1,2) = -(1.0_dp+out%xi)*out%beta/real(nu,dp)
            out%varcov(2,1) = out%varcov(1,2)
            out%information = 'expected'
        else
            call numerical_hessian(gpd_objective, ctx, xbest, hess, cov, hstatus)
            if (hstatus == 0) then
                out%varcov = cov
            else
                out%status = evir_singular_hessian
                out%varcov = safe_nan()
            end if
        end if
        out%se = sqrt(max([out%varcov(1,1),out%varcov(2,2)],0.0_dp))
    end function fit_gpd

    function fit_pot(data, times, threshold, nextremes, run) result(out)
        real(dp), intent(in) :: data(:)
        real(dp), intent(in), optional :: times(:), threshold, run
        integer, intent(in), optional :: nextremes
        type(pot_fit_result) :: out
        type(pot_context) :: ctx
        real(dp), allocatable :: alltimes(:), exceed(:), exctimes(:)
        real(dp) :: u, span, xbar, s2, shape0, betahat, extra, scale0
        real(dp) :: x0(3), xbest(3), fbest, hess(3,3), cov(3,3)
        integer :: i, nu, hstatus, istat
        logical :: conv
        type(decluster_result) :: dc

        if (size(data) < 3 .or. (present(threshold) .eqv. present(nextremes))) then
            out%status = evir_invalid_input
            return
        end if
        allocate(alltimes(size(data)))
        if (present(times)) then
            if (size(times) /= size(data)) then
                out%status = evir_invalid_input
                return
            end if
            alltimes = times
        else
            alltimes = [(real(i,dp), i=1,size(data))]
        end if
        if (present(threshold)) then
            u = threshold
        else
            u = find_threshold(data,nextremes,istat)
            if (istat /= evir_ok) then
                out%status = istat
                return
            end if
        end if
        nu = count(data > u)
        if (nu < 2) then
            out%status = evir_no_exceedances
            return
        end if
        exceed = pack(data,data > u)
        exctimes = pack(alltimes,data > u)
        if (present(run)) then
            dc = decluster(exceed,exctimes,run)
            if (dc%status /= evir_ok) then
                out%status = dc%status
                return
            end if
            exceed = dc%values
            exctimes = dc%times
            nu = size(exceed)
            out%run = run
        end if
        span = alltimes(size(alltimes))-alltimes(1)
        if (span <= 0.0_dp) then
            out%status = evir_invalid_input
            return
        end if
        xbar = mean_value(exceed)-u
        s2 = sample_variance(exceed)
        if (s2 > epsilon(1.0_dp)) then
            shape0 = -0.5_dp*((xbar*xbar)/s2-1.0_dp)
            betahat = 0.5_dp*xbar*((xbar*xbar)/s2+1.0_dp)
        else
            shape0 = 0.1_dp
            betahat = max(xbar,1.0e-3_dp)
        end if
        if (abs(shape0) <= 1.0e-8_dp) then
            extra = -log(real(nu,dp)/span)
        else
            extra = ((real(nu,dp)/span)**(-shape0)-1.0_dp)/shape0
        end if
        scale0 = betahat/(1.0_dp+shape0*extra)
        if (scale0 <= 0.0_dp) scale0 = max(betahat,1.0e-3_dp)
        x0 = [shape0,scale0,0.0_dp]
        ctx%data = exceed
        ctx%threshold = u
        ctx%span = span
        call nelder_mead(pot_objective,ctx,x0,xbest,fbest,conv,initial_step=0.08_dp)
        out%n = size(data)
        out%n_exceed = nu
        out%period = [alltimes(1),alltimes(size(alltimes))]
        out%span = span
        out%threshold = u
        out%p_less_threshold = 1.0_dp-real(count(data>u),dp)/real(size(data),dp)
        out%intensity = real(nu,dp)/span
        out%xi = xbest(1)
        out%sigma = xbest(2)
        out%mu = xbest(3)
        out%beta = out%sigma+out%xi*(u-out%mu)
        out%nllh = fbest
        out%converged = conv
        out%exceedances = exceed
        out%times = exctimes
        if (.not. conv) out%status = evir_optimization_failed
        call numerical_hessian(pot_objective,ctx,xbest,hess,cov,hstatus)
        if (hstatus == 0) then
            out%varcov = cov
            out%se = sqrt(max(diagonal3(cov),0.0_dp))
        else
            out%status = evir_singular_hessian
            out%varcov = safe_nan()
            out%se = safe_nan()
        end if
    end function fit_pot

    pure real(dp) function gev_return_level(fit, k_blocks) result(level)
        type(gev_fit_result), intent(in) :: fit
        real(dp), intent(in) :: k_blocks
        if (k_blocks <= 1.0_dp .or. fit%sigma <= 0.0_dp) then
            level = safe_nan()
        else
            level = qgev(1.0_dp-1.0_dp/k_blocks,fit%xi,fit%mu,fit%sigma)
        end if
    end function gev_return_level

    function gev_return_level_profile(fit,k_blocks,ci_p) result(out)
        type(gev_fit_result), intent(in) :: fit
        real(dp), intent(in) :: k_blocks
        real(dp), intent(in), optional :: ci_p
        type(profile_result) :: out
        real(dp), parameter :: mult(21) = [0.5_dp,0.6_dp,0.7_dp,0.8_dp,0.85_dp,0.9_dp,0.95_dp,1.0_dp, &
            1.1_dp,1.2_dp,1.25_dp,1.5_dp,1.75_dp,2.0_dp,2.25_dp,2.5_dp,2.75_dp,3.0_dp,3.25_dp,3.5_dp,4.5_dp]
        type(gev_profile_context) :: ctx
        real(dp) :: ci, estimate, x0(2), xbest(2), fbest, cutoff
        logical :: conv
        integer :: i
        ci = 0.95_dp
        if (present(ci_p)) ci = ci_p
        if (fit%status /= evir_ok .or. k_blocks <= 1.0_dp) then
            out%status = evir_invalid_input
            return
        end if
        estimate = gev_return_level(fit,k_blocks)
        out%estimate = estimate
        allocate(out%x(size(mult)),out%loglik(size(mult)))
        ctx%data = fit%data
        ctx%pp = 1.0_dp/k_blocks
        x0 = [merge(fit%xi,0.01_dp,abs(fit%xi)>1.0e-5_dp),fit%sigma]
        do i = 1,size(mult)
            out%x(i) = estimate*mult(i)
            ctx%return_level = out%x(i)
            call nelder_mead(gev_profile_objective,ctx,x0,xbest,fbest,conv,initial_step=0.08_dp)
            out%loglik(i) = -fbest
            if (conv) x0 = xbest
        end do
        cutoff = -fit%nllh-0.5_dp*chi_square_quantile_df1(ci)
        call profile_interval(out%x,out%loglik,cutoff,out%lower,out%upper)
        if (ieee_is_nan(out%lower) .or. ieee_is_nan(out%upper)) out%status = evir_domain_error
    end function gev_return_level_profile

    pure real(dp) function gpd_quantile_estimate(fit,p,tail_scale) result(q)
        type(gpd_fit_result), intent(in) :: fit
        real(dp), intent(in) :: p
        logical, intent(in), optional :: tail_scale
        real(dp) :: lambda,a
        logical :: tail
        tail = .true.
        if (present(tail_scale)) tail = tail_scale
        lambda = 1.0_dp
        if (tail) lambda = 1.0_dp/(1.0_dp-fit%p_less_threshold)
        a = lambda*(1.0_dp-p)
        if (a <= 0.0_dp .or. fit%beta <= 0.0_dp) then
            q = safe_nan()
        else if (abs(fit%xi) <= 1.0e-8_dp) then
            q = fit%threshold-fit%beta*log(a)
        else
            q = fit%threshold+fit%beta*(a**(-fit%xi)-1.0_dp)/fit%xi
        end if
    end function gpd_quantile_estimate

    subroutine gpd_quantile_wald(fit,p,lower,estimate,se,upper,ci_p,tail_scale,status)
        type(gpd_fit_result), intent(in) :: fit
        real(dp), intent(in) :: p
        real(dp), intent(out) :: lower,estimate,se,upper
        real(dp), intent(in), optional :: ci_p
        logical, intent(in), optional :: tail_scale
        integer, intent(out), optional :: status
        real(dp) :: ci,lambda,a,g,gp,varq,z
        logical :: tail
        ci = 0.95_dp
        if (present(ci_p)) ci = ci_p
        tail = .true.
        if (present(tail_scale)) tail = tail_scale
        lambda = merge(1.0_dp/(1.0_dp-fit%p_less_threshold),1.0_dp,tail)
        a = lambda*(1.0_dp-p)
        estimate = gpd_quantile_estimate(fit,p,tail)
        if (abs(fit%xi) <= 1.0e-7_dp) then
            g = -log(a)
            gp = 0.5_dp*log(a)**2
        else
            g = (a**(-fit%xi)-1.0_dp)/fit%xi
            gp = (-(a**(-fit%xi)-1.0_dp)/fit%xi-a**(-fit%xi)*log(a))/fit%xi
        end if
        varq = g*g*fit%varcov(2,2)+(fit%beta*gp)**2*fit%varcov(1,1)+ &
            2.0_dp*g*fit%beta*gp*fit%varcov(1,2)
        if (varq < 0.0_dp) then
            se = safe_nan(); lower=se; upper=se
            if (present(status)) status=evir_domain_error
            return
        end if
        se = sqrt(varq)
        z = normal_quantile(1.0_dp-(1.0_dp-ci)/2.0_dp)
        lower = estimate-z*se
        upper = estimate+z*se
        if (present(status)) status=evir_ok
    end subroutine gpd_quantile_wald

    function gpd_quantile_profile(fit,p,plotmax,ci_p,n_grid,tail_scale) result(out)
        type(gpd_fit_result), intent(in) :: fit
        real(dp), intent(in) :: p,plotmax
        real(dp), intent(in), optional :: ci_p
        integer, intent(in), optional :: n_grid
        logical, intent(in), optional :: tail_scale
        type(profile_result) :: out
        type(gpd_profile_context) :: ctx
        integer :: ng,i
        real(dp) :: ci,lambda,xlo,x0(1),xbest(1),fbest,cutoff
        logical :: tail,conv
        ci=0.95_dp; if(present(ci_p)) ci=ci_p
        ng=50; if(present(n_grid)) ng=n_grid
        tail=.true.; if(present(tail_scale)) tail=tail_scale
        lambda=merge(1.0_dp/(1.0_dp-fit%p_less_threshold),1.0_dp,tail)
        ctx%a=lambda*(1.0_dp-p)
        ctx%threshold=fit%threshold
        ctx%excess=fit%exceedances-fit%threshold
        ctx%shortfall=.false.
        out%estimate=gpd_quantile_estimate(fit,p,tail)
        if(plotmax<=fit%threshold .or. ng<5) then
            out%status=evir_invalid_input; return
        end if
        allocate(out%x(ng),out%loglik(ng))
        xlo=max(fit%threshold+epsilon(1.0_dp),min(out%estimate,fit%threshold*1.0001_dp))
        if(xlo<=0.0_dp) then
            do i=1,ng
                out%x(i)=xlo+(plotmax-xlo)*real(i-1,dp)/real(ng-1,dp)
            end do
        else
            do i=1,ng
                out%x(i)=exp(log(xlo)+(log(plotmax)-log(xlo))*real(i-1,dp)/real(ng-1,dp))
            end do
        end if
        x0(1)=fit%xi
        do i=1,ng
            ctx%target=out%x(i)
            call nelder_mead(gpd_profile_objective,ctx,x0,xbest,fbest,conv,initial_step=0.05_dp)
            out%loglik(i)=-fbest
            if(conv) x0=xbest
        end do
        cutoff=-fit%nllh-0.5_dp*chi_square_quantile_df1(ci)
        call profile_interval(out%x,out%loglik,cutoff,out%lower,out%upper)
    end function gpd_quantile_profile

    pure real(dp) function gpd_shortfall_estimate(fit,p,tail_scale) result(es)
        type(gpd_fit_result), intent(in) :: fit
        real(dp), intent(in) :: p
        logical, intent(in), optional :: tail_scale
        real(dp) :: q
        q=gpd_quantile_estimate(fit,p,tail_scale)
        if(fit%xi>=1.0_dp) then
            es=huge(1.0_dp)
        else
            es=q+(fit%beta+fit%xi*(q-fit%threshold))/(1.0_dp-fit%xi)
        end if
    end function gpd_shortfall_estimate

    function gpd_shortfall_profile(fit,p,plotmax,ci_p,n_grid,tail_scale) result(out)
        type(gpd_fit_result), intent(in) :: fit
        real(dp), intent(in) :: p,plotmax
        real(dp), intent(in), optional :: ci_p
        integer, intent(in), optional :: n_grid
        logical, intent(in), optional :: tail_scale
        type(profile_result) :: out
        type(gpd_profile_context) :: ctx
        integer :: ng,i
        real(dp) :: ci,lambda,xlo,x0(1),xbest(1),fbest,cutoff
        logical :: tail,conv
        ci=0.95_dp; if(present(ci_p)) ci=ci_p
        ng=50; if(present(n_grid)) ng=n_grid
        tail=.true.; if(present(tail_scale)) tail=tail_scale
        lambda=merge(1.0_dp/(1.0_dp-fit%p_less_threshold),1.0_dp,tail)
        ctx%a=lambda*(1.0_dp-p); ctx%threshold=fit%threshold
        ctx%excess=fit%exceedances-fit%threshold; ctx%shortfall=.true.
        out%estimate=gpd_shortfall_estimate(fit,p,tail)
        if(plotmax<=fit%threshold .or. ng<5) then
            out%status=evir_invalid_input; return
        end if
        allocate(out%x(ng),out%loglik(ng))
        xlo=max(fit%threshold+epsilon(1.0_dp),min(out%estimate,fit%threshold*1.0001_dp))
        if(xlo<=0.0_dp) then
            do i=1,ng; out%x(i)=xlo+(plotmax-xlo)*real(i-1,dp)/real(ng-1,dp); end do
        else
            do i=1,ng; out%x(i)=exp(log(xlo)+(log(plotmax)-log(xlo))*real(i-1,dp)/real(ng-1,dp)); end do
        end if
        x0(1)=min(fit%xi,0.95_dp)
        do i=1,ng
            ctx%target=out%x(i)
            call nelder_mead(gpd_profile_objective,ctx,x0,xbest,fbest,conv,initial_step=0.05_dp)
            out%loglik(i)=-fbest
            if(conv) x0=xbest
        end do
        cutoff=-fit%nllh-0.5_dp*chi_square_quantile_df1(ci)
        call profile_interval(out%x,out%loglik,cutoff,out%lower,out%upper)
    end function gpd_shortfall_profile

    subroutine risk_measures(fit,p,quantile,shortfall)
        type(gpd_fit_result), intent(in) :: fit
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: quantile(size(p)),shortfall(size(p))
        integer :: i
        do i=1,size(p)
            quantile(i)=gpd_quantile_estimate(fit,p(i),.true.)
            shortfall(i)=gpd_shortfall_estimate(fit,p(i),.true.)
        end do
    end subroutine risk_measures

    function gev_objective(x,context) result(f)
        real(dp),intent(in)::x(:)
        class(*),intent(in)::context
        real(dp)::f
        select type(context)
        type is(vector_context); f=gev_negloglik_value(x,context%data)
        class default; f=huge(1.0_dp)/100.0_dp
        end select
    end function gev_objective

    function gumbel_objective(x,context) result(f)
        real(dp),intent(in)::x(:)
        class(*),intent(in)::context
        real(dp)::f
        select type(context)
        type is(vector_context); f=gumbel_negloglik_value(x,context%data)
        class default; f=huge(1.0_dp)/100.0_dp
        end select
    end function gumbel_objective

    function gpd_objective(x,context) result(f)
        real(dp),intent(in)::x(:)
        class(*),intent(in)::context
        real(dp)::f
        select type(context)
        type is(vector_context); f=gpd_negloglik_value(x,context%data)
        class default; f=huge(1.0_dp)/100.0_dp
        end select
    end function gpd_objective

    function pot_objective(x,context) result(f)
        real(dp),intent(in)::x(:)
        class(*),intent(in)::context
        real(dp)::f
        select type(context)
        type is(pot_context); f=pot_negloglik_value(x,context%data,context%threshold,context%span)
        class default; f=huge(1.0_dp)/100.0_dp
        end select
    end function pot_objective

    function gev_profile_objective(x,context) result(f)
        real(dp),intent(in)::x(:)
        class(*),intent(in)::context
        real(dp)::f,xi,sigma,mu,a
        real(dp)::theta(3)
        select type(context)
        type is(gev_profile_context)
            xi=x(1); sigma=x(2)
            if(sigma<=0.0_dp) then; f=huge(1.0_dp)/100.0_dp; return; end if
            a=-log(1.0_dp-context%pp)
            if(abs(xi)<=1.0e-7_dp) then
                mu=context%return_level+sigma*log(a)
            else
                mu=context%return_level-sigma*(a**(-xi)-1.0_dp)/xi
            end if
            theta=[xi,sigma,mu]
            f=gev_negloglik_value(theta,context%data)
        class default
            f=huge(1.0_dp)/100.0_dp
        end select
    end function gev_profile_objective

    function gpd_profile_objective(x,context) result(f)
        real(dp),intent(in)::x(:)
        class(*),intent(in)::context
        real(dp)::f,xi,beta,denom
        real(dp)::theta(2)
        select type(context)
        type is(gpd_profile_context)
            xi=x(1)
            if(context%shortfall) then
                if(xi>=1.0_dp) then; f=huge(1.0_dp)/100.0_dp; return; end if
                if(abs(xi)<=1.0e-7_dp) then
                    beta=(context%target-context%threshold)/(1.0_dp-log(context%a))
                else
                    denom=((context%a**(-xi)-1.0_dp)/xi)+1.0_dp
                    beta=(1.0_dp-xi)*(context%target-context%threshold)/denom
                end if
            else
                if(abs(xi)<=1.0e-7_dp) then
                    beta=(context%target-context%threshold)/(-log(context%a))
                else
                    beta=xi*(context%target-context%threshold)/(context%a**(-xi)-1.0_dp)
                end if
            end if
            theta=[xi,beta]
            f=gpd_negloglik_value(theta,context%excess)
        class default
            f=huge(1.0_dp)/100.0_dp
        end select
    end function gpd_profile_objective

    pure function diagonal3(a) result(d)
        real(dp),intent(in)::a(3,3)
        real(dp)::d(3)
        d=[a(1,1),a(2,2),a(3,3)]
    end function diagonal3

    pure character(len=8) function lowercase(text) result(out)
        character(len=*),intent(in)::text
        integer::i,c,n
        out=' '
        n=min(len_trim(text),len(out))
        do i=1,n
            c=iachar(text(i:i))
            if(c>=iachar('A').and.c<=iachar('Z')) c=c+32
            out(i:i)=achar(c)
        end do
    end function lowercase

    subroutine profile_interval(x,ll,cutoff,lower,upper)
        real(dp),intent(in)::x(:),ll(:),cutoff
        real(dp),intent(out)::lower,upper
        integer::i,first,last
        real(dp)::w
        first=0; last=0
        do i=1,size(x)
            if(ll(i)>=cutoff) then
                if(first==0) first=i
                last=i
            end if
        end do
        if(first==0) then
            lower=safe_nan(); upper=safe_nan(); return
        end if
        lower=x(first); upper=x(last)
        if(first>1.and.abs(ll(first)-ll(first-1))>epsilon(1.0_dp)) then
            w=(cutoff-ll(first-1))/(ll(first)-ll(first-1))
            lower=x(first-1)+w*(x(first)-x(first-1))
        end if
        if(last<size(x).and.abs(ll(last+1)-ll(last))>epsilon(1.0_dp)) then
            w=(cutoff-ll(last))/(ll(last+1)-ll(last))
            upper=x(last)+w*(x(last+1)-x(last))
        end if
    end subroutine profile_interval

end module evir_fitting
