! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_math
    use r_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: log_two_pi = 1.8378770664093454835606594728112353_dp

    public :: normal_logpdf
    public :: normal_cdf
    public :: normal_quantile
    public :: logistic
    public :: log1pexp
    public :: logsumexp
    public :: poisson_logpmf
    public :: geometric_logpmf
    public :: exponential_logpdf
    public :: weibull_logpdf
    public :: family_loglik

contains

    pure elemental real(dp) function normal_logpdf(x, mean_value, sd_value) result(value)
        real(dp), intent(in) :: x !! Point at which the normal log-density is evaluated.
        real(dp), intent(in) :: mean_value !! Normal distribution mean.
        real(dp), intent(in) :: sd_value !! Positive normal standard deviation.
        real(dp) :: z

        if (sd_value <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
        end if
        z = (x - mean_value) / sd_value
        value = -0.5_dp * (log_two_pi + 2.0_dp * log(sd_value) + z * z)
    end function normal_logpdf

    pure elemental real(dp) function normal_cdf(x) result(value)
        real(dp), intent(in) :: x !! Standard-normal quantile at which the cumulative probability is evaluated.

        value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
    end function normal_cdf

    pure elemental real(dp) function normal_quantile(probability) result(value)
        real(dp), intent(in) :: probability !! Standard-normal cumulative probability in the closed interval zero to one.
        real(dp) :: lower
        real(dp) :: midpoint
        real(dp) :: upper
        integer :: iteration

        if (probability <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
        end if
        if (probability >= 1.0_dp) then
            value = huge(1.0_dp)
            return
        end if
        lower = -1.0_dp
        upper = 1.0_dp
        do while (normal_cdf(lower) > probability)
            lower = 2.0_dp * lower
        end do
        do while (normal_cdf(upper) < probability)
            upper = 2.0_dp * upper
        end do
        do iteration = 1, 120
            midpoint = 0.5_dp * (lower + upper)
            if (normal_cdf(midpoint) < probability) then
                lower = midpoint
            else
                upper = midpoint
            end if
        end do
        value = 0.5_dp * (lower + upper)
    end function normal_quantile

    pure elemental real(dp) function logistic(x) result(value)
        real(dp), intent(in) :: x !! Real-valued logit.

        if (x >= 0.0_dp) then
            value = 1.0_dp / (1.0_dp + exp(-x))
        else
            value = exp(x) / (1.0_dp + exp(x))
        end if
    end function logistic

    pure elemental real(dp) function log1pexp(x) result(value)
        real(dp), intent(in) :: x !! Real argument in log(1 + exp(x)).

        if (x > 0.0_dp) then
            value = x + log(1.0_dp + exp(-x))
        else
            value = log(1.0_dp + exp(x))
        end if
    end function log1pexp

    pure real(dp) function logsumexp(x) result(value)
        real(dp), intent(in) :: x(:) !! Vector whose exponentials are summed in logarithmic space.
        real(dp) :: maximum_value

        if (size(x) == 0) then
            value = -huge(1.0_dp)
            return
        end if
        maximum_value = maxval(x)
        value = maximum_value + log(sum(exp(x - maximum_value)))
    end function logsumexp

    pure elemental real(dp) function poisson_logpmf(count, lambda) result(value)
        integer, intent(in) :: count !! Nonnegative Poisson count.
        real(dp), intent(in) :: lambda !! Positive Poisson mean.

        if (count < 0 .or. lambda <= 0.0_dp) then
            value = -huge(1.0_dp)
        else
            value = real(count, dp) * log(lambda) - lambda - log_gamma(real(count + 1, dp))
        end if
    end function poisson_logpmf

    pure elemental real(dp) function geometric_logpmf(count, mean_value) result(value)
        integer, intent(in) :: count !! Nonnegative geometric count using the R geometric convention.
        real(dp), intent(in) :: mean_value !! Positive mean count, so success probability is 1/(1+mean).
        real(dp) :: probability

        if (count < 0 .or. mean_value < 0.0_dp) then
            value = -huge(1.0_dp)
            return
        end if
        probability = 1.0_dp / (1.0_dp + mean_value)
        value = log(probability) + real(count, dp) * log(1.0_dp - probability)
    end function geometric_logpmf

    pure elemental real(dp) function exponential_logpdf(x, rate) result(value)
        real(dp), intent(in) :: x !! Nonnegative exponential observation.
        real(dp), intent(in) :: rate !! Positive exponential rate.

        if (x < 0.0_dp .or. rate <= 0.0_dp) then
            value = -huge(1.0_dp)
        else
            value = log(rate) - rate * x
        end if
    end function exponential_logpdf

    pure elemental real(dp) function weibull_logpdf(x, shape, scale) result(value)
        real(dp), intent(in) :: x !! Positive Weibull observation.
        real(dp), intent(in) :: shape !! Positive Weibull shape.
        real(dp), intent(in) :: scale !! Positive Weibull scale.

        if (x <= 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
            value = -huge(1.0_dp)
        else
            value = log(shape) - log(scale) + (shape - 1.0_dp) * log(x / scale) - (x / scale) ** shape
        end if
    end function weibull_logpdf

    pure real(dp) function family_loglik(family, y, eta, extra) result(value)
        integer, intent(in) :: family !! MCMCglmm native family code for the supported scalar likelihood subset.
        real(dp), intent(in) :: y !! Observed scalar response on the natural response scale.
        real(dp), intent(in) :: eta !! Linear predictor; interpretation follows the selected family.
        real(dp), intent(in), optional :: extra !! Optional positive dispersion, shape, or scale parameter used by some families.
        real(dp) :: parameter
        real(dp) :: probability

        parameter = 1.0_dp
        if (present(extra)) parameter = extra

        select case (family)
        case (1)
            value = normal_logpdf(y, eta, sqrt(max(parameter, tiny(1.0_dp))))
        case (2)
            value = poisson_logpmf(nint(y), exp(eta))
        case (4)
            value = weibull_logpdf(y, max(parameter, tiny(1.0_dp)), exp(eta))
        case (5)
            value = exponential_logpdf(y, exp(eta))
        case (17)
            value = geometric_logpmf(nint(y), exp(-eta))
        case (19)
            probability = logistic(eta)
            if (nint(y) == 0) then
                value = log(max(1.0_dp - probability, tiny(1.0_dp)))
            else if (nint(y) == 1) then
                value = log(max(probability, tiny(1.0_dp)))
            else
                value = -huge(1.0_dp)
            end if
        case default
            value = -huge(1.0_dp)
        end select
    end function family_loglik

end module mcmcglmm_math
