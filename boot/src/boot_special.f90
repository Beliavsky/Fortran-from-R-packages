module boot_special
    use boot_kinds, only : dp
    implicit none
    private
    public :: normal_cdf, normal_pdf, normal_quantile, logistic_cdf, logit_fn

contains

    elemental real(dp) function normal_pdf(x) result(y)
        real(dp), intent(in) :: x
        y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*acos(-1.0_dp))
    end function normal_pdf

    elemental real(dp) function normal_cdf(x) result(p)
        real(dp), intent(in) :: x
        p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
    end function normal_cdf

    elemental real(dp) function logistic_cdf(x) result(p)
        real(dp), intent(in) :: x
        if (x >= 0.0_dp) then
            p = 1.0_dp/(1.0_dp + exp(-x))
        else
            p = exp(x)/(1.0_dp + exp(x))
        end if
    end function logistic_cdf

    elemental real(dp) function logit_fn(p) result(x)
        real(dp), intent(in) :: p
        if (p <= 0.0_dp) then
            x = -huge(1.0_dp)
        else if (p >= 1.0_dp) then
            x = huge(1.0_dp)
        else
            x = log(p) - log(1.0_dp-p)
        end if
    end function logit_fn

    elemental real(dp) function normal_quantile(p) result(x)
        real(dp), intent(in) :: p
        real(dp) :: q, r
        real(dp), parameter :: a(6) = [ &
            -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
            -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
            -3.066479806614716e+01_dp, 2.506628277459239e+00_dp ]
        real(dp), parameter :: b(5) = [ &
            -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
            -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
            -1.328068155288572e+01_dp ]
        real(dp), parameter :: c(6) = [ &
            -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
            -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
             4.374664141464968e+00_dp,  2.938163982698783e+00_dp ]
        real(dp), parameter :: d(4) = [ &
             7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
             2.445134137142996e+00_dp, 3.754408661907416e+00_dp ]
        real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
        if (p <= 0.0_dp) then
            x = -huge(1.0_dp)
        else if (p >= 1.0_dp) then
            x = huge(1.0_dp)
        else if (p < plow) then
            q = sqrt(-2.0_dp*log(p))
            x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
                ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
        else if (p <= phigh) then
            q = p-0.5_dp
            r = q*q
            x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
                (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
        else
            q = sqrt(-2.0_dp*log(1.0_dp-p))
            x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
                 ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
        end if
        if (abs(x) < huge(1.0_dp)/10.0_dp) then
            x = x - (normal_cdf(x)-p)/normal_pdf(x)
        end if
    end function normal_quantile
end module boot_special
