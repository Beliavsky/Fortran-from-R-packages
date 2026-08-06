module lme4_family
   use lme4_kinds, only : dp, pi
   implicit none
   private

   abstract interface
      function link_scalar(mu) result(eta)
         import :: dp
         real(dp), intent(in) :: mu
         real(dp) :: eta
      end function link_scalar

      function inverse_link_scalar(eta) result(mu)
         import :: dp
         real(dp), intent(in) :: eta
         real(dp) :: mu
      end function inverse_link_scalar

      function derivative_scalar(eta) result(value)
         import :: dp
         real(dp), intent(in) :: eta
         real(dp) :: value
      end function derivative_scalar

      function variance_scalar(mu, dispersion) result(value)
         import :: dp
         real(dp), intent(in) :: mu, dispersion
         real(dp) :: value
      end function variance_scalar

      function loglik_scalar(y, mu, weight, dispersion) result(value)
         import :: dp
         real(dp), intent(in) :: y, mu, weight, dispersion
         real(dp) :: value
      end function loglik_scalar

      function valid_scalar(y) result(ok)
         import :: dp
         real(dp), intent(in) :: y
         logical :: ok
      end function valid_scalar
   end interface

   type, public :: family_spec_t
      character(len=:), allocatable :: name
      procedure(link_scalar), pointer, nopass :: link => null()
      procedure(inverse_link_scalar), pointer, nopass :: inverse_link => null()
      procedure(derivative_scalar), pointer, nopass :: dmu_deta => null()
      procedure(variance_scalar), pointer, nopass :: variance => null()
      procedure(loglik_scalar), pointer, nopass :: log_likelihood => null()
      procedure(valid_scalar), pointer, nopass :: valid_response => null()
   contains
      procedure :: validate => validate_family_spec
   end type family_spec_t

   public :: gaussian_identity_family, binomial_probit_family
   public :: binomial_cloglog_family, quasipoisson_log_family
   public :: logistic_scalar, normal_cdf_scalar, normal_pdf_scalar, inverse_normal_cdf

contains

   subroutine validate_family_spec(self, ok, message)
      class(family_spec_t), intent(in) :: self
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: message

      ok = associated(self%link) .and. associated(self%inverse_link) .and. &
         associated(self%dmu_deta) .and. associated(self%variance) .and. &
         associated(self%log_likelihood) .and. associated(self%valid_response)
      if (ok) then
         message = 'ok'
      else
         message = 'custom family requires link, inverse-link, derivative, variance, '// &
            'log-likelihood, and response-validation callbacks'
      end if
   end subroutine validate_family_spec

   function gaussian_identity_family() result(family)
      type(family_spec_t) :: family
      family%name = 'gaussian(identity)'
      family%link => identity_link
      family%inverse_link => identity_inverse
      family%dmu_deta => identity_derivative
      family%variance => gaussian_variance
      family%log_likelihood => gaussian_loglik
      family%valid_response => finite_response
   end function gaussian_identity_family

   function binomial_probit_family() result(family)
      type(family_spec_t) :: family
      family%name = 'binomial(probit)'
      family%link => probit_link
      family%inverse_link => probit_inverse
      family%dmu_deta => probit_derivative
      family%variance => binomial_variance
      family%log_likelihood => binomial_loglik
      family%valid_response => binomial_response
   end function binomial_probit_family

   function binomial_cloglog_family() result(family)
      type(family_spec_t) :: family
      family%name = 'binomial(cloglog)'
      family%link => cloglog_link
      family%inverse_link => cloglog_inverse
      family%dmu_deta => cloglog_derivative
      family%variance => binomial_variance
      family%log_likelihood => binomial_loglik
      family%valid_response => binomial_response
   end function binomial_cloglog_family

   function quasipoisson_log_family() result(family)
      type(family_spec_t) :: family
      family%name = 'quasipoisson(log)'
      family%link => log_link
      family%inverse_link => log_inverse
      family%dmu_deta => log_derivative
      family%variance => poisson_variance
      family%log_likelihood => poisson_loglik
      family%valid_response => count_response
   end function quasipoisson_log_family

   pure real(dp) function identity_link(mu) result(eta)
      real(dp), intent(in) :: mu
      eta = mu
   end function identity_link

   pure real(dp) function identity_inverse(eta) result(mu)
      real(dp), intent(in) :: eta
      mu = eta
   end function identity_inverse

   pure real(dp) function identity_derivative(eta) result(value)
      real(dp), intent(in) :: eta
      value = 1.0_dp + 0.0_dp*eta
   end function identity_derivative

   pure real(dp) function gaussian_variance(mu, dispersion) result(value)
      real(dp), intent(in) :: mu, dispersion
      value = max(tiny(1.0_dp),dispersion) + 0.0_dp*mu
   end function gaussian_variance

   pure real(dp) function gaussian_loglik(y, mu, weight, dispersion) result(value)
      real(dp), intent(in) :: y, mu, weight, dispersion
      real(dp) :: variance
      variance = max(tiny(1.0_dp),dispersion/weight)
      value = -0.5_dp*(log(2.0_dp*pi*variance)+(y-mu)**2/variance)
   end function gaussian_loglik

   pure logical function finite_response(y) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: y
      ok = ieee_is_finite(y)
   end function finite_response

   pure real(dp) function probit_link(mu) result(eta)
      real(dp), intent(in) :: mu
      eta = inverse_normal_cdf(min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,mu)))
   end function probit_link

   pure real(dp) function probit_inverse(eta) result(mu)
      real(dp), intent(in) :: eta
      mu = min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,normal_cdf_scalar(eta)))
   end function probit_inverse

   pure real(dp) function probit_derivative(eta) result(value)
      real(dp), intent(in) :: eta
      value = max(1.0e-12_dp,normal_pdf_scalar(eta))
   end function probit_derivative

   pure real(dp) function cloglog_link(mu) result(eta)
      real(dp), intent(in) :: mu
      eta = log(-log(max(1.0e-12_dp,1.0_dp-min(1.0_dp-1.0e-12_dp,mu))))
   end function cloglog_link

   pure real(dp) function cloglog_inverse(eta) result(mu)
      real(dp), intent(in) :: eta
      real(dp) :: e
      e = exp(min(40.0_dp,max(-40.0_dp,eta)))
      mu = min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,1.0_dp-exp(-e)))
   end function cloglog_inverse

   pure real(dp) function cloglog_derivative(eta) result(value)
      real(dp), intent(in) :: eta
      real(dp) :: e
      e = exp(min(40.0_dp,max(-40.0_dp,eta)))
      value = max(1.0e-12_dp,e*exp(-e))
   end function cloglog_derivative

   pure real(dp) function log_link(mu) result(eta)
      real(dp), intent(in) :: mu
      eta = log(max(1.0e-12_dp,mu))
   end function log_link

   pure real(dp) function log_inverse(eta) result(mu)
      real(dp), intent(in) :: eta
      mu = exp(min(30.0_dp,max(-30.0_dp,eta)))
   end function log_inverse

   pure real(dp) function log_derivative(eta) result(value)
      real(dp), intent(in) :: eta
      value = log_inverse(eta)
   end function log_derivative

   pure real(dp) function binomial_variance(mu, dispersion) result(value)
      real(dp), intent(in) :: mu, dispersion
      value = max(1.0e-12_dp,mu*(1.0_dp-mu))*max(1.0_dp,dispersion)
   end function binomial_variance

   pure real(dp) function poisson_variance(mu, dispersion) result(value)
      real(dp), intent(in) :: mu, dispersion
      value = max(1.0e-12_dp,mu)*max(1.0_dp,dispersion)
   end function poisson_variance

   pure real(dp) function binomial_loglik(y, mu, weight, dispersion) result(value)
      real(dp), intent(in) :: y, mu, weight, dispersion
      value = weight*(y*log(max(1.0e-12_dp,mu)) + &
         (1.0_dp-y)*log(max(1.0e-12_dp,1.0_dp-mu)))/max(1.0_dp,dispersion)
   end function binomial_loglik

   pure real(dp) function poisson_loglik(y, mu, weight, dispersion) result(value)
      real(dp), intent(in) :: y, mu, weight, dispersion
      value = weight*(y*log(max(1.0e-12_dp,mu))-mu-log_gamma(y+1.0_dp))/ &
         max(1.0_dp,dispersion)
   end function poisson_loglik

   pure logical function binomial_response(y) result(ok)
      real(dp), intent(in) :: y
      ok = y >= 0.0_dp .and. y <= 1.0_dp
   end function binomial_response

   pure logical function count_response(y) result(ok)
      real(dp), intent(in) :: y
      ok = y >= 0.0_dp
   end function count_response

   pure real(dp) function logistic_scalar(x) result(value)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp+exp(-min(700.0_dp,x)))
      else
         value = exp(max(-700.0_dp,x))/(1.0_dp+exp(max(-700.0_dp,x)))
      end if
   end function logistic_scalar

   pure real(dp) function normal_cdf_scalar(x) result(value)
      real(dp), intent(in) :: x
      value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf_scalar

   pure real(dp) function normal_pdf_scalar(x) result(value)
      real(dp), intent(in) :: x
      value = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf_scalar

   pure real(dp) function inverse_normal_cdf(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
         4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
         7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
         2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
      real(dp) :: q, r

      if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p-0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
            (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
      x = x-(normal_cdf_scalar(x)-p)/max(1.0e-15_dp,normal_pdf_scalar(x))
   end function inverse_normal_cdf

end module lme4_family
