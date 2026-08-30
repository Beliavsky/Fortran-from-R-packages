module learnbayes_special
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   use learnbayes_kinds, only: dp
   implicit none
   private

   real(dp), parameter :: pi_dp = acos(-1.0_dp)
   real(dp), parameter :: sqrt2_dp = sqrt(2.0_dp)

   public :: beta_log
   public :: log1pexp
   public :: logistic
   public :: normal_cdf
   public :: normal_quantile
   public :: regularized_beta
   public :: quantile_type7

contains

   pure function beta_log(a, b) result(value)
      real(dp), intent(in) :: a !! First positive beta-function shape parameter.
      real(dp), intent(in) :: b !! Second positive beta-function shape parameter.
      real(dp) :: value

      value = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
   end function beta_log

   pure function log1pexp(x) result(value)
      real(dp), intent(in) :: x !! Real argument for a stable evaluation of log(1 + exp(x)).
      real(dp) :: value

      if (x > 0.0_dp) then
         value = x + log(1.0_dp + exp(-x))
      else
         value = log(1.0_dp + exp(x))
      end if
   end function log1pexp

   pure function logistic(x) result(value)
      real(dp), intent(in) :: x !! Real-valued log-odds transformed to a probability in [0, 1].
      real(dp) :: value

      if (x >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp + exp(-x))
      else
         value = exp(x)/(1.0_dp + exp(x))
      end if
   end function logistic

   pure function normal_cdf(x, mean, sd) result(value)
      real(dp), intent(in) :: x !! Point at which the Gaussian cumulative probability is evaluated.
      real(dp), intent(in), optional :: mean !! Optional Gaussian mean; defaults to zero.
      real(dp), intent(in), optional :: sd !! Optional positive Gaussian standard deviation; defaults to one.
      real(dp) :: value
      real(dp) :: mu
      real(dp) :: sigma

      mu = 0.0_dp
      if (present(mean)) mu = mean
      sigma = 1.0_dp
      if (present(sd)) sigma = sd

      if (sigma <= 0.0_dp) then
         value = merge(0.0_dp, 1.0_dp, x < mu)
         return
      end if
      if (x <= -0.25_dp*huge(1.0_dp)) then
         value = 0.0_dp
         return
      else if (x >= 0.25_dp*huge(1.0_dp)) then
         value = 1.0_dp
         return
      end if
      value = 0.5_dp*erfc(-(x - mu)/(sigma*sqrt2_dp))
   end function normal_cdf

   pure function normal_quantile(p, mean, sd) result(value)
      real(dp), intent(in) :: p !! Cumulative probability in [0, 1]; endpoints map to infinities.
      real(dp), intent(in), optional :: mean !! Optional Gaussian mean; defaults to zero.
      real(dp), intent(in), optional :: sd !! Optional nonnegative Gaussian standard deviation; defaults to one.
      real(dp) :: value
      real(dp) :: mu
      real(dp) :: sigma
      real(dp) :: q
      real(dp) :: r
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
         -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
         -3.066479806614716e+01_dp, 2.506628277459239e+00_dp]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
         -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
         -1.328068155288572e+01_dp]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
         -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
         4.374664141464968e+00_dp, 2.938163982698783e+00_dp]
      real(dp), parameter :: d(4) = [ &
         7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
         2.445134137142996e+00_dp, 3.754408661907416e+00_dp]
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow

      mu = 0.0_dp
      if (present(mean)) mu = mean
      sigma = 1.0_dp
      if (present(sd)) sigma = sd

      if (p <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if

      if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         r = (((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6))/ &
            ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         r = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6))*q/ &
            (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + 1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp - p))
         r = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6))/ &
            ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0_dp)
      end if

      ! One Halley correction brings the rational approximation close to full double precision.
      q = 0.5_dp*erfc(-r/sqrt2_dp) - p
      r = r - q*sqrt(2.0_dp*pi_dp)*exp(0.5_dp*r*r)/(1.0_dp + 0.5_dp*r*q*sqrt(2.0_dp*pi_dp)*exp(0.5_dp*r*r))
      value = mu + sigma*r
   end function normal_quantile

   pure function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x !! Evaluation point for the beta CDF; values outside [0, 1] saturate.
      real(dp), intent(in) :: a !! Positive first beta shape parameter.
      real(dp), intent(in) :: b !! Positive second beta shape parameter.
      real(dp) :: value
      real(dp) :: bt

      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      end if
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a*log(x) + b*log(1.0_dp - x))
      if (x < (a + 1.0_dp)/(a + b + 2.0_dp)) then
         value = bt*beta_cont_frac(a, b, x)/a
      else
         value = 1.0_dp - bt*beta_cont_frac(b, a, 1.0_dp - x)/b
      end if
      value = max(0.0_dp, min(1.0_dp, value))
   end function regularized_beta

   pure function beta_cont_frac(a, b, x) result(value)
      real(dp), intent(in) :: a !! Positive first beta shape used in the continued fraction.
      real(dp), intent(in) :: b !! Positive second beta shape used in the continued fraction.
      real(dp), intent(in) :: x !! Fraction argument strictly between zero and one.
      real(dp) :: value
      integer :: m
      integer :: m2
      real(dp) :: aa
      real(dp) :: c
      real(dp) :: d
      real(dp) :: del
      real(dp) :: h
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps

      c = 1.0_dp
      d = 1.0_dp - (a + b)*x/(a + 1.0_dp)
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, 200
         m2 = 2*m
         aa = real(m, dp)*(b - real(m, dp))*x/((a + real(m2 - 1, dp))*(a + real(m2, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a + real(m, dp))*(a + b + real(m, dp))*x/ &
            ((a + real(m2, dp))*(a + real(m2 + 1, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del - 1.0_dp) <= eps) exit
      end do
      value = h
   end function beta_cont_frac

   function quantile_type7(x, p) result(value)
      real(dp), intent(in) :: x(:) !! Sample values; finite values are sorted internally.
      real(dp), intent(in) :: p !! Requested probability in [0, 1], using R's default type-7 interpolation.
      real(dp) :: value
      real(dp), allocatable :: work(:)
      real(dp) :: h
      real(dp) :: frac
      integer :: j
      integer :: n

      n = size(x)
      if (n == 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(work(n))
      work = x
      call insertion_sort(work)
      if (p <= 0.0_dp) then
         value = work(1)
      else if (p >= 1.0_dp) then
         value = work(n)
      else
         h = 1.0_dp + real(n - 1, dp)*p
         j = int(floor(h))
         frac = h - real(j, dp)
         if (j >= n) then
            value = work(n)
         else
            value = (1.0_dp - frac)*work(j) + frac*work(j + 1)
         end if
      end if
   end function quantile_type7

   pure subroutine insertion_sort(x)
      real(dp), intent(inout) :: x(:) !! Array sorted into nondecreasing order in place.
      integer :: i
      integer :: j
      real(dp) :: key

      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine insertion_sort

end module learnbayes_special
