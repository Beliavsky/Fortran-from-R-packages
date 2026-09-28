module psych_basic
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use psych_distributions, only : normal_cdf, normal_quantile, student_t_cdf
   use psych_stats, only : finite_mean, finite_sd, finite_variance, finite_quantile, finite_median
   implicit none
   private
   public :: alpha_to_r, r_to_alpha, fisher_z, fisher_r, r_to_d, d_to_r, t_to_r, r_to_t
   public :: g_to_r, chi2_to_r, r_to_chi2, r_to_covariance, logistic_value, logit_value
   public :: graded_response_probability
   public :: d_to_overlap, d_to_overlap2, d_to_common_language, d_to_u3
   public :: geometric_mean, harmonic_mean, skewness, excess_kurtosis, sd_value
   public :: mssd_value, rmssd_value, shannon_entropy, winsorize, winsor_mean, winsor_variance, winsor_sd
   public :: phi_coefficient, yule_q, yule_y, yule_inverse, yule_to_phi
   public :: correlation_ci, correlation_p_value, auc_from_table, rmsea_value
   public :: manhattan_distance, minkowski_distance, trace_value
   public :: d_to_t, t_to_d, means_to_d, means_to_t

contains

   pure elemental real(dp) function alpha_to_r(alpha, n_var) result(rbar)
      real(dp), intent(in) :: alpha !! Cronbach alpha coefficient.
      integer, intent(in) :: n_var !! Number of items entering alpha; must exceed one for usual interpretation.
      rbar = alpha/(real(n_var,dp)-alpha*real(n_var-1,dp))
   end function alpha_to_r

   pure elemental real(dp) function r_to_alpha(r, n_var) result(alpha)
      real(dp), intent(in) :: r !! Average inter-item correlation.
      integer, intent(in) :: n_var !! Number of items entering the scale.
      alpha = real(n_var,dp)*r/(1.0_dp+real(n_var-1,dp)*r)
   end function r_to_alpha

   pure elemental real(dp) function fisher_z(rho) result(z)
      real(dp), intent(in) :: rho !! Correlation in the open interval `(-1,1)`.
      z = 0.5_dp*log((1.0_dp+rho)/(1.0_dp-rho))
   end function fisher_z

   pure elemental real(dp) function fisher_r(z) result(rho)
      real(dp), intent(in) :: z !! Fisher z-transformed correlation.
      rho = tanh(z)
   end function fisher_r

   pure elemental real(dp) function r_to_d(rho) result(d)
      real(dp), intent(in) :: rho !! Correlation in `(-1,1)`.
      d = 2.0_dp*rho/sqrt(max(tiny(1.0_dp),1.0_dp-rho*rho))
   end function r_to_d

   pure elemental real(dp) function d_to_r(d) result(rho)
      real(dp), intent(in) :: d !! Standardized mean difference.
      rho = d/sqrt(d*d+4.0_dp)
   end function d_to_r

   pure elemental real(dp) function t_to_r(t, df) result(rho)
      real(dp), intent(in) :: t !! Student-t statistic.
      real(dp), intent(in) :: df !! Positive degrees of freedom.
      rho = sign(1.0_dp,t)*sqrt(t*t/(t*t+df))
   end function t_to_r

   pure elemental real(dp) function r_to_t(rho, n) result(t)
      real(dp), intent(in) :: rho !! Correlation in `(-1,1)`.
      integer, intent(in) :: n !! Sample size; values greater than two are required.
      t = rho*sqrt(real(n-2,dp)/max(tiny(1.0_dp),1.0_dp-rho*rho))
   end function r_to_t

   pure elemental real(dp) function g_to_r(g, df, n) result(rho)
      real(dp), intent(in) :: g !! Standardized effect statistic used by `psych::g2r`.
      real(dp), intent(in) :: df !! Associated degrees of freedom.
      integer, intent(in) :: n !! Sample size.
      rho = sign(1.0_dp,g)*abs(g)/sqrt(g*g+4.0_dp*df/real(n,dp))
   end function g_to_r

   pure elemental real(dp) function chi2_to_r(chi2, n) result(rho)
      real(dp), intent(in) :: chi2 !! Nonnegative chi-square statistic.
      integer, intent(in) :: n !! Sample size used in the phi-equivalent transformation.
      rho = sqrt(max(0.0_dp,chi2)/real(n,dp))
   end function chi2_to_r

   pure elemental real(dp) function r_to_chi2(rho, n) result(chi2)
      real(dp), intent(in) :: rho !! Correlation coefficient.
      integer, intent(in) :: n !! Sample size.
      chi2 = rho*rho*real(n,dp)
   end function r_to_chi2

   pure subroutine r_to_covariance(rho, sigma, covariance)
      real(dp), intent(in) :: rho(:, :) !! Correlation matrix with shape `(p,p)`.
      real(dp), intent(in) :: sigma(:) !! Marginal standard deviations with length `p`.
      real(dp), intent(out) :: covariance(:, :) !! Covariance matrix `diag(sigma)*rho*diag(sigma)`.
      integer :: i, j
      do j = 1, size(rho,2)
         do i = 1, size(rho,1)
            covariance(i,j) = sigma(i)*rho(i,j)*sigma(j)
         end do
      end do
   end subroutine r_to_covariance

   pure elemental real(dp) function logistic_value(x, difficulty, discrimination, lower, upper) result(p)
      real(dp), intent(in) :: x !! Latent score or predictor value.
      real(dp), intent(in) :: difficulty !! Logistic midpoint parameter.
      real(dp), intent(in) :: discrimination !! Logistic slope parameter.
      real(dp), intent(in) :: lower !! Lower asymptote.
      real(dp), intent(in) :: upper !! Upper asymptote.
      p = lower+(upper-lower)/(1.0_dp+exp(discrimination*(difficulty-x)))
   end function logistic_value

   pure real(dp) function graded_response_probability(x, category, thresholds, discrimination, lower, upper) result(p)
      real(dp), intent(in) :: x !! Latent trait value at which the graded-response category probability is evaluated.
      integer, intent(in) :: category !! One-based response category from 1 through `size(thresholds) + 1`.
      real(dp), intent(in) :: thresholds(:) !! Ordered category thresholds separating adjacent graded-response categories.
      real(dp), intent(in), optional :: discrimination !! Logistic slope; defaults to the upstream value `1.5`.
      real(dp), intent(in), optional :: lower !! Lower logistic asymptote; defaults to zero.
      real(dp), intent(in), optional :: upper !! Upper logistic asymptote; defaults to one.
      real(dp) :: a, c, z
      integer :: n_thresholds

      a = 1.5_dp
      c = 0.0_dp
      z = 1.0_dp
      if (present(discrimination)) a = discrimination
      if (present(lower)) c = lower
      if (present(upper)) z = upper
      n_thresholds = size(thresholds)

      if (n_thresholds == 0 .or. category < 1 .or. category > n_thresholds + 1) then
         p = 0.0_dp
      else if (category == 1) then
         p = 1.0_dp - logistic_value(x, thresholds(1), a, c, z)
      else if (category == n_thresholds + 1) then
         p = logistic_value(x, thresholds(n_thresholds), a, c, z)
      else
         p = logistic_value(x, thresholds(category - 1), a, c, z) - &
            logistic_value(x, thresholds(category), a, c, z)
      end if
   end function graded_response_probability

   pure elemental real(dp) function logit_value(p) result(value)
      real(dp), intent(in) :: p !! Probability in the open interval `(0,1)`.
      value = log(p/(1.0_dp-p))
   end function logit_value

   pure elemental real(dp) function d_to_overlap(d) result(value)
      real(dp), intent(in) :: d !! Standardized mean difference.
      value = 2.0_dp*normal_cdf(-abs(d)/2.0_dp)
   end function d_to_overlap

   pure elemental real(dp) function d_to_overlap2(d) result(value)
      real(dp), intent(in) :: d !! Standardized mean difference.
      real(dp) :: overlap
      overlap = d_to_overlap(d)
      value = overlap/(2.0_dp-overlap)
   end function d_to_overlap2

   pure elemental real(dp) function d_to_common_language(d) result(value)
      real(dp), intent(in) :: d !! Standardized mean difference.
      value = normal_cdf(abs(d)/sqrt(2.0_dp))
   end function d_to_common_language

   pure elemental real(dp) function d_to_u3(d) result(value)
      real(dp), intent(in) :: d !! Standardized mean difference.
      value = normal_cdf(abs(d))
   end function d_to_u3

   pure real(dp) function geometric_mean(x) result(value)
      real(dp), intent(in) :: x(:) !! Positive observations; NaNs are omitted.
      real(dp), allocatable :: logs(:)
      integer :: i, n
      n = 0
      allocate(logs(size(x)))
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i)) .and. x(i) > 0.0_dp) then
            n = n + 1
            logs(n) = log(x(i))
         end if
      end do
      if (n == 0) then
         value = 0.0_dp
      else
         value = exp(sum(logs(1:n))/real(n,dp))
      end if
   end function geometric_mean

   pure real(dp) function harmonic_mean(x, keep_zero) result(value)
      real(dp), intent(in) :: x(:) !! Numeric observations; NaNs are omitted.
      logical, intent(in), optional :: keep_zero !! Include zero values when true; default follows upstream and includes them.
      logical :: include_zero
      real(dp) :: sum_inverse
      integer :: i, n
      include_zero = .true.
      if (present(keep_zero)) include_zero = keep_zero
      sum_inverse = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) cycle
         if (abs(x(i)) <= tiny(1.0_dp)) then
            if (include_zero) then
               value = 0.0_dp
               return
            end if
            cycle
         end if
         sum_inverse = sum_inverse+1.0_dp/x(i)
         n = n+1
      end do
      if (n > 0 .and. abs(sum_inverse) > tiny(1.0_dp)) then
         value = real(n,dp)/sum_inverse
      else
         value = 0.0_dp
      end if
   end function harmonic_mean

   pure real(dp) function skewness(x, kind) result(value)
      real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted.
      integer, intent(in), optional :: kind !! Upstream estimator type 1, 2, or 3; default is 3.
      real(dp) :: m, s2, s3, sdv
      integer :: i, n, mode
      mode = 3
      if (present(kind)) mode = kind
      m = finite_mean(x)
      s2 = 0.0_dp
      s3 = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            s2 = s2+(x(i)-m)**2
            s3 = s3+(x(i)-m)**3
            n = n+1
         end if
      end do
      if (n < 3 .or. s2 <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      select case (mode)
      case (1)
         value = sqrt(real(n,dp))*s3/(s2**1.5_dp)
      case (2)
         value = real(n,dp)*sqrt(real(n-1,dp))*s3/(real(n-2,dp)*s2**1.5_dp)
      case default
         sdv = sqrt(s2/real(n-1,dp))
         value = s3/(real(n,dp)*sdv**3)
      end select
   end function skewness

   pure real(dp) function excess_kurtosis(x, kind) result(value)
      real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted.
      integer, intent(in), optional :: kind !! Upstream estimator type 1, 2, or 3; default is 3.
      real(dp) :: m, s2, s4, sdv
      integer :: i, n, mode
      mode = 3
      if (present(kind)) mode = kind
      m = finite_mean(x)
      s2 = 0.0_dp
      s4 = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            s2 = s2+(x(i)-m)**2
            s4 = s4+(x(i)-m)**4
            n = n+1
         end if
      end do
      if (n < 4 .or. s2 <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      select case (mode)
      case (1)
         value = real(n,dp)*s4/(s2*s2)-3.0_dp
      case (2)
         value = real(n*(n+1),dp)*s4/(real((n-1)*(n-2)*(n-3),dp)*(s2/real(n-1,dp))**2) &
            - 3.0_dp*real((n-1)*(n-1),dp)/real((n-2)*(n-3),dp)
      case default
         sdv = sqrt(s2/real(n-1,dp))
         value = s4/(real(n,dp)*sdv**4)-3.0_dp
      end select
   end function excess_kurtosis

   pure real(dp) function sd_value(x) result(value)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaNs are omitted.
      value = finite_sd(x)
   end function sd_value

   pure real(dp) function mssd_value(x, lag) result(value)
      real(dp), intent(in) :: x(:) !! Ordered observations; NaNs in either endpoint omit a difference.
      integer, intent(in), optional :: lag !! Positive successive-difference lag; default is one.
      integer :: i, k, n
      real(dp) :: s
      k = 1
      if (present(lag)) k = max(1,lag)
      s = 0.0_dp
      n = 0
      do i = 1, size(x)-k
         if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(x(i+k))) then
            s = s+(x(i+k)-x(i))**2
            n = n+1
         end if
      end do
      if (n > 0) then
         value = s/real(n,dp)
      else
         value = 0.0_dp
      end if
   end function mssd_value

   pure real(dp) function rmssd_value(x, lag) result(value)
      real(dp), intent(in) :: x(:) !! Ordered observations; NaNs in either endpoint omit a difference.
      integer, intent(in), optional :: lag !! Positive successive-difference lag; default is one.
      if (present(lag)) then
         value = sqrt(max(0.0_dp,mssd_value(x,lag)))
      else
         value = sqrt(max(0.0_dp,mssd_value(x)))
      end if
   end function rmssd_value

   pure real(dp) function shannon_entropy(probabilities, base, normalize) result(value)
      real(dp), intent(in) :: probabilities(:) !! Category counts or nonnegative weights; normalized internally.
      real(dp), intent(in), optional :: base !! Logarithm base; default is two.
      logical, intent(in), optional :: normalize !! Divide by maximum entropy when true.
      real(dp) :: b, total, p, maximum
      logical :: corrected
      integer :: i, k
      b = 2.0_dp
      if (present(base)) b = base
      corrected = .false.
      if (present(normalize)) corrected = normalize
      total = sum(max(probabilities,0.0_dp))
      value = 0.0_dp
      k = 0
      if (total <= 0.0_dp) return
      do i = 1, size(probabilities)
         if (probabilities(i) > 0.0_dp) then
            p = probabilities(i)/total
            value = value-p*log(p)/log(b)
            k = k+1
         end if
      end do
      if (corrected .and. k > 1) then
         maximum = log(real(k,dp))/log(b)
         value = value/maximum
      end if
   end function shannon_entropy

   pure subroutine winsorize(x, trim, result)
      real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are preserved.
      real(dp), intent(in) :: trim !! Symmetric winsorization fraction in `[0,0.5]`.
      real(dp), intent(out) :: result(:) !! Winsorized sample in original order.
      real(dp) :: lower, upper, median
      integer :: i
      lower = finite_quantile(x,max(0.0_dp,min(0.5_dp,trim)))
      upper = finite_quantile(x,1.0_dp-max(0.0_dp,min(0.5_dp,trim)))
      median = finite_median(x)
      result = x
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) cycle
         if (trim >= 0.5_dp) then
            result(i) = median
         else
            result(i) = min(upper,max(lower,x(i)))
         end if
      end do
   end subroutine winsorize

   pure real(dp) function winsor_mean(x, trim) result(value)
      real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted.
      real(dp), intent(in), optional :: trim !! Symmetric winsorization fraction; default `0.2`.
      real(dp), allocatable :: w(:)
      real(dp) :: fraction
      fraction = 0.2_dp
      if (present(trim)) fraction = trim
      allocate(w(size(x)))
      call winsorize(x,fraction,w)
      value = finite_mean(w)
   end function winsor_mean

   pure real(dp) function winsor_variance(x, trim) result(value)
      real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted.
      real(dp), intent(in), optional :: trim !! Symmetric winsorization fraction; default `0.2`.
      real(dp), allocatable :: w(:)
      real(dp) :: fraction
      fraction = 0.2_dp
      if (present(trim)) fraction = trim
      allocate(w(size(x)))
      call winsorize(x,fraction,w)
      value = finite_variance(w)
   end function winsor_variance

   pure real(dp) function winsor_sd(x, trim) result(value)
      real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted.
      real(dp), intent(in), optional :: trim !! Symmetric winsorization fraction; default `0.2`.
      if (present(trim)) then
         value = sqrt(max(0.0_dp,winsor_variance(x,trim)))
      else
         value = sqrt(max(0.0_dp,winsor_variance(x)))
      end if
   end function winsor_sd

   pure real(dp) function phi_coefficient(table2x2) result(phi)
      real(dp), intent(in) :: table2x2(2,2) !! Two-by-two contingency table of nonnegative counts or proportions.
      real(dp) :: total, r1, r2, c1, c2, denom
      total = sum(table2x2)
      if (total <= 0.0_dp) then
         phi = 0.0_dp
         return
      end if
      r1 = sum(table2x2(1,:))/total
      r2 = sum(table2x2(2,:))/total
      c1 = sum(table2x2(:,1))/total
      c2 = sum(table2x2(:,2))/total
      denom = sqrt(max(0.0_dp,r1*r2*c1*c2))
      if (denom > 0.0_dp) then
         phi = (table2x2(1,1)/total-c1*r1)/denom
      else
         phi = 0.0_dp
      end if
   end function phi_coefficient

   pure real(dp) function yule_q(table2x2) result(q)
      real(dp), intent(in) :: table2x2(2,2) !! Two-by-two contingency table.
      real(dp) :: ad, bc
      ad = table2x2(1,1)*table2x2(2,2)
      bc = table2x2(1,2)*table2x2(2,1)
      if (ad+bc > 0.0_dp) then
         q = (ad-bc)/(ad+bc)
      else
         q = 0.0_dp
      end if
   end function yule_q

   pure real(dp) function yule_y(table2x2) result(y)
      real(dp), intent(in) :: table2x2(2,2) !! Two-by-two contingency table.
      real(dp) :: ad, bc, sad, sbc
      ad = max(0.0_dp,table2x2(1,1)*table2x2(2,2))
      bc = max(0.0_dp,table2x2(1,2)*table2x2(2,1))
      sad = sqrt(ad)
      sbc = sqrt(bc)
      if (sad+sbc > 0.0_dp) then
         y = (sad-sbc)/(sad+sbc)
      else
         y = 0.0_dp
      end if
   end function yule_y

   pure subroutine yule_inverse(q, margins, table2x2)
      real(dp), intent(in) :: q !! Requested Yule-Q value in `[-1,1]`.
      real(dp), intent(in) :: margins(2) !! First-row and first-column proportions or counts on a common total of one.
      real(dp), intent(out) :: table2x2(2,2) !! Reconstructed table with the requested fixed margins.
      real(dp) :: lo, hi, mid, fmid, r1, c1
      integer :: iter
      r1 = margins(1)
      c1 = margins(2)
      lo = max(0.0_dp,r1+c1-1.0_dp)
      hi = min(r1,c1)
      do iter = 1, 100
         mid = 0.5_dp*(lo+hi)
         table2x2(1,1) = mid
         table2x2(1,2) = r1-mid
         table2x2(2,1) = c1-mid
         table2x2(2,2) = 1.0_dp-r1-c1+mid
         fmid = yule_q(table2x2)-q
         if (abs(fmid) < 1.0e-12_dp) exit
         if (fmid > 0.0_dp) then
            hi = mid
         else
            lo = mid
         end if
      end do
   end subroutine yule_inverse

   pure real(dp) function yule_to_phi(q, margins) result(phi)
      real(dp), intent(in) :: q !! Requested Yule-Q value in `[-1,1]`.
      real(dp), intent(in) :: margins(2) !! First-row and first-column proportions.
      real(dp) :: table2x2(2,2)
      call yule_inverse(q,margins,table2x2)
      phi = phi_coefficient(table2x2)
   end function yule_to_phi

   pure subroutine correlation_ci(rho, n, confidence, lower, upper)
      real(dp), intent(in) :: rho !! Observed correlation in `(-1,1)`.
      integer, intent(in) :: n !! Sample size; must exceed three.
      real(dp), intent(in), optional :: confidence !! Central confidence probability; default `0.95`.
      real(dp), intent(out) :: lower !! Lower confidence limit on the correlation scale.
      real(dp), intent(out) :: upper !! Upper confidence limit on the correlation scale.
      real(dp) :: conf, z, se, critical
      conf = 0.95_dp
      if (present(confidence)) conf = confidence
      z = fisher_z(rho)
      se = 1.0_dp/sqrt(real(max(1,n-3),dp))
      critical = normal_quantile(0.5_dp+0.5_dp*conf)
      lower = fisher_r(z-critical*se)
      upper = fisher_r(z+critical*se)
   end subroutine correlation_ci

   pure real(dp) function correlation_p_value(rho, n, two_tailed) result(p)
      real(dp), intent(in) :: rho !! Observed correlation in `(-1,1)`.
      integer, intent(in) :: n !! Sample size; must exceed two.
      logical, intent(in), optional :: two_tailed !! Return a two-sided p-value when true; default true.
      real(dp) :: t, tail
      logical :: two
      two = .true.
      if (present(two_tailed)) two = two_tailed
      t = r_to_t(rho,n)
      tail = 1.0_dp-student_t_cdf(abs(t),real(n-2,dp))
      p = merge(min(1.0_dp,2.0_dp*tail),tail,two)
   end function correlation_p_value

   pure real(dp) function auc_from_table(table2x2) result(auc)
      real(dp), intent(in) :: table2x2(2,2) !! Confusion table with true classes in rows and predictions in columns.
      real(dp) :: sensitivity, specificity, false_positive, dprime
      if (sum(table2x2(1,:)) <= 0.0_dp .or. sum(table2x2(2,:)) <= 0.0_dp) then
         auc = 0.5_dp
         return
      end if
      sensitivity = table2x2(1,1)/sum(table2x2(1,:))
      specificity = table2x2(2,2)/sum(table2x2(2,:))
      sensitivity = min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,sensitivity))
      false_positive = min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,1.0_dp-specificity))
      dprime = normal_quantile(sensitivity)-normal_quantile(false_positive)
      auc = normal_cdf(dprime/sqrt(2.0_dp))
   end function auc_from_table

   pure real(dp) function rmsea_value(chi_square, df, n_obs) result(value)
      real(dp), intent(in) :: chi_square !! Model chi-square statistic.
      real(dp), intent(in) :: df !! Positive model degrees of freedom.
      integer, intent(in) :: n_obs !! Sample size greater than one.
      value = sqrt(max(chi_square/(df*real(n_obs,dp))-1.0_dp/real(n_obs-1,dp),0.0_dp))
   end function rmsea_value

   pure real(dp) function manhattan_distance(x, y) result(value)
      real(dp), intent(in) :: x(:) !! First point.
      real(dp), intent(in) :: y(:) !! Second point with the same dimension as `x`.
      value = sum(abs(x-y))
   end function manhattan_distance

   pure real(dp) function minkowski_distance(x, y, power) result(value)
      real(dp), intent(in) :: x(:) !! First point.
      real(dp), intent(in) :: y(:) !! Second point with the same dimension as `x`.
      real(dp), intent(in), optional :: power !! Positive Minkowski order; default two.
      real(dp) :: p
      p = 2.0_dp
      if (present(power)) p = power
      value = sum(abs(x-y)**p)**(1.0_dp/p)
   end function minkowski_distance

   pure real(dp) function trace_value(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Matrix whose diagonal is summed.
      integer :: i
      value = 0.0_dp
      do i = 1, min(size(a,1),size(a,2))
         value = value+a(i,i)
      end do
   end function trace_value

   pure real(dp) function d_to_t(d, n, n1, n2) result(value)
      real(dp), intent(in) :: d !! Standardized mean difference.
      integer, intent(in), optional :: n !! Optional total sample size for equal-size two-group conversion.
      integer, intent(in), optional :: n1 !! Optional first-group sample size.
      integer, intent(in), optional :: n2 !! Optional second-group sample size.
      if (present(n1)) then
         if (present(n2)) then
            value = d / sqrt(1.0_dp / real(n1, dp) + 1.0_dp / real(n2, dp))
         else
            value = d * sqrt(real(n1, dp))
         end if
      else if (present(n)) then
         value = d * sqrt(real(n, dp)) / 2.0_dp
      else
         value = 0.0_dp
      end if
   end function d_to_t

   pure real(dp) function t_to_d(t, n, n1, n2) result(value)
      real(dp), intent(in) :: t !! Student t statistic.
      integer, intent(in), optional :: n !! Optional total sample size for equal-size two-group conversion.
      integer, intent(in), optional :: n1 !! Optional first-group sample size.
      integer, intent(in), optional :: n2 !! Optional second-group sample size.
      if (present(n1)) then
         if (present(n2)) then
            value = t * sqrt(1.0_dp / real(n1, dp) + 1.0_dp / real(n2, dp))
         else
            value = t / sqrt(real(n1, dp))
         end if
      else if (present(n)) then
         value = 2.0_dp * t / sqrt(real(n, dp))
      else
         value = 0.0_dp
      end if
   end function t_to_d

   pure real(dp) function means_to_d(mean1, mean2, sd1, sd2, n1, n2, pooled) result(value)
      real(dp), intent(in) :: mean1 !! Mean of the first group.
      real(dp), intent(in) :: mean2 !! Mean of the second group.
      real(dp), intent(in) :: sd1 !! Standard deviation of the first group.
      real(dp), intent(in) :: sd2 !! Standard deviation of the second group.
      integer, intent(in), optional :: n1 !! Optional first-group sample size for pooled variance.
      integer, intent(in), optional :: n2 !! Optional second-group sample size for pooled variance.
      logical, intent(in), optional :: pooled !! Use a pooled standard deviation when true; default true.
      real(dp) :: denominator
      logical :: use_pooled
      use_pooled = .true.
      if (present(pooled)) use_pooled = pooled
      if (use_pooled) then
         if (present(n1) .and. present(n2) .and. n1 + n2 > 2) then
            denominator = sqrt((real(n1 - 1, dp) * sd1 * sd1 + real(n2 - 1, dp) * sd2 * sd2) / &
               real(n1 + n2 - 2, dp))
         else
            denominator = sqrt(0.5_dp * (sd1 * sd1 + sd2 * sd2))
         end if
      else
         denominator = sd1
      end if
      if (denominator > tiny(1.0_dp)) then
         value = (mean1 - mean2) / denominator
      else
         value = 0.0_dp
      end if
   end function means_to_d

   pure subroutine means_to_t(mean1, mean2, sd1, sd2, n1, n2, t, df, d, pooled)
      real(dp), intent(in) :: mean1 !! Mean of the first group.
      real(dp), intent(in) :: mean2 !! Mean of the second group.
      real(dp), intent(in) :: sd1 !! Standard deviation of the first group.
      real(dp), intent(in) :: sd2 !! Standard deviation of the second group.
      integer, intent(in) :: n1 !! First-group sample size.
      integer, intent(in) :: n2 !! Second-group sample size.
      real(dp), intent(out) :: t !! Two-group t statistic.
      real(dp), intent(out) :: df !! Degrees of freedom, pooled or Welch-Satterthwaite.
      real(dp), intent(out) :: d !! Standardized effect corresponding to t.
      logical, intent(in), optional :: pooled !! Use equal-variance pooled standard error when true; default true.
      real(dp) :: se, vp, a, b
      logical :: use_pooled
      use_pooled = .true.
      if (present(pooled)) use_pooled = pooled
      if (use_pooled) then
         vp = (real(n1 - 1, dp) * sd1 * sd1 + real(n2 - 1, dp) * sd2 * sd2) / real(max(n1 + n2 - 2, 1), dp)
         se = sqrt(vp * (1.0_dp / real(n1, dp) + 1.0_dp / real(n2, dp)))
         df = real(n1 + n2 - 2, dp)
      else
         a = sd1 * sd1 / real(n1, dp)
         b = sd2 * sd2 / real(n2, dp)
         se = sqrt(a + b)
         df = (a + b) ** 2 / max(a * a / real(max(n1 - 1, 1), dp) + b * b / real(max(n2 - 1, 1), dp), tiny(1.0_dp))
      end if
      if (se > tiny(1.0_dp)) then
         t = (mean1 - mean2) / se
      else
         t = 0.0_dp
      end if
      d = t_to_d(t, n1=n1, n2=n2)
   end subroutine means_to_t

end module psych_basic
