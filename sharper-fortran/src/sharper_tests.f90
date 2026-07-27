! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on SharpeR, copyright 2012-2025 Steven E. Pav.
module sharper_tests
   use sharper_kinds, only: dp, huge_dp
   use sharper_types, only: sr_result, sropt_result, test_result, moment_vcov_result
   use sharper_math, only: normal_cdf, normal_quantile, student_t_cdf, student_t_quantile
   use sharper_math, only: chisq_cdf, chisq_quantile, binomial_coefficient
   use sharper_distributions, only: prt, psropt, qsropt
   use sharper_estimation, only: fit_sr, sr_standard_error, sr_vcov
   use sharper_linalg, only: solve_linear, quadratic_form
   implicit none
   private

   public :: sr_test, paired_sr_test, unpaired_sr_test
   public :: sr_equality_test, sropt_test
   public :: sr_max_test, sr_conditional_test
   public :: power_sr_test, required_n_sr_test
   public :: power_sropt_test, required_df2_sropt_test

contains

   pure function lowercase(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, code
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            out(i:i) = achar(code+32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lowercase

   pure elemental function two_sided_probability(p) result(v)
      real(dp), intent(in) :: p
      real(dp) :: v
      v = min(1.0_dp,2.0_dp*min(p,1.0_dp-p))
   end function two_sided_probability

   function sr_test(z, zeta, alternative, method, conf_level) result(out)
      type(sr_result), intent(in) :: z
      real(dp), intent(in), optional :: zeta, conf_level
      character(len=*), intent(in), optional :: alternative, method
      type(test_result) :: out
      real(dp) :: null_value, level, cdf, tstat, midpoint, se, alpha, q
      real(dp), allocatable :: se_values(:)
      character(len=24) :: alt, local_method
      null_value = 0.0_dp
      level = 0.95_dp
      alt = 'two.sided'
      local_method = 'exact'
      if (present(zeta)) null_value = zeta
      if (present(conf_level)) level = conf_level
      if (present(alternative)) alt = lowercase(trim(alternative))
      if (present(method)) local_method = lowercase(trim(method))
      tstat = z%value(1)/(z%rescal(1)*sqrt(z%ope))
      out%statistic = tstat
      out%estimate = z%value(1)
      out%null_value = null_value
      out%df1 = z%df(1)
      out%alternative = alt
      out%method = 'one-sample Sharpe ratio test'
      alpha = 1.0_dp-level
      if (trim(local_method) == 'exact') then
         cdf = prt(z%value(1),real(z%df(1),dp),z%rescal(1)*sqrt(z%ope),null_value)
         select case (trim(alt))
         case ('less')
            out%p_value = cdf
            out%conf_low = -huge_dp
            out%conf_high = qlambdap_wrapper(level,z,tstat)
         case ('greater')
            out%p_value = 1.0_dp-cdf
            out%conf_low = qlambdap_wrapper(alpha,z,tstat)
            out%conf_high = huge_dp
         case default
            out%p_value = two_sided_probability(cdf)
            out%conf_low = qlambdap_wrapper(0.5_dp*alpha,z,tstat)
            out%conf_high = qlambdap_wrapper(1.0_dp-0.5_dp*alpha,z,tstat)
         end select
      else
         midpoint = z%value(1)
         if (trim(local_method) == 'z') then
            midpoint = midpoint*(1.0_dp-1.0_dp/(4.0_dp*real(z%df(1),dp)))
         end if
         se_values = sr_standard_error(z,'t')
         se = se_values(1)
         q = (midpoint-null_value)/se
         select case (trim(alt))
         case ('less')
            out%p_value = normal_cdf(q)
            out%conf_low = -huge_dp
            out%conf_high = midpoint+normal_quantile(level)*se
         case ('greater')
            out%p_value = 1.0_dp-normal_cdf(q)
            out%conf_low = midpoint+normal_quantile(alpha)*se
            out%conf_high = huge_dp
         case default
            out%p_value = two_sided_probability(normal_cdf(q))
            out%conf_low = midpoint+normal_quantile(0.5_dp*alpha)*se
            out%conf_high = midpoint+normal_quantile(1.0_dp-0.5_dp*alpha)*se
         end select
      end if
   end function sr_test

   function qlambdap_wrapper(p, z, tstat) result(value)
      use sharper_distributions, only: qlambdap
      real(dp), intent(in) :: p, tstat
      type(sr_result), intent(in) :: z
      real(dp) :: value
      value = qlambdap(p,real(z%df(1),dp),tstat)*z%rescal(1)*sqrt(z%ope)
   end function qlambdap_wrapper

   function paired_sr_test(x, y, ope, alternative) result(out)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(in), optional :: ope
      character(len=*), intent(in), optional :: alternative
      type(test_result) :: out
      real(dp), allocatable :: xy(:, :), contrast(:, :)
      real(dp) :: local_ope
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      if (size(x) /= size(y)) then
         out%status = 1
         return
      end if
      allocate(xy(size(x),2),contrast(1,2))
      xy(:,1) = x
      xy(:,2) = y
      contrast(1,:) = [1.0_dp,-1.0_dp]
      out = sr_equality_test(xy,'t',alternative,contrast,local_ope)
      out%method = 'paired Sharpe ratio test'
   end function paired_sr_test

   function unpaired_sr_test(srs, contrasts, null_value, alternative, ope, conf_level) result(out)
      type(sr_result), intent(in) :: srs(:)
      real(dp), intent(in), optional :: contrasts(:), null_value, ope, conf_level
      character(len=*), intent(in), optional :: alternative
      type(test_result) :: out
      real(dp), allocatable :: c(:), se(:), one_se(:)
      real(dp) :: nullv, local_ope, level, estimate_native, variance_native, statistic
      real(dp) :: qlo, qhi
      integer :: i, min_df
      character(len=24) :: alt
      allocate(c(size(srs)),se(size(srs)))
      do i = 1, size(srs)
         if (present(contrasts)) then
            c(i) = contrasts(i)
         else
            c(i) = merge(1.0_dp,-1.0_dp,mod(i,2) == 1)
         end if
         one_se = sr_standard_error(srs(i),'t')
         se(i) = one_se(1)/sqrt(srs(i)%ope)
      end do
      nullv = 0.0_dp
      local_ope = srs(1)%ope
      level = 0.95_dp
      alt = 'two.sided'
      if (present(null_value)) nullv = null_value
      if (present(ope)) local_ope = ope
      if (present(conf_level)) level = conf_level
      if (present(alternative)) alt = lowercase(trim(alternative))
      estimate_native = 0.0_dp
      variance_native = 0.0_dp
      min_df = huge(1)
      do i = 1, size(srs)
         estimate_native = estimate_native+c(i)*srs(i)%value(1)/sqrt(srs(i)%ope)
         variance_native = variance_native+(c(i)*se(i))**2
         min_df = min(min_df,srs(i)%df(1))
      end do
      statistic = (estimate_native-nullv/sqrt(local_ope))/sqrt(variance_native)
      out%statistic = statistic
      out%estimate = estimate_native*sqrt(local_ope)
      out%null_value = nullv
      out%df1 = min_df
      out%alternative = alt
      out%method = 'unpaired k-sample Sharpe ratio test'
      select case (trim(alt))
      case ('less')
         out%p_value = student_t_cdf(statistic,real(min_df,dp))
         qlo = -huge_dp
         qhi = student_t_quantile(level,real(min_df,dp))
      case ('greater')
         out%p_value = 1.0_dp-student_t_cdf(statistic,real(min_df,dp))
         qlo = student_t_quantile(1.0_dp-level,real(min_df,dp))
         qhi = huge_dp
      case default
         out%p_value = two_sided_probability(student_t_cdf(statistic,real(min_df,dp)))
         qlo = student_t_quantile(0.5_dp*(1.0_dp-level),real(min_df,dp))
         qhi = student_t_quantile(0.5_dp*(1.0_dp+level),real(min_df,dp))
      end select
      if (qlo <= -0.5_dp*huge_dp) then
         out%conf_low = -huge_dp
      else
         out%conf_low = (nullv/sqrt(local_ope)+qlo*sqrt(variance_native))*sqrt(local_ope)
      end if
      if (qhi >= 0.5_dp*huge_dp) then
         out%conf_high = huge_dp
      else
         out%conf_high = (nullv/sqrt(local_ope)+qhi*sqrt(variance_native))*sqrt(local_ope)
      end if
   end function unpaired_sr_test

   function sr_equality_test(x, test_type, alternative, contrasts, ope) result(out)
      real(dp), intent(in) :: x(:, :)
      character(len=*), intent(in), optional :: test_type, alternative
      real(dp), intent(in), optional :: contrasts(:, :), ope
      type(test_result) :: out
      real(dp), allocatable :: c(:, :), esr(:), coc(:, :), sol(:)
      type(moment_vcov_result) :: moments
      real(dp) :: local_ope, statistic
      integer :: p, k, i, status
      character(len=24) :: kind, alt
      p = size(x,2)
      local_ope = 1.0_dp
      kind = 'chisq'
      alt = 'two.sided'
      if (present(ope)) local_ope = ope
      if (present(test_type)) kind = lowercase(trim(test_type))
      if (present(alternative)) alt = lowercase(trim(alternative))
      if (present(contrasts)) then
         allocate(c(size(contrasts,1),size(contrasts,2)))
         c = contrasts
      else
         allocate(c(p-1,p))
         c = 0.0_dp
         do i = 1, p-1
            c(i,i) = 1.0_dp
            c(i,i+1) = -1.0_dp
         end do
      end if
      k = size(c,1)
      moments = sr_vcov(x,local_ope)
      allocate(esr(k),coc(k,k),sol(k))
      esr = matmul(c,moments%mean)
      coc = matmul(c,matmul(moments%covariance,transpose(c)))
      out%alternative = alt
      out%parameter = k
      out%df1 = k
      out%df2 = size(x,1)-k
      if (trim(kind) == 't') then
         if (k /= 1) then
            out%status = 1
            return
         end if
         statistic = esr(1)/sqrt(coc(1,1))
         out%statistic = statistic
         select case (trim(alt))
         case ('less')
            out%p_value = student_t_cdf(statistic,real(size(x,1)-1,dp))
         case ('greater')
            out%p_value = 1.0_dp-student_t_cdf(statistic,real(size(x,1)-1,dp))
         case default
            out%p_value = two_sided_probability(student_t_cdf(statistic,real(size(x,1)-1,dp)))
         end select
      else
         call solve_linear(coc,esr,sol,status)
         if (status /= 0) then
            out%status = status
            return
         end if
         statistic = dot_product(esr,sol)
         out%statistic = statistic
         if (trim(kind) == 'f') then
            statistic = real(size(x,1)-k,dp)*statistic / &
                        (real(size(x,1)-1,dp)*real(k,dp))
            out%p_value = 1.0_dp-f_cdf_wrapper(statistic,real(k,dp),real(size(x,1)-k,dp))
         else
            out%p_value = 1.0_dp-chisq_cdf(statistic,real(k,dp))
         end if
         out%alternative = 'two.sided'
      end if
      out%method = 'test for equality of Sharpe ratios'
   end function sr_equality_test

   function f_cdf_wrapper(x, df1, df2) result(v)
      use sharper_math, only: f_cdf
      real(dp), intent(in) :: x, df1, df2
      real(dp) :: v
      v = f_cdf(x,df1,df2)
   end function f_cdf_wrapper

   function sropt_test(z, zeta_s, alternative) result(out)
      type(sropt_result), intent(in) :: z
      real(dp), intent(in), optional :: zeta_s
      character(len=*), intent(in), optional :: alternative
      type(test_result) :: out
      real(dp) :: nullv, cdf
      character(len=24) :: alt
      nullv = 0.0_dp
      alt = 'greater'
      if (present(zeta_s)) nullv = zeta_s
      if (present(alternative)) alt = lowercase(trim(alternative))
      cdf = psropt(z%value,real(z%df1,dp),real(z%df2,dp),nullv,z%ope,z%drag)
      select case (trim(alt))
      case ('less')
         out%p_value = cdf
      case ('two.sided')
         out%p_value = two_sided_probability(cdf)
      case default
         out%p_value = 1.0_dp-cdf
      end select
      out%statistic = z%t2
      out%estimate = z%value
      out%null_value = nullv
      out%df1 = z%df1
      out%df2 = z%df2
      out%alternative = alt
      out%method = 'one-sample optimal Sharpe ratio test'
   end function sropt_test

   function sr_max_test(srs, df, ope, kappa, rho, zeta0, conf_level, method, loglog) result(out)
      real(dp), intent(in) :: srs(:)
      integer, intent(in) :: df
      real(dp), intent(in), optional :: ope, kappa, rho, zeta0, conf_level
      character(len=*), intent(in), optional :: method
      logical, intent(in), optional :: loglog
      type(test_result) :: out
      real(dp) :: local_ope, local_kappa, local_rho, z0, level, alpha
      real(dp) :: native(size(srs)), root_xi(size(srs)), meansr, a0, a2, b0, b2, cval
      real(dp) :: statistic, g0, g1, chi_value, p_value, weight
      integer :: nday, nstrat, nabove, j
      logical :: use_loglog
      character(len=24) :: local_method
      local_ope = 1.0_dp
      local_kappa = 1.0_dp
      local_rho = 0.0_dp
      z0 = 0.0_dp
      level = 0.95_dp
      local_method = 'bonferroni'
      use_loglog = .true.
      if (present(ope)) local_ope = ope
      if (present(kappa)) local_kappa = kappa
      if (present(rho)) local_rho = rho
      if (present(zeta0)) z0 = zeta0/sqrt(local_ope)
      if (present(conf_level)) level = conf_level
      if (present(method)) local_method = lowercase(trim(method))
      if (present(loglog)) use_loglog = loglog
      nday = df+1
      nstrat = size(srs)
      native = srs/sqrt(local_ope)
      meansr = sum(native)/real(nstrat,dp)
      a0 = (1.0_dp-local_rho)+0.5_dp*local_kappa*z0*z0*(1.0_dp-local_rho**2)
      a2 = local_rho+0.25_dp*(local_kappa-1.0_dp)*z0*z0 + &
           0.5_dp*local_kappa*z0*z0*local_rho**2
      b0 = 1.0_dp/sqrt(a0)
      b2 = -1.0_dp/(real(nstrat,dp)*sqrt(a0)) + &
           1.0_dp/(real(nstrat,dp)*sqrt(a0+real(nstrat,dp)*a2))
      cval = b0+real(nstrat,dp)*b2
      root_xi = sqrt(real(nday,dp))*(b0*(native-meansr)+cval*(meansr-z0))
      nabove = nstrat
      if (use_loglog .and. nday > 2) then
         nabove = count(root_xi > sqrt(real(nday,dp))*cval*z0 - &
                        sqrt(2.0_dp*log(log(real(nday,dp)))))
      end if
      out%parameter = nabove
      out%alternative = 'greater'
      out%null_value = z0*sqrt(local_ope)
      out%method = local_method
      if (nabove < 1) then
         out%statistic = huge_dp
         out%p_value = 1.0_dp
         out%conf_low = -huge_dp
         out%conf_high = huge_dp
         return
      end if
      alpha = 1.0_dp-level
      select case (trim(local_method))
      case ('chi-bar-square','chi_bar_square','chibarsquare')
         statistic = sum(max(root_xi,0.0_dp)**2)
         p_value = 0.0_dp
         do j = 0, nabove
            weight = binomial_coefficient(nabove,j)/2.0_dp**nabove
            if (j == 0) then
               if (statistic <= tiny(1.0_dp)) p_value = p_value+weight
            else
               p_value = p_value+weight*(1.0_dp-chisq_cdf(statistic,real(j,dp)))
            end if
         end do
         out%statistic = statistic
         out%p_value = min(1.0_dp,p_value)
         out%conf_low = maxval(native)*sqrt(local_ope) - &
                        normal_quantile(1.0_dp-alpha)*sqrt(local_ope/real(nday,dp))
      case ('follman')
         g0 = real(nday,dp)*b0*b0*sum((native-meansr)**2)
         g1 = real(nday,dp)*real(nstrat,dp)*(cval*(meansr-z0))**2
         statistic = g0+g1
         chi_value = 1.0_dp-chisq_cdf(statistic,real(nabove,dp))
         if (meansr > z0) then
            p_value = 0.5_dp*chi_value
         else
            p_value = 1.0_dp-0.5_dp*chi_value
         end if
         out%statistic = statistic
         out%p_value = p_value
         out%conf_low = (meansr-sqrt(max(0.0_dp, &
            chisq_quantile(1.0_dp-2.0_dp*alpha,real(nabove,dp))/real(nday,dp)-g0) / &
            real(nstrat,dp))/cval)*sqrt(local_ope)
      case default
         statistic = maxval(root_xi)
         out%statistic = statistic
         out%p_value = min(1.0_dp,real(nabove,dp)*(1.0_dp-normal_cdf(statistic)))
         out%conf_low = (meansr+(b0*(maxval(native)-meansr) - &
            normal_quantile(1.0_dp-alpha/real(nabove,dp))/sqrt(real(nday,dp)))/cval)*sqrt(local_ope)
      end select
      out%conf_high = huge_dp
   end function sr_max_test

   pure function truncated_normal_cdf(x, mu, sigma, lower, upper) result(v)
      real(dp), intent(in) :: x, mu, sigma, lower, upper
      real(dp) :: v, flo, fhi, fx
      flo = merge(0.0_dp,normal_cdf((lower-mu)/sigma),lower <= -0.5_dp*huge_dp)
      fhi = merge(1.0_dp,normal_cdf((upper-mu)/sigma),upper >= 0.5_dp*huge_dp)
      fx = normal_cdf((x-mu)/sigma)
      v = (fx-flo)/max(fhi-flo,tiny(1.0_dp))
      v = min(1.0_dp,max(0.0_dp,v))
   end function truncated_normal_cdf

   function sr_conditional_test(srs, df, covariance, zeta0, ope, alternative, conf_level) result(out)
      real(dp), intent(in) :: srs(:), covariance(:, :)
      integer, intent(in) :: df
      real(dp), intent(in), optional :: zeta0, ope, conf_level
      character(len=*), intent(in), optional :: alternative
      type(test_result) :: out
      real(dp) :: local_ope, nullv, level, native(size(srs)), eta(size(srs))
      real(dp) :: cvec(size(srs)), zvec(size(srs)), avec(size(srs))
      real(dp) :: variance, sigma, observed, alpha_coeff, rhs, lower, upper, cdf
      real(dp) :: alpha, lo_mu, hi_mu, mid_mu, target
      integer :: kmax, j, iter
      character(len=24) :: alt
      local_ope = 1.0_dp
      nullv = 0.0_dp
      level = 0.95_dp
      alt = 'two.sided'
      if (present(ope)) local_ope = ope
      if (present(zeta0)) nullv = zeta0/sqrt(local_ope)
      if (present(conf_level)) level = conf_level
      if (present(alternative)) alt = lowercase(trim(alternative))
      native = srs/sqrt(local_ope)
      kmax = maxloc(native,dim=1)
      eta = 0.0_dp
      eta(kmax) = 1.0_dp
      variance = dot_product(eta,matmul(covariance,eta))
      if (variance <= tiny(1.0_dp)) then
         out%status = 1
         return
      end if
      sigma = sqrt(variance)
      cvec = matmul(covariance,eta)/variance
      observed = native(kmax)
      zvec = native-cvec*observed
      lower = -huge_dp
      upper = huge_dp
      do j = 1, size(srs)
         if (j == kmax) cycle
         avec = 0.0_dp
         avec(j) = 1.0_dp
         avec(kmax) = -1.0_dp
         alpha_coeff = dot_product(avec,cvec)
         rhs = -dot_product(avec,zvec)
         if (alpha_coeff > tiny(1.0_dp)) then
            upper = min(upper,rhs/alpha_coeff)
         else if (alpha_coeff < -tiny(1.0_dp)) then
            lower = max(lower,rhs/alpha_coeff)
         end if
      end do
      cdf = truncated_normal_cdf(observed,nullv,sigma,lower,upper)
      select case (trim(alt))
      case ('less')
         out%p_value = cdf
      case ('greater')
         out%p_value = 1.0_dp-cdf
      case default
         out%p_value = two_sided_probability(cdf)
      end select
      alpha = 1.0_dp-level
      lo_mu = observed-20.0_dp*sigma-10.0_dp
      hi_mu = observed+20.0_dp*sigma+10.0_dp
      target = 1.0_dp-0.5_dp*alpha
      do iter = 1, 120
         mid_mu = 0.5_dp*(lo_mu+hi_mu)
         if (truncated_normal_cdf(observed,mid_mu,sigma,lower,upper) > target) then
            lo_mu = mid_mu
         else
            hi_mu = mid_mu
         end if
      end do
      out%conf_low = 0.5_dp*(lo_mu+hi_mu)*sqrt(local_ope)
      lo_mu = observed-20.0_dp*sigma-10.0_dp
      hi_mu = observed+20.0_dp*sigma+10.0_dp
      target = 0.5_dp*alpha
      do iter = 1, 120
         mid_mu = 0.5_dp*(lo_mu+hi_mu)
         if (truncated_normal_cdf(observed,mid_mu,sigma,lower,upper) > target) then
            lo_mu = mid_mu
         else
            hi_mu = mid_mu
         end if
      end do
      out%conf_high = 0.5_dp*(lo_mu+hi_mu)*sqrt(local_ope)
      out%statistic = observed
      out%estimate = srs(kmax)
      out%null_value = nullv*sqrt(local_ope)
      out%parameter = df
      out%alternative = alt
      out%method = 'polyhedral conditional max test'
   end function sr_conditional_test

   function power_sr_test(n, zeta, sig_level, two_sided, ope) result(value)
      integer, intent(in) :: n
      real(dp), intent(in) :: zeta, sig_level
      logical, intent(in), optional :: two_sided
      real(dp), intent(in), optional :: ope
      real(dp) :: value, native_zeta, delta, critical
      logical :: use_two
      use_two = .false.
      if (present(two_sided)) use_two = two_sided
      native_zeta = zeta
      if (present(ope)) native_zeta = zeta/sqrt(ope)
      delta = native_zeta*sqrt(real(n,dp))
      if (use_two) then
         critical = student_t_quantile(1.0_dp-0.5_dp*sig_level,real(n-1,dp))
         value = nct_tail_power(critical,real(n-1,dp),delta,.true.)
      else
         critical = student_t_quantile(1.0_dp-sig_level,real(n-1,dp))
         value = nct_tail_power(critical,real(n-1,dp),delta,.false.)
      end if
   end function power_sr_test

   function nct_tail_power(critical, df, delta, two_sided) result(value)
      use sharper_math, only: nct_cdf
      real(dp), intent(in) :: critical, df, delta
      logical, intent(in) :: two_sided
      real(dp) :: value
      value = 1.0_dp-nct_cdf(critical,df,delta)
      if (two_sided) value = value+nct_cdf(-critical,df,delta)
      value = min(1.0_dp,max(0.0_dp,value))
   end function nct_tail_power

   function required_n_sr_test(zeta, target_power, sig_level, two_sided, ope) result(n)
      real(dp), intent(in) :: zeta, target_power, sig_level
      logical, intent(in), optional :: two_sided
      real(dp), intent(in), optional :: ope
      integer :: n, lo, hi, mid
      lo = 2
      hi = 4
      do while (power_sr_test(hi,zeta,sig_level,two_sided,ope) < target_power .and. hi < 100000000)
         hi = 2*hi
      end do
      do while (lo+1 < hi)
         mid = (lo+hi)/2
         if (power_sr_test(mid,zeta,sig_level,two_sided,ope) < target_power) then
            lo = mid
         else
            hi = mid
         end if
      end do
      n = hi
   end function required_n_sr_test

   function power_sropt_test(df1, df2, zeta_s, sig_level, ope) result(value)
      integer, intent(in) :: df1, df2
      real(dp), intent(in) :: zeta_s, sig_level
      real(dp), intent(in), optional :: ope
      real(dp) :: value, local_ope, critical
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      critical = qsropt(sig_level,real(df1,dp),real(df2,dp),ope=local_ope,lower_tail=.false.)
      value = psropt(critical,real(df1,dp),real(df2,dp),zeta_s,local_ope,lower_tail=.false.)
   end function power_sropt_test

   function required_df2_sropt_test(df1, zeta_s, target_power, sig_level, ope) result(df2)
      integer, intent(in) :: df1
      real(dp), intent(in) :: zeta_s, target_power, sig_level
      real(dp), intent(in), optional :: ope
      integer :: df2, lo, hi, mid
      lo = df1+1
      hi = max(df1+2,2*df1)
      do while (power_sropt_test(df1,hi,zeta_s,sig_level,ope) < target_power .and. hi < 100000000)
         hi = 2*hi
      end do
      do while (lo+1 < hi)
         mid = (lo+hi)/2
         if (power_sropt_test(df1,mid,zeta_s,sig_level,ope) < target_power) then
            lo = mid
         else
            hi = mid
         end if
      end do
      df2 = hi
   end function required_df2_sropt_test

end module sharper_tests
