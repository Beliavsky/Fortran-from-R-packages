program test_reference
   use matrixtests
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   implicit none
   real(dp) :: x(5, 12), y(5, 12), t(12), na, expected, actual, column_value, mu, conf, tolerance
   integer :: g(12), i, unit, ios, mode, row, old_mode, checks, failures
   character(64) :: name, old_name, field, token
   character(9) :: alt
   type(test_result), allocatable :: r(:), c(:)
   na = ieee_value(0.0_dp, ieee_quiet_nan)
   x(1, :) = [1.1_dp, 3.2_dp, 2.5_dp, 5.8_dp, 4.7_dp, 7.4_dp, 6.6_dp, 8.3_dp, 10.2_dp, 9.5_dp, 12.4_dp, 11.6_dp]
   x(2, :) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 6.0_dp]
   x(3, :) = [na, 2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp, na, 2.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, na]
   x(4, :) = 1.0e7_dp + [1.0_dp, 4.0_dp, 2.0_dp, 8.0_dp, 3.0_dp, 9.0_dp, 5.0_dp, 7.0_dp, 12.0_dp, 6.0_dp, 10.0_dp, 11.0_dp]
   x(5, :) = 2.0_dp
   y(1, :) = [1.5_dp, 2.8_dp, 3.5_dp, 4.1_dp, 5.9_dp, 6.3_dp, 8.2_dp, 7.7_dp, 9.8_dp, 12.2_dp, 10.3_dp, 11.1_dp]
   y(2, :) = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 5.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp]
   y(3, :) = [1.0_dp, 2.0_dp, na, 4.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, na, 7.0_dp, 9.0_dp, 10.0_dp, 11.0_dp]
   y(4, :) = 1.0e7_dp + [2.0_dp, 5.0_dp, 1.0_dp, 7.0_dp, 4.0_dp, 10.0_dp, 6.0_dp, 8.0_dp, 11.0_dp, 7.0_dp, 9.0_dp, 13.0_dp]
   y(5, :) = 3.0_dp
   g = [(1 + (i - 1) / 4, i = 1, 12)]
   t = [(real(2 * (i - 1), dp), i = 1, 12)]
   open(newunit = unit, file = 'test/reference.txt', status = 'old')
   old_name = ''
   old_mode = -1
   checks = 0
   failures = 0
   do
      read(unit, *, iostat = ios) name, mode, row, field, token
      if (ios /= 0) exit
      if (name /= old_name .or. mode /= old_mode) then
         alt = 'two.sided'
         if (mode == 2) alt = 'less'
         if (mode == 3) alt = 'greater'
         mu = 0.0_dp
         if (index(name, 'f_var') > 0) mu = 1.0_dp
         conf = 0.95_dp
         if (mode == 2 .or. mode == 3) then
            mu = 0.4_dp
            if (index(name, 'f_var') > 0) mu = 1.3_dp
            conf = 0.9_dp
         end if
         select case (name)
         case ('row_t_equalvar')
            r = row_t_equalvar(x, y, null = mu, alternative = alt, conf_level = conf)
            c = col_t_equalvar(transpose(x), transpose(y), null = mu, alternative = alt, conf_level = conf)
         case ('row_t_paired')
            r = row_t_paired(x, y, null = mu, alternative = alt, conf_level = conf)
            c = col_t_paired(transpose(x), transpose(y), null = mu, alternative = alt, conf_level = conf)
         case ('row_t_onesample')
            r = row_t_onesample(x, null = mu, alternative = alt, conf_level = conf)
            c = col_t_onesample(transpose(x), null = mu, alternative = alt, conf_level = conf)
         case ('row_t_welch')
            r = row_t_welch(x, y, null = mu, alternative = alt, conf_level = conf)
            c = col_t_welch(transpose(x), transpose(y), null = mu, alternative = alt, conf_level = conf)
         case ('row_cor_pearson')
            r = row_cor_pearson(x, y, alternative = alt, conf_level = conf)
            c = col_cor_pearson(transpose(x), transpose(y), alternative = alt, conf_level = conf)
         case ('row_f_var')
            r = row_f_var(x, y, null = mu, alternative = alt, conf_level = conf)
            c = col_f_var(transpose(x), transpose(y), null = mu, alternative = alt, conf_level = conf)
         case ('row_oneway_equalvar')
            r = row_oneway_equalvar(x, g)
            c = col_oneway_equalvar(transpose(x), g)
         case ('row_oneway_welch')
            r = row_oneway_welch(x, g)
            c = col_oneway_welch(transpose(x), g)
         case ('row_kruskalwallis')
            r = row_kruskalwallis(x, g)
            c = col_kruskalwallis(transpose(x), g)
         case ('row_waerden')
            r = row_waerden(x, g)
            c = col_waerden(transpose(x), g)
         case ('row_bartlett')
            r = row_bartlett(x, g)
            c = col_bartlett(transpose(x), g)
         case ('row_jarquebera')
            r = row_jarquebera(x)
            c = col_jarquebera(transpose(x))
         case ('row_flignerkilleen')
            r = row_flignerkilleen(x, g)
            c = col_flignerkilleen(transpose(x), g)
         case ('row_levene')
            r = row_levene(x, g)
            c = col_levene(transpose(x), g)
         case ('row_brownforsythe')
            r = row_brownforsythe(x, g)
            c = col_brownforsythe(transpose(x), g)
         case ('row_wilcoxon_onesample')
            if (mode == 4) then
            r = row_wilcoxon_onesample(x, null = mu, alternative = alt, exact = .false., correct = .false.)
            c = col_wilcoxon_onesample(transpose(x), null = mu, alternative = alt, exact = .false., correct = .false.)
            else
            r = row_wilcoxon_onesample(x, null = mu, alternative = alt)
            c = col_wilcoxon_onesample(transpose(x), null = mu, alternative = alt)
            end if
         case ('row_wilcoxon_paired')
            if (mode == 4) then
            r = row_wilcoxon_paired(x, y, null = mu, alternative = alt, exact = .false., correct = .false.)
            c = col_wilcoxon_paired(transpose(x), transpose(y), null = mu, alternative = alt, &
               exact = .false., correct = .false.)
            else
            r = row_wilcoxon_paired(x, y, null = mu, alternative = alt)
            c = col_wilcoxon_paired(transpose(x), transpose(y), null = mu, alternative = alt)
            end if
         case ('row_wilcoxon_twosample')
            if (mode == 4) then
            r = row_wilcoxon_twosample(x, y, null = mu, alternative = alt, exact = .false., correct = .false.)
            c = col_wilcoxon_twosample(transpose(x), transpose(y), null = mu, alternative = alt, &
               exact = .false., correct = .false.)
            else
            r = row_wilcoxon_twosample(x, y, null = mu, alternative = alt)
            c = col_wilcoxon_twosample(transpose(x), transpose(y), null = mu, alternative = alt)
            end if
         case ('row_cosinor')
            r = row_cosinor(x, t)
            c = col_cosinor(transpose(x), t)
         case ('row_andersondarling')
            r = row_andersondarling(x)
            c = col_andersondarling(transpose(x))
         case ('row_kolmogorovsmirnov_twosample')
            if (mode == 4) then
            r = row_kolmogorovsmirnov_twosample(x, y, alternative = alt, exact = .false.)
            c = col_kolmogorovsmirnov_twosample(transpose(x), transpose(y), alternative = alt, exact = .false.)
            else
            r = row_kolmogorovsmirnov_twosample(x, y, alternative = alt)
            c = col_kolmogorovsmirnov_twosample(transpose(x), transpose(y), alternative = alt)
            end if
         case default
            error stop 'unknown fixture family'
         end select
         old_name = name
         old_mode = mode
      end if
      if (token == 'NA' .or. token == 'NaN') then
         expected = na
      else if (token == 'Inf') then
         expected = ieee_value(0.0_dp, ieee_positive_inf)
      else if (token == '-Inf') then
         expected = ieee_value(0.0_dp, ieee_negative_inf)
      else
         read(token, *) expected
      end if
      actual = value(r(row), field)
      column_value = value(c(row), field)
      tolerance = 2.0e-8_dp
      if (name == 'row_kolmogorovsmirnov_twosample' .and. mode == 4) tolerance = 5.0e-7_dp
      if (.not. same(actual, expected, tolerance)) then
         failures = failures + 1
         print *, trim(name), mode, row, trim(field), ' got/expected ', actual, expected
      end if
      if (.not. same(actual, column_value, 1.0e-13_dp)) error stop "row/column disagreement"
      checks = checks + 1
   end do
   close(unit)
   print *, checks, "actual-R fields checked, failures:", failures
   if (failures /= 0) error stop "R reference mismatch"
contains
   pure logical function same(a, b, tolerance) result(equal)
      real(dp), intent(in) :: a !! Observed Fortran result.
      real(dp), intent(in) :: b !! Expected R result, including missing and infinite diagnostics.
      real(dp), intent(in) :: tolerance !! Relative/absolute comparison tolerance.
      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         equal = ieee_is_nan(a) .and. ieee_is_nan(b)
      else if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
         equal = a == b
      else
         equal = abs(a - b) <= tolerance * max(1.0_dp, abs(b))
      end if
   end function same
   pure real(dp) function value(r, field) result(v)
      type(test_result), intent(in) :: r !! One test result containing the union of numerical fields.
      character(*), intent(in) :: field !! Exact upstream numerical column name to compare.
      select case (field)
      case ('obs')
         v = real(r%obs, dp)
      case ('obs.x')
         v = real(r%obs_x, dp)
      case ('obs.y')
         v = real(r%obs_y, dp)
      case ('obs.tot')
         v = real(r%obs_tot, dp)
      case ('obs.paired')
         v = real(r%obs_paired, dp)
      case ('obs.groups')
         v = real(r%obs_groups, dp)
      case ('mean')
         v = r%mean
      case ('mean.x')
         v = r%mean_x
      case ('mean.y')
         v = r%mean_y
      case ('mean.diff')
         v = r%mean_diff
      case ('var')
         v = r%var
      case ('var.x')
         v = r%var_x
      case ('var.y')
         v = r%var_y
      case ('var.diff')
         v = r%var_diff
      case ('var.pooled')
         v = r%var_pooled
      case ('var.ratio')
         v = r%var_ratio
      case ('stderr')
         v = r%stderr
      case ('df')
         v = r%df
      case ('df.num')
         v = r%df_num
      case ('df.denom')
         v = r%df_denom
      case ('df.between')
         v = r%df_num
      case ('df.within')
         v = r%df_denom
      case ('df.model')
         v = r%df_num
      case ('df.residual')
         v = r%df_denom
      case ('statistic')
         v = r%statistic
      case ('pvalue')
         v = r%pvalue
      case ('conf.low')
         v = r%conf_low
      case ('conf.high')
         v = r%conf_high
      case ('cor')
         v = r%cor
      case ('skewness')
         v = r%skewness
      case ('kurtosis')
         v = r%kurtosis
      case ('sumsq.between')
         v = r%sumsq_between
      case ('sumsq.within')
         v = r%sumsq_within
      case ('meansq.between')
         v = r%meansq_between
      case ('meansq.within')
         v = r%meansq_within
      case ('mesor')
         v = r%mesor
      case ('amplitude')
         v = r%amplitude
      case ('acrophase')
         v = r%acrophase
      case ('rsquared')
         v = r%rsquared
      case ('period')
         v = r%period
      case ('exact')
         v = merge(1.0_dp, 0.0_dp, r%exact)
      case ('corrected')
         v = merge(1.0_dp, 0.0_dp, r%corrected)
      case default
         v = ieee_value(0.0_dp, ieee_quiet_nan)
      end select
   end function value
end program test_reference
