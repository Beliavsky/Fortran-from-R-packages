"""Regenerate thin, documented row/column adapters (no numerical algorithms)."""
from pathlib import Path
root = Path(__file__).resolve().parents[1]
names = ['t_equalvar', 't_paired', 't_onesample', 't_welch', 'cor_pearson', 'f_var',
 'oneway_equalvar', 'oneway_welch', 'kruskalwallis', 'waerden', 'bartlett', 'jarquebera',
 'flignerkilleen', 'levene', 'brownforsythe', 'wilcoxon_onesample', 'wilcoxon_paired',
 'wilcoxon_twosample', 'cosinor', 'andersondarling', 'kolmogorovsmirnov_twosample']
s = ['! SPDX-License-Identifier: GPL-2.0-only', '! Row/column numerical adapters for matrixTests 0.2.3.1.',
'module matrixtests', '   use matrixtests_types', '   use matrixtests_kernels',
'   use matrixtests_ks, only : ks_test', '   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan',
'   implicit none', '   private', '   public :: dp, test_result, mt_ok, mt_invalid, mt_insufficient, mt_degenerate, mt_limit, mt_nonfinite']
for name in names: s += [f'   public :: row_{name}, col_{name}']
s += ['contains']
groups = {'oneway_equalvar':1,'oneway_welch':2,'kruskalwallis':3,'waerden':4,'bartlett':5,
          'flignerkilleen':6,'levene':7,'brownforsythe':8}
for name in names:
    pair = name in ['t_equalvar','t_paired','t_welch','cor_pearson','f_var','wilcoxon_paired','wilcoxon_twosample','kolmogorovsmirnov_twosample']
    is_t = name.startswith('t_')
    is_w = name.startswith('wilcoxon_')
    has_null = is_t or is_w or name == 'f_var'
    has_alt = is_t or is_w or name in ['cor_pearson','f_var','kolmogorovsmirnov_twosample']
    has_conf = is_t or name in ['cor_pearson','f_var']
    has_exact = is_w or name == 'kolmogorovsmirnov_twosample'
    args = ['x'] + (['y'] if pair else []) + (['g'] if name in groups else []) + (['t','period'] if name == 'cosinor' else [])
    if has_null: args += ['null']
    if has_alt: args += ['alternative']
    if has_conf: args += ['conf_level']
    if has_exact: args += ['exact']
    if is_w: args += ['correct']
    pure = '' if name == 'jarquebera' else 'pure '
    for axis in ['row','col']:
        prefix=f'{axis}_{name}'
        s += ['',f'   {pure}function {prefix}({", ".join(args)}) result(r)',
              f'      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by {axis}; IEEE NaNs are missing.']
        if pair: s += [f'      real(dp), intent(in) :: y(:, :) !! Second matrix; one {axis} broadcasts over x, otherwise test counts must agree.']
        if name in groups: s += ['      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.']
        if name == 'cosinor': s += ['      real(dp), intent(in) :: t(:) !! Times along the observation dimension; NaN times are omitted.',
             '      real(dp), intent(in), optional :: period !! Positive period in time units, default 24.']
        if has_null: s += ['      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).']
        if has_alt: s += ['      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.']
        if has_conf: s += ['      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.']
        if has_exact: s += ['      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.']
        if is_w: s += ['      logical, intent(in), optional :: correct !! Apply continuity correction to approximate Wilcoxon p-values, default true.']
        s += ['      type(test_result), allocatable :: r(:)']
        if axis == 'col':
            actual=['transpose(x)']+(['transpose(y)'] if pair else [])+(['g'] if name in groups else [])+(['t','period'] if name=='cosinor' else [])
            actual += [a for a in args if a in ['null','alternative','conf_level','exact','correct']]
            s += [f'      r = row_{name}({", ".join(actual)})', f'   end function {prefix}']
            continue
        s += ['      integer :: i']
        if pair: s += ['      integer :: j']
        if name=='kolmogorovsmirnov_twosample': s += ['      real(dp), allocatable :: a(:), b(:)', '      integer :: status']
        if has_null: s += ['      real(dp) :: mu']
        if has_alt: s += ['      character(9) :: alt']
        if has_conf: s += ['      real(dp) :: conf']
        if name=='cosinor': s += ['      real(dp) :: cycle']
        s += ['      allocate(r(size(x, 1)))', '      r = empty_result()']
        if pair:
            s += ['      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return']
            if name in ['t_paired','cor_pearson','wilcoxon_paired']: s += ['      if (size(x, 2) /= size(y, 2)) return']
        if has_null: s += [f'      mu = {"1.0_dp" if name=="f_var" else "0.0_dp"}', '      if (present(null)) mu = null']
        if has_alt: s += ["      alt = 'two.sided'", '      if (present(alternative)) then',
                          '         if (.not. valid_alternative(alternative)) return', '         alt = alternative', '      end if']
        if has_conf: s += ['      conf = 0.95_dp','      if (present(conf_level)) conf = conf_level']
        if name=='cosinor': s += ['      cycle = 24.0_dp', '      if (present(period)) cycle = period']
        s += ['      do i = 1, size(x, 1)']
        if pair: s += ['         j = i', '         if (size(y, 1) == 1) j = 1']
        b='y(j, :)' if pair else '[real(dp) ::]'
        if is_t: s += [f'         r(i) = t_kernel(x(i, :), {b}, {dict(t_onesample=1,t_equalvar=2,t_welch=3,t_paired=4)[name]}, mu, alt, conf)']
        elif name=='f_var': s += ['         r(i) = f_kernel(x(i, :), y(j, :), mu, alt, conf)']
        elif name=='cor_pearson': s += ['         r(i) = correlation_kernel(x(i, :), y(j, :), alt, conf)']
        elif name in groups: s += [f'         r(i) = group_kernel(x(i, :), g, {groups[name]})']
        elif is_w: s += [f'         r(i) = wilcoxon_kernel(x(i, :), {b}, {dict(wilcoxon_onesample=1,wilcoxon_twosample=2,wilcoxon_paired=3)[name]}, &',
                         '            mu, alt, exact, correct)']
        elif name=='andersondarling': s += ['         r(i) = ad_kernel(x(i, :))']
        elif name=='jarquebera': s += ['         r(i) = jarque_kernel(x(i, :))']
        elif name=='cosinor': s += ['         r(i) = cosinor_kernel(x(i, :), t, cycle)']
        else:
            s += ['         a = pack(x(i, :), .not. ieee_is_nan(x(i, :)))', '         b = pack(y(j, :), .not. ieee_is_nan(y(j, :)))',
                  '         r(i)%obs_x = size(a)', '         r(i)%obs_y = size(b)', '         r(i)%obs_tot = size(a) + size(b)',
                  '         r(i)%alternative = alt',
                  '         call ks_test(a, b, r(i)%statistic, r(i)%pvalue, status, alt, exact, r(i)%exact)',
                  '         r(i)%status = mt_ok', '         if (status /= 0) r(i)%status = mt_invalid',
                  '         if (status == 3) r(i)%status = mt_limit']
        s += ['      end do', f'   end function {prefix}']
s += ['end module matrixtests']
(root/'src/matrixtests.f90').write_text('\n'.join(s)+'\n')
