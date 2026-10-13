"""Generate exhaustive row/column comparisons against committed actual-R records."""
from pathlib import Path
import re
root=Path(__file__).resolve().parents[1]
exports=re.findall(r'export\((row_[^)]+)\)',(root/'upstream/NAMESPACE').read_text())
s=['program test_reference', '   use matrixtests', '   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, &',
'      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf', '   implicit none',
'   real(dp) :: x(5, 12), y(5, 12), t(12), na, expected, actual, column_value, mu, conf, tolerance',
'   integer :: g(12), i, unit, ios, mode, row, old_mode, checks, failures',
'   character(64) :: name, old_name, field, token', "   character(9) :: alt", '   type(test_result), allocatable :: r(:), c(:)',
'   na = ieee_value(0.0_dp, ieee_quiet_nan)',
'   x(1, :) = [1.1_dp, 3.2_dp, 2.5_dp, 5.8_dp, 4.7_dp, 7.4_dp, 6.6_dp, 8.3_dp, 10.2_dp, 9.5_dp, 12.4_dp, 11.6_dp]',
'   x(2, :) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 6.0_dp]',
'   x(3, :) = [na, 2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp, na, 2.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, na]',
'   x(4, :) = 1.0e7_dp + [1.0_dp, 4.0_dp, 2.0_dp, 8.0_dp, 3.0_dp, 9.0_dp, 5.0_dp, 7.0_dp, 12.0_dp, 6.0_dp, 10.0_dp, 11.0_dp]',
'   x(5, :) = 2.0_dp',
'   y(1, :) = [1.5_dp, 2.8_dp, 3.5_dp, 4.1_dp, 5.9_dp, 6.3_dp, 8.2_dp, 7.7_dp, 9.8_dp, 12.2_dp, 10.3_dp, 11.1_dp]',
'   y(2, :) = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 5.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp]',
'   y(3, :) = [1.0_dp, 2.0_dp, na, 4.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, na, 7.0_dp, 9.0_dp, 10.0_dp, 11.0_dp]',
'   y(4, :) = 1.0e7_dp + [2.0_dp, 5.0_dp, 1.0_dp, 7.0_dp, 4.0_dp, 10.0_dp, 6.0_dp, 8.0_dp, 11.0_dp, 7.0_dp, 9.0_dp, 13.0_dp]',
'   y(5, :) = 3.0_dp', '   g = [(1 + (i - 1) / 4, i = 1, 12)]', '   t = [(real(2 * (i - 1), dp), i = 1, 12)]',
"   open(newunit = unit, file = 'test/reference.txt', status = 'old')", "   old_name = ''", '   old_mode = -1',
'   checks = 0', '   failures = 0', '   do',
'      read(unit, *, iostat = ios) name, mode, row, field, token', '      if (ios /= 0) exit',
'      if (name /= old_name .or. mode /= old_mode) then',
"         alt = 'two.sided'", "         if (mode == 2) alt = 'less'", "         if (mode == 3) alt = 'greater'",
'         mu = 0.0_dp', "         if (index(name, 'f_var') > 0) mu = 1.0_dp", '         conf = 0.95_dp',
'         if (mode == 2 .or. mode == 3) then', '            mu = 0.4_dp', "            if (index(name, 'f_var') > 0) mu = 1.3_dp", '            conf = 0.9_dp', '         end if',
'         select case (name)']
groups=['oneway_equalvar','oneway_welch','kruskalwallis','waerden','bartlett','flignerkilleen','levene','brownforsythe']
for name in exports:
    n=name[4:]; pair=n in ['t_equalvar','t_paired','t_welch','cor_pearson','f_var','wilcoxon_paired','wilcoxon_twosample','kolmogorovsmirnov_twosample']
    args=['x']+(['y'] if pair else [])+(['g'] if n in groups else [])+(['t'] if n=='cosinor' else [])
    if n.startswith('t_') or n.startswith('wilcoxon') or n=='f_var': args+=['null = mu']
    if n.startswith('t_') or n.startswith('wilcoxon') or n in ['f_var','cor_pearson','kolmogorovsmirnov_twosample']: args+=['alternative = alt']
    if n.startswith('t_') or n in ['f_var','cor_pearson']: args+=['conf_level = conf']
    s += [f"         case ('{name}')"]
    exact=n.startswith('wilcoxon') or n=='kolmogorovsmirnov_twosample'
    if exact: s += ['            if (mode == 4) then']
    for variant in ([True,False] if exact else [False]):
        more=['exact = .false.'] if variant else []
        if variant and n.startswith('wilcoxon'): more+=['correct = .false.']
        for axis,var in [('row','r'),('col','c')]:
            aa=[('transpose('+a+')' if a in ['x','y'] and axis=='col' else a) for a in args]+more
            line=f'            {var} = {axis}_{n}('+', '.join(aa)+')'
            if len(line)>125:
                sep=line.rfind(', ',0,110);line=line[:sep+1]+' &\n               '+line[sep+2:]
            s += [line]
        if variant: s += ['            else']
    if exact: s += ['            end if']
s += ["         case default", "            error stop 'unknown fixture family'", '         end select',
'         old_name = name', '         old_mode = mode', '      end if',
"      if (token == 'NA' .or. token == 'NaN') then", '         expected = na', "      else if (token == 'Inf') then",
'         expected = ieee_value(0.0_dp, ieee_positive_inf)', "      else if (token == '-Inf') then", '         expected = ieee_value(0.0_dp, ieee_negative_inf)',
'      else', '         read(token, *) expected', '      end if',
'      actual = value(r(row), field)', '      column_value = value(c(row), field)',
'      tolerance = 2.0e-8_dp',
"      if (name == 'row_kolmogorovsmirnov_twosample' .and. mode == 4) tolerance = 5.0e-7_dp",
'      if (.not. same(actual, expected, tolerance)) then',
'         failures = failures + 1', "         print *, trim(name), mode, row, trim(field), ' got/expected ', actual, expected", '      end if',
'      if (.not. same(actual, column_value, 1.0e-13_dp)) error stop "row/column disagreement"',
'      checks = checks + 1', '   end do', '   close(unit)', '   print *, checks, "actual-R fields checked, failures:", failures',
'   if (failures /= 0) error stop "R reference mismatch"', 'contains',
'   pure logical function same(a, b, tolerance) result(equal)', '      real(dp), intent(in) :: a !! Observed Fortran result.',
'      real(dp), intent(in) :: b !! Expected R result, including missing and infinite diagnostics.',
'      real(dp), intent(in) :: tolerance !! Relative/absolute comparison tolerance.',
'      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then', '         equal = ieee_is_nan(a) .and. ieee_is_nan(b)',
'      else if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then', '         equal = a == b',
'      else', '         equal = abs(a - b) <= tolerance * max(1.0_dp, abs(b))', '      end if', '   end function same',
'   pure real(dp) function value(r, field) result(v)',
'      type(test_result), intent(in) :: r !! One test result containing the union of numerical fields.',
'      character(*), intent(in) :: field !! Exact upstream numerical column name to compare.', '      select case (field)']
fields=['obs','obs.x','obs.y','obs.tot','obs.paired','obs.groups','mean','mean.x','mean.y','mean.diff','var','var.x','var.y','var.diff','var.pooled','var.ratio','stderr','df','df.num','df.denom','df.between','df.within','df.model','df.residual','statistic','pvalue','conf.low','conf.high','cor','skewness','kurtosis','sumsq.between','sumsq.within','meansq.between','meansq.within','mesor','amplitude','acrophase','rsquared','period','exact','corrected']
for f in fields:
    dest={'df.between':'df_num','df.model':'df_num','df.within':'df_denom','df.residual':'df_denom'}.get(f,f.replace('.','_'))
    vv='r%'+dest
    if f.startswith('obs'):vv='real('+vv+', dp)'
    if f in ['exact','corrected']:vv='merge(1.0_dp, 0.0_dp, '+vv+')'
    s += [f"      case ('{f}')",f'         v = {vv}']
s += ['      case default', '         v = ieee_value(0.0_dp, ieee_quiet_nan)', '      end select', '   end function value', 'end program test_reference']
(root/'test/test_reference.f90').write_text('\n'.join(s)+'\n')
