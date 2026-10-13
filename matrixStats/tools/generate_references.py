#!/usr/bin/env python3
"""Regenerate deterministic matrixStats 1.5.0 references using Rscript."""
from pathlib import Path
import subprocess, json
root=Path(__file__).resolve().parents[1]
cases=[]
def case(label,r,f,kind='vector'): cases.append(dict(label=label,r=r,f=f,kind=kind))
# All row/column numerical reduction mappings.
reducers={'Sums2':'sums2','Means2':'means2','Vars':'vars','Sds':'sds','Medians':'medians','Mads':'mads','IQRs':'iqrs','Mins':'mins','Maxs':'maxs','LogSumExps':'log_sum_exps','Counts':'counts','Alls':'alls','Anys':'anys','Prods':'prods'}
for p in ('row','col'):
 for r,f in reducers.items():
  case(p+r,f'{p}{r}(x)',f'{p}_{f}(x)')
  case(p+r+' selected',f'{p}{r}(xn,rows=c(4,1,4),cols=c(5,2,3),na.rm=TRUE)',f'{p}_{f}(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.)')
  case(p+r+' missing',f'{p}{r}(xn)',f'{p}_{f}(xn)')
 for r,f in {'Ranges':'ranges','Quantiles':'quantiles','Ranks':'ranks','Cumsums':'cumsums','Cumprods':'cumprods','Cummins':'cummins','Cummaxs':'cummaxs','Diffs':'diffs'}.items():
  ra=',probs=c(0,.25,.5,.75,1)' if r=='Quantiles' else (',ties.method="average"' if r=='Ranks' else '')
  fa=', probs = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]' if r=='Quantiles' else (", ties_method = 'average'" if r=='Ranks' else '')
  case(p+r,f'{p}{r}(x{ra})',f'{p}_{f}(x{fa})','matrix')
  case(p+r+' missing',f'{p}{r}(xn{ra})',f'{p}_{f}(xn{fa})','matrix') if r not in ('Quantiles',) else None
 case(p+'OrderStats',f'{p}OrderStats(x,which=2)',f'{p}_order_stats(x, 2)')
 for r,f in {'AnyMissings':'any_missings','AnyNAs':'any_nas'}.items(): case(p+r,f'{p}{r}(xn)',f'merge(1.0_dp, 0.0_dp, {p}_{f}(xn))')
 for r,f in {'Var':'var','Sd':'sd','Mad':'mad','IQR':'iqr'}.items():
  case(p+r+'Diffs',f'{p}{r}Diffs(x,diff=2)',f'{p}_{f}_diffs(x, diff = 2)')
  case(p+r+'Diffs missing',f'{p}{r}Diffs(xn,na.rm=TRUE,trim=.2)',f'{p}_{f}_diffs(xn, na_rm = .true., trim = 0.2_dp)')
 for r,f in {'Means':'means','Medians':'medians','Vars':'vars','Sds':'sds','Mads':'mads'}.items():
  weight='wr' if p=='row' else 'wc'
  case(p+'Weighted'+r+' no weights',f'{p}Weighted{r}(x)',f'{p}_weighted_{f}(x)')
  case(p+'Weighted'+r,f'{p}Weighted{r}(x,w={weight})',f'{p}_weighted_{f}(x, {weight})')
  case(p+'Weighted'+r+' missing',f'{p}Weighted{r}(xn,w={weight},na.rm=TRUE,rows=c(4,2),cols=c(5,2,1))',f'{p}_weighted_{f}(xn, {weight}, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.)')
 case(p+'Tabulates',f'{p}Tabulates(matrix(as.integer(xn),4,5),values=c(1,4,NA))',f'real({p}_tabulates(xn, [1.0_dp, 4.0_dp, nan]), dp)','matrix')
 case(p+'Collapse',f'{p}Collapse(x,idxs=c(3,1,2))',f'{p}_collapse(x, [3, 1, 2])')
for r,f in {'sum2':'sum2','mean2':'mean2','product':'product','allValue':'all_value','anyValue':'any_value','count':'count_value','logSumExp':'log_sum_exp','iqr':'iqr'}.items():
 case(r,f'{r}(v)',f'{f}(v)','scalar')
 case(r+' subset',f'{r}(vn,idxs=c(6,1,4,2),na.rm=TRUE)',f'{f}(vn, idxs = [6, 1, 4, 2], na_rm = .true.)','scalar')
case('anyMissing','anyMissing(vn)','merge(1.0_dp, 0.0_dp, any_missing(vn))','scalar')
case('diff2','diff2(v,lag=2,differences=2)','diff2(v, lag = 2, differences = 2)')
for r,f in {'varDiff':'var_diff','sdDiff':'sd_diff','madDiff':'mad_diff','iqrDiff':'iqr_diff'}.items():
 case(r,f'{r}(v,diff=2,trim=.2)',f'{f}(v, diff = 2, trim = 0.2_dp)','scalar')
for r,f in {'Mean':'mean','Median':'median','Var':'var','Sd':'sd','Mad':'mad'}.items():
 case('weighted'+r,f'weighted{r}(v,w)',f'weighted_{f}(v, w)','scalar')
 case('weighted'+r+' missing subset',f'weighted{r}(vn,w,idxs=c(6,3,2,1),na.rm=TRUE)',f'weighted_{f}(vn, w, idxs = [6, 3, 2, 1], na_rm = .true.)','scalar')
case('binCounts','binCounts(v,bx=c(0,2,4,8,10))','real(bin_counts(v, [0.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 10.0_dp]), dp)')
case('binCounts right','binCounts(v,bx=c(0,2,4,8,10),right=TRUE)','real(bin_counts(v, [0.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 10.0_dp], right = .true.), dp)')
case('binMeans','binMeans(vn,v,bx=c(0,2,4,8,10))','bm')
case('signTabulate','signTabulate(c(-Inf,-2,0,NA,NaN,Inf,1))','real(sign_tabulate([ninf, -2.0_dp, 0.0_dp, nan, nan, inf, 1.0_dp]), dp)')
case('indexByRow','indexByRow(c(4,5),idxs=c(1,2,5,6,20))','real(index_by_row([4, 5], [1, 2, 5, 6, 20]), dp)')
case('rowAvgsPerColSet','rowAvgsPerColSet(x,S=matrix(c(1,3,5,2,4,NA),3,2))','row_avgs_per_col_set(x, reshape([1, 3, 5, 2, 4, 0], [3, 2]))','matrix')
case('colAvgsPerRowSet','colAvgsPerRowSet(x,S=matrix(c(1,3,4,2,4,NA),3,2))','col_avgs_per_row_set(x, reshape([1, 3, 4, 2, 4, 0], [3, 2]))','matrix')
for r,f in {'x_OP_y':'x_op_y','t_tx_OP_y':'t_tx_op_y'}.items():
 for op in '+-*/':
  case(r+op,f'{r}(xn,c(2,NA,4),OP="{op}",commute=TRUE,na.rm=TRUE)',f"{f}(xn, [2.0_dp, nan, 4.0_dp], '{op}', commute = .true., na_rm = .true.)",'matrix')
# Stable special values, selection, and parameter edge contracts.
for p in ['row','col']:
 case(p+'Means centered',f'{p}Vars(x,center={p}Means(x))',f'{p}_vars(x, center = {p}_means2(x))')
 case(p+'MAD centered',f'{p}Mads(x,center={p}Medians(x),constant=2)',f'{p}_mads(x, center = {p}_medians(x), constant = 2.0_dp)')
 case(p+'negative subset',f'{p}Sums2(x,rows=c(-1,-3),cols=c(0,2,5))',f'{p}_sums2(x, rows = [-1, -3], cols = [0, 2, 5])')
for r,f in {'sum2':'sum2','mean2':'mean2','product':'product','logSumExp':'log_sum_exp','iqr':'iqr'}.items():
 case(r+' empty',f'{r}(numeric())',f'{f}([real(dp) ::])','scalar')
 case(r+' infinities',f'{r}(c(-Inf,Inf))',f'{f}([ninf, inf])','scalar')
case('all missing decisive','allValue(c(NA,0),value=1)','all_value([nan, 0.0_dp], value = 1.0_dp)','scalar')
case('any missing decisive','anyValue(c(NA,1),value=1)','any_value([nan, 1.0_dp], value = 1.0_dp)','scalar')
case('count NA','count(c(NA,1,NaN),value=NA_real_)','count_value([nan, 1.0_dp, nan], value = nan)','scalar')
case('colRanks preserveShape','colRanks(x,ties.method="dense",preserveShape=TRUE)',"col_ranks(x, ties_method = 'dense', preserve_shape = .true.)",'matrix')
# Generate an auditable R script. Fail rather than silently use a different package.
r='''library(matrixStats)
stopifnot(as.character(packageVersion("matrixStats")) == "1.5.0")
x <- matrix(c(1,4,-2,9,2,4,-1,8,3,8,0,7,4,8,1,6,5,9,2,5),4,5)
xn <- x; xn[1,2] <- NA; xn[3,4] <- NaN
v <- c(4,1,7,1,9,2); vn <- v; vn[2] <- NA
w <- c(1,2,0,3,4,2); wr <- c(1,0,2,3,1); wc <- c(1,2,0,3)
output <- commandArgs(trailingOnly=TRUE)[1]
con <- file(output, "w")
write_ref <- function(x) {
 x <- as.numeric(x)
 text <- vapply(x, function(v) if (is.na(v)) "NaN" else if (v == Inf) "Inf" else if (v == -Inf) "-Inf" else sprintf("%.17g",v), "")
 writeLines(paste(text,collapse=","), con)
}
'''
for c in cases:r+=f'# {c["label"]}\nwrite_ref({c["r"]})\n'
r+='close(con)\n'
(root/'test/generate_reference.R').write_text(r)
out=root/'test/reference_values.txt'
subprocess.run(['Rscript',str(root/'test/generate_reference.R'),str(out)],check=True)
lines=out.read_text().splitlines()
assert len(lines)==len(cases)
(root/'test/reference_cases.json').write_text(json.dumps(cases,indent=2)+'\n')
def literal(v):
 if v=='NaN':return 'nan'
 if v=='Inf':return 'inf'
 if v=='-Inf':return 'ninf'
 if '.' not in v and 'e' not in v:v+='.0'
 return v+'_dp'
f='''! SPDX-License-Identifier: Artistic-2.0
! Generated deterministic references: R 4.5.0, matrixStats 1.5.0.
program test_reference
   use matrixstats
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
   implicit none
   real(dp) :: x(4, 5), xn(4, 5), v(6), vn(6), w(6), wr(5), wc(4)
   real(dp) :: nan, inf, ninf
   real(dp), allocatable :: bm(:)
   integer :: checks, failures
   checks = 0
   failures = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   ninf = -inf
   x = reshape([1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, 2.0_dp, 4.0_dp, -1.0_dp, 8.0_dp, &
      3.0_dp, 8.0_dp, 0.0_dp, 7.0_dp, 4.0_dp, 8.0_dp, 1.0_dp, 6.0_dp, 5.0_dp, 9.0_dp, 2.0_dp, 5.0_dp], [4, 5])
   xn = x
   xn(1, 2) = nan
   xn(3, 4) = nan
   v = [4.0_dp, 1.0_dp, 7.0_dp, 1.0_dp, 9.0_dp, 2.0_dp]
   vn = v
   vn(2) = nan
   w = [1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 4.0_dp, 2.0_dp]
   wr = [1.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 1.0_dp]
   wc = [1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp]
   call bin_means(vn, v, [0.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 10.0_dp], bm)
'''
for c,values in zip(cases,lines):
 vals=[literal(x) for x in values.split(',') if x]
 expected='[real(dp) :: '+', &\n      '.join(', '.join(vals[j:j+4]) for j in range(0,len(vals),4))+']'
 expr=c['f']
 if c['kind']=='scalar':expr='['+expr+']'
 routine='check_matrix' if c['kind']=='matrix' else 'check'
 # avoid long lines in generated Fortran expressions
 import textwrap
 call=f"   call {routine}('{c['label']}', {expr}, &\n      {expected})\n"
 # Break long function argument lists only at comma-space boundaries.
 ls=[]
 for line in call.splitlines():
  while len(line)>124:
   pos=line.rfind(', ',0,115)
   if pos<0:raise ValueError(line)
   ls.append(line[:pos+1]+' &');line='      '+line[pos+2:]
  ls.append(line)
 f+='\n'.join(ls)+'\n'
f+='''   if (failures > 0) error stop 1
   print '(a,i0,a)', 'matrixStats R references: ', checks, ' passed'
contains
   subroutine check(label, actual, expected)
      character(*), intent(in) :: label !! Description of the upstream R reference case.
      real(dp), intent(in) :: actual(:) !! Computed values, flattened in Fortran/R column-major order.
      real(dp), intent(in) :: expected(:) !! Independently generated matrixStats 1.5.0 reference values.
      integer :: i
      logical :: ok
      if (size(actual) /= size(expected)) then
         print *, label, ': shape mismatch ', size(actual), size(expected)
         error stop 1
      end if
      do i = 1, size(actual)
         if (ieee_is_nan(expected(i))) then
            ok = ieee_is_nan(actual(i))
         else if (abs(expected(i)) == inf) then
            ok = actual(i) == expected(i)
         else
            ok = .not. ieee_is_nan(actual(i))
            if (ok) ok = abs(actual(i) - expected(i)) <= 5.0e-11_dp * max(1.0_dp, abs(expected(i)))
         end if
         if (.not. ok) then
            print *, label, ': value ', i, 'actual ', actual(i), ' expected ', expected(i)
            failures = failures + 1
         end if
      end do
      checks = checks + 1
   end subroutine check
   subroutine check_matrix(label, actual, expected)
      character(*), intent(in) :: label !! Description of the matrix-valued reference case.
      real(dp), intent(in) :: actual(:, :) !! Computed matrix in the documented output orientation.
      real(dp), intent(in) :: expected(:) !! Column-major flattened upstream reference matrix.
      call check(label, reshape(actual, [size(actual)]), expected)
   end subroutine check_matrix
end program test_reference
'''
(root/'test/test_reference.f90').write_text(f)
print('Generated',len(cases),'reference comparisons')
