! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on SharpeR, copyright 2012-2025 Steven E. Pav.
module sharper_distributions
   use sharper_kinds, only: dp, huge_dp
   use sharper_math, only: student_t_pdf, student_t_cdf, student_t_quantile
   use sharper_math, only: nct_pdf, nct_cdf, nct_quantile
   use sharper_math, only: f_pdf, f_cdf, f_quantile, ncf_pdf, ncf_cdf, ncf_quantile
   use sharper_math, only: random_student_t, random_nct, random_f, random_ncf
   use sharper_math, only: random_normal, random_chisq, clamp_probability
   implicit none
   private

   public :: drt, prt, qrt, rrt
   public :: dsr, psr, qsr, rsr
   public :: dt2, pt2, qt2, rt2
   public :: dsropt, psropt, qsropt, rsropt
   public :: plambdap, qlambdap, rlambdap
   public :: pco_f, qco_f, pco_sropt, qco_sropt
   public :: t2_to_f, f_to_t2, t2_to_sropt, sropt_to_t2

contains

   pure function use_lower_tail(value, lower_tail) result(v)
      real(dp), intent(in) :: value
      logical, intent(in), optional :: lower_tail
      real(dp) :: v
      if (present(lower_tail)) then
         if (lower_tail) then
            v = value
         else
            v = 1.0_dp-value
         end if
      else
         v = value
      end if
      v = clamp_probability(v)
   end function use_lower_tail

   pure elemental function t2_to_f(t2, df1, df2) result(fv)
      real(dp), intent(in) :: t2, df1, df2
      real(dp) :: fv
      fv = t2*(df2-df1)/(df1*(df2-1.0_dp))
   end function t2_to_f

   pure elemental function f_to_t2(fv, df1, df2) result(t2)
      real(dp), intent(in) :: fv, df1, df2
      real(dp) :: t2
      t2 = fv*(df1*(df2-1.0_dp))/(df2-df1)
   end function f_to_t2

   pure elemental function t2_to_sropt(t2, n) result(z)
      real(dp), intent(in) :: t2, n
      real(dp) :: z
      z = sqrt(max(0.0_dp,t2)/n)
   end function t2_to_sropt

   pure elemental function sropt_to_t2(z, n) result(t2)
      real(dp), intent(in) :: z, n
      real(dp) :: t2
      t2 = n*max(0.0_dp,z)**2
   end function sropt_to_t2

   function drt(x, df, k, rho, log_density) result(v)
      real(dp), intent(in) :: x, df, k
      real(dp), intent(in), optional :: rho
      logical, intent(in), optional :: log_density
      real(dp) :: v, delta
      logical :: want_log
      want_log = .false.
      if (present(log_density)) want_log = log_density
      delta = 0.0_dp
      if (present(rho)) delta = rho/k
      if (abs(delta) <= 1.0e-14_dp) then
         v = student_t_pdf(x/k,df)/k
      else
         v = nct_pdf(x/k,df,delta)/k
      end if
      if (want_log) v = log(max(v,tiny(1.0_dp)))
   end function drt

   function prt(q, df, k, rho, lower_tail) result(v)
      real(dp), intent(in) :: q, df, k
      real(dp), intent(in), optional :: rho
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, delta
      delta = 0.0_dp
      if (present(rho)) delta = rho/k
      if (abs(delta) <= 1.0e-14_dp) then
         v = student_t_cdf(q/k,df)
      else
         v = nct_cdf(q/k,df,delta)
      end if
      v = use_lower_tail(v,lower_tail)
   end function prt

   function qrt(p, df, k, rho, lower_tail) result(v)
      real(dp), intent(in) :: p, df, k
      real(dp), intent(in), optional :: rho
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, pp, delta
      pp = p
      if (present(lower_tail)) then
         if (.not. lower_tail) pp = 1.0_dp-p
      end if
      delta = 0.0_dp
      if (present(rho)) delta = rho/k
      if (abs(delta) <= 1.0e-14_dp) then
         v = k*student_t_quantile(pp,df)
      else
         v = k*nct_quantile(pp,df,delta)
      end if
   end function qrt

   function rrt(df, k, rho) result(v)
      real(dp), intent(in) :: df, k
      real(dp), intent(in), optional :: rho
      real(dp) :: v, delta
      delta = 0.0_dp
      if (present(rho)) delta = rho/k
      if (abs(delta) <= 1.0e-14_dp) then
         v = k*random_student_t(df)
      else
         v = k*random_nct(df,delta)
      end if
   end function rrt

   function dsr(x, df, zeta, ope, log_density) result(v)
      real(dp), intent(in) :: x, df
      real(dp), intent(in), optional :: zeta, ope
      logical, intent(in), optional :: log_density
      real(dp) :: v, k, local_ope
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      k = sqrt(local_ope/df)
      if (present(zeta)) then
         v = drt(x,df-1.0_dp,k,zeta,log_density)
      else
         v = drt(x,df-1.0_dp,k,log_density=log_density)
      end if
   end function dsr

   function psr(q, df, zeta, ope, lower_tail) result(v)
      real(dp), intent(in) :: q, df
      real(dp), intent(in), optional :: zeta, ope
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, k, local_ope
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      k = sqrt(local_ope/df)
      if (present(zeta)) then
         v = prt(q,df-1.0_dp,k,zeta,lower_tail)
      else
         v = prt(q,df-1.0_dp,k,lower_tail=lower_tail)
      end if
   end function psr

   function qsr(p, df, zeta, ope, lower_tail) result(v)
      real(dp), intent(in) :: p, df
      real(dp), intent(in), optional :: zeta, ope
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, k, local_ope
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      k = sqrt(local_ope/df)
      if (present(zeta)) then
         v = qrt(p,df-1.0_dp,k,zeta,lower_tail)
      else
         v = qrt(p,df-1.0_dp,k,lower_tail=lower_tail)
      end if
   end function qsr

   function rsr(df, zeta, ope) result(v)
      real(dp), intent(in) :: df
      real(dp), intent(in), optional :: zeta, ope
      real(dp) :: v, k, local_ope
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      k = sqrt(local_ope/df)
      if (present(zeta)) then
         v = rrt(df-1.0_dp,k,zeta)
      else
         v = rrt(df-1.0_dp,k)
      end if
   end function rsr

   function dt2(x, df1, df2, delta2, log_density) result(v)
      real(dp), intent(in) :: x, df1, df2
      real(dp), intent(in), optional :: delta2
      logical, intent(in), optional :: log_density
      real(dp) :: v, fv, jac
      logical :: want_log
      want_log = .false.
      if (present(log_density)) want_log = log_density
      fv = t2_to_f(x,df1,df2)
      jac = (df2-df1)/(df1*(df2-1.0_dp))
      if (present(delta2)) then
         v = ncf_pdf(fv,df1,df2-df1,delta2)*jac
      else
         v = f_pdf(fv,df1,df2-df1)*jac
      end if
      if (want_log) v = log(max(v,tiny(1.0_dp)))
   end function dt2

   function pt2(q, df1, df2, delta2, lower_tail) result(v)
      real(dp), intent(in) :: q, df1, df2
      real(dp), intent(in), optional :: delta2
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, fv
      fv = t2_to_f(q,df1,df2)
      if (present(delta2)) then
         v = ncf_cdf(fv,df1,df2-df1,delta2)
      else
         v = f_cdf(fv,df1,df2-df1)
      end if
      v = use_lower_tail(v,lower_tail)
   end function pt2

   function qt2(p, df1, df2, delta2, lower_tail) result(v)
      real(dp), intent(in) :: p, df1, df2
      real(dp), intent(in), optional :: delta2
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, pp, fq
      pp = p
      if (present(lower_tail)) then
         if (.not. lower_tail) pp = 1.0_dp-p
      end if
      if (present(delta2)) then
         fq = ncf_quantile(pp,df1,df2-df1,delta2)
      else
         fq = f_quantile(pp,df1,df2-df1)
      end if
      v = f_to_t2(fq,df1,df2)
   end function qt2

   function rt2(df1, df2, delta2) result(v)
      real(dp), intent(in) :: df1, df2
      real(dp), intent(in), optional :: delta2
      real(dp) :: v, fv
      if (present(delta2)) then
         fv = random_ncf(df1,df2-df1,delta2)
      else
         fv = random_f(df1,df2-df1)
      end if
      v = f_to_t2(fv,df1,df2)
   end function rt2

   function dsropt(x, df1, df2, zeta_s, ope, drag, log_density) result(v)
      real(dp), intent(in) :: x, df1, df2
      real(dp), intent(in), optional :: zeta_s, ope, drag
      logical, intent(in), optional :: log_density
      real(dp) :: v, xx, zz, local_ope, local_drag, t2, delta2, jac
      logical :: want_log
      local_ope = 1.0_dp
      local_drag = 0.0_dp
      want_log = .false.
      if (present(ope)) local_ope = ope
      if (present(drag)) local_drag = drag
      if (present(log_density)) want_log = log_density
      xx = (x+local_drag)/sqrt(local_ope)
      if (xx <= 0.0_dp) then
         if (want_log) then
            v = log(tiny(1.0_dp))
         else
            v = 0.0_dp
         end if
         return
      end if
      delta2 = 0.0_dp
      if (present(zeta_s)) then
         zz = zeta_s/sqrt(local_ope)
         delta2 = sropt_to_t2(zz,df2)
      end if
      t2 = sropt_to_t2(xx,df2)
      jac = 2.0_dp*df2*xx/sqrt(local_ope)
      if (want_log) then
         v = dt2(t2,df1,df2,delta2,.true.)+log(jac)
      else
         v = dt2(t2,df1,df2,delta2)*jac
      end if
   end function dsropt

   function psropt(q, df1, df2, zeta_s, ope, drag, lower_tail) result(v)
      real(dp), intent(in) :: q, df1, df2
      real(dp), intent(in), optional :: zeta_s, ope, drag
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, qq, zz, local_ope, local_drag, delta2
      local_ope = 1.0_dp
      local_drag = 0.0_dp
      if (present(ope)) local_ope = ope
      if (present(drag)) local_drag = drag
      qq = (q+local_drag)/sqrt(local_ope)
      delta2 = 0.0_dp
      if (present(zeta_s)) then
         zz = zeta_s/sqrt(local_ope)
         delta2 = sropt_to_t2(zz,df2)
      end if
      v = pt2(sropt_to_t2(qq,df2),df1,df2,delta2,lower_tail)
   end function psropt

   function qsropt(p, df1, df2, zeta_s, ope, drag, lower_tail) result(v)
      real(dp), intent(in) :: p, df1, df2
      real(dp), intent(in), optional :: zeta_s, ope, drag
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, zz, local_ope, local_drag, delta2
      local_ope = 1.0_dp
      local_drag = 0.0_dp
      if (present(ope)) local_ope = ope
      if (present(drag)) local_drag = drag
      delta2 = 0.0_dp
      if (present(zeta_s)) then
         zz = zeta_s/sqrt(local_ope)
         delta2 = sropt_to_t2(zz,df2)
      end if
      v = sqrt(local_ope)*t2_to_sropt(qt2(p,df1,df2,delta2,lower_tail),df2)-local_drag
   end function qsropt

   function rsropt(df1, df2, zeta_s, ope, drag) result(v)
      real(dp), intent(in) :: df1, df2
      real(dp), intent(in), optional :: zeta_s, ope, drag
      real(dp) :: v, zz, local_ope, local_drag, delta2
      local_ope = 1.0_dp
      local_drag = 0.0_dp
      if (present(ope)) local_ope = ope
      if (present(drag)) local_drag = drag
      delta2 = 0.0_dp
      if (present(zeta_s)) then
         zz = zeta_s/sqrt(local_ope)
         delta2 = sropt_to_t2(zz,df2)
      end if
      v = sqrt(local_ope)*t2_to_sropt(rt2(df1,df2,delta2),df2)-local_drag
   end function rsropt

   function plambdap(q, df, tstat, lower_tail) result(v)
      real(dp), intent(in) :: q, df, tstat
      logical, intent(in), optional :: lower_tail
      real(dp) :: v
      v = 1.0_dp-nct_cdf(tstat,df,q)
      if (present(lower_tail)) then
         if (.not. lower_tail) v = 1.0_dp-v
      end if
      v = clamp_probability(v)
   end function plambdap

   function qlambdap(p, df, tstat, lower_tail) result(v)
      real(dp), intent(in) :: p, df, tstat
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, target, lo, hi, mid, pmid
      integer :: iter
      logical :: lower
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (p <= 0.0_dp) then
         if (lower) then
            v = -huge_dp
         else
            v = huge_dp
         end if
         return
      else if (p >= 1.0_dp) then
         if (lower) then
            v = huge_dp
         else
            v = -huge_dp
         end if
         return
      end if
      target = p
      lo = min(-1.0_dp,tstat-12.0_dp)
      hi = max(1.0_dp,tstat+12.0_dp)
      do iter = 1, 160
         mid = 0.5_dp*(lo+hi)
         pmid = plambdap(mid,df,tstat,lower)
         if (pmid < target) then
            lo = mid
         else
            hi = mid
         end if
      end do
      v = 0.5_dp*(lo+hi)
   end function qlambdap

   function rlambdap(df, tstat) result(v)
      real(dp), intent(in) :: df, tstat
      real(dp) :: v
      v = random_normal()+tstat*sqrt(random_chisq(df)/df)
   end function rlambdap

   function pco_f(q, df1, df2, x, lower_tail) result(v)
      real(dp), intent(in) :: q, df1, df2, x
      logical, intent(in), optional :: lower_tail
      real(dp) :: v
      v = 1.0_dp-ncf_cdf(x,df1,df2,q)
      if (present(lower_tail)) then
         if (.not. lower_tail) v = 1.0_dp-v
      end if
      v = clamp_probability(v)
   end function pco_f

   function qco_f(p, df1, df2, x, lower_tail, lower_bound, upper_bound) result(v)
      real(dp), intent(in) :: p, df1, df2, x
      logical, intent(in), optional :: lower_tail
      real(dp), intent(in), optional :: lower_bound, upper_bound
      real(dp) :: v, lo, hi, mid, pmid, target
      logical :: lower
      integer :: iter
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      lo = 0.0_dp
      if (present(lower_bound)) lo = max(0.0_dp,lower_bound)
      hi = max(1.0_dp,lo+1.0_dp)
      if (present(upper_bound)) hi = upper_bound
      target = p
      if (.not. present(upper_bound)) then
         do while (pco_f(hi,df1,df2,x,lower) < target .and. hi < huge_dp/4.0_dp)
            hi = 2.0_dp*hi
         end do
      end if
      do iter = 1, 140
         mid = 0.5_dp*(lo+hi)
         pmid = pco_f(mid,df1,df2,x,lower)
         if (pmid < target) then
            lo = mid
         else
            hi = mid
         end if
      end do
      v = 0.5_dp*(lo+hi)
   end function qco_f

   function pco_sropt(q, df1, df2, z_s, ope, lower_tail) result(v)
      real(dp), intent(in) :: q, df1, df2, z_s
      real(dp), intent(in), optional :: ope
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, local_ope
      logical :: lower
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      v = psropt(z_s,df1,df2,q,local_ope,lower_tail=.not. lower)
   end function pco_sropt

   function qco_sropt(p, df1, df2, z_s, ope, lower_tail, lower_bound, upper_bound) result(v)
      real(dp), intent(in) :: p, df1, df2, z_s
      real(dp), intent(in), optional :: ope, lower_bound, upper_bound
      logical, intent(in), optional :: lower_tail
      real(dp) :: v, lo, hi, mid, pmid, local_ope
      logical :: lower
      integer :: iter
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      lo = 0.0_dp
      if (present(lower_bound)) lo = max(0.0_dp,lower_bound)
      hi = max(1.0_dp,lo+1.0_dp)
      if (present(upper_bound)) hi = upper_bound
      if (.not. present(upper_bound)) then
         do while (pco_sropt(hi,df1,df2,z_s,local_ope,lower) < p .and. hi < huge_dp/4.0_dp)
            hi = 2.0_dp*hi
         end do
      end if
      do iter = 1, 140
         mid = 0.5_dp*(lo+hi)
         pmid = pco_sropt(mid,df1,df2,z_s,local_ope,lower)
         if (pmid < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      v = 0.5_dp*(lo+hi)
   end function qco_sropt

end module sharper_distributions
