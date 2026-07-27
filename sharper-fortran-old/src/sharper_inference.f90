! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on SharpeR, copyright 2012-2025 Steven E. Pav.
module sharper_inference
   use sharper_kinds, only: dp
   use sharper_types, only: sropt_result, del_sropt_result
   use sharper_math, only: ncf_pdf, ncchisq_quantile
   use sharper_distributions, only: t2_to_f, qco_sropt
   implicit none
   private

   type, public :: spanning_f_result
      real(dp) :: f_value = 0.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
      real(dp) :: r1 = 0.0_dp
      real(dp) :: effective_scale = 0.0_dp
   end type spanning_f_result

   public :: f_ncp_unbiased, f_ncp_krs, f_ncp_mle
   public :: f_inference, t2_inference
   public :: infer_sropt, infer_del_sropt
   public :: sropt_confint, achieved_snr_confint
   public :: achieved_delta_snr_confint, del_sropt_as_f

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

   pure elemental function f_ncp_unbiased(f_value, df1, df2) result(value)
      real(dp), intent(in) :: f_value, df1, df2
      real(dp) :: value
      value = f_value*(df2-2.0_dp)*df1/df2-df1
   end function f_ncp_unbiased

   pure elemental function f_ncp_krs(f_value, df1, df2) result(value)
      real(dp), intent(in) :: f_value, df1, df2
      real(dp) :: value, xbs, delta0, phi2
      xbs = f_value*df1/df2
      delta0 = (df2-2.0_dp)*xbs-df1
      phi2 = 2.0_dp*xbs*(df2-2.0_dp)/(df1+2.0_dp)
      value = max(delta0,phi2)
   end function f_ncp_krs

   function f_ncp_mle(f_value, df1, df2) result(value)
      real(dp), intent(in) :: f_value, df1, df2
      real(dp) :: value, lo, hi, x1, x2, f1, f2, ratio
      integer :: iter
      if (f_value <= 1.0_dp) then
         value = 0.0_dp
         return
      end if
      lo = 0.0_dp
      hi = 1.0_dp
      f1 = log(max(ncf_pdf(f_value,df1,df2,hi),tiny(1.0_dp)))
      do
         f2 = log(max(ncf_pdf(f_value,df1,df2,2.0_dp*hi),tiny(1.0_dp)))
         if (f2 <= f1 .or. hi > 1.0e8_dp) exit
         hi = 2.0_dp*hi
         f1 = f2
      end do
      hi = 2.0_dp*hi
      ratio = 0.5_dp*(sqrt(5.0_dp)-1.0_dp)
      x1 = hi-ratio*(hi-lo)
      x2 = lo+ratio*(hi-lo)
      f1 = log(max(ncf_pdf(f_value,df1,df2,x1),tiny(1.0_dp)))
      f2 = log(max(ncf_pdf(f_value,df1,df2,x2),tiny(1.0_dp)))
      do iter = 1, 120
         if (f1 < f2) then
            lo = x1
            x1 = x2
            f1 = f2
            x2 = lo+ratio*(hi-lo)
            f2 = log(max(ncf_pdf(f_value,df1,df2,x2),tiny(1.0_dp)))
         else
            hi = x2
            x2 = x1
            f2 = f1
            x1 = hi-ratio*(hi-lo)
            f1 = log(max(ncf_pdf(f_value,df1,df2,x1),tiny(1.0_dp)))
         end if
      end do
      value = 0.5_dp*(lo+hi)
   end function f_ncp_mle

   function f_inference(f_value, df1, df2, method) result(value)
      real(dp), intent(in) :: f_value, df1, df2
      character(len=*), intent(in), optional :: method
      real(dp) :: value
      character(len=24) :: local_method
      local_method = 'krs'
      if (present(method)) local_method = lowercase(trim(method))
      select case (trim(local_method))
      case ('mle')
         value = f_ncp_mle(f_value,df1,df2)
      case ('unbiased')
         value = f_ncp_unbiased(f_value,df1,df2)
      case default
         value = f_ncp_krs(f_value,df1,df2)
      end select
   end function f_inference

   function t2_inference(t2, df1, df2, method) result(value)
      real(dp), intent(in) :: t2, df1, df2
      character(len=*), intent(in), optional :: method
      real(dp) :: value, f_value
      f_value = t2_to_f(t2,df1,df2)
      value = f_inference(f_value,df1,df2-df1,method)
   end function t2_inference

   function infer_sropt(z, method) result(value)
      type(sropt_result), intent(in) :: z
      character(len=*), intent(in), optional :: method
      real(dp) :: value, delta2
      delta2 = t2_inference(z%t2,real(z%df1,dp),real(z%df2,dp),method)
      value = sqrt(z%ope)*sqrt(max(0.0_dp,delta2)/real(z%df2,dp))-z%drag
   end function infer_sropt

   pure function del_sropt_as_f(z) result(out)
      type(del_sropt_result), intent(in) :: z
      type(spanning_f_result) :: out
      real(dp) :: nmin, znum, zdenom, beta_value
      nmin = real(z%df2-1,dp)
      znum = 1.0_dp+z%t2_sub/nmin
      zdenom = 1.0_dp+z%t2/nmin
      out%r1 = 1.0_dp-1.0_dp/znum
      beta_value = znum/zdenom
      out%df1 = real(z%df1-z%df1_sub,dp)
      out%df2 = real(z%df2-z%df1,dp)
      out%f_value = out%df2*(1.0_dp-beta_value)/(out%df1*beta_value)
      out%effective_scale = real(z%df2,dp)*(1.0_dp-out%r1)
   end function del_sropt_as_f

   function infer_del_sropt(z, method) result(value)
      type(del_sropt_result), intent(in) :: z
      character(len=*), intent(in), optional :: method
      real(dp) :: value, delta2
      type(spanning_f_result) :: fdata
      fdata = del_sropt_as_f(z)
      delta2 = f_inference(fdata%f_value,fdata%df1,fdata%df2,method)
      value = sqrt(z%ope)*sqrt(max(0.0_dp,delta2/fdata%effective_scale))
   end function infer_del_sropt

   function sropt_confint(z, level) result(ci)
      type(sropt_result), intent(in) :: z
      real(dp), intent(in), optional :: level
      real(dp) :: ci(2), local_level, alpha
      local_level = 0.95_dp
      if (present(level)) local_level = level
      alpha = 0.5_dp*(1.0_dp-local_level)
      ci(2) = qco_sropt(1.0_dp-alpha,real(z%df1,dp),real(z%df2,dp),z%value,z%ope)
      ci(1) = qco_sropt(alpha,real(z%df1,dp),real(z%df2,dp),z%value,z%ope, &
                         upper_bound=ci(2))
   end function sropt_confint

   function achieved_snr_confint(z, level) result(ci)
      type(sropt_result), intent(in) :: z
      real(dp), intent(in), optional :: level
      real(dp) :: ci(2), local_level, lower_p, inferred, lambda_quarter
      real(dp) :: native_observed, adjustment
      local_level = 0.95_dp
      if (present(level)) local_level = level
      lower_p = 1.0_dp-local_level
      inferred = max(0.0_dp,infer_sropt(z,'krs')/sqrt(z%ope))
      lambda_quarter = real(z%df2,dp)*inferred*inferred/4.0_dp
      native_observed = sqrt(max(0.0_dp,z%t2/real(z%df2,dp)))
      adjustment = ncchisq_quantile(1.0_dp-lower_p,real(z%df1,dp),lambda_quarter)-lambda_quarter
      if (native_observed <= tiny(1.0_dp)) then
         ci(1) = -huge(1.0_dp)
      else
         ci(1) = sqrt(z%ope)*(native_observed-adjustment/(real(z%df2,dp)*native_observed))
      end if
      ci(2) = huge(1.0_dp)
   end function achieved_snr_confint

   function achieved_delta_snr_confint(z, level) result(ci)
      type(del_sropt_result), intent(in) :: z
      real(dp), intent(in), optional :: level
      real(dp) :: ci(2), local_level, lower_p, inferred, lambda_quarter
      real(dp) :: native_observed, adjustment, nstrat
      local_level = 0.95_dp
      if (present(level)) local_level = level
      lower_p = 1.0_dp-local_level
      inferred = max(0.0_dp,infer_del_sropt(z,'krs')/sqrt(z%ope))
      lambda_quarter = real(z%df2,dp)*inferred*inferred/4.0_dp
      nstrat = real(z%df1-z%df1_sub,dp)
      native_observed = sqrt(max(0.0_dp,z%t2_delta/real(z%df2,dp)))
      adjustment = ncchisq_quantile(1.0_dp-lower_p,nstrat,lambda_quarter)-lambda_quarter
      if (native_observed <= tiny(1.0_dp)) then
         ci(1) = -huge(1.0_dp)
      else
         ci(1) = sqrt(z%ope)*(native_observed-adjustment/(real(z%df2,dp)*native_observed))
      end if
      ci(2) = huge(1.0_dp)
   end function achieved_delta_snr_confint

end module sharper_inference
