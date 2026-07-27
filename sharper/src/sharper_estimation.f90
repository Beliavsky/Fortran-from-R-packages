! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on SharpeR, copyright 2012-2025 Steven E. Pav.
module sharper_estimation
   use sharper_kinds, only: dp
   use sharper_types, only: sr_result, sropt_result, del_sropt_result
   use sharper_types, only: moment_vcov_result
   use sharper_math, only: normal_quantile
   use sharper_linalg, only: column_mean, sample_covariance, covariance_of_mean
   use sharper_linalg, only: solve_linear, invert_matrix, quadratic_form
   use sharper_linalg, only: symmetric_vech, symmetric_ivech
   use sharper_distributions, only: qlambdap, pt2
   implicit none
   private

   public :: compute_sr, fit_sr, fit_sr_matrix, reannualize_sr
   public :: sample_k_statistics, sample_standardized_cumulants
   public :: sr_bias, sr_variance, sr_standard_error, sr_confint, predint
   public :: sr_vcov
   public :: markowitz_from_moments, fit_sropt, make_sropt
   public :: reannualize_sropt, sric
   public :: make_del_sropt, fit_del_sropt
   public :: sm_vcov, ism_vcov

   interface fit_sr
      module procedure fit_sr_vector
      module procedure fit_sr_matrix
   end interface fit_sr

contains

   pure elemental function compute_sr(mu, c0, sigma, ope) result(value)
      real(dp), intent(in) :: mu, c0, sigma
      real(dp), intent(in), optional :: ope
      real(dp) :: value
      value = (mu-c0)/sigma
      if (present(ope)) value = value*sqrt(ope)
   end function compute_sr

   pure function lower_string(text) result(out)
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
   end function lower_string

   function fit_sr_vector(x, c0, ope, higher_order) result(z)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: c0, ope
      logical, intent(in), optional :: higher_order
      type(sr_result) :: z
      real(dp) :: mu, sigma, local_c0, local_ope
      logical :: use_higher
      integer :: n
      n = size(x)
      local_c0 = 0.0_dp
      local_ope = 1.0_dp
      use_higher = .false.
      if (present(c0)) local_c0 = c0
      if (present(ope)) local_ope = ope
      if (present(higher_order)) use_higher = higher_order
      allocate(z%value(1),z%df(1),z%rescal(1))
      if (n < 2) then
         z%value = 0.0_dp
         z%df = 0
         z%rescal = 0.0_dp
         return
      end if
      mu = sum(x)/real(n,dp)
      sigma = sqrt(sum((x-mu)**2)/real(n-1,dp))
      z%value(1) = compute_sr(mu,local_c0,sigma,local_ope)
      z%df(1) = n-1
      z%rescal(1) = 1.0_dp/sqrt(real(n,dp))
      z%c0 = local_c0
      z%ope = local_ope
      if (use_higher) then
         allocate(z%cumulants(4,1))
         z%cumulants(:,1) = sample_standardized_cumulants(x)
      end if
   end function fit_sr_vector

   function fit_sr_matrix(x, c0, ope, higher_order) result(z)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: c0, ope
      logical, intent(in), optional :: higher_order
      type(sr_result) :: z
      real(dp) :: mu(size(x,2)), sigma(size(x,2)), local_c0, local_ope
      logical :: use_higher
      integer :: n, p, j
      n = size(x,1)
      p = size(x,2)
      local_c0 = 0.0_dp
      local_ope = 1.0_dp
      use_higher = .false.
      if (present(c0)) local_c0 = c0
      if (present(ope)) local_ope = ope
      if (present(higher_order)) use_higher = higher_order
      allocate(z%value(p),z%df(p),z%rescal(p))
      if (n < 2) then
         z%value = 0.0_dp
         z%df = 0
         z%rescal = 0.0_dp
         return
      end if
      mu = column_mean(x)
      do j = 1, p
         sigma(j) = sqrt(sum((x(:,j)-mu(j))**2)/real(n-1,dp))
      end do
      z%value = compute_sr(mu,local_c0,sigma,local_ope)
      z%df = n-1
      z%rescal = 1.0_dp/sqrt(real(n,dp))
      z%c0 = local_c0
      z%ope = local_ope
      if (use_higher) then
         allocate(z%cumulants(4,p))
         do j = 1, p
            z%cumulants(:,j) = sample_standardized_cumulants(x(:,j))
         end do
      end if
   end function fit_sr_matrix

   function reannualize_sr(z, new_ope) result(out)
      type(sr_result), intent(in) :: z
      real(dp), intent(in) :: new_ope
      type(sr_result) :: out
      out = z
      out%value = z%value*sqrt(new_ope/z%ope)
      out%ope = new_ope
   end function reannualize_sr

   function sample_k_statistics(x) result(k)
      real(dp), intent(in) :: x(:)
      real(dp) :: k(6)
      real(dp) :: s(6), nn(6), nd(6), xr(size(x))
      real(dp) :: nreal
      integer :: i, n
      n = size(x)
      k = 0.0_dp
      if (n < 6) return
      xr = x
      do i = 1, 6
         s(i) = sum(xr**i)
         nn(i) = real(n,dp)**i
         nd(i) = exp(log_gamma(real(n+1,dp))-log_gamma(real(n-i+1,dp)))
      end do
      nreal = real(n,dp)
      k(1) = s(1)/nreal
      k(2) = (nreal*s(2)-s(1)**2)/nd(2)
      k(3) = (nn(2)*s(3)-3.0_dp*nreal*s(2)*s(1)+2.0_dp*s(1)**3)/nd(3)
      k(4) = ((nn(3)+nn(2))*s(4)-4.0_dp*(nn(2)+nreal)*s(3)*s(1) - &
              3.0_dp*(nn(2)-nreal)*s(2)**2+12.0_dp*nreal*s(2)*s(1)**2 - &
              6.0_dp*s(1)**4)/nd(4)
      k(5) = ((nn(4)+5.0_dp*nn(3))*s(5)-5.0_dp*(nn(3)+5.0_dp*nn(2))*s(4)*s(1) - &
              10.0_dp*(nn(3)-nn(2))*s(3)*s(2)+20.0_dp*(nn(2)+2.0_dp*nreal)*s(3)*s(1)**2 + &
              30.0_dp*(nn(2)-nreal)*s(2)**2*s(1)-60.0_dp*nreal*s(2)*s(1)**3 + &
              24.0_dp*s(1)**5)/nd(5)
      k(6) = ((nn(5)+16.0_dp*nn(4)+11.0_dp*nn(3)-4.0_dp*nn(2))*s(6) - &
              6.0_dp*(nn(4)+16.0_dp*nn(3)+11.0_dp*nn(2)-4.0_dp*nreal)*s(5)*s(1) - &
              15.0_dp*nreal*(nreal-1.0_dp)**2*(nreal+4.0_dp)*s(4)*s(2) - &
              10.0_dp*(nn(4)-2.0_dp*nn(3)+5.0_dp*nn(2)-4.0_dp*nreal)*s(3)**2 + &
              30.0_dp*(nn(3)+9.0_dp*nn(2)+2.0_dp*nreal)*s(4)*s(1)**2 + &
              120.0_dp*(nn(3)-nreal)*s(3)*s(2)*s(1) + &
              30.0_dp*(nn(3)-3.0_dp*nn(2)+2.0_dp*nreal)*s(2)**3 - &
              120.0_dp*(nn(2)+3.0_dp*nreal)*s(3)*s(1)**3 - &
              270.0_dp*(nn(2)-nreal)*s(2)**2*s(1)**2 + &
              360.0_dp*nreal*s(2)*s(1)**4-120.0_dp*s(1)**6)/nd(6)
   end function sample_k_statistics

   function sample_standardized_cumulants(x) result(gamma_values)
      real(dp), intent(in) :: x(:)
      real(dp) :: gamma_values(4), k(6), centered(size(x)), mu
      mu = sum(x)/real(size(x),dp)
      centered = x-mu
      k = sample_k_statistics(centered)
      if (k(2) <= tiny(1.0_dp)) then
         gamma_values = 0.0_dp
      else
         gamma_values(1) = k(3)/k(2)**1.5_dp
         gamma_values(2) = k(4)/k(2)**2
         gamma_values(3) = k(5)/k(2)**2.5_dp
         gamma_values(4) = k(6)/k(2)**3
      end if
   end function sample_standardized_cumulants

   pure function sr_bias(snr, n, cumulants, method) result(value)
      real(dp), intent(in) :: snr
      integer, intent(in) :: n
      real(dp), intent(in) :: cumulants(4)
      character(len=*), intent(in), optional :: method
      real(dp) :: value, tt, s
      character(len=32) :: local_method
      tt = real(n,dp)
      s = snr
      local_method = 'simple'
      if (present(method)) local_method = lower_string(trim(method))
      if (trim(local_method) == 'second_order') then
         value = 3.0_dp*s/(4.0_dp*tt)+49.0_dp*s/(32.0_dp*tt**2) - &
            cumulants(1)*(1.0_dp/(2.0_dp*tt)+3.0_dp/(8.0_dp*tt**2)) + &
            s*cumulants(2)*(3.0_dp/(8.0_dp*tt)-15.0_dp/(32.0_dp*tt**2)) + &
            3.0_dp*cumulants(3)/(8.0_dp*tt**2)-5.0_dp*s*cumulants(4)/(16.0_dp*tt**2) - &
            5.0_dp*s*cumulants(1)**2/(4.0_dp*tt**2) + &
            105.0_dp*s*cumulants(2)**2/(128.0_dp*tt**2) - &
            15.0_dp*cumulants(1)*cumulants(2)/(16.0_dp*tt**2)
      else
         value = (2.0_dp+cumulants(2))*s*3.0_dp/(8.0_dp*tt) - &
                 cumulants(1)/(2.0_dp*tt)
      end if
   end function sr_bias

   pure function sr_variance(snr, n, cumulants) result(value)
      real(dp), intent(in) :: snr
      integer, intent(in) :: n
      real(dp), intent(in) :: cumulants(4)
      real(dp) :: value, tt, s
      tt = real(n,dp)
      s = snr
      value = (1.0_dp+s*s/2.0_dp)/tt+(19.0_dp*s*s/8.0_dp+2.0_dp)/tt**2 - &
         cumulants(1)*s*(1.0_dp/tt+2.5_dp/tt**2) + &
         cumulants(2)*s*s*(1.0_dp/(4.0_dp*tt)+3.0_dp/(8.0_dp*tt**2)) + &
         5.0_dp*cumulants(3)*s/(4.0_dp*tt**2)-3.0_dp*cumulants(4)*s*s/(8.0_dp*tt**2) + &
         cumulants(1)**2*(7.0_dp/(4.0_dp*tt**2)-3.0_dp*s*s/(2.0_dp*tt**2)) - &
         15.0_dp*cumulants(1)*cumulants(2)*s/(4.0_dp*tt**2) + &
         39.0_dp*cumulants(2)**2*s*s/(32.0_dp*tt**2)
   end function sr_variance

   function sr_standard_error(z, method) result(se)
      type(sr_result), intent(in) :: z
      character(len=*), intent(in), optional :: method
      real(dp), allocatable :: se(:)
      real(dp) :: tstat(size(z%value)), t_se(size(z%value)), s
      integer :: j
      character(len=32) :: local_method
      allocate(se(size(z%value)))
      local_method = 't'
      if (present(method)) local_method = lower_string(trim(method))
      tstat = z%value/(z%rescal*sqrt(z%ope))
      do j = 1, size(z%value)
         select case (trim(local_method))
         case ('mertens')
            if (.not. allocated(z%cumulants)) then
               se(j) = huge(1.0_dp)
               cycle
            end if
            s = 1.0_dp-z%cumulants(1,j)*tstat(j)/sqrt(real(z%df(j)+1,dp)) + &
                (z%cumulants(2,j)+2.0_dp)*tstat(j)**2/(4.0_dp*real(z%df(j)+1,dp))
            t_se(j) = sqrt(max(0.0_dp,s))
         case ('bao')
            if (.not. allocated(z%cumulants)) then
               se(j) = huge(1.0_dp)
               cycle
            end if
            s = tstat(j)/sqrt(real(z%df(j)+1,dp))
            t_se(j) = sqrt(real(z%df(j)+1,dp)*sr_variance(s,z%df(j)+1,z%cumulants(:,j)))
         case default
            t_se(j) = sqrt(1.0_dp+tstat(j)**2/(2.0_dp*real(z%df(j),dp)))
         end select
         se(j) = t_se(j)*z%rescal(j)*sqrt(z%ope)
      end do
   end function sr_standard_error

   function sr_confint(z, level, method, inflate_by) result(ci)
      type(sr_result), intent(in) :: z
      real(dp), intent(in), optional :: level, inflate_by
      character(len=*), intent(in), optional :: method
      real(dp), allocatable :: ci(:, :)
      real(dp) :: local_level, alpha, qlo, qhi, inflation
      real(dp) :: tstat(size(z%value)), tse(size(z%value)), midpoint(size(z%value))
      real(dp), allocatable :: srse(:)
      integer :: j
      character(len=32) :: local_method
      local_level = 0.95_dp
      inflation = 1.0_dp
      local_method = 'exact'
      if (present(level)) local_level = level
      if (present(inflate_by)) inflation = inflate_by
      if (present(method)) local_method = lower_string(trim(method))
      alpha = 0.5_dp*(1.0_dp-local_level)
      allocate(ci(size(z%value),2))
      tstat = z%value/(z%rescal*sqrt(z%ope))
      if (trim(local_method) == 'exact') then
         do j = 1, size(z%value)
            ci(j,1) = qlambdap(alpha,real(z%df(j),dp),tstat(j))*z%rescal(j)*sqrt(z%ope)
            ci(j,2) = qlambdap(1.0_dp-alpha,real(z%df(j),dp),tstat(j))*z%rescal(j)*sqrt(z%ope)
         end do
      else
         srse = sr_standard_error(z,local_method)
         tse = srse/(z%rescal*sqrt(z%ope))
         midpoint = tstat
         if (trim(local_method) == 'z') midpoint = tstat*(1.0_dp-1.0_dp/(4.0_dp*real(z%df,dp)))
         if (trim(local_method) == 'bao' .and. allocated(z%cumulants)) then
            do j = 1, size(z%value)
               midpoint(j) = tstat(j)-sqrt(real(z%df(j)+1,dp))* &
                  sr_bias(tstat(j)/sqrt(real(z%df(j)+1,dp)),z%df(j)+1, &
                          z%cumulants(:,j),'second_order')
            end do
         end if
         qlo = normal_quantile(alpha)
         qhi = normal_quantile(1.0_dp-alpha)
         ci(:,1) = (midpoint+qlo*inflation*tse)*z%rescal*sqrt(z%ope)
         ci(:,2) = (midpoint+qhi*inflation*tse)*z%rescal*sqrt(z%ope)
      end if
   end function sr_confint

   function predint(z, oosdf, level, method, oosrescal) result(ci)
      type(sr_result), intent(in) :: z
      integer, intent(in) :: oosdf
      real(dp), intent(in), optional :: level, oosrescal
      character(len=*), intent(in), optional :: method
      real(dp), allocatable :: ci(:, :)
      real(dp) :: rescal_value, inflation
      rescal_value = 1.0_dp/sqrt(real(oosdf+1,dp))
      if (present(oosrescal)) rescal_value = oosrescal
      inflation = sqrt(1.0_dp+real(z%df(1)+1,dp)*rescal_value**2)
      ci = sr_confint(z,level,method,inflation)
   end function predint

   function sr_vcov(x, ope) result(out)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: ope
      type(moment_vcov_result) :: out
      real(dp) :: local_ope, mu(size(x,2)), m2(size(x,2)), sig2(size(x,2))
      real(dp), allocatable :: y(:, :), shat(:, :), dmat(:, :)
      integer :: n, p, j
      n = size(x,1)
      p = size(x,2)
      local_ope = 1.0_dp
      if (present(ope)) local_ope = ope
      out%n = n
      out%p = p
      allocate(out%mean(p),out%covariance(p,p),y(n,2*p),shat(2*p,2*p),dmat(2*p,p))
      y(:,1:p) = x
      y(:,p+1:2*p) = x*x
      mu = column_mean(x)
      m2 = column_mean(x*x)
      sig2 = m2-mu*mu
      out%mean = mu/sqrt(sig2)*sqrt(local_ope)
      shat = covariance_of_mean(y)
      dmat = 0.0_dp
      do j = 1, p
         dmat(j,j) = m2(j)/sig2(j)**1.5_dp
         dmat(p+j,j) = -mu(j)/(2.0_dp*sig2(j)**1.5_dp)
      end do
      out%covariance = local_ope*matmul(transpose(dmat),matmul(shat,dmat))
   end function sr_vcov

   function markowitz_from_moments(mu, sigma, df2, ope, drag, h) result(z)
      real(dp), intent(in) :: mu(:), sigma(:, :)
      integer, intent(in) :: df2
      real(dp), intent(in), optional :: ope, drag
      real(dp), intent(in), optional :: h(:, :)
      type(sropt_result) :: z
      real(dp) :: local_ope, local_drag, zeta_sq
      real(dp), allocatable :: work_mu(:), work_sigma(:, :), work_w(:), mapped_w(:)
      integer :: p, status
      local_ope = 1.0_dp
      local_drag = 0.0_dp
      if (present(ope)) local_ope = ope
      if (present(drag)) local_drag = drag
      p = size(mu)
      if (present(h)) then
         allocate(work_mu(size(h,1)),work_sigma(size(h,1),size(h,1)),work_w(size(h,1)))
         work_mu = matmul(h,mu)
         work_sigma = matmul(h,matmul(sigma,transpose(h)))
         call solve_linear(work_sigma,work_mu,work_w,status)
         allocate(mapped_w(p))
         mapped_w = matmul(transpose(h),work_w)
         z%df1 = size(h,1)
      else
         allocate(work_mu(p),work_sigma(p,p),work_w(p),mapped_w(p))
         work_mu = mu
         work_sigma = sigma
         call solve_linear(work_sigma,work_mu,work_w,status)
         mapped_w = work_w
         z%df1 = p
      end if
      if (status /= 0) then
         z%value = 0.0_dp
         z%t2 = 0.0_dp
         return
      end if
      zeta_sq = dot_product(work_mu,work_w)
      z%df2 = df2
      z%ope = local_ope
      z%drag = local_drag
      z%t2 = real(df2,dp)*zeta_sq
      z%value = sqrt(max(0.0_dp,zeta_sq)*local_ope)-local_drag
      z%p_value = pt2(z%t2,real(z%df1,dp),real(z%df2,dp),lower_tail=.false.)
      allocate(z%weights(p),z%mean(p),z%covariance(p,p))
      z%weights = mapped_w
      z%mean = mu
      z%covariance = sigma
   end function markowitz_from_moments

   function fit_sropt(x, ope, drag) result(z)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: ope, drag
      type(sropt_result) :: z
      z = markowitz_from_moments(column_mean(x),sample_covariance(x),size(x,1),ope,drag)
   end function fit_sropt

   function make_sropt(value, df1, df2, ope, drag) result(z)
      real(dp), intent(in) :: value
      integer, intent(in) :: df1, df2
      real(dp), intent(in), optional :: ope, drag
      type(sropt_result) :: z
      real(dp) :: local_ope, local_drag, native
      local_ope = 1.0_dp
      local_drag = 0.0_dp
      if (present(ope)) local_ope = ope
      if (present(drag)) local_drag = drag
      z%value = value
      z%df1 = df1
      z%df2 = df2
      z%ope = local_ope
      z%drag = local_drag
      native = (value+local_drag)/sqrt(local_ope)
      z%t2 = real(df2,dp)*native*native
      z%p_value = pt2(z%t2,real(df1,dp),real(df2,dp),lower_tail=.false.)
   end function make_sropt

   function reannualize_sropt(z, new_ope) result(out)
      type(sropt_result), intent(in) :: z
      real(dp), intent(in) :: new_ope
      type(sropt_result) :: out
      out = z
      out%value = z%value*sqrt(new_ope/z%ope)
      out%ope = new_ope
   end function reannualize_sropt

   pure function sric(z) result(value)
      type(sropt_result), intent(in) :: z
      real(dp) :: value, native
      native = z%value/sqrt(z%ope)
      if (abs(native) <= tiny(1.0_dp)) then
         value = -huge(1.0_dp)
      else
         value = sqrt(z%ope)*(native-real(z%df1-1,dp)/(native*real(z%df2,dp)))
      end if
   end function sric

   function make_del_sropt(value, sub_value, df1, df2, df1_sub, ope, drag) result(z)
      real(dp), intent(in) :: value, sub_value
      integer, intent(in) :: df1, df2, df1_sub
      real(dp), intent(in), optional :: ope, drag
      type(del_sropt_result) :: z
      real(dp) :: local_ope, local_drag
      local_ope = 1.0_dp
      local_drag = 0.0_dp
      if (present(ope)) local_ope = ope
      if (present(drag)) local_drag = drag
      z%value = value
      z%sub_value = sub_value
      z%delta_value = sqrt(max(0.0_dp,value*value-sub_value*sub_value))
      z%df1 = df1
      z%df2 = df2
      z%df1_sub = df1_sub
      z%ope = local_ope
      z%drag = local_drag
      z%t2 = real(df2,dp)*((value+local_drag)/sqrt(local_ope))**2
      z%t2_sub = real(df2,dp)*((sub_value+local_drag)/sqrt(local_ope))**2
      z%t2_delta = z%t2-z%t2_sub
   end function make_del_sropt

   function fit_del_sropt(x, h, ope, drag) result(z)
      real(dp), intent(in) :: x(:, :), h(:, :)
      real(dp), intent(in), optional :: ope, drag
      type(del_sropt_result) :: z
      type(sropt_result) :: full, sub
      real(dp) :: local_ope, local_drag
      local_ope = 1.0_dp
      local_drag = 0.0_dp
      if (present(ope)) local_ope = ope
      if (present(drag)) local_drag = drag
      full = fit_sropt(x,local_ope,local_drag)
      sub = markowitz_from_moments(column_mean(x),sample_covariance(x),size(x,1), &
                                   local_ope,local_drag,h)
      z = make_del_sropt(full%value,sub%value,full%df1,full%df2,sub%df1,local_ope,local_drag)
   end function fit_del_sropt

   function sm_vcov(x, normal_model, fit_intercept) result(out)
      real(dp), intent(in) :: x(:, :)
      logical, intent(in), optional :: normal_model, fit_intercept
      type(moment_vcov_result) :: out
      logical :: use_normal, use_intercept
      integer :: n, p, q, i, j, k, a, b, d
      real(dp), allocatable :: y(:, :), sigma(:, :)
      real(dp) :: mu(size(x,2))
      n = size(x,1)
      p = size(x,2)
      use_normal = .false.
      use_intercept = .true.
      if (present(normal_model)) use_normal = normal_model
      if (present(fit_intercept)) use_intercept = fit_intercept
      q = p*(p+1)/2
      if (use_intercept) q = q+p
      out%n = n
      out%p = p
      allocate(out%mean(q),out%covariance(q,q),y(n,q))
      k = 0
      if (use_intercept) then
         y(:,1:p) = x
         k = p
      end if
      do i = 1, p
         do j = i, p
            k = k+1
            y(:,k) = x(:,i)*x(:,j)
         end do
      end do
      out%mean = column_mean(y)
      if (.not. use_normal) then
         out%covariance = covariance_of_mean(y)
         return
      end if
      mu = column_mean(x)
      sigma = sample_covariance(x)
      out%covariance = 0.0_dp
      if (use_intercept) out%covariance(1:p,1:p) = sigma/real(n,dp)
      k = merge(p,0,use_intercept)
      do i = 1, p
         do j = i, p
            k = k+1
            if (use_intercept) then
               do a = 1, p
                  out%covariance(a,k) = (mu(i)*sigma(a,j)+mu(j)*sigma(a,i))/real(n,dp)
                  out%covariance(k,a) = out%covariance(a,k)
               end do
            end if
            d = merge(p,0,use_intercept)
            do a = 1, p
               do b = a, p
                  d = d+1
                  out%covariance(k,d) = (mu(i)*mu(a)*sigma(j,b)+mu(i)*mu(b)*sigma(j,a) + &
                     mu(j)*mu(a)*sigma(i,b)+mu(j)*mu(b)*sigma(i,a) + &
                     sigma(i,a)*sigma(j,b)+sigma(i,b)*sigma(j,a))/real(n,dp)
                  out%covariance(d,k) = out%covariance(k,d)
               end do
            end do
         end do
      end do
   end function sm_vcov

   function inverse_unified_vector(theta_values, fit_intercept, status) result(v)
      real(dp), intent(in) :: theta_values(:)
      logical, intent(in) :: fit_intercept
      integer, intent(out) :: status
      real(dp), allocatable :: v(:), packed(:), theta(:, :), itheta(:, :)
      if (fit_intercept) then
         allocate(packed(size(theta_values)+1))
         packed(1) = 1.0_dp
         packed(2:) = theta_values
      else
         allocate(packed(size(theta_values)))
         packed = theta_values
      end if
      theta = symmetric_ivech(packed,status)
      if (status /= 0) then
         allocate(v(0))
         return
      end if
      allocate(itheta(size(theta,1),size(theta,2)))
      call invert_matrix(theta,itheta,status)
      if (status /= 0) then
         allocate(v(0))
         return
      end if
      packed = symmetric_vech(itheta)
      if (fit_intercept) then
         allocate(v(size(packed)-1))
         v = packed(2:)
      else
         allocate(v(size(packed)))
         v = packed
      end if
   end function inverse_unified_vector

   function ism_vcov(x, normal_model, fit_intercept) result(out)
      real(dp), intent(in) :: x(:, :)
      logical, intent(in), optional :: normal_model, fit_intercept
      type(moment_vcov_result) :: out
      type(moment_vcov_result) :: sm
      logical :: use_intercept
      real(dp), allocatable :: base(:), plus(:), minus(:), jac(:, :), shifted(:)
      real(dp) :: hstep, sym_value
      integer :: i, j, q, status
      use_intercept = .true.
      if (present(fit_intercept)) use_intercept = fit_intercept
      sm = sm_vcov(x,normal_model,use_intercept)
      base = inverse_unified_vector(sm%mean,use_intercept,status)
      out%n = sm%n
      out%p = sm%p
      out%status = status
      if (status /= 0) return
      q = size(sm%mean)
      allocate(out%mean(size(base)),out%covariance(size(base),size(base)))
      allocate(jac(size(base),q),shifted(q))
      out%mean = base
      do i = 1, q
         shifted = sm%mean
         hstep = sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(sm%mean(i)))
         shifted(i) = shifted(i)+hstep
         plus = inverse_unified_vector(shifted,use_intercept,status)
         shifted(i) = sm%mean(i)-hstep
         minus = inverse_unified_vector(shifted,use_intercept,status)
         if (status /= 0) then
            out%status = status
            return
         end if
         jac(:,i) = (plus-minus)/(2.0_dp*hstep)
      end do
      out%covariance = matmul(jac,matmul(sm%covariance,transpose(jac)))

      ! The delta-method covariance is symmetric analytically, but separate
      ! matrix products can leave platform-dependent roundoff differences in
      ! opposite triangles. Copy a common average to both entries so callers
      ! receive an exactly symmetric covariance matrix on every compiler.
      do i = 1, size(out%covariance,1)
         do j = i+1, size(out%covariance,2)
            sym_value = 0.5_dp*(out%covariance(i,j)+out%covariance(j,i))
            out%covariance(i,j) = sym_value
            out%covariance(j,i) = sym_value
         end do
      end do
   end function ism_vcov

end module sharper_estimation
