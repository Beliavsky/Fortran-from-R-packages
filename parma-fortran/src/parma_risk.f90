! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_risk
   use parma_kinds, only: dp, parma_small, parma_big, pi
   use parma_types, only: risk_mad, risk_ev, risk_minimax, risk_cvar, risk_cdar, &
      risk_lpm, risk_upm, risk_rachev
   implicit none
   private
   public :: sentropy, centropy, smooth_max, smooth_max_gradient
   public :: smooth_abs, smooth_abs_gradient, smooth_abs2, smooth_abs2_gradient
   public :: portfolio_returns, mad_risk, variance_risk, quadratic_variance
   public :: minimax_risk, lpm_risk, upm_risk, cvar_risk, cdar_risk
   public :: cumulative_return, drawdown_series, rachev_ratio, risk_value, riskfun
   public :: benchmark_variance
   public :: scenario_mad, scenario_variance, scenario_lpm, scenario_cvar
   public :: empirical_quantile, sort_real

contains

   function sentropy(w, info) result(value)
      real(dp), intent(in) :: w(:)
      integer, intent(out), optional :: info
      real(dp) :: value
      real(dp), allocatable :: x(:)
      integer :: ierr

      ierr = 0
      allocate(x(size(w)))
      x = w
      where (x < 0.0_dp .and. abs(x) < 1.0e-4_dp) x = abs(x)
      if (any(x < 0.0_dp)) then
         ierr = 1
         value = huge(1.0_dp)
      else
         value = -sum(merge(x * log(max(x,tiny(1.0_dp))), 0.0_dp, x > 0.0_dp))
      end if
      if (present(info)) info = ierr
   end function sentropy

   function centropy(w, q, info) result(value)
      real(dp), intent(in) :: w(:), q(:)
      integer, intent(out), optional :: info
      real(dp) :: value
      real(dp), allocatable :: xq(:)
      integer :: ierr

      ierr = 0
      if (size(w) /= size(q)) then
         ierr = -1
         value = huge(1.0_dp)
      else
         allocate(xq(size(q)))
         xq = q
         where (xq < 0.0_dp .and. abs(xq) < 1.0e-4_dp) xq = abs(xq)
         if (any(xq <= 0.0_dp)) then
            ierr = 1
            value = huge(1.0_dp)
         else
            value = sum((w - xq)**2 / xq)
         end if
      end if
      if (present(info)) info = ierr
   end function centropy

   elemental function smooth_max(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = 0.5_dp * (sqrt(x*x + parma_small) + x)
   end function smooth_max

   elemental function smooth_max_gradient(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = 0.5_dp + 0.5_dp*x/sqrt(x*x + parma_small)
   end function smooth_max_gradient

   elemental function smooth_abs(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = sqrt((x + parma_small)*(x + parma_small))
   end function smooth_abs

   elemental function smooth_abs_gradient(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = (x + parma_small) / sqrt((x + parma_small)*(x + parma_small))
   end function smooth_abs_gradient

   elemental function smooth_abs2(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = (2.0_dp*x/pi) * atan(parma_big*x)
   end function smooth_abs2

   elemental function smooth_abs2_gradient(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = (2.0_dp/pi)*atan(parma_big*x) + &
         (2.0_dp*x/pi)*parma_big/(1.0_dp + (parma_big*x)**2)
   end function smooth_abs2_gradient

   function portfolio_returns(weights, data, benchmark) result(r)
      real(dp), intent(in) :: weights(:), data(:,:)
      real(dp), intent(in), optional :: benchmark(:)
      real(dp) :: r(size(data,1))

      r = matmul(data, weights)
      if (present(benchmark)) r = r - benchmark
   end function portfolio_returns

   function mad_risk(weights, data, benchmark) result(value)
      real(dp), intent(in) :: weights(:), data(:,:)
      real(dp), intent(in), optional :: benchmark(:)
      real(dp) :: value
      real(dp) :: r(size(data,1))

      r = portfolio_returns(weights, data)
      if (present(benchmark)) then
         r = r - benchmark
      else
         r = r - sum(r)/real(size(r),dp)
      end if
      value = sum(abs(r))/real(size(r),dp)
   end function mad_risk

   function variance_risk(weights, data, benchmark) result(value)
      real(dp), intent(in) :: weights(:), data(:,:)
      real(dp), intent(in), optional :: benchmark(:)
      real(dp) :: value
      real(dp) :: r(size(data,1))

      r = portfolio_returns(weights, data)
      if (present(benchmark)) then
         r = r - benchmark
      else
         r = r - sum(r)/real(size(r),dp)
      end if
      value = sum(r*r)/real(size(r),dp)
   end function variance_risk

   function quadratic_variance(weights, covariance) result(value)
      real(dp), intent(in) :: weights(:), covariance(:,:)
      real(dp) :: value
      value = dot_product(weights, matmul(covariance, weights))
   end function quadratic_variance

   function minimax_risk(weights, data, benchmark) result(value)
      real(dp), intent(in) :: weights(:), data(:,:)
      real(dp), intent(in), optional :: benchmark(:)
      real(dp) :: value
      real(dp) :: r(size(data,1))

      r = portfolio_returns(weights, data)
      if (present(benchmark)) r = r - benchmark
      value = -minval(r)
   end function minimax_risk

   function lpm_risk(weights, data, threshold, moment, legacy) result(value)
      real(dp), intent(in) :: weights(:), data(:,:), threshold, moment
      logical, intent(in), optional :: legacy
      real(dp) :: value, thresholdx
      real(dp) :: r(size(data,1)), p(size(data,1))
      logical :: use_legacy

      use_legacy = .false.
      if (present(legacy)) use_legacy = legacy
      r = portfolio_returns(weights, data)
      thresholdx = threshold
      if (abs(threshold-999.0_dp) <= epsilon(1.0_dp)) then
         r = r-sum(r)/real(size(r),dp)
         thresholdx = 0.0_dp
      end if
      if (use_legacy) then
         p = max(r - thresholdx, 0.0_dp)
      else
         p = max(thresholdx - r, 0.0_dp)
      end if
      if (moment <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = (sum(p**moment)/real(size(p),dp))**(1.0_dp/moment)
      end if
   end function lpm_risk

   function upm_risk(weights, data, threshold, moment) result(value)
      real(dp), intent(in) :: weights(:), data(:,:), threshold, moment
      real(dp) :: value, thresholdx
      real(dp) :: r(size(data,1)), p(size(data,1))

      r = portfolio_returns(weights, data)
      thresholdx = threshold
      if (abs(threshold-999.0_dp) <= epsilon(1.0_dp)) then
         r = r-sum(r)/real(size(r),dp)
         thresholdx = 0.0_dp
      end if
      p = max(r - thresholdx, 0.0_dp)
      if (moment <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = (sum(p**moment)/real(size(p),dp))**(1.0_dp/moment)
      end if
   end function upm_risk

   function empirical_quantile(x, probability) result(value)
      real(dp), intent(in) :: x(:), probability
      real(dp) :: value
      real(dp), allocatable :: y(:)
      real(dp) :: pos, frac
      integer :: i, n

      n = size(x)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      allocate(y(n))
      y = x
      call sort_real(y)
      if (probability <= 0.0_dp) then
         value = y(1)
      else if (probability >= 1.0_dp) then
         value = y(n)
      else
         pos = 1.0_dp + probability*real(n-1,dp)
         i = floor(pos)
         frac = pos - real(i,dp)
         if (i >= n) then
            value = y(n)
         else
            value = (1.0_dp-frac)*y(i) + frac*y(i+1)
         end if
      end if
   end function empirical_quantile

   function cvar_risk(weights, data, alpha, benchmark, var_level) result(value)
      real(dp), intent(in) :: weights(:), data(:,:), alpha
      real(dp), intent(in), optional :: benchmark(:)
      real(dp), intent(out), optional :: var_level
      real(dp) :: value, q
      real(dp) :: r(size(data,1))

      r = portfolio_returns(weights, data)
      if (present(benchmark)) r = r - benchmark
      if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
         value = huge(1.0_dp)
         if (present(var_level)) var_level = 0.0_dp
         return
      end if
      q = empirical_quantile(r, alpha)
      value = -q + sum(max(q-r,0.0_dp))/real(size(r),dp)/alpha
      if (present(var_level)) var_level = q
   end function cvar_risk

   function cumulative_return(weights, data, benchmark) result(cumret)
      real(dp), intent(in) :: weights(:), data(:,:)
      real(dp), intent(in), optional :: benchmark(:)
      real(dp) :: cumret(size(data,1))
      real(dp) :: r(size(data,1))
      integer :: i

      r = portfolio_returns(weights, data)
      if (present(benchmark)) r = r - benchmark
      cumret(1) = r(1)
      do i = 2, size(r)
         cumret(i) = cumret(i-1) + r(i)
      end do
   end function cumulative_return

   function drawdown_series(weights, data, lower, benchmark) result(dd)
      real(dp), intent(in) :: weights(:), data(:,:)
      logical, intent(in), optional :: lower
      real(dp), intent(in), optional :: benchmark(:)
      real(dp) :: dd(size(data,1))
      real(dp) :: c(size(data,1)), running
      logical :: is_lower
      integer :: i

      is_lower = .true.
      if (present(lower)) is_lower = lower
      if (present(benchmark)) then
         c = cumulative_return(weights, data, benchmark)
      else
         c = cumulative_return(weights, data)
      end if
      running = c(1)
      dd(1) = 0.0_dp
      do i = 2, size(c)
         if (is_lower) then
            running = max(running,c(i))
            dd(i) = running - c(i)
         else
            running = min(running,c(i))
            dd(i) = c(i) - running
         end if
      end do
   end function drawdown_series

   function cdar_risk(weights, data, alpha, lower, benchmark, dar_level) result(value)
      real(dp), intent(in) :: weights(:), data(:,:), alpha
      logical, intent(in), optional :: lower
      real(dp), intent(in), optional :: benchmark(:)
      real(dp), intent(out), optional :: dar_level
      real(dp) :: value, q
      real(dp) :: dd(size(data,1))

      if (present(benchmark)) then
         if (present(lower)) then
            dd = drawdown_series(weights,data,lower,benchmark)
         else
            dd = drawdown_series(weights,data,benchmark=benchmark)
         end if
      else if (present(lower)) then
         dd = drawdown_series(weights,data,lower)
      else
         dd = drawdown_series(weights,data)
      end if
      if (alpha <= 0.0_dp .or. alpha > 1.0_dp) then
         value = huge(1.0_dp)
         if (present(dar_level)) dar_level = 0.0_dp
         return
      end if
      q = empirical_quantile(dd, 1.0_dp-alpha)
      value = q + sum(max(dd-q,0.0_dp))/real(size(dd),dp)/alpha
      if (present(dar_level)) dar_level = q
   end function cdar_risk

   function rachev_ratio(weights, data, alpha_down, alpha_up) result(value)
      real(dp), intent(in) :: weights(:), data(:,:), alpha_down, alpha_up
      real(dp) :: value, downside, upside

      downside = cvar_risk(weights,data,alpha_down)
      upside = cvar_risk(weights,-data,alpha_up)
      if (abs(upside) <= tiny(1.0_dp)) then
         value = huge(1.0_dp)
      else
         value = downside/upside
      end if
   end function rachev_ratio

   function risk_value(weights, data, risk, alpha, moment, threshold, benchmark, &
      covariance, lpm_legacy, var_level, dar_level) result(value)
      real(dp), intent(in) :: weights(:)
      real(dp), intent(in), optional :: data(:,:)
      integer, intent(in) :: risk
      real(dp), intent(in), optional :: alpha, moment, threshold
      real(dp), intent(in), optional :: benchmark(:), covariance(:,:)
      logical, intent(in), optional :: lpm_legacy
      real(dp), intent(out), optional :: var_level, dar_level
      real(dp) :: value, a, p, t, v, d
      logical :: legacy

      a = 0.05_dp
      p = 1.0_dp
      t = 0.0_dp
      legacy = .false.
      if (present(alpha)) a = alpha
      if (present(moment)) p = moment
      if (present(threshold)) t = threshold
      if (present(lpm_legacy)) legacy = lpm_legacy
      v = 0.0_dp
      d = 0.0_dp
      select case (risk)
      case (risk_mad)
         if (.not. present(data)) then
            value = huge(1.0_dp)
         else if (present(benchmark)) then
            value = mad_risk(weights,data,benchmark)
         else
            value = mad_risk(weights,data)
         end if
      case (risk_ev)
         if (present(covariance)) then
            value = quadratic_variance(weights,covariance)
         else if (.not. present(data)) then
            value = huge(1.0_dp)
         else if (present(benchmark)) then
            value = variance_risk(weights,data,benchmark)
         else
            value = variance_risk(weights,data)
         end if
      case (risk_minimax)
         if (.not. present(data)) then
            value = huge(1.0_dp)
         else if (present(benchmark)) then
            value = minimax_risk(weights,data,benchmark)
         else
            value = minimax_risk(weights,data)
         end if
      case (risk_cvar)
         if (.not. present(data)) then
            value = huge(1.0_dp)
         else if (present(benchmark)) then
            value = cvar_risk(weights,data,a,benchmark,v)
         else
            value = cvar_risk(weights,data,a,var_level=v)
         end if
      case (risk_cdar)
         if (.not. present(data)) then
            value = huge(1.0_dp)
         else if (present(benchmark)) then
            value = cdar_risk(weights,data,a,benchmark=benchmark,dar_level=d)
         else
            value = cdar_risk(weights,data,a,dar_level=d)
         end if
      case (risk_lpm)
         if (present(data)) then
            value = lpm_risk(weights,data,t,p,legacy)
         else
            value = huge(1.0_dp)
         end if
      case (risk_upm)
         if (present(data)) then
            value = upm_risk(weights,data,t,p)
         else
            value = huge(1.0_dp)
         end if
      case (risk_rachev)
         if (present(data)) then
            value = rachev_ratio(weights,data,a,a)
         else
            value = huge(1.0_dp)
         end if
      case default
         value = huge(1.0_dp)
      end select
      value = abs(value)
      if (present(var_level)) var_level = v
      if (present(dar_level)) dar_level = d
   end function risk_value

   function scenario_mad(data) result(values)
      real(dp), intent(in) :: data(:,:)
      real(dp) :: values(size(data,1),size(data,2))
      real(dp) :: mu(size(data,2))
      mu = sum(data,dim=1)/real(size(data,1),dp)
      values = abs(data-spread(mu,1,size(data,1)))
   end function scenario_mad

   function scenario_variance(data) result(values)
      real(dp), intent(in) :: data(:,:)
      real(dp) :: values(size(data,1),size(data,2))
      real(dp) :: mu(size(data,2))
      mu = sum(data,dim=1)/real(size(data,1),dp)
      values = (data-spread(mu,1,size(data,1)))**2
   end function scenario_variance

   function scenario_lpm(data, threshold, moment, legacy) result(values)
      real(dp), intent(in) :: data(:,:), threshold, moment
      logical, intent(in), optional :: legacy
      real(dp) :: values(size(data,1),size(data,2))
      real(dp) :: work(size(data,1),size(data,2)), thresholdx
      logical :: use_legacy
      use_legacy = .false.
      if (present(legacy)) use_legacy = legacy
      work = data
      thresholdx = threshold
      if (abs(threshold-999.0_dp) <= epsilon(1.0_dp)) then
         work = work-spread(sum(work,dim=1)/real(size(work,1),dp),1,size(work,1))
         thresholdx = 0.0_dp
      end if
      if (use_legacy) then
         values = max(work-thresholdx,0.0_dp)**moment
      else
         values = max(thresholdx-work,0.0_dp)**moment
      end if
   end function scenario_lpm

   function scenario_cvar(data, alpha) result(values)
      real(dp), intent(in) :: data(:,:), alpha
      real(dp) :: values(size(data,1),size(data,2))
      real(dp) :: q
      integer :: j
      do j = 1, size(data,2)
         q = empirical_quantile(data(:,j),alpha)
         values(:,j) = -q + smooth_max(q-data(:,j))
      end do
   end function scenario_cvar

   function benchmark_variance(weights,covariance,benchmark_covariance) result(value)
      real(dp), intent(in) :: weights(:),covariance(:,:),benchmark_covariance(:)
      real(dp) :: value
      if (size(benchmark_covariance) /= size(weights)+1) then
         value = huge(1.0_dp)
      else
         value = benchmark_covariance(1)+quadratic_variance(weights,covariance)- &
            2.0_dp*dot_product(weights,benchmark_covariance(2:))
      end if
   end function benchmark_variance

   function riskfun(weights,data,risk,benchmark,alpha,moment,threshold,var_level,dar_level,legacy_lpm) result(value)
      real(dp), intent(in) :: weights(:),data(:,:)
      character(len=*), intent(in) :: risk
      real(dp), intent(in), optional :: benchmark(:),alpha,moment,threshold
      real(dp), intent(out), optional :: var_level,dar_level
      logical, intent(in), optional :: legacy_lpm
      real(dp) :: value
      integer :: code
      character(len=:), allocatable :: name

      name = lowercase(trim(adjustl(risk)))
      select case (name)
      case ('mad')
         code = risk_mad
      case ('ev','variance')
         code = risk_ev
      case ('minimax','minmax')
         code = risk_minimax
      case ('cvar')
         code = risk_cvar
      case ('cdar')
         code = risk_cdar
      case ('lpm')
         code = risk_lpm
      case ('upm')
         code = risk_upm
      case ('rachev')
         code = risk_rachev
      case default
         value = huge(1.0_dp)
         if (present(var_level)) var_level = 0.0_dp
         if (present(dar_level)) dar_level = 0.0_dp
         return
      end select
      if (present(benchmark)) then
         value = risk_value(weights,data,code,alpha,moment,threshold,benchmark=benchmark, &
            lpm_legacy=legacy_lpm,var_level=var_level,dar_level=dar_level)
      else
         value = risk_value(weights,data,code,alpha,moment,threshold,lpm_legacy=legacy_lpm, &
            var_level=var_level,dar_level=dar_level)
      end if
   end function riskfun

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i,code
      do i = 1,len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            lower(i:i) = achar(code+32)
         else
            lower(i:i) = text(i:i)
         end if
      end do
   end function lowercase

   recursive subroutine quicksort(a, left, right)
      real(dp), intent(inout) :: a(:)
      integer, intent(in) :: left, right
      integer :: i, j
      real(dp) :: pivot, temp
      if (left >= right) return
      i = left
      j = right
      pivot = a((left+right)/2)
      do
         do while (a(i) < pivot)
            i = i + 1
         end do
         do while (a(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            temp = a(i)
            a(i) = a(j)
            a(j) = temp
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (left < j) call quicksort(a,left,j)
      if (i < right) call quicksort(a,i,right)
   end subroutine quicksort

   subroutine sort_real(a)
      real(dp), intent(inout) :: a(:)
      if (size(a) > 1) call quicksort(a,1,size(a))
   end subroutine sort_real

end module parma_risk
