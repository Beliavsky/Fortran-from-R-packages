! SPDX-License-Identifier: GPL-3.0-only
module portvine_stats
   use portvine_kinds, only : dp
   use portvine_types, only : risk_var, risk_es_mean, risk_es_median, risk_es_mc
   implicit none
   private
   public :: empirical_quantile, est_var, est_es, estimate_risk_measure
   public :: sample_mean, sample_sd, sort_real

contains

   pure real(dp) function sample_mean(x) result(value)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = sum(x)/real(size(x),dp)
      end if
   end function sample_mean

   pure real(dp) function sample_sd(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) < 2) then
         value = 0.0_dp
      else
         m = sample_mean(x)
         value = sqrt(max(0.0_dp,sum((x-m)**2)/real(size(x)-1,dp)))
      end if
   end function sample_sd

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   real(dp) function empirical_quantile(x, probability) result(value)
      real(dp), intent(in) :: x(:), probability
      real(dp), allocatable :: work(:)
      real(dp) :: h, frac, p
      integer :: lo, hi, n

      n = size(x)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      allocate(work(n))
      work = x
      call sort_real(work)
      p = max(0.0_dp,min(1.0_dp,probability))
      ! R's default type-7 quantile: h = 1 + (n-1)p.
      h = 1.0_dp + real(n-1,dp)*p
      lo = max(1,min(n,int(floor(h))))
      hi = max(1,min(n,lo+1))
      frac = h-real(lo,dp)
      value = (1.0_dp-frac)*work(lo)+frac*work(hi)
   end function empirical_quantile

   subroutine est_var(sample, alpha, value)
      real(dp), intent(in) :: sample(:), alpha(:)
      real(dp), intent(out) :: value(size(alpha))
      integer :: i
      do i = 1, size(alpha)
         value(i) = empirical_quantile(sample,alpha(i))
      end do
   end subroutine est_var

   subroutine est_es(sample, alpha, value, method, mc_samples)
      real(dp), intent(in) :: sample(:), alpha(:)
      real(dp), intent(out) :: value(size(alpha))
      integer, intent(in), optional :: method, mc_samples
      integer :: i, j, m, n_tail, imethod
      real(dp) :: q, u, accum
      real(dp), allocatable :: tail(:)

      imethod = risk_es_mean
      if (present(method)) imethod = method
      m = 100
      if (present(mc_samples)) m = max(1,mc_samples)
      do i = 1, size(alpha)
         q = empirical_quantile(sample,alpha(i))
         if (imethod == risk_es_mc) then
            accum = 0.0_dp
            do j = 1, m
               call random_number(u)
               accum = accum+empirical_quantile(sample,u*alpha(i))
            end do
            value(i) = accum/real(m,dp)
         else
            n_tail = count(sample <= q)
            if (n_tail <= 0) then
               value(i) = q
            else
               allocate(tail(n_tail))
               tail = pack(sample,sample <= q)
               if (imethod == risk_es_median) then
                  value(i) = empirical_quantile(tail,0.5_dp)
               else
                  value(i) = sum(tail)/real(n_tail,dp)
               end if
               deallocate(tail)
            end if
         end if
      end do
   end subroutine est_es

   subroutine estimate_risk_measure(sample, alpha, measure, value, mc_samples)
      real(dp), intent(in) :: sample(:), alpha(:)
      integer, intent(in) :: measure
      real(dp), intent(out) :: value(size(alpha))
      integer, intent(in), optional :: mc_samples
      select case (measure)
      case (risk_var)
         call est_var(sample,alpha,value)
      case (risk_es_mean, risk_es_median, risk_es_mc)
         if (present(mc_samples)) then
            call est_es(sample,alpha,value,measure,mc_samples)
         else
            call est_es(sample,alpha,value,measure)
         end if
      case default
         value = huge(1.0_dp)
      end select
   end subroutine estimate_risk_measure

end module portvine_stats
