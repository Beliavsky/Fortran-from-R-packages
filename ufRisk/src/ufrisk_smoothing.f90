! SPDX-License-Identifier: GPL-3.0-only
module ufrisk_smoothing
   use kind_mod, only : dp
   use ufrisk_types, only : long_memory_smooth_result, ufrisk_ok, &
      ufrisk_invalid_input, ufrisk_smoothing_failed
   use smoots_smoothing, only : local_polynomial_smooth, trend_kernel_constants
   use smoots_stats, only : factorial_real
   use smoots_status, only : sm_ok
   use fracdiff_model_api, only : fracdiff_fit
   use fracdiff_types, only : fracdiff_model
   use fracdiff_status, only : fd_ok, fd_iteration_limit
   implicit none
   private
   public :: long_memory_smooth
contains
   subroutine long_memory_smooth(y, result, p, mu, b_start, iterations, p_min, p_max, &
      q_min, q_max, m_terms)
      real(dp), intent(in) :: y(:)
      type(long_memory_smooth_result), intent(out) :: result
      integer, intent(in), optional :: p, mu, iterations, p_min, p_max, q_min, q_max, m_terms
      real(dp), intent(in), optional :: b_start
      integer :: pp, mm, max_iterations, ar_min, ar_max, ma_min, ma_max, terms
      integer :: n, iteration, k, derivative_degree, boundary, status
      integer :: best_p, best_q
      real(dp) :: bandwidth, old_bandwidth, inflated, d, cf0, rp, muk
      real(dp) :: curvature, c1, c2, exponent, candidate, relative_change
      real(dp), allocatable :: estimate(:), derivative(:), weights(:,:)
      type(fracdiff_model) :: fit

      pp = 3
      if (present(p)) pp = p
      mm = 1
      if (present(mu)) mm = mu
      bandwidth = 0.15_dp
      if (present(b_start)) bandwidth = b_start
      max_iterations = 40
      if (present(iterations)) max_iterations = max(1,iterations)
      ar_min = 0
      if (present(p_min)) ar_min = max(0,p_min)
      ar_max = ar_min
      if (present(p_max)) ar_max = max(ar_min,p_max)
      ma_min = 0
      if (present(q_min)) ma_min = max(0,q_min)
      ma_max = ma_min
      if (present(q_max)) ma_max = max(ma_min,q_max)
      terms = 100
      if (present(m_terms)) terms = max(1,m_terms)
      n = size(y)
      if (n < 20 .or. (pp /= 1 .and. pp /= 3) .or. mm < 0 .or. &
         bandwidth <= 0.0_dp .or. bandwidth >= 0.5_dp) then
         result%status = ufrisk_invalid_input
         return
      end if
      allocate(result%bandwidth_path(max_iterations)); result%bandwidth_path = 0.0_dp
      call trend_kernel_constants(pp,mm,0.05_dp,rp,muk,status)
      if (status /= sm_ok .or. abs(muk) <= tiny(1.0_dp)) then
         result%status = ufrisk_smoothing_failed
         return
      end if
      k = pp+1
      derivative_degree = pp+2
      boundary = int(0.05_dp*real(n,dp))
      bandwidth = min(0.49_dp,max(real(n,dp)**(-0.45_dp),bandwidth))
      old_bandwidth = bandwidth
      do iteration = 1, max_iterations
         call local_polynomial_smooth(y,0,pp,mm,bandwidth,1,estimate,weights,status)
         if (status /= sm_ok) then
            result%status = ufrisk_smoothing_failed
            return
         end if
         call select_fractional_model(y-estimate,ar_min,ar_max,ma_min,ma_max,terms, &
            fit,best_p,best_q,status)
         if (status /= ufrisk_ok) then
            result%status = status
            return
         end if
         d = max(-0.49_dp,min(0.49_dp,fit%d))
         cf0 = fractional_spectrum_zero(fit)
         inflated = min(0.49_dp,max(real(n,dp)**(-0.45_dp), &
            inflation_bandwidth(pp,bandwidth,d)))
         call local_polynomial_smooth(y,k,derivative_degree,mm,inflated,1, &
            derivative,weights,status)
         if (status /= sm_ok) then
            result%status = ufrisk_smoothing_failed
            return
         end if
         if (n-2*boundary > 0) then
            curvature = sum(derivative(boundary+1:n-boundary)**2)/ &
               real(n-2*boundary,dp)
         else
            curvature = sum(derivative**2)/real(n,dp)
         end if
         c1 = factorial_real(k)**2/(2.0_dp*real(k,dp))
         c2 = 0.9_dp*rp/(muk*muk)
         exponent = 1.0_dp/(2.0_dp*real(k,dp)+1.0_dp-2.0_dp*d)
         candidate = (c1*c2*max(2.0_dp*acos(-1.0_dp)*cf0,epsilon(1.0_dp))/ &
            max(curvature,epsilon(1.0_dp)))**exponent * &
            real(n,dp)**(-(1.0_dp-2.0_dp*d)*exponent)
         candidate = min(0.49_dp,max(real(n,dp)**(-0.45_dp),candidate))
         result%bandwidth_path(iteration) = candidate
         relative_change = abs(candidate-bandwidth)/max(candidate,epsilon(1.0_dp))
         old_bandwidth = bandwidth
         bandwidth = candidate
         if (iteration >= 3 .and. relative_change < 1.0_dp/real(n,dp)) exit
         if (iteration >= 4 .and. abs(candidate-old_bandwidth)/ &
            max(candidate,epsilon(1.0_dp)) < 1.0_dp/real(n,dp)) exit
      end do
      call local_polynomial_smooth(y,0,pp,mm,bandwidth,1,result%estimate,weights,status)
      if (status /= sm_ok) then
         result%status = ufrisk_smoothing_failed
         return
      end if
      allocate(result%residuals(n)); result%residuals = y-result%estimate
      result%iterations = min(iteration,max_iterations)
      if (result%iterations < max_iterations) &
         result%bandwidth_path = result%bandwidth_path(1:result%iterations)
      result%bandwidth = bandwidth
      result%d = fit%d
      result%cf0 = cf0
      result%ar_order = best_p
      result%ma_order = best_q
      result%fractional_fit = fit
      result%status = ufrisk_ok
   end subroutine long_memory_smooth

   subroutine select_fractional_model(x,p_min,p_max,q_min,q_max,m_terms,best,best_p,best_q,status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: p_min,p_max,q_min,q_max,m_terms
      type(fracdiff_model), intent(out) :: best
      integer, intent(out) :: best_p,best_q,status
      type(fracdiff_model) :: candidate
      real(dp) :: bic, best_bic
      integer :: p,q,npar
      best_bic = huge(1.0_dp); best_p = p_min; best_q = q_min
      status = ufrisk_smoothing_failed
      do p = p_min,p_max
         do q = q_min,q_max
            candidate = fracdiff_fit(x,nar=p,nma=q,m_terms=m_terms)
            if (candidate%status == fd_ok .or. candidate%status == fd_iteration_limit) then
               npar = p+q+3
               bic = -2.0_dp*candidate%log_likelihood+log(real(size(x),dp))*real(npar,dp)
               if (bic < best_bic) then
                  best_bic = bic
                  best = candidate
                  best_p = p
                  best_q = q
                  status = ufrisk_ok
               end if
            end if
         end do
      end do
   end subroutine select_fractional_model

   pure real(dp) function fractional_spectrum_zero(model) result(value)
      type(fracdiff_model), intent(in) :: model
      real(dp) :: numerator, denominator
      numerator = 1.0_dp
      denominator = 1.0_dp
      if (allocated(model%ma)) numerator = 1.0_dp-sum(model%ma)
      if (allocated(model%ar)) denominator = 1.0_dp-sum(model%ar)
      if (abs(denominator) <= sqrt(epsilon(1.0_dp))) then
         value = huge(1.0_dp)
      else
         value = (numerator/denominator)**2*model%sigma**2/(2.0_dp*acos(-1.0_dp))
      end if
   end function fractional_spectrum_zero

   pure real(dp) function inflation_bandwidth(p,bandwidth,d) result(value)
      integer, intent(in) :: p
      real(dp), intent(in) :: bandwidth,d
      real(dp) :: exponent
      if (p == 1) then
         exponent = (5.0_dp-2.0_dp*d)/(7.0_dp-2.0_dp*d)
      else
         exponent = (9.0_dp-2.0_dp*d)/(11.0_dp-2.0_dp*d)
      end if
      value = bandwidth**exponent
   end function inflation_bandwidth
end module ufrisk_smoothing
