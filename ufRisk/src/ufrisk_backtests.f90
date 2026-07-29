! SPDX-License-Identifier: GPL-3.0-only
module ufrisk_backtests
   use kind_mod, only : dp
   use ufrisk_types
   use ufrisk_math, only : normal_cdf, student_cdf, binomial_cdf, chi_square_cdf, &
      log_probability
   implicit none
   private
   public :: lossfunc, covtest, trafftest
contains
   pure function lossfunc(loss, expected_shortfall, beta) result(out)
      real(dp), intent(in) :: loss(:), expected_shortfall(:)
      real(dp), intent(in), optional :: beta
      type(ufrisk_loss_result) :: out
      real(dp) :: opportunity, difference
      integer :: i
      opportunity = 1.0e-4_dp
      if (present(beta)) opportunity = beta
      if (size(loss) <= 1 .or. size(expected_shortfall) /= size(loss)) then
         out%status = ufrisk_invalid_input
         return
      end if
      do i = 1, size(loss)
         difference = loss(i)-expected_shortfall(i)
         if (loss(i) > expected_shortfall(i)) then
            out%regulatory = out%regulatory+difference*difference
            out%firm = out%firm+difference*difference
            out%abad = out%abad+difference*difference
            out%feng = out%feng+difference*difference
         else
            out%firm = out%firm+opportunity*expected_shortfall(i)
            out%abad = out%abad+opportunity*abs(difference)
            if (loss(i) >= 0.0_dp) then
               out%feng = out%feng+opportunity*abs(difference)
            else
               out%feng = out%feng+opportunity*expected_shortfall(i)
            end if
         end if
      end do
      out%regulatory = 10000.0_dp*out%regulatory
      out%firm = 10000.0_dp*out%firm
      out%abad = 10000.0_dp*out%abad
      out%feng = 10000.0_dp*out%feng
   end function lossfunc

   pure function covtest(loss, value_at_risk, tail_probability) result(out)
      real(dp), intent(in) :: loss(:), value_at_risk(:), tail_probability
      type(ufrisk_coverage_result) :: out
      logical, allocatable :: hit(:)
      integer :: i, n, n0, n1
      real(dp) :: phat, p01, p11, log_null, log_unconditional, log_markov
      n = size(loss)
      out%tail_probability = tail_probability
      if (n <= 1 .or. size(value_at_risk) /= n .or. tail_probability <= 0.0_dp .or. &
         tail_probability >= 1.0_dp) then
         out%status = ufrisk_invalid_input
         return
      end if
      allocate(hit(n)); hit = loss > value_at_risk
      n1 = count(hit); n0 = n-n1
      if (n1 == 0) then
         out%status = ufrisk_no_violations
         return
      end if
      do i = 2, n
         if (.not.hit(i-1) .and. .not.hit(i)) out%n00 = out%n00+1
         if (.not.hit(i-1) .and. hit(i)) out%n01 = out%n01+1
         if (hit(i-1) .and. .not.hit(i)) out%n10 = out%n10+1
         if (hit(i-1) .and. hit(i)) out%n11 = out%n11+1
      end do
      phat = real(n1,dp)/real(n,dp)
      p01 = 0.0_dp
      if (out%n00+out%n01 > 0) p01 = real(out%n01,dp)/real(out%n00+out%n01,dp)
      p11 = 0.0_dp
      if (out%n10+out%n11 > 0) p11 = real(out%n11,dp)/real(out%n10+out%n11,dp)
      log_null = log_probability(1.0_dp-tail_probability,n0) + &
         log_probability(tail_probability,n1)
      log_unconditional = log_probability(1.0_dp-phat,n0) + log_probability(phat,n1)
      log_markov = log_probability(1.0_dp-p01,out%n00) + log_probability(p01,out%n01) + &
         log_probability(1.0_dp-p11,out%n10) + log_probability(p11,out%n11)
      out%lr_unconditional = max(0.0_dp,-2.0_dp*(log_null-log_unconditional))
      out%lr_independence = max(0.0_dp,-2.0_dp*(log_unconditional-log_markov))
      out%lr_conditional = max(0.0_dp,-2.0_dp*(log_null-log_markov))
      out%p_unconditional = 1.0_dp-chi_square_cdf(out%lr_unconditional,1)
      out%p_independence = 1.0_dp-chi_square_cdf(out%lr_independence,1)
      out%p_conditional = 1.0_dp-chi_square_cdf(out%lr_conditional,2)
   end function covtest

   pure function trafftest(result) result(out)
      type(ufrisk_result), intent(in) :: result
      type(ufrisk_traffic_result) :: out
      real(dp), allocatable :: loss(:)
      real(dp) :: sdev, standardized_loss, breach, mean_e, mean_v, mean_b, sd_b
      integer :: i, n
      if (.not.allocated(result%returns_out) .or. .not.allocated(result%sigma_forecast) .or. &
         .not.allocated(result%var_es_level) .or. .not.allocated(result%var_var_level) .or. &
         .not.allocated(result%expected_shortfall)) then
         out%status = ufrisk_invalid_input
         return
      end if
      n = size(result%returns_out)
      if (n < 1) then
         out%status = ufrisk_invalid_input
         return
      end if
      allocate(loss(n)); loss = -result%returns_out
      sdev = 1.0_dp
      if (result%distribution == ufrisk_distribution_student) &
         sdev = sqrt(result%degrees_freedom/(result%degrees_freedom-2.0_dp))
      do i = 1, n
         if (result%var_es_level(i) < loss(i)) then
            out%violations_es_var = out%violations_es_var+1
            standardized_loss = -(result%returns_out(i)-result%mean_return)/ &
               result%sigma_forecast(i)*sdev
            if (result%distribution == ufrisk_distribution_student) then
               breach = 1.0_dp-(1.0_dp-student_cdf(standardized_loss, &
                  result%degrees_freedom))/result%es_tail_probability
            else
               breach = 1.0_dp-(1.0_dp-normal_cdf(standardized_loss))/ &
                  result%es_tail_probability
            end if
            out%breach_sum = out%breach_sum+breach
         end if
         if (result%var_var_level(i) < loss(i)) out%violations_var = out%violations_var+1
         if (result%expected_shortfall(i) < loss(i)) out%violations_es = out%violations_es+1
      end do
      out%p_es_var = binomial_cdf(out%violations_es_var,n,result%es_tail_probability)
      out%p_var = binomial_cdf(out%violations_var,n,result%var_tail_probability)
      mean_e = real(n,dp)*result%es_tail_probability
      mean_v = real(n,dp)*result%var_tail_probability
      mean_b = 0.5_dp*result%es_tail_probability*real(n,dp)
      sd_b = sqrt(result%es_tail_probability*(4.0_dp-3.0_dp*result%es_tail_probability)/ &
         12.0_dp*real(n,dp))
      if (sd_b > 0.0_dp) out%p_es = normal_cdf((out%breach_sum-mean_b)/sd_b)
      if (mean_e > 0.0_dp .and. mean_v > 0.0_dp .and. mean_b > 0.0_dp) then
         out%weighted_absolute_deviation = abs(real(out%violations_es_var,dp)-mean_e)/mean_e + &
            abs(real(out%violations_var,dp)-mean_v)/mean_v + &
            abs(out%breach_sum-mean_b)/mean_b
      end if
   end function trafftest
end module ufrisk_backtests
