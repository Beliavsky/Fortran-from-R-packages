module psych_reliability
   use r_kinds, only : dp
      use psych_stats, only : covariance_matrix, correlation_matrix, trace_matrix
   use psych_correlations, only : squared_multiple_correlations
   use psych_types, only : alpha_result
   implicit none
   private
   public :: cronbach_alpha_from_matrices, cronbach_alpha_data, guttman_lambdas
   public :: split_half_extremes, spearman_brown, reliability_summary

contains

   subroutine cronbach_alpha_from_matrices(covariance, correlation, result, info)
      real(dp), intent(in) :: covariance(:, :) !! Item covariance matrix.
      real(dp), intent(in) :: correlation(:, :) !! Item correlation matrix.
      type(alpha_result), intent(out) :: result !! Raw/standardized alpha, G6, mean correlation, and signal/noise ratio.
      integer, intent(out) :: info !! Zero when squared multiple correlations are computed successfully.
      real(dp), allocatable :: smc(:)
      real(dp) :: sum_r
      integer :: p
      p = size(correlation,1)
      allocate(smc(p))
      if (p <= 1) then
         result = alpha_result()
         info = 0
         return
      end if
      if (abs(sum(covariance)) > tiny(1.0_dp)) then
         result%raw_alpha = (1.0_dp-trace_matrix(covariance)/sum(covariance))*real(p,dp)/real(p-1,dp)
      else
         result%raw_alpha = 0.0_dp
      end if
      sum_r = sum(correlation)
      if (abs(sum_r) > tiny(1.0_dp)) then
         result%standardized_alpha = (1.0_dp-real(p,dp)/sum_r)*real(p,dp)/real(p-1,dp)
      else
         result%standardized_alpha = 0.0_dp
      end if
      result%average_r = (sum_r-real(p,dp))/real(p*(p-1),dp)
      if (abs(1.0_dp - result%average_r) > tiny(1.0_dp)) then
         result%signal_noise = real(p,dp)*result%average_r/(1.0_dp-result%average_r)
      else
         result%signal_noise = huge(1.0_dp)
      end if
      call squared_multiple_correlations(correlation,smc,info)
      if (abs(sum_r) > tiny(1.0_dp)) then
         result%g6 = 1.0_dp-(real(p,dp)-sum(smc))/sum_r
      else
         result%g6 = 0.0_dp
      end if
   end subroutine cronbach_alpha_from_matrices

   subroutine cronbach_alpha_data(x, result, info)
      real(dp), intent(in) :: x(:, :) !! Observation-by-item data matrix.
      type(alpha_result), intent(out) :: result !! Raw/standardized alpha and associated reliability summaries.
      integer, intent(out) :: info !! Zero on successful matrix calculations.
      real(dp), allocatable :: covariance(:, :), correlation(:, :)
      integer :: p
      p = size(x,2)
      allocate(covariance(p,p),correlation(p,p))
      call covariance_matrix(x,covariance)
      call correlation_matrix(x,correlation)
      call cronbach_alpha_from_matrices(covariance,correlation,result,info)
   end subroutine cronbach_alpha_data

   subroutine guttman_lambdas(r, lambdas, info)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix.
      real(dp), intent(out) :: lambdas(6) !! Guttman lambda 1, 2, 3, approximate 4/5, and lambda 6.
      integer, intent(out) :: info !! Zero when squared multiple correlations are available.
      real(dp), allocatable :: smc(:)
      real(dp) :: total, offsum, offsq, lambda1, cmax
      integer :: i, p
      p = size(r,1)
      allocate(smc(p))
      total = sum(r)
      if (abs(total) <= tiny(1.0_dp) .or. p <= 1) then
         lambdas = 0.0_dp
         info = 0
         return
      end if
      lambda1 = 1.0_dp-trace_matrix(r)/total
      offsum = total-trace_matrix(r)
      offsq = 0.0_dp
      cmax = 0.0_dp
      do i = 1, p
         offsq = offsq+sum(r(:,i)*r(:,i))-r(i,i)*r(i,i)
         cmax = max(cmax,sum(r(:,i)*r(:,i))-r(i,i)*r(i,i))
      end do
      lambdas(1) = lambda1
      lambdas(2) = (offsum+sqrt(max(0.0_dp,offsq*real(p,dp)/real(p-1,dp))))/total
      lambdas(3) = real(p,dp)*lambda1/real(p-1,dp)
      lambdas(4) = lambdas(3)
      lambdas(5) = lambda1+2.0_dp*sqrt(max(0.0_dp,cmax))/total
      call squared_multiple_correlations(r,smc,info)
      lambdas(6) = (total+sum(smc)-trace_matrix(r))/total
   end subroutine guttman_lambdas

   pure real(dp) function spearman_brown(r_half) result(reliability)
      real(dp), intent(in) :: r_half !! Correlation between two parallel test halves.
      reliability = 2.0_dp*r_half/(1.0_dp+r_half)
   end function spearman_brown

   subroutine split_half_extremes(r, minimum_reliability, maximum_reliability)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix; exhaustive balanced splits are used up to 20 items.
      real(dp), intent(out) :: minimum_reliability !! Smallest Spearman-Brown reliability among examined balanced splits.
      real(dp), intent(out) :: maximum_reliability !! Largest Spearman-Brown reliability among examined balanced splits.
      integer :: p, mask, maxmask, i, n1, n2
      real(dp) :: v1, v2, covariance12, rhalf, rel
      logical, allocatable :: first(:)
      p = size(r,1)
      minimum_reliability = huge(1.0_dp)
      maximum_reliability = -huge(1.0_dp)
      if (p < 2 .or. p > 20) then
         minimum_reliability = 0.0_dp
         maximum_reliability = 0.0_dp
         return
      end if
      allocate(first(p))
      maxmask = 2**(p-1)-1
      do mask = 0, maxmask
         first = .false.
         first(1) = .true.
         do i = 2, p
            first(i) = btest(mask,i-2)
         end do
         n1 = count(first)
         n2 = p-n1
         if (abs(n1-n2) > 1 .or. n1 == 0 .or. n2 == 0) cycle
         call split_variances(r,first,v1,v2,covariance12)
         if (v1 <= 0.0_dp .or. v2 <= 0.0_dp) cycle
         rhalf = covariance12/sqrt(v1*v2)
         rel = spearman_brown(rhalf)
         minimum_reliability = min(minimum_reliability,rel)
         maximum_reliability = max(maximum_reliability,rel)
      end do
      if (minimum_reliability > 0.5_dp * huge(1.0_dp)) minimum_reliability = 0.0_dp
      if (maximum_reliability < -0.5_dp * huge(1.0_dp)) maximum_reliability = 0.0_dp
   end subroutine split_half_extremes

   pure subroutine split_variances(r, first, v1, v2, covariance12)
      real(dp), intent(in) :: r(:, :) !! Item covariance/correlation matrix.
      logical, intent(in) :: first(:) !! True for items in the first half and false for the second half.
      real(dp), intent(out) :: v1 !! Variance of the first-half unit-weighted score.
      real(dp), intent(out) :: v2 !! Variance of the second-half unit-weighted score.
      real(dp), intent(out) :: covariance12 !! Covariance between the two unit-weighted half scores.
      integer :: i, j
      v1 = 0.0_dp
      v2 = 0.0_dp
      covariance12 = 0.0_dp
      do j = 1, size(r,2)
         do i = 1, size(r,1)
            if (first(i) .and. first(j)) v1 = v1+r(i,j)
            if (.not. first(i) .and. .not. first(j)) v2 = v2+r(i,j)
            if (first(i) .and. .not. first(j)) covariance12 = covariance12+r(i,j)
         end do
      end do
   end subroutine split_variances

   subroutine reliability_summary(r, alpha_value, g6_value, min_split, max_split, info)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix.
      real(dp), intent(out) :: alpha_value !! Standardized Cronbach alpha.
      real(dp), intent(out) :: g6_value !! Guttman lambda-6 reliability.
      real(dp), intent(out) :: min_split !! Minimum balanced split-half reliability.
      real(dp), intent(out) :: max_split !! Maximum balanced split-half reliability.
      integer, intent(out) :: info !! Zero when SMC calculations succeed.
      type(alpha_result) :: ar
      real(dp), allocatable :: dummy_cov(:, :)
      allocate(dummy_cov(size(r,1),size(r,2)))
      dummy_cov = r
      call cronbach_alpha_from_matrices(dummy_cov,r,ar,info)
      alpha_value = ar%standardized_alpha
      g6_value = ar%g6
      call split_half_extremes(r,min_split,max_split)
   end subroutine reliability_summary

end module psych_reliability
