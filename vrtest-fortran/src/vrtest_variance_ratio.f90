! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
!
! Variance-ratio and rank/sign tests translated from vrtest 1.2.
module vrtest_variance_ratio
   use vrtest_kinds, only : dp, pi
   use vrtest_types
   use vrtest_utils, only : mean_value, variance_value, autocorrelation, &
      normal_quantile, chi_square_quantile, sample_quantile, average_ranks, &
      solve_linear, wild_weights, random_normal_vector, random_permutation
   implicit none
   private

   public :: ar1_fit, adjust_thin, abel_bandwidth, quadratic_spectral_kernel
   public :: fast_variance_ratio, lm_statistic, lmcd_statistics
   public :: automatic_variance_ratio, automatic_vr_bootstrap
   public :: lo_mackinlay, chow_denning, variance_ratio_bootstrap
   public :: wald_test, wright_tests, joint_wright_tests
   public :: wright_critical_values, joint_wright_critical_values
   public :: subsample_variance_ratio, panel_variance_ratio
   public :: variance_ratio_minus_one, variance_ratio_curve

contains

   pure function ar1_fit(x, legacy_indexing) result(fit)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: legacy_indexing
      type(ar1_result) :: fit
      integer :: n, last_pair
      real(dp) :: den
      logical :: legacy
      real(dp), allocatable :: residual(:)

      n = size(x)
      if (n < 3) return
      legacy = .false.
      if (present(legacy_indexing)) legacy = legacy_indexing
      last_pair = n - 1
      ! vrtest::AR1 omits the final adjacent pair because of an R indexing
      ! error.  The optional compatibility mode reproduces that behavior.
      if (legacy) last_pair = n - 2
      if (last_pair < 2) return
      den = sum(x(1:last_pair)**2)
      if (den <= tiny(1.0_dp)) return
      fit%coefficient = dot_product(x(1:last_pair),x(2:last_pair+1))/den
      allocate(residual(last_pair))
      residual = x(2:last_pair+1) - fit%coefficient*x(1:last_pair)
      fit%innovation_variance = sum(residual**2)/real(max(1,last_pair-1),dp)
      fit%standard_error = sqrt(max(0.0_dp,fit%innovation_variance/den))
   end function ar1_fit

   function adjust_thin(y) result(adjusted)
      real(dp), intent(in) :: y(:)
      real(dp), allocatable :: adjusted(:)
      integer :: n
      real(dp) :: mu, b, den

      n = size(y)
      if (n < 2) then
         allocate(adjusted(0))
         return
      end if
      mu = mean_value(y)
      den = sum((y(1:n-1)-mu)**2)
      if (den <= tiny(1.0_dp)) then
         b = 0.0_dp
      else
         b = dot_product(y(1:n-1)-mu,y(2:n)-mu)/den
      end if
      allocate(adjusted(n-1))
      adjusted = ((y(2:n)-mu)-b*(y(1:n-1)-mu))/max(1.0e-12_dp,1.0_dp-b)
   end function adjust_thin

   pure real(dp) function abel_bandwidth(n, coefficient) result(bandwidth)
      integer, intent(in) :: n
      real(dp), intent(in) :: coefficient
      real(dp) :: alpha, c
      c = min(1.0_dp-1.0e-10_dp,max(-1.0_dp+1.0e-10_dp,coefficient))
      alpha = 4.0_dp*c*c/(1.0_dp-c)**4
      bandwidth = 1.3221_dp*(max(alpha*real(n,dp),tiny(1.0_dp)))**0.2_dp
      bandwidth = max(bandwidth,1.0e-8_dp)
   end function abel_bandwidth

   pure real(dp) function quadratic_spectral_kernel(x) result(weight)
      real(dp), intent(in) :: x
      real(dp) :: z
      if (abs(x) < sqrt(epsilon(1.0_dp))) then
         weight = 1.0_dp
      else
         z = 6.0_dp*pi*x/5.0_dp
         weight = 25.0_dp/(12.0_dp*pi*pi*x*x)*(sin(z)/z-cos(z))
      end if
   end function quadratic_spectral_kernel

   pure real(dp) function rolling_variance_ratio(y, k) result(vr)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: k
      integer :: i, n
      real(dp) :: mu, summ, window, vr1

      n = size(y)
      if (k < 1 .or. k > n .or. n < 2) then
         vr = 0.0_dp
         return
      end if
      mu = mean_value(y)
      vr1 = sum((y-mu)**2)/real(n,dp)
      if (vr1 <= tiny(1.0_dp)) then
         vr = 0.0_dp
         return
      end if
      window = sum(y(1:k))
      summ = (window-real(k,dp)*mu)**2
      do i = 2, n-k+1
         window = window-y(i-1)+y(i+k-1)
         summ = summ+(window-real(k,dp)*mu)**2
      end do
      vr = (summ/(real(n*k,dp)))/vr1
   end function rolling_variance_ratio

   pure function fast_variance_ratio(x, kvec) result(vr)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: kvec(:)
      real(dp) :: vr(size(kvec))
      integer :: i, j, k
      real(dp) :: total

      do i = 1, size(kvec)
         k = kvec(i)
         if (k <= 1 .or. k > size(x)) then
            vr(i) = 1.0_dp
            cycle
         end if
         total = 0.0_dp
         do j = 1, k-1
            total = total+(1.0_dp-real(j,dp)/real(k,dp))*autocorrelation(x,j)
         end do
         vr(i) = 1.0_dp+2.0_dp*total
      end do
   end function fast_variance_ratio

   pure function lm_statistic(y, k, fast) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: k
      logical, intent(in), optional :: fast
      type(lm_result) :: ans
      integer :: j, n
      real(dp) :: tem1, tem2, summ, denom, fast_one(1)
      real(dp), allocatable :: y1(:)
      logical :: use_fast

      n = size(y)
      if (k <= 1 .or. k > n) return
      use_fast = .false.
      if (present(fast)) use_fast = fast
      if (use_fast) then
         fast_one = fast_variance_ratio(y,[k])
         ans%variance_ratio = fast_one(1)
      else
         ans%variance_ratio = rolling_variance_ratio(y,k)
      end if
      tem1 = 2.0_dp*real(2*k-1,dp)*real(k-1,dp)
      tem2 = 3.0_dp*real(k,dp)
      ans%homoskedastic = sqrt(real(n,dp))*(ans%variance_ratio-1.0_dp)/sqrt(tem1/tem2)
      allocate(y1(n))
      y1 = (y-mean_value(y))**2
      denom = sum(y1)**2
      if (denom <= tiny(1.0_dp)) return
      summ = 0.0_dp
      do j = 1, k-1
         summ = summ+4.0_dp*(1.0_dp-real(j,dp)/real(k,dp))**2 * &
            sum(y1(j+1:n)*y1(1:n-j))/denom
      end do
      if (summ > tiny(1.0_dp)) then
         ans%heteroskedastic = sqrt(real(n,dp))*(ans%variance_ratio-1.0_dp)/sqrt(real(n,dp)*summ)
      end if
   end function lm_statistic

   pure function lmcd_statistics(y, kvec, fast) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      logical, intent(in), optional :: fast
      type(lmcd_result) :: ans
      type(lm_result) :: lm
      integer :: i
      logical :: use_fast

      use_fast = .false.
      if (present(fast)) use_fast = fast
      allocate(ans%homoskedastic(size(kvec)),ans%heteroskedastic(size(kvec)))
      do i = 1, size(kvec)
         lm = lm_statistic(y,kvec(i),use_fast)
         ans%homoskedastic(i) = lm%homoskedastic
         ans%heteroskedastic(i) = lm%heteroskedastic
      end do
      if (size(kvec) > 0) then
         ans%cd_homoskedastic = maxval(abs(ans%homoskedastic))
         ans%cd_heteroskedastic = maxval(abs(ans%heteroskedastic))
      end if
   end function lmcd_statistics

   pure function automatic_variance_ratio(y, legacy_ar1) result(ans)
      real(dp), intent(in) :: y(:)
      logical, intent(in), optional :: legacy_ar1
      type(auto_vr_result) :: ans
      type(ar1_result) :: fit
      integer :: i, n
      real(dp) :: den, cross
      logical :: legacy

      n = size(y)
      if (n < 3) return
      legacy = .false.
      if (present(legacy_ar1)) legacy = legacy_ar1
      fit = ar1_fit(y,legacy)
      ans%ar1_coefficient = fit%coefficient
      ans%bandwidth = abel_bandwidth(n,fit%coefficient)
      den = sum(y*y)
      if (den <= tiny(1.0_dp)) return
      ans%variance_ratio_sum = 1.0_dp
      do i = 1, n-1
         cross = dot_product(y(1:n-i),y(i+1:n))/den
         ans%variance_ratio_sum = ans%variance_ratio_sum + &
            2.0_dp*quadratic_spectral_kernel(real(i,dp)/ans%bandwidth)*cross
      end do
      ans%statistic = sqrt(real(n,dp)/ans%bandwidth)* &
         (ans%variance_ratio_sum-1.0_dp)/sqrt(2.0_dp)
   end function automatic_variance_ratio

   function automatic_vr_bootstrap(y, nboot, wild, probabilities) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: nboot
      character(len=*), intent(in) :: wild
      real(dp), intent(in), optional :: probabilities(:)
      type(auto_bootstrap_result) :: ans
      type(auto_vr_result) :: observed, simulated
      real(dp), allocatable :: weights(:), ys(:), stats(:), sums(:), probs(:)
      integer :: i, j, exceed

      observed = automatic_variance_ratio(y)
      ans%test_statistic = observed%statistic
      ans%variance_ratio_sum = observed%variance_ratio_sum
      if (present(probabilities)) then
         allocate(probs(size(probabilities)))
         probs = probabilities
      else
         allocate(probs(2))
         probs = [0.025_dp,0.975_dp]
      end if
      allocate(ans%statistic_interval(size(probs)),ans%vr_sum_interval(size(probs)))
      if (nboot <= 0) return
      allocate(weights(size(y)),ys(size(y)),stats(nboot),sums(nboot))
      exceed = 0
      do i = 1, nboot
         call wild_weights(wild,weights)
         ys = y*weights
         simulated = automatic_variance_ratio(ys)
         stats(i) = simulated%statistic
         sums(i) = simulated%variance_ratio_sum
         if (abs(stats(i)) > abs(observed%statistic)) exceed = exceed+1
      end do
      ans%p_value = real(exceed,dp)/real(nboot,dp)
      do j = 1, size(probs)
         ans%statistic_interval(j) = sample_quantile(stats,probs(j))
         ans%vr_sum_interval(j) = sample_quantile(sums,probs(j))
      end do
   end function automatic_vr_bootstrap

   pure function lo_mackinlay(y, kvec) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      type(lmcd_result) :: ans
      ans = lmcd_statistics(y,kvec,.false.)
   end function lo_mackinlay

   pure function chow_denning(y, kvec) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      type(chow_denning_result) :: ans
      type(lmcd_result) :: lmcd
      real(dp), parameter :: alpha(3) = [0.1_dp,0.05_dp,0.01_dp]
      real(dp) :: per
      integer :: i

      allocate(ans%holding_periods(size(kvec)))
      ans%holding_periods = kvec
      lmcd = lmcd_statistics(y,kvec,.false.)
      ans%cd_homoskedastic = lmcd%cd_homoskedastic
      ans%cd_heteroskedastic = lmcd%cd_heteroskedastic
      do i = 1, 3
         per = 0.5_dp*(1.0_dp-(1.0_dp-alpha(i))**(1.0_dp/real(size(kvec),dp)))
         ans%critical_values(i) = normal_quantile(1.0_dp-per)
      end do
   end function chow_denning

   function variance_ratio_bootstrap(y, kvec, nboot, wild, probabilities) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:), nboot
      character(len=*), intent(in) :: wild
      real(dp), intent(in), optional :: probabilities(:)
      type(bootstrap_result) :: ans
      type(lmcd_result) :: obs, sim
      real(dp), allocatable :: probs(:), weights(:), ys(:), statmat(:,:), observed(:)
      integer, allocatable :: indices(:)
      integer :: i, j, exceed
      logical :: ordinary

      allocate(ans%holding_periods(size(kvec)),ans%lm_p_values(size(kvec)))
      ans%holding_periods = kvec
      if (present(probabilities)) then
         allocate(probs(size(probabilities)))
         probs = probabilities
      else
         allocate(probs(2))
         probs = [0.025_dp,0.975_dp]
      end if
      allocate(ans%confidence_intervals(size(kvec),size(probs)))
      ordinary = trim(adjustl(wild)) == 'No' .or. trim(adjustl(wild)) == 'no'
      obs = lmcd_statistics(y,kvec,.true.)
      allocate(observed(size(kvec)+1),statmat(max(0,nboot),size(kvec)+1))
      if (ordinary) then
         observed = [obs%homoskedastic,obs%cd_homoskedastic]
      else
         observed = [obs%heteroskedastic,obs%cd_heteroskedastic]
      end if
      if (nboot <= 0) return
      allocate(weights(size(y)),ys(size(y)),indices(size(y)))
      do i = 1, nboot
         if (ordinary) then
            call bootstrap_indices(size(y),indices)
            ys = y(indices)
         else
            call wild_weights(wild,weights)
            ys = y*weights
         end if
         sim = lmcd_statistics(ys,kvec,.true.)
         if (ordinary) then
            statmat(i,:) = [sim%homoskedastic,sim%cd_homoskedastic]
         else
            statmat(i,:) = [sim%heteroskedastic,sim%cd_heteroskedastic]
         end if
      end do
      do j = 1, size(kvec)+1
         exceed = count(abs(statmat(:,j)) > abs(observed(j)))
         if (j <= size(kvec)) then
            ans%lm_p_values(j) = real(exceed,dp)/real(nboot,dp)
            do i = 1, size(probs)
               ans%confidence_intervals(j,i) = sample_quantile(statmat(:,j),probs(i))
            end do
         else
            ans%cd_p_value = real(exceed,dp)/real(nboot,dp)
         end if
      end do
   end function variance_ratio_bootstrap

   subroutine bootstrap_indices(n, idx)
      integer, intent(in) :: n
      integer, intent(out) :: idx(n)
      real(dp) :: u(n)
      call random_number(u)
      idx = min(n,1+int(u*real(n,dp)))
   end subroutine bootstrap_indices

   pure function covariance_matrix(kvec) result(mat)
      integer, intent(in) :: kvec(:)
      real(dp) :: mat(size(kvec),size(kvec))
      integer :: i, j, kmin, kmax
      do i = 1, size(kvec)
         do j = 1, size(kvec)
            kmin = min(kvec(i),kvec(j))
            kmax = max(kvec(i),kvec(j))
            mat(i,j) = 2.0_dp*real((3*kmax-kmin-1)*(kmin-1),dp)/real(3*kmax,dp)
         end do
      end do
   end function covariance_matrix

   pure function wald_test(y, kvec) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      type(wald_result) :: ans
      real(dp), allocatable :: vr(:), solution(:), mat(:,:)
      real(dp), parameter :: alpha(3) = [0.1_dp,0.05_dp,0.01_dp]
      integer :: i

      allocate(ans%holding_periods(size(kvec)),vr(size(kvec)),solution(size(kvec)), &
         mat(size(kvec),size(kvec)))
      ans%holding_periods = kvec
      do i = 1, size(kvec)
         vr(i) = rolling_variance_ratio(y,kvec(i))-1.0_dp
      end do
      mat = covariance_matrix(kvec)
      call solve_linear(mat,vr,solution,ans%solve_info)
      if (ans%solve_info == 0) ans%statistic = real(size(y),dp)*dot_product(vr,solution)
      do i = 1, 3
         ans%critical_values(i) = chi_square_quantile(1.0_dp-alpha(i),real(size(kvec),dp))
      end do
   end function wald_test

   pure real(dp) function standardized_vr_statistic(x, k) result(ans)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      integer :: i, n
      real(dp) :: window, summ, vr1, vr2, tem1, tem2

      n = size(x)
      if (k <= 1 .or. k > n) then
         ans = 0.0_dp
         return
      end if
      window = sum(x(1:k))
      summ = window*window
      do i = 2, n-k+1
         window = window-x(i-1)+x(i+k-1)
         summ = summ+window*window
      end do
      vr1 = sum(x*x)/real(n,dp)
      if (vr1 <= tiny(1.0_dp)) then
         ans = 0.0_dp
         return
      end if
      vr2 = summ/real(n*k,dp)
      tem1 = 2.0_dp*real(2*k-1,dp)*real(k-1,dp)
      tem2 = 3.0_dp*real(k*n,dp)
      ans = (vr2/vr1-1.0_dp)/sqrt(tem1/tem2)
   end function standardized_vr_statistic

   pure function wright_statistics_one(y, k) result(stats)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: k
      real(dp) :: stats(3)
      real(dp), allocatable :: ranks(:), r1(:), r2(:), signs(:)
      integer :: n

      n = size(y)
      allocate(ranks(n),r1(n),r2(n),signs(n))
      call average_ranks(y,ranks)
      r1 = (ranks-0.5_dp*real(n+1,dp))/sqrt(real((n-1)*(n+1),dp)/12.0_dp)
      r2 = normal_quantile_array(ranks/real(n+1,dp))
      signs = merge(1.0_dp,-1.0_dp,y > 0.0_dp)
      stats(1) = standardized_vr_statistic(r1,k)
      stats(2) = standardized_vr_statistic(r2,k)
      stats(3) = standardized_vr_statistic(signs,k)
   end function wright_statistics_one

   pure function normal_quantile_array(p) result(x)
      real(dp), intent(in) :: p(:)
      real(dp) :: x(size(p))
      integer :: i
      do i = 1, size(p)
         x(i) = normal_quantile(p(i))
      end do
   end function normal_quantile_array

   pure function wright_tests(y, kvec) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      type(wright_result) :: ans
      integer :: i
      allocate(ans%statistics(size(kvec),3))
      do i = 1, size(kvec)
         ans%statistics(i,:) = wright_statistics_one(y,kvec(i))
      end do
   end function wright_tests

   pure function joint_wright_tests(y, kvec) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      type(joint_wright_result) :: ans
      type(wright_result) :: all_stats
      allocate(ans%holding_periods(size(kvec)))
      ans%holding_periods = kvec
      all_stats = wright_tests(y,kvec)
      if (size(kvec) > 0) then
         ans%rank_uniform = maxval(abs(all_stats%statistics(:,1)))
         ans%rank_normal = maxval(abs(all_stats%statistics(:,2)))
         ans%sign = maxval(abs(all_stats%statistics(:,3)))
      end if
   end function joint_wright_tests

   function wright_critical_values(n, k, niter, probabilities) result(ans)
      integer, intent(in) :: n, k, niter
      real(dp), intent(in), optional :: probabilities(:)
      type(wright_critical_result) :: ans
      real(dp), allocatable :: probs(:), mat(:,:), r1(:), r2(:), signs(:), z(:)
      integer, allocatable :: perm(:)
      integer :: i, j

      allocate(ans%holding_periods(1))
      ans%holding_periods = k
      if (present(probabilities)) then
         allocate(probs(size(probabilities)))
         probs = probabilities
      else
         allocate(probs(6))
         probs = [0.005_dp,0.025_dp,0.05_dp,0.95_dp,0.975_dp,0.995_dp]
      end if
      allocate(ans%critical_values(size(probs),3),mat(max(0,niter),3), &
         r1(n),r2(n),signs(n),z(n),perm(n))
      do i = 1, niter
         call random_permutation(n,perm)
         r1 = (real(perm,dp)-0.5_dp*real(n+1,dp))/sqrt(real((n-1)*(n+1),dp)/12.0_dp)
         r2 = normal_quantile_array(real(perm,dp)/real(n+1,dp))
         call random_normal_vector(z)
         signs = merge(1.0_dp,-1.0_dp,z > 0.0_dp)
         mat(i,:) = [standardized_vr_statistic(r1,k),standardized_vr_statistic(r2,k), &
            standardized_vr_statistic(signs,k)]
      end do
      do j = 1, 3
         do i = 1, size(probs)
            ans%critical_values(i,j) = sample_quantile(mat(:,j),probs(i))
         end do
      end do
   end function wright_critical_values

   function joint_wright_critical_values(n, kvec, niter, probabilities) result(ans)
      integer, intent(in) :: n, kvec(:), niter
      real(dp), intent(in), optional :: probabilities(:)
      type(wright_critical_result) :: ans
      real(dp), allocatable :: probs(:), mat(:,:), r1(:), r2(:), signs(:), z(:), row(:)
      integer, allocatable :: perm(:)
      integer :: i, j, q

      allocate(ans%holding_periods(size(kvec)))
      ans%holding_periods = kvec
      if (present(probabilities)) then
         allocate(probs(size(probabilities)))
         probs = probabilities
      else
         allocate(probs(3))
         probs = [0.90_dp,0.95_dp,0.99_dp]
      end if
      allocate(ans%critical_values(size(probs),3),mat(max(0,niter),3), &
         r1(n),r2(n),signs(n),z(n),perm(n),row(size(kvec)))
      do i = 1, niter
         call random_permutation(n,perm)
         r1 = (real(perm,dp)-0.5_dp*real(n+1,dp))/sqrt(real((n-1)*(n+1),dp)/12.0_dp)
         r2 = normal_quantile_array(real(perm,dp)/real(n+1,dp))
         call random_normal_vector(z)
         signs = merge(1.0_dp,-1.0_dp,z > 0.0_dp)
         do q = 1, size(kvec)
            row(q) = standardized_vr_statistic(r1,kvec(q))
         end do
         mat(i,1) = maxval(abs(row))
         do q = 1, size(kvec)
            row(q) = standardized_vr_statistic(r2,kvec(q))
         end do
         mat(i,2) = maxval(abs(row))
         do q = 1, size(kvec)
            row(q) = standardized_vr_statistic(signs,kvec(q))
         end do
         mat(i,3) = maxval(abs(row))
      end do
      do j = 1, 3
         do i = 1, size(probs)
            ans%critical_values(i,j) = sample_quantile(mat(:,j),probs(i))
         end do
      end do
   end function joint_wright_critical_values

   pure real(dp) function wk_statistic(y, k) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: k
      ans = sqrt(real(size(y),dp))*(rolling_variance_ratio(y,k)-1.0_dp)
   end function wk_statistic

   pure real(dp) function wk_joint_statistic(y, kvec) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      real(dp) :: stats(size(kvec))
      integer :: i
      do i = 1, size(kvec)
         stats(i) = wk_statistic(y,kvec(i))
      end do
      if (size(stats) == 0) then
         ans = 0.0_dp
      else
         ans = maxval(abs(stats))
      end if
   end function wk_joint_statistic

   pure function subsample_variance_ratio(y, kvec) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      type(subsample_result) :: ans
      integer :: n, b_low, b_high, step, i, j, b, nblocks, exceed
      real(dp) :: observed

      n = size(y)
      allocate(ans%holding_periods(size(kvec)),ans%block_lengths(6),ans%p_values(6))
      ans%holding_periods = kvec
      b_low = int(2.5_dp*real(n,dp)**0.3_dp)
      b_high = int(3.5_dp*real(n,dp)**0.6_dp)
      step = max(1,(b_high-b_low)/7)
      do i = 1, 6
         ans%block_lengths(i) = min(n,b_low+i*step)
      end do
      observed = wk_joint_statistic(y,kvec)
      do i = 1, 6
         b = ans%block_lengths(i)
         nblocks = n-b+1
         exceed = 0
         do j = 1, nblocks
            if (wk_joint_statistic(y(j:j+b-1),kvec) > observed) exceed = exceed+1
         end do
         ans%p_values(i) = real(exceed,dp)/real(max(1,nblocks),dp)
      end do
   end function subsample_variance_ratio

   function panel_variance_ratio(data, nboot) result(ans)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in) :: nboot
      type(panel_vr_result) :: ans
      real(dp), allocatable :: stats(:), boot_stats(:), weights(:), ys(:,:)
      integer :: i, j, k, exceed1, exceed2, exceed3
      type(auto_vr_result) :: av

      k = size(data,2)
      allocate(stats(k),boot_stats(k),weights(size(data,1)),ys(size(data,1),k))
      do j = 1, k
         av = automatic_variance_ratio(data(:,j))
         stats(j) = av%statistic
      end do
      ans%max_absolute_statistic = maxval(abs(stats))
      ans%sum_square_statistic = sum(stats**2)
      ans%mean_statistic = sqrt(real(k,dp))*mean_value(stats)
      exceed1 = 0
      exceed2 = 0
      exceed3 = 0
      do i = 1, nboot
         call wild_weights('Mammen',weights)
         do j = 1, k
            ys(:,j) = data(:,j)*weights
            av = automatic_variance_ratio(ys(:,j))
            boot_stats(j) = av%statistic
         end do
         if (maxval(abs(boot_stats)) > ans%max_absolute_statistic) exceed1 = exceed1+1
         if (sum(boot_stats**2) > ans%sum_square_statistic) exceed2 = exceed2+1
         if (abs(sqrt(real(k,dp))*mean_value(boot_stats)) > abs(ans%mean_statistic)) exceed3 = exceed3+1
      end do
      if (nboot > 0) then
         ans%max_absolute_p_value = real(exceed1,dp)/real(nboot,dp)
         ans%sum_square_p_value = real(exceed2,dp)/real(nboot,dp)
         ans%mean_p_value = real(exceed3,dp)/real(nboot,dp)
      end if
   end function panel_variance_ratio

   pure function variance_ratio_minus_one(y, kvec) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kvec(:)
      type(vr_minus_one_result) :: ans
      type(auto_vr_result) :: av
      integer :: i
      allocate(ans%holding_periods(size(kvec)),ans%values(size(kvec)))
      ans%holding_periods = kvec
      av = automatic_variance_ratio(y)
      ans%automatic = av%variance_ratio_sum-1.0_dp
      do i = 1, size(kvec)
         ans%values(i) = rolling_variance_ratio(y,kvec(i))-1.0_dp
      end do
   end function variance_ratio_minus_one

   pure function variance_ratio_curve(y, max_holding_period) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: max_holding_period
      type(vr_curve_result) :: ans
      integer :: i, k, nout, n
      real(dp) :: tem1, tem2, se

      n = size(y)
      nout = max(0,min(max_holding_period,n)-1)
      allocate(ans%holding_periods(nout),ans%variance_ratios(nout), &
         ans%lower_95(nout),ans%upper_95(nout))
      do i = 1, nout
         k = i+1
         ans%holding_periods(i) = k
         ans%variance_ratios(i) = rolling_variance_ratio(y,k)
         tem1 = 2.0_dp*real(2*k-1,dp)*real(k-1,dp)
         tem2 = 3.0_dp*real(k*n,dp)
         se = sqrt(tem1/tem2)
         ans%lower_95(i) = 1.0_dp-1.96_dp*se
         ans%upper_95(i) = 1.0_dp+1.96_dp*se
      end do
   end function variance_ratio_curve

end module vrtest_variance_ratio
