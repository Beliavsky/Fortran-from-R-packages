! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
!
! Portmanteau, spectral, and martingale-difference tests from vrtest 1.2.
module vrtest_spectral
   use vrtest_kinds, only : dp, pi
   use vrtest_types
   use vrtest_utils, only : mean_value, variance_value, chi_square_cdf, &
      chi_square_quantile, sample_quantile, wild_weights, log_mean_exp, &
      simpson_integral, solve_linear
   implicit none
   private

   public :: automatic_portmanteau, average_exponential_test
   public :: spectral_shape_test, generalized_spectral_test
   public :: dominguez_lobato_statistic, dominguez_lobato_test
   public :: chen_deo_test

contains

   pure function automatic_portmanteau(y, lags) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: lags
      type(auto_q_result) :: ans
      real(dp), allocatable :: data(:), acov(:), fourth(:), ratio(:), bp(:), criterion(:)
      real(dp) :: penalty, maxrho, threshold
      integer :: i, n

      n = size(y)
      if (n < 2 .or. lags < 1 .or. lags >= n) return
      allocate(data(n),acov(lags),fourth(lags),ratio(lags),bp(lags),criterion(lags))
      data = y-mean_value(y)
      do i = 1, lags
         acov(i) = sum(data(1:n-i)*data(i+1:n))/real(n,dp)
         fourth(i) = sum(data(1:n-i)**2*data(i+1:n)**2)/real(n-i,dp)
         if (fourth(i) > tiny(1.0_dp)) then
            ratio(i) = acov(i)**2/fourth(i)
         else
            ratio(i) = 0.0_dp
         end if
      end do
      bp(1) = real(n,dp)*ratio(1)
      do i = 2, lags
         bp(i) = bp(i-1)+real(n,dp)*ratio(i)
      end do
      maxrho = sqrt(real(n,dp))*sqrt(maxval(ratio))
      threshold = sqrt(2.4_dp*log(real(n,dp)))
      penalty = 2.0_dp
      if (maxrho <= threshold) penalty = log(real(n,dp))
      do i = 1, lags
         criterion(i) = bp(i)-real(i,dp)*penalty
      end do
      ans%selected_lag = maxloc(criterion,dim=1)
      ans%statistic = bp(ans%selected_lag)
      ans%p_value = 1.0_dp-chi_square_cdf(ans%statistic,1.0_dp)
   end function automatic_portmanteau

   pure function average_exponential_test(y) result(ans)
      real(dp), intent(in) :: y(:)
      type(average_exponential_result) :: ans
      integer, parameter :: ngrid = 161
      real(dp), allocatable :: b(:), half_lm(:), half_lr(:)
      real(dp) :: sigma_t, sigma_h, c, cross, bsum, lm, lr
      integer :: j, k, n

      n = size(y)
      if (n < 3) return
      allocate(b(n-1),half_lm(ngrid),half_lr(ngrid))
      sigma_t = sum(y*y)/real(n-1,dp)
      if (sigma_t <= tiny(1.0_dp)) return
      do j = 1, ngrid
         c = -0.8_dp+0.01_dp*real(j-1,dp)
         b(1) = y(1)
         do k = 2, n-1
            b(k) = y(k)+c*b(k-1)
         end do
         cross = dot_product(y(2:n),b)
         bsum = sum(b*b)
         sigma_h = sigma_t-cross*cross/(real(n,dp)*max(bsum,tiny(1.0_dp)))
         lm = cross*cross/real(n,dp)*(1.0_dp-c*c)/(sigma_t*sigma_t)
         if (sigma_h > tiny(1.0_dp)) then
            lr = real(n,dp)*log(sigma_t/sigma_h)
         else
            lr = huge(1.0_dp)/100.0_dp
         end if
         half_lm(j) = 0.5_dp*lm
         half_lr(j) = 0.5_dp*lr
      end do
      ans%exponential_lm = log_mean_exp(half_lm)
      ans%exponential_lr = log_mean_exp(half_lr)
   end function average_exponential_test

   pure function integrated_acm(r, residuals) result(values)
      real(dp), intent(in) :: r, residuals(:)
      real(dp) :: values(3)
      integer :: n, s, t, upper
      real(dp) :: gamma0, fhat, lambda, cosine_sum, sine_sum, periodogram
      real(dp) :: u, us

      n = size(residuals)
      values = 0.0_dp
      if (n < 2) return
      gamma0 = sum(residuals**2)/real(n,dp)
      if (gamma0 <= tiny(1.0_dp)) return
      fhat = 0.0_dp
      if (r > 0.0_dp) then
         upper = int(real(n,dp)*r/2.0_dp)
         do s = 1, upper
            lambda = 2.0_dp*pi*real(s,dp)/real(n,dp)
            cosine_sum = 0.0_dp
            sine_sum = 0.0_dp
            do t = 0, n-1
               cosine_sum = cosine_sum+cos(lambda*real(t,dp))*residuals(t+1)
               sine_sum = sine_sum+sin(lambda*real(t,dp))*residuals(t+1)
            end do
            periodogram = (cosine_sum**2+sine_sum**2)/(2.0_dp*pi*real(n,dp))
            fhat = fhat+periodogram
         end do
      end if
      fhat = 2.0_dp*pi*fhat/real(n,dp)
      u = fhat-gamma0*r/2.0_dp
      us = sqrt(2.0_dp*real(n,dp))*u/gamma0
      if (r > 0.0_dp .and. r < 1.0_dp) values(1) = us*us/(r*(1.0_dp-r))
      values(2) = us*us
      values(3) = abs(us)
   end function integrated_acm

   pure function spectral_shape_test(x) result(ans)
      real(dp), intent(in) :: x(:)
      type(spectral_shape_result) :: ans
      real(dp) :: ad(11), cvm(11), mx(11), values(3), r
      integer :: j

      do j = 1, 11
         r = 0.1_dp*real(j-1,dp)
         values = integrated_acm(r,x)
         ad(j) = values(1)
         cvm(j) = values(2)
         mx(j) = values(3)
      end do
      ans%anderson_darling = simpson_integral(0.0_dp,1.0_dp,ad)
      ans%cramer_von_mises = simpson_integral(0.0_dp,1.0_dp,cvm)
      ans%maximum = simpson_integral(0.0_dp,1.0_dp,mx)
   end function spectral_shape_test

   pure real(dp) function generalized_spectral_statistic_weighted(residuals, weight) result(statistic)
      real(dp), intent(in) :: residuals(:), weight(:,:)
      real(dp) :: v, term
      integer :: n, i, j, lag, m

      n = size(residuals)
      statistic = 0.0_dp
      if (n < 3 .or. size(weight,1) /= n-1 .or. size(weight,2) /= n-1) return
      v = variance_value(residuals,.true.)
      if (v <= tiny(1.0_dp)) return
      do lag = 1, n-1
         m = n-lag
         term = 0.0_dp
         do i = 1, m
            do j = 1, m
               term = term+residuals(i+lag)*weight(i,j)*residuals(j+lag)
            end do
         end do
         statistic = statistic+term/(real(lag,dp)*pi)**2/real(n-lag+1,dp)
      end do
      statistic = statistic/v
   end function generalized_spectral_statistic_weighted

   pure subroutine generalized_spectral_weights(y, weight)
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: weight(size(y)-1,size(y)-1)
      real(dp) :: diff
      integer :: i, j
      do i = 1, size(y)-1
         do j = 1, size(y)-1
            diff = y(i)-y(j)
            weight(i,j) = exp(-0.5_dp*diff*diff)
         end do
      end do
   end subroutine generalized_spectral_weights

   function generalized_spectral_test(y, nboot) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: nboot
      type(generalized_spectral_result) :: ans
      real(dp), allocatable :: residuals(:), weights(:), boot(:), ys(:), kernel(:,:)
      integer :: i, exceed, n

      n = size(y)
      if (n < 3) return
      allocate(residuals(n),kernel(n-1,n-1))
      residuals = y-mean_value(y)
      call generalized_spectral_weights(y,kernel)
      ans%statistic = generalized_spectral_statistic_weighted(residuals,kernel)
      if (nboot <= 0) return
      allocate(weights(n),boot(nboot),ys(n))
      exceed = 0
      do i = 1, nboot
         call wild_weights('Mammen',weights)
         ys = residuals*weights
         ys = ys-mean_value(ys)
         boot(i) = generalized_spectral_statistic_weighted(ys,kernel)
         if (boot(i) > ans%statistic) exceed = exceed+1
      end do
      ans%p_value = real(exceed,dp)/real(nboot,dp)
      ans%bootstrap_critical_values = [sample_quantile(boot,0.90_dp), &
         sample_quantile(boot,0.95_dp),sample_quantile(boot,0.99_dp)]
   end function generalized_spectral_test

   pure function dominguez_lobato_statistic(y, p) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: p
      type(dl_result) :: ans
      real(dp), allocatable :: centered(:), point_stats(:)
      real(dp) :: variance, inner
      integer :: n, i, j, q, m
      logical :: indicator

      n = size(y)
      if (p < 1 .or. p >= n) return
      m = n-p
      allocate(centered(n),point_stats(m))
      centered = y-mean_value(y)
      variance = sum(centered**2)/real(m,dp)
      if (variance <= tiny(1.0_dp)) return
      ans%cp_statistic = 0.0_dp
      do j = p+1, n
         inner = 0.0_dp
         do i = p+1, n
            indicator = .true.
            do q = 1, p
               if (centered(i-q) > centered(j-q)) then
                  indicator = .false.
                  exit
               end if
            end do
            if (indicator) inner = inner+centered(i)
         end do
         ans%cp_statistic = ans%cp_statistic+inner*inner
         point_stats(j-p) = abs(inner/sqrt(real(m,dp)))
      end do
      ans%cp_statistic = ans%cp_statistic/(variance*real(m*m,dp))
      ans%kp_statistic = maxval(point_stats)/sqrt(variance)
   end function dominguez_lobato_statistic

   function dominguez_lobato_test(y, nboot, p) result(ans)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: nboot, p
      type(dl_result) :: ans
      type(dl_result) :: boot_result
      real(dp), allocatable :: centered(:), weights(:), ys(:)
      integer :: i, cp_exceed, kp_exceed, n

      n = size(y)
      ans = dominguez_lobato_statistic(y,p)
      if (nboot <= 0) return
      allocate(centered(n),weights(n),ys(n))
      centered = y-mean_value(y)
      cp_exceed = 0
      kp_exceed = 0
      do i = 1, nboot
         call wild_weights('Mammen',weights)
         weights = weights-mean_value(weights)
         ys = centered*weights
         boot_result = dominguez_lobato_statistic(ys,p)
         if (abs(boot_result%cp_statistic) > abs(ans%cp_statistic)) cp_exceed = cp_exceed+1
         if (abs(boot_result%kp_statistic) > abs(ans%kp_statistic)) kp_exceed = kp_exceed+1
      end do
      ans%cp_p_value = real(cp_exceed,dp)/real(nboot,dp)
      ans%kp_p_value = real(kp_exceed,dp)/real(nboot,dp)
   end function dominguez_lobato_test

   pure function chen_deo_test(x, kvec) result(ans)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: kvec(:)
      type(chen_deo_result) :: ans
      real(dp), allocatable :: z(:), lambda(:), wmat(:,:), w1(:), w2(:), w3(:)
      real(dp), allocatable :: beta(:), periodogram(:), sums(:), vp(:), tau(:)
      real(dp), allocatable :: lmat(:,:), bvec(:), avec(:), sigma_core(:,:), sigma(:,:)
      real(dp), allocatable :: mu_beta(:), sigma_beta(:,:), delta(:), solved(:)
      complex(dp) :: dft, phase
      real(dp) :: sigma2, cnk, denom
      integer :: n, nfreq, d, ks, i, j, t
      real(dp), parameter :: upper_probs(5) = [0.01_dp,0.02_dp,0.05_dp,0.10_dp,0.20_dp]

      n = size(x)
      d = size(kvec)
      allocate(ans%holding_periods(d))
      ans%holding_periods = kvec
      if (n < 8 .or. d == 0 .or. minval(kvec) < 2 .or. maxval(kvec) >= n) return
      ks = maxval(kvec)
      nfreq = (n-1)/2
      sigma2 = variance_value(x,.true.)
      if (sigma2 <= tiny(1.0_dp)) return
      allocate(z(n),lambda(nfreq),wmat(nfreq,d),w1(d),w2(d),w3(d),beta(d), &
         periodogram(nfreq),sums(d),vp(d),tau(ks),lmat(ks+1,d),bvec(ks), &
         avec(ks),sigma_core(ks+1,ks+1),sigma(d,d),mu_beta(d), &
         sigma_beta(d,d),delta(d),solved(d))
      z = x-mean_value(x)
      do i = 1, nfreq
         lambda(i) = 2.0_dp*pi*real(i,dp)/real(n,dp)
         do j = 1, d
            wmat(i,j) = (sin(real(kvec(j),dp)*lambda(i)/2.0_dp)/ &
               sin(lambda(i)/2.0_dp))**2/real(kvec(j),dp)
         end do
         dft = (0.0_dp,0.0_dp)
         do t = 1, n
            phase = exp(cmplx(0.0_dp,-lambda(i)*real(t,dp),kind=dp))
            dft = dft+cmplx(z(t),0.0_dp,kind=dp)*phase
         end do
         periodogram(i) = abs(dft)**2/(2.0_dp*pi*real(n,dp))
      end do
      w1 = sum(wmat,dim=1)
      w2 = sum(wmat*wmat,dim=1)
      w3 = sum(wmat*wmat*wmat,dim=1)
      beta = 1.0_dp-(2.0_dp/3.0_dp)*(w1*w3)/(w2*w2)
      do j = 1, d
         sums(j) = sum(periodogram*wmat(:,j)) / (1.0_dp-real(kvec(j),dp)/real(n,dp)) * &
            4.0_dp*pi/(real(n,dp)*sigma2)
         vp(j) = exp(beta(j)*log(max(sums(j),tiny(1.0_dp))))
      end do
      do j = 1, ks
         denom = real(n-j-4,dp)
         if (abs(denom) <= tiny(1.0_dp)) denom = sign(1.0_dp,denom)*tiny(1.0_dp)
         tau(j) = sum(z(j+1:n)**2*z(1:n-j)**2)/(sigma2*sigma2*denom)
      end do
      lmat = 0.0_dp
      do i = 1, d
         cnk = real(n,dp)/real(n-kvec(i),dp)
         do j = 1, kvec(i)-1
            lmat(j,i) = 2.0_dp*cnk*(1.0_dp-real(j,dp)/real(kvec(i),dp))
         end do
         lmat(ks+1,i) = -(real(kvec(i),dp)*cnk-real(n,dp)/real(n-1,dp))
      end do
      do j = 1, ks
         bvec(j) = 2.0_dp*real(n-j,dp)*tau(j)/real(n,dp)**3 + &
            2.0_dp*real(j,dp)/real(n,dp)**3
         avec(j) = real(n-j,dp)*tau(j)/real(n,dp)**2 + real(j,dp)/real(n,dp)**2
      end do
      sigma_core = 0.0_dp
      do j = 1, ks
         sigma_core(j,j) = avec(j)
         sigma_core(j,ks+1) = bvec(j)
         sigma_core(ks+1,j) = bvec(j)
      end do
      sigma_core(ks+1,ks+1) = 2.0_dp/real(n,dp)**2
      sigma = matmul(transpose(lmat),matmul(sigma_core,lmat))
      do i = 1, d
         mu_beta(i) = 1.0_dp+0.5_dp*beta(i)*(beta(i)-1.0_dp)*sigma(i,i)
         do j = 1, d
            sigma_beta(i,j) = beta(i)*beta(j)*sigma(i,j)
         end do
      end do
      ans%variance_ratio_sum = sum(vp-1.0_dp)
      delta = vp-mu_beta
      call solve_linear(sigma_beta,delta,solved,ans%solve_info)
      if (ans%solve_info == 0) ans%qp_statistic = dot_product(delta,solved)
      do i = 1, 5
         ans%chi_square_upper_quantiles(i) = chi_square_quantile(1.0_dp-upper_probs(i),real(d,dp))
      end do
   end function chen_deo_test

end module vrtest_spectral
