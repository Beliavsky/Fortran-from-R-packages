program test_v02
   use vgam_kinds, only : dp
   use vgam_vglm, only : family_gaussian, family_poisson
   use vgam_distributions, only : dpois_v
   use vgam_reduced_rank, only : rrvglm_result_t, fit_rrvglm
   use vgam_count_models, only : dztpois_v, dhurdlepois_v, dzinb_v, &
      zero_truncated_poisson_result_t, hurdle_count_result_t, zinb_result_t, &
      fit_zero_truncated_poisson, fit_hurdle_poisson, &
      fit_hurdle_negative_binomial, fit_zero_inflated_negative_binomial
   use vgam_gaitd, only : gaitd_distribution_t, gaitd_poisson
   use vgam_qreg, only : yeo_johnson, yeo_johnson_inverse, fit_yj_normal, yj_normal_result_t
   use vgam_timeseries, only : ar1_result_t, fit_ar1
   implicit none

   call test_rrvglm
   call test_count_models
   call test_gaitd
   call test_yj
   call test_ar1
   print '(a)', 'test_v02: PASS'

contains

   subroutine assert_true(ok, message)
      logical, intent(in) :: ok
      character(*), intent(in) :: message
      if (.not. ok) then
         print '(a)', 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(x, y, tol, message)
      real(dp), intent(in) :: x, y, tol
      character(*), intent(in) :: message
      call assert_true(abs(x - y) <= tol, message)
   end subroutine assert_close

   subroutine test_rrvglm
      integer, parameter :: n = 80, p = 3, m = 3
      real(dp) :: x(n, p), y(n, m), eta_true(n, m), btrue(p, m), t, rmse
      integer :: i
      integer :: fam(m)
      type(rrvglm_result_t) :: fit

      btrue(1, :) = [0.3_dp, -0.2_dp, 0.5_dp]
      btrue(2, :) = [0.8_dp, -0.4_dp, 1.2_dp]
      btrue(3, :) = [-0.4_dp, 0.2_dp, -0.6_dp]
      do i = 1, n
         t = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
         x(i, :) = [1.0_dp, t, sin(2.0_dp*t)]
      end do
      eta_true = matmul(x, btrue)
      y(:, 1:2) = eta_true(:, 1:2)
      y(:, 3) = exp(eta_true(:, 3))
      do i = 1, n
         y(i, 1:2) = y(i, 1:2) + 0.01_dp*sin(real(i, dp))*[1.0_dp, -0.5_dp]
      end do
      fam = [family_gaussian, family_gaussian, family_poisson]
      call fit_rrvglm(y, x, 1, fam, fit, max_iter=100, tol=1.0e-8_dp)
      call assert_true(allocated(fit%coefficients), 'rrvglm produced coefficients')
      call assert_true(fit%effective_rank == 1, 'rrvglm reduced block rank')
      rmse = sqrt(sum((fit%fitted - y)**2)/real(n*m, dp))
      call assert_true(rmse < 0.03_dp, 'rrvglm Gaussian fit accuracy')
   end subroutine test_rrvglm

   subroutine test_count_models
      real(dp), parameter :: lambda = 1.7_dp
      real(dp) :: s, expected
      integer :: k
      integer, parameter :: n = 20
      integer :: ypos(n), yh(n)
      real(dp) :: x(n, 1)
      type(zero_truncated_poisson_result_t) :: ztp
      type(hurdle_count_result_t) :: hp, hnb
      type(zinb_result_t) :: zinb

      expected = dpois_v(1, lambda)/(1.0_dp - exp(-lambda))
      call assert_close(dztpois_v(1, lambda), expected, 1.0e-13_dp, 'zero-truncated Poisson pmf')
      s = 0.0_dp
      do k = 0, 80
         s = s + dhurdlepois_v(k, lambda, 0.3_dp)
      end do
      call assert_close(s, 1.0_dp, 1.0e-11_dp, 'hurdle Poisson normalization')
      s = 0.0_dp
      do k = 0, 150
         s = s + dzinb_v(k, 2.0_dp, 3.0_dp, 0.25_dp)
      end do
      call assert_close(s, 1.0_dp, 1.0e-10_dp, 'ZINB normalization')

      x(:, 1) = 1.0_dp
      ypos = [1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4,4,5]
      call fit_zero_truncated_poisson(ypos, x, ztp, max_iter=200)
      call assert_true(allocated(ztp%coefficients), 'zero-truncated Poisson fit')
      call assert_true(ztp%fitted_mean(1) > 1.0_dp, 'zero-truncated Poisson mean')

      yh = [0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,5]
      call fit_hurdle_poisson(yh, x, x, hp, max_iter=250)
      call assert_true(allocated(hp%count_coefficients), 'hurdle Poisson fit')
      call assert_true(hp%zero_probability(1) > 0.1_dp .and. hp%zero_probability(1) < 0.6_dp, &
                       'hurdle zero probability')
      call fit_hurdle_negative_binomial(yh, x, x, hnb, max_iter=300)
      call assert_true(allocated(hnb%count_coefficients) .and. hnb%size > 0.0_dp, &
                       'hurdle negative-binomial fit')
      call fit_zero_inflated_negative_binomial(yh, x, x, zinb, max_iter=350)
      call assert_true(allocated(zinb%count_coefficients) .and. zinb%size > 0.0_dp, &
                       'zero-inflated negative-binomial fit')
   end subroutine test_count_models

   subroutine test_gaitd
      type(gaitd_distribution_t) :: d
      integer :: trunc(1), ap(1), ip(1), dpnt(1)
      real(dp) :: av(1), iv(1), df(1)
      trunc = [0]
      ap = [1]
      av = [0.20_dp]
      ip = [2]
      iv = [0.10_dp]
      dpnt = [3]
      df = [0.5_dp]
      call gaitd_poisson(2.0_dp, 40, d, truncate=trunc, altered_points=ap, &
                         altered_probabilities=av, inflated_points=ip, &
                         inflation_probabilities=iv, deflated_points=dpnt, &
                         deflation_factors=df)
      call assert_true(d%status == 0, 'GAITD Poisson construction')
      call assert_close(sum(d%pmf), 1.0_dp, 1.0e-13_dp, 'GAITD normalization')
      call assert_close(d%probability(0), 0.0_dp, 1.0e-15_dp, 'GAITD truncation')
      call assert_close(d%probability(1), 0.20_dp, 1.0e-13_dp, 'GAITD altered mass')
      call assert_true(d%variance > 0.0_dp, 'GAITD moments')
   end subroutine test_gaitd

   subroutine test_yj
      integer, parameter :: n = 41
      real(dp) :: x(n, 2), y(n), value, back, t
      integer :: i
      type(yj_normal_result_t) :: fit

      do i = 1, n
         t = -2.0_dp + 4.0_dp*real(i - 1, dp)/real(n - 1, dp)
         value = yeo_johnson(t, 0.4_dp)
         back = yeo_johnson_inverse(value, 0.4_dp)
         call assert_close(back, t, 2.0e-12_dp, 'Yeo-Johnson inverse')
         call assert_close(yeo_johnson(t, 1.0_dp), t, 2.0e-12_dp, 'Yeo-Johnson lambda=1')
         x(i, :) = [1.0_dp, t/2.0_dp]
         y(i) = 0.4_dp + 0.8_dp*x(i, 2) + 0.15_dp*sin(1.7_dp*real(i, dp))
      end do
      call fit_yj_normal(y, x, fit, max_iter=250)
      call assert_true(allocated(fit%coefficients), 'Yeo-Johnson normal fit')
      call assert_true(fit%sigma > 0.0_dp .and. fit%sigma < 1.0_dp, 'Yeo-Johnson sigma')
   end subroutine test_yj

   subroutine test_ar1
      integer, parameter :: n = 300
      real(dp) :: x(n)
      real(dp), allocatable :: fmean(:), fvar(:)
      integer :: i
      type(ar1_result_t) :: fit
      x(1) = 0.2_dp
      do i = 2, n
         x(i) = 0.15_dp + 0.65_dp*x(i - 1) + 0.15_dp*(2.0_dp* &
            modulo(sin(12.9898_dp*real(i, dp))*43758.5453_dp, 1.0_dp) - 1.0_dp)
      end do
      call fit_ar1(x, fit, max_iter=250)
      call assert_true(fit%innovation_sd > 0.0_dp, 'AR1 positive innovation sd')
      call assert_true(fit%rho > 0.3_dp .and. fit%rho < 0.95_dp, 'AR1 rho estimate')
      call fit%forecast(x(n), 3, fmean, fvar)
      call assert_true(size(fmean) == 3 .and. all(fvar > 0.0_dp), 'AR1 forecast')
   end subroutine test_ar1

end program test_v02
