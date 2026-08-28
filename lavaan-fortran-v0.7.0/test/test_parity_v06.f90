program test_parity_v06
   use lavaan_kinds, only : dp
   use lavaan_muthen1984, only : muthen1984_result, muthen1984_ordinal
   use lavaan_miiv_markers, only : miiv_marker, miiv_marker_equation, ram_miiv_marker_equations
   use lavaan_ram, only : ram_model
   use lavaan_mml_adaptive, only : mml_mixed_loglik_adaptive
   use lavaan_multilevel_random, only : random_coefficient_loglik, random_coefficient_result, fit_random_coefficient_ml
   use lavaan_multilevel_random, only : random_effects_result, random_effects_eb
   use lavaan_sam_gamma, only : sam_continuous_gamma, sam_browne_unbiased_gamma
   use lavaan_browne, only : browne_test_result, browne_residual_nt
   implicit none
   integer :: fails
   fails = 0
   call test_muthen(fails)
   call test_miiv_markers(fails)
   call test_adaptive_mml(fails)
   call test_random_coefficients(fails)
   call test_sam_gamma(fails)
   call test_browne(fails)
   if (fails /= 0) then
      write(*,'(a,i0)') 'test_parity_v06: FAIL ', fails
      error stop 1
   end if
   write(*,'(a)') 'test_parity_v06: PASS'
contains
   subroutine check(ok, msg, f)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: msg
      integer, intent(inout) :: f
      if (.not.ok) then
      f = f + 1
      write(*,'(a)') 'FAIL: '//trim(msg)
      end if
   end subroutine check

   subroutine test_muthen(f)
      integer, intent(inout) :: f
      integer :: data(600,2)
      type(muthen1984_result) :: r
      real(dp) :: pi
      pi = acos(-1.0_dp)
      data(1:200, :) = spread([1,1], 1, 200)
      data(201:300, :) = spread([1,2], 1, 100)
      data(301:400, :) = spread([2,1], 1, 100)
      data(401:600, :) = spread([2,2], 1, 200)
      call muthen1984_ordinal(data, r)
      call check(r%status == 0, 'Muthen 1984 status', f)
      call check(abs(r%categorical%correlation(1,2) - 0.5_dp) < 2.0e-6_dp, 'Muthen polychoric rho', f)
      call check(abs(r%categorical%thresholds(1)) < 1.0e-10_dp .and. &
                 abs(r%categorical%thresholds(2)) < 1.0e-10_dp, 'Muthen binary thresholds', f)
      call check(abs(r%categorical%gamma(1,1) - 0.5_dp*pi) < 2.0e-3_dp, 'Muthen threshold NACOV 1', f)
      call check(abs(r%categorical%gamma(2,2) - 0.5_dp*pi) < 2.0e-3_dp, 'Muthen threshold NACOV 2', f)
      call check(abs(r%categorical%gamma(3,3) - pi*pi/6.0_dp) < 2.0e-3_dp, 'Muthen binary rho NACOV', f)
      call check(maxval(abs(sum(r%score, dim=1))) < 2.0e-4_dp, 'Muthen stacked score equations', f)
      call check(maxval(abs(r%categorical%gamma - transpose(r%categorical%gamma))) < 1.0e-12_dp, &
                 'Muthen Gamma symmetry', f)
   end subroutine test_muthen

   subroutine test_miiv_markers(f)
      integer, intent(inout) :: f
      type(ram_model) :: m
      type(miiv_marker) :: marker(2)
      type(miiv_marker_equation), allocatable :: eq(:)
      integer :: status, k, hit
      allocate(m%a(6,6), m%s(6,6), m%m(6), m%observed(4))
      m%a = 0.0_dp
      m%s = 0.0_dp
      m%m = 0.0_dp
      m%a(2,1) = 0.5_dp
      m%a(3,1) = 1.0_dp
      m%a(4,2) = 1.0_dp
      m%a(5,1) = 0.8_dp
      m%a(6,2) = 0.9_dp
      m%s(1,1) = 1.0_dp
      m%s(2,2) = 0.75_dp
      do k = 3,6
      m%s(k,k) = 0.2_dp
      end do
      m%observed = [3,4,5,6]
      marker(1)%latent_node = 1
      marker(1)%marker_node = 3
      marker(2)%latent_node = 2
      marker(2)%marker_node = 4
      call ram_miiv_marker_equations(m, marker, eq, status)
      call check(status == 0, 'marker MIIV status', f)
      hit = 0
      do k = 1, size(eq)
         if (eq(k)%outcome_node == 2) then
            hit = k
            exit
         end if
      end do
      call check(hit > 0, 'marker MIIV structural equation found', f)
      if (hit > 0) then
         call check(eq(hit)%proxy_outcome_node == 4, 'marker MIIV outcome rewrite', f)
         call check(size(eq(hit)%proxy_predictor_nodes) == 1 .and. eq(hit)%proxy_predictor_nodes(1) == 3, &
                    'marker MIIV predictor rewrite', f)
         call check(abs(eq(hit)%proxy_coefficients(1) - 0.5_dp) < 1.0e-12_dp, 'marker MIIV coefficient scaling', f)
         call check(any(eq(hit)%instrument_nodes == 5), 'marker MIIV valid alternate indicator', f)
         call check(.not.any(eq(hit)%instrument_nodes == 6), 'marker MIIV descendant exclusion', f)
      end if
   end subroutine test_miiv_markers

   subroutine test_adaptive_mml(f)
      integer, intent(inout) :: f
      real(dp) :: data(1,1), load(1,1), intercept(1), rsd(1), thr(1,1), lmean(1), lcov(1,1), ll, exact
      logical :: ord(1)
      integer :: ncat(1)
      data = 0.0_dp
      load = 1.0_dp
      intercept = 0.0_dp
      rsd = 1.0_dp
      thr = 0.0_dp
      lmean = 0.0_dp
      lcov = 1.0_dp
      ord = .false.
      ncat = 0
      ll = mml_mixed_loglik_adaptive(data, ord, load, intercept, rsd, thr, ncat, lmean, lcov, 5)
      exact = -0.5_dp * log(4.0_dp * acos(-1.0_dp))
      call check(abs(ll - exact) < 2.0e-6_dp, 'adaptive GH exact Gaussian marginal', f)
   end subroutine test_adaptive_mml

   subroutine test_random_coefficients(f)
      integer, intent(inout) :: f
      real(dp) :: y(4,1), x(4,1), beta(1,1), z(4,1), g(1,1), r(1,1), ll, exact, pi2
      integer :: cl(4)
      real(dp) :: yf(24,1), xf(24,1), zf(24,1)
      integer :: clf(24), i, c, pos
      type(random_coefficient_result) :: fit
      type(random_effects_result) :: eb
      y(:,1) = [1.0_dp,-1.0_dp,0.5_dp,-0.5_dp]
      x = 1.0_dp
      beta = 0.0_dp
      z = 1.0_dp
      g = 0.5_dp
      r = 1.0_dp
      cl = [1,1,2,2]
      ll = random_coefficient_loglik(y, cl, x, beta, z, g, r)
      pi2 = 2.0_dp * acos(-1.0_dp)
      exact = -0.5_dp*(2.0_dp*log(pi2)+log(2.0_dp)+2.0_dp) - &
              0.5_dp*(2.0_dp*log(pi2)+log(2.0_dp)+0.5_dp)
      call check(abs(ll-exact) < 1.0e-11_dp, 'random-intercept exact likelihood', f)
      call random_effects_eb(y, cl, x, beta, z, g, r, eb)
      call check(eb%status == 0 .and. maxval(abs(eb%mean(:,1))) < 1.0e-12_dp, 'random-effect EB mean', f)
      call check(maxval(abs(eb%vcov(1,1,:) - 0.25_dp)) < 1.0e-12_dp, 'random-effect EB covariance', f)
      pos = 0
      do c = 1, 8
         do i = 1, 3
            pos = pos + 1
            clf(pos) = c
            xf(pos,1) = 1.0_dp
            zf(pos,1) = 1.0_dp
            yf(pos,1) = 1.2_dp + 0.35_dp*sin(real(c,dp)) + 0.15_dp*cos(real(3*i+c,dp))
         end do
      end do
      call fit_random_coefficient_ml(yf, clf, xf, zf, fit, compute_se=.false.)
      call check(fit%status == 0, 'random-coefficient fit status', f)
      call check(abs(fit%beta(1,1) - 1.2_dp) < 0.25_dp, 'random-coefficient fixed intercept', f)
      call check(fit%random_cov(1,1) > 0.0_dp .and. fit%residual_cov(1,1) > 0.0_dp, &
                 'random-coefficient positive covariance', f)
   end subroutine test_random_coefficients

   subroutine test_sam_gamma(f)
      integer, intent(inout) :: f
      real(dp) :: x(8,2)
      real(dp), allocatable :: gb(:,:), gu(:,:)
      integer :: i, status
      do i = 1, 8
         x(i,1) = real(i-4,dp)
         x(i,2) = 0.5_dp*x(i,1) + real(mod(i,3)-1,dp)
      end do
      call sam_continuous_gamma(x, gb, unbiased=.false., status=status)
      call check(status == 0 .and. size(gb,1) == 5, 'SAM biased ADF Gamma shape', f)
      call sam_browne_unbiased_gamma(x, gu, status)
      call check(status == 0 .and. size(gu,1) == 5, 'SAM Browne Gamma status', f)
      call check(maxval(abs(gu-transpose(gu))) < 1.0e-12_dp, 'SAM Browne Gamma symmetry', f)
      call check(gu(1,1) > 0.0_dp .and. gu(2,2) > 0.0_dp, 'SAM Browne mean block positive', f)
      call check(maxval(abs(gu-gb)) > 1.0e-6_dp, 'SAM Browne finite-sample correction active', f)
   end subroutine test_sam_gamma
   subroutine test_browne(f)
      integer, intent(inout) :: f
      real(dp) :: s(2,2), m(2,2)
      type(browne_test_result) :: r
      s = reshape([1.2_dp,0.3_dp,0.3_dp,0.8_dp],[2,2])
      m = s
      call browne_residual_nt(s,m,200,1,r)
      call check(r%status == 0 .and. abs(r%statistic) < 1.0e-14_dp, 'Browne NT exact-fit statistic', f)
      m(2,1) = 0.15_dp
      m(1,2) = 0.15_dp
      call browne_residual_nt(s,m,200,1,r)
      call check(r%status == 0 .and. r%statistic > 0.0_dp .and. abs(r%df-2.0_dp)<1.0e-12_dp, &
                 'Browne NT positive residual test', f)
   end subroutine test_browne
end program test_parity_v06
