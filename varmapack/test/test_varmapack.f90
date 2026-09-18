program test_varmapack
   use varmapack, only : dp, varmapack_autocov, varmapack_cov2corr, varmapack_model
   use varmapack, only : varmapack_model_type, varmapack_testcase, varmapack_testcases
   use randompack, only : randompack_rng, randompack_rng_type
   implicit none

   call test_analysis()
   call test_sample_statistics()
   call test_testcases()
   call test_simulation()
   call test_varmax()
   print '(a)', 'All varmapack tests passed.'

contains

   subroutine require(condition, message)
      logical, intent(in) :: condition !! Assertion condition that must be true.
      character(len=*), intent(in) :: message !! Diagnostic printed before stopping when `condition` is false.
      if (.not. condition) then
         print '(a)', 'FAILED: ' // trim(message)
         error stop 1
      end if
   end subroutine require

   subroutine test_analysis()
      type(varmapack_model_type) :: model
      real(dp), allocatable :: a(:, :, :), gamma(:, :, :), psi(:, :, :), theta(:, :, :)
      real(dp) :: rho
      integer :: info

      allocate(a(1, 1, 1))
      a = 0.5_dp
      model = varmapack_model(reshape([1.0_dp], [1, 1]), a=a, info=info)
      call require(info == 0, 'construct AR(1)')
      call model%specrad(rho, info)
      call require(info == 0 .and. abs(rho - 0.5_dp) < 1.0e-12_dp, 'AR spectral radius')
      call model%ma_specrad(rho, info)
      call require(info == 0 .and. rho == 0.0_dp, 'zero MA spectral radius')
      call model%acvf(1, gamma, info)
      call require(info == 0, 'AR autocovariance')
      call require(abs(gamma(1, 1, 1) - 4.0_dp / 3.0_dp) < 2.0e-12_dp, 'AR variance')
      call require(abs(gamma(1, 1, 2) - 2.0_dp / 3.0_dp) < 2.0e-12_dp, 'AR lag-one covariance')
      call model%psi(2, psi, info)
      call require(info == 0, 'AR psi')
      call require(maxval(abs(reshape(psi, [3]) - [1.0_dp, 0.5_dp, 0.25_dp])) < 1.0e-14_dp, 'AR psi values')
      call model%irf(2, theta, info)
      call require(info == 0, 'AR irf')
      call require(maxval(abs(theta - psi)) < 1.0e-14_dp, 'AR unit-covariance irf')

      deallocate(a)
      allocate(a(2, 2, 1))
      a(:, :, 1) = reshape([0.1_dp, 0.3_dp, 0.2_dp, 0.4_dp], [2, 2])
      model = varmapack_model(reshape([1.0_dp, 0.25_dp, 0.25_dp, 2.0_dp], [2, 2]), a=a, &
         b=reshape([0.5_dp, 0.7_dp, 0.6_dp, 0.8_dp], [2, 2, 1]), info=info)
      call require(info == 0, 'construct multivariate ARMA')
      call model%psi(2, psi, info)
      call require(info == 0, 'multivariate psi')
      call require(maxval(abs(psi - reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
         0.6_dp, 1.0_dp, 0.8_dp, 1.2_dp, 0.26_dp, 0.58_dp, 0.32_dp, 0.72_dp], [2, 2, 3]))) < &
         1.0e-12_dp, 'multivariate psi orientation')
   end subroutine test_analysis

   subroutine test_sample_statistics()
      real(dp), allocatable :: cov(:, :, :), corr(:, :, :), matrix_corr(:, :)
      real(dp) :: x(2, 8), expected(2, 2), input_cov(2, 2, 2)
      integer :: info, i

      do i = 1, 8
         x(1, i) = real(i, dp)
         x(2, i) = 2.0_dp * real(i, dp)
      end do
      call varmapack_autocov(x, 1, cov, info)
      call require(info == 0, 'sample autocovariance')
      call require(abs(cov(1, 2, 1) - 2.0_dp * cov(1, 1, 1)) < 1.0e-13_dp, 'sample orientation')
      expected = matmul(x(:, 2:8) - spread(sum(x, dim=2) / 8.0_dp, 2, 7), &
         transpose(x(:, 1:7) - spread(sum(x, dim=2) / 8.0_dp, 2, 7))) / 8.0_dp
      call require(maxval(abs(cov(:, :, 2) - expected)) < 1.0e-12_dp, 'lag-one ML normalization')

      input_cov = reshape([4.0_dp, 3.0_dp, 3.0_dp, 9.0_dp, 2.0_dp, 6.0_dp, -3.0_dp, 1.5_dp], [2, 2, 2])
      call varmapack_cov2corr(input_cov, corr, info)
      call require(info == 0, 'covariance sequence to correlation')
      call require(abs(corr(1, 1, 1) - 1.0_dp) < 1.0e-15_dp, 'correlation diagonal')
      call require(abs(corr(1, 2, 1) - 0.5_dp) < 1.0e-15_dp, 'correlation off diagonal')
      call varmapack_cov2corr(input_cov(:, :, 1), matrix_corr, info)
      call require(info == 0 .and. maxval(abs(matrix_corr - corr(:, :, 1))) < 1.0e-15_dp, 'matrix cov2corr')

      x = 0.0_dp
      x(1, 1:4) = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp]
      call varmapack_autocov(x(1:1, 1:4), 2, cov, info, corrected=.true.)
      call require(info == 0, 'corrected sample autocovariance')
      call require(abs(cov(1, 1, 3) + 3.0_dp) < 1.0e-14_dp, 'corrected lag-two normalization')
   end subroutine test_sample_statistics

   subroutine test_testcases()
      character(len=16) :: names(16)
      integer :: p(16), q(16), r(16), info
      type(varmapack_model_type) :: named, indexed, rho_model, random1, random2
      type(randompack_rng_type) :: rng1, rng2
      real(dp) :: radius

      call varmapack_testcases(names, p, q, r)
      call require(trim(names(1)) == 'tinyAR', 'first testcase name')
      call require(trim(names(16)) == 'largeARMA', 'last testcase name')
      call require(all([p(15), q(15), r(15)] == [5, 0, 7]), 'largeAR dimensions')
      named = varmapack_testcase('smallARMA1', info=info)
      call require(info == 0, 'named testcase')
      indexed = varmapack_testcase(8, info=info)
      call require(info == 0, 'indexed testcase')
      call require(maxval(abs(named%a - indexed%a)) == 0.0_dp, 'named/indexed A parity')
      call require(maxval(abs(named%b - indexed%b)) == 0.0_dp, 'named/indexed B parity')
      call require(maxval(abs(named%sigma - indexed%sigma)) == 0.0_dp, 'named/indexed sigma parity')

      rho_model = varmapack_testcase('rho', p=2, q=1, r=2, rho=0.8_dp, info=info)
      call require(info == 0, 'rho testcase')
      call rho_model%specrad(radius, info)
      call require(info == 0 .and. abs(radius - 0.8_dp) < 2.0e-6_dp, 'rho testcase radius')

      rng1 = randompack_rng(seed=321)
      rng2 = randompack_rng(seed=321)
      random1 = varmapack_testcase('random', p=2, q=1, r=2, rng=rng1, info=info)
      call require(info == 0, 'random testcase 1')
      random2 = varmapack_testcase('random', p=2, q=1, r=2, rng=rng2, info=info)
      call require(info == 0, 'random testcase 2')
      call require(maxval(abs(random1%a - random2%a)) == 0.0_dp, 'random testcase reproducibility')
   end subroutine test_testcases

   subroutine test_simulation()
      type(varmapack_model_type) :: white, arma, mean_model, singular_model, unstable
      type(randompack_rng_type) :: rng1, rng2
      real(dp), allocatable :: a(:, :, :), b(:, :, :), x1(:, :, :), x2(:, :, :), e1(:, :, :), e2(:, :, :)
      real(dp), allocatable :: x0(:, :, :), mu(:, :)
      real(dp) :: sigma2(2, 2)
      integer :: info, t

      white = varmapack_model(reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), info=info)
      rng1 = randompack_rng(seed=17)
      call white%sim(10, x1, e1, info, nrep=2, rng=rng1)
      call require(info == 0, 'white-noise simulation')
      call require(maxval(abs(x1 - e1)) == 0.0_dp, 'white-noise X equals shocks')

      allocate(a(1, 1, 1), b(1, 1, 1))
      a = 0.4_dp
      b = 0.2_dp
      arma = varmapack_model(reshape([1.0_dp], [1, 1]), a=a, b=b, info=info)
      rng1 = randompack_rng(seed=123)
      rng2 = randompack_rng(seed=123)
      call arma%sim(20, x1, e1, info, nrep=3, rng=rng1)
      call require(info == 0, 'stationary ARMA simulation 1')
      call arma%sim(20, x2, e2, info, nrep=3, rng=rng2)
      call require(info == 0, 'stationary ARMA simulation 2')
      call require(maxval(abs(x1 - x2)) == 0.0_dp .and. maxval(abs(e1 - e2)) == 0.0_dp, 'ARMA RNG reproducibility')

      allocate(mu(2, 3))
      mu = reshape([10.0_dp, 20.0_dp, 11.0_dp, 21.0_dp, 12.0_dp, 22.0_dp], [2, 3])
      mean_model = varmapack_model(reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), mu=mu, info=info)
      rng1 = randompack_rng(seed=22)
      call mean_model%sim(5, x1, e1, info, nrep=2, rng=rng1)
      call require(info == 0, 'mean-path simulation')
      do t = 1, 5
         call require(maxval(abs(x1(:, t, 1) - e1(:, t, 1) - mu(:, min(t, 3)))) < 1.0e-13_dp, 'mean path')
      end do

      sigma2 = 1.0_dp
      singular_model = varmapack_model(sigma2, info=info)
      rng1 = randompack_rng(seed=23)
      call singular_model%sim(6, x1, e1, info, nrep=2, rng=rng1)
      call require(info == 0, 'singular covariance simulation')
      call require(maxval(abs(x1(1, :, :) - x1(2, :, :))) < 1.0e-13_dp, 'singular covariance equality')

      a = 1.25_dp
      b = 0.5_dp
      unstable = varmapack_model(reshape([1.0_dp], [1, 1]), a=a, b=b, info=info)
      allocate(x0(1, 3, 1))
      x0(1, :, 1) = [2.0_dp, 3.0_dp, 4.0_dp]
      rng1 = randompack_rng(seed=111)
      call unstable%sim(6, x1, e1, info, x0=x0, rng=rng1)
      call require(info == 0, 'nonstationary supplied-start simulation')
      call require(maxval(abs(x1(1, 1:3, 1) - x0(1, :, 1))) < 1.0e-13_dp, 'startup observations preserved')
      do t = 4, 6
         call require(abs(x1(1, t, 1) - (1.25_dp * x1(1, t - 1, 1) + e1(1, t, 1) + &
            0.5_dp * e1(1, t - 1, 1))) < 1.0e-12_dp, 'nonstationary recursion')
      end do
   end subroutine test_simulation

   subroutine test_varmax()
      type(varmapack_model_type) :: model
      type(randompack_rng_type) :: rng
      real(dp), allocatable :: a(:, :, :), b(:, :, :), c(:, :, :), x0(:, :, :), z(:, :, :), x(:, :, :), e(:, :, :)
      real(dp) :: expected
      integer :: info, j, t

      allocate(a(1, 1, 1), b(1, 1, 1), c(1, 1, 2), x0(1, 1, 1), z(1, 6, 1))
      a = 0.6_dp
      b = 0.2_dp
      c(1, 1, :) = [0.4_dp, -0.1_dp]
      x0 = 0.5_dp
      z(1, :, 1) = [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp, 3.0_dp, -3.0_dp]
      model = varmapack_model(reshape([1.0_dp], [1, 1]), a=a, b=b, c=c, info=info)
      call require(info == 0, 'construct VARMAX')
      rng = randompack_rng(seed=109)
      call model%sim(6, x, e, info, nrep=2, x0=x0, z=z, rng=rng)
      call require(info == 0, 'VARMAX simulation')
      do j = 1, 2
         do t = 2, 6
            expected = 0.6_dp * x(1, t - 1, j) + e(1, t, j) + 0.2_dp * e(1, t - 1, j) + &
               0.4_dp * z(1, t, 1) - 0.1_dp * z(1, t - 1, 1)
            call require(abs(x(1, t, j) - expected) < 1.0e-12_dp, 'VARMAX recurrence')
         end do
      end do

      deallocate(a, b, c, x0, z)
      allocate(c(1, 1, 1), z(1, 5, 1))
      c = 0.5_dp
      z(1, :, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
      model = varmapack_model(reshape([1.0_dp], [1, 1]), c=c, info=info)
      call require(info == 0, 'construct zero-start VARMAX')
      rng = randompack_rng(seed=110)
      call model%sim(5, x, e, info, nrep=2, z=z, rng=rng)
      call require(info == 0, 'zero-start VARMAX simulation')
      do j = 1, 2
         call require(maxval(abs(x(1, :, j) - e(1, :, j) - 0.5_dp * z(1, :, 1))) < 1.0e-12_dp, &
            'zero-start VARMAX recurrence')
      end do
   end subroutine test_varmax

end program test_varmapack
