program test_polca
   use polca, only : dp, polca_model, polca_simulation, polca_fit, polca_update_prior, &
        polca_item_likelihood, polca_posterior, polca_predcell, polca_entropy, &
        polca_reorder_probs, polca_table_1d, polca_table_2d, polca_rmulti, &
        polca_simdata, polca_coef, polca_vcov
   implicit none

   call test_rmulti()
   call test_prior()
   call test_one_class()
   call test_two_class()
   call test_regression()
   call test_simulation()
   print '(a)', 'All poLCA tests passed.'

contains

   subroutine test_rmulti()
      real(dp) :: p(3, 3), u(3)
      integer :: cat(3)

      p = reshape([0.2_dp, 0.1_dp, 0.7_dp, &
                   0.3_dp, 0.8_dp, 0.1_dp, &
                   0.5_dp, 0.1_dp, 0.2_dp], [3, 3], order=[2, 1])
      u = [0.10_dp, 0.50_dp, 0.95_dp]
      call polca_rmulti(p, u, cat)
      call assert_true(all(cat == [1, 2, 3]), 'rmulti categories')
   end subroutine test_rmulti

   subroutine test_prior()
      real(dp) :: coeff(2, 2), x(2, 2), prior(2, 3)
      real(dp) :: e2, e3, den

      coeff(:, 1) = [log(2.0_dp), 0.5_dp]
      coeff(:, 2) = [log(3.0_dp), -0.25_dp]
      x = reshape([1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp], [2, 2], order=[2, 1])
      call polca_update_prior(coeff, x, prior)
      e2 = 2.0_dp
      e3 = 3.0_dp
      den = 1.0_dp + e2 + e3
      call assert_close(prior(1, 1), 1.0_dp / den, 1.0e-12_dp, 'prior reference')
      call assert_close(prior(1, 2), e2 / den, 1.0e-12_dp, 'prior class 2')
      call assert_close(prior(1, 3), e3 / den, 1.0e-12_dp, 'prior class 3')
      call assert_close(sum(prior(2, :)), 1.0_dp, 1.0e-12_dp, 'prior row sum')
   end subroutine test_prior

   subroutine test_one_class()
      integer :: y(6, 2), k(2), cells(2, 2)
      real(dp) :: start(1, 3, 2), prob(2), post(6, 1), tab1(2), tab2(2, 3)
      real(dp), allocatable :: coef(:, :), cov(:, :)
      type(polca_model) :: model

      y = reshape([1, 1, 2, 2, 1, 0, &
                   1, 2, 3, 2, 3, 1], [6, 2])
      k = [2, 3]
      start = 0.0_dp
      start(1, 1:2, 1) = [0.5_dp, 0.5_dp]
      start(1, 1:3, 2) = [0.3_dp, 0.4_dp, 0.3_dp]
      call polca_fit(y, k, 1, model, probs_start=start, calc_se=.true.)
      call assert_close(model%probs(1, 1, 1), 3.0_dp / 5.0_dp, 1.0e-12_dp, 'one class item 1')
      call assert_close(model%probs(1, 2, 1), 2.0_dp / 5.0_dp, 1.0e-12_dp, 'one class item 1 category 2')
      call assert_close(model%probs(1, 1, 2), 2.0_dp / 6.0_dp, 1.0e-12_dp, 'one class item 2')
      call assert_true(model%has_se, 'one class standard errors')
      call assert_true(model%nobs_complete == 5, 'complete observation count')
      call polca_posterior(model, y, post)
      call assert_close(maxval(abs(post - 1.0_dp)), 0.0_dp, 1.0e-12_dp, 'one class posterior')
      cells = reshape([1, 2, 1, 3], [2, 2])
      call polca_predcell(model, cells, prob)
      call assert_close(prob(1), model%probs(1, 1, 1) * model%probs(1, 1, 2), 1.0e-12_dp, 'predcell 1')
      call assert_true(polca_entropy(model) > 0.0_dp, 'positive entropy')
      call polca_table_1d(model, 1, tab1)
      call assert_close(sum(tab1), real(model%n, dp), 1.0e-10_dp, 'one dimensional table total')
      call polca_table_2d(model, 1, 2, tab2)
      call assert_close(sum(tab2), real(model%n, dp), 1.0e-10_dp, 'two dimensional table total')
      call polca_coef(model, coef)
      call polca_vcov(model, cov)
      call assert_true(size(coef) == 0, 'one class coefficient accessor')
      call assert_true(size(cov) == 0, 'one class covariance accessor')
   end subroutine test_one_class

   subroutine test_two_class()
      integer, parameter :: n = 2400
      integer :: y(n, 3), k(3), cell(3), i, idx, rep
      real(dp) :: start(2, 2, 3), prob, share(2), reordered(2, 2, 3)
      type(polca_model) :: model

      k = 2
      start = 0.0_dp
      start(1, :, 1) = [0.90_dp, 0.10_dp]
      start(1, :, 2) = [0.80_dp, 0.20_dp]
      start(1, :, 3) = [0.85_dp, 0.15_dp]
      start(2, :, 1) = [0.15_dp, 0.85_dp]
      start(2, :, 2) = [0.25_dp, 0.75_dp]
      start(2, :, 3) = [0.20_dp, 0.80_dp]
      share = [0.4_dp, 0.6_dp]

      idx = 0
      do i = 0, 7
         cell = [1 + iand(i, 1), 1 + iand(ishft(i, -1), 1), 1 + iand(ishft(i, -2), 1)]
         prob = mixture_cell(start, share, cell)
         do rep = 1, nint(prob * real(n, dp))
            if (idx < n) then
               idx = idx + 1
               y(idx, :) = cell
            end if
         end do
      end do
      do while (idx < n)
         idx = idx + 1
         y(idx, :) = [2, 2, 2]
      end do

      call polca_fit(y, k, 2, model, probs_start=start, maxiter=1000, tol=1.0e-12_dp, calc_se=.true.)
      call assert_true(model%converged, 'two class convergence')
      call assert_true(model%loglik < 0.0_dp, 'two class finite likelihood')
      call assert_close(sum(model%class_share), 1.0_dp, 1.0e-12_dp, 'two class shares sum')
      call assert_true(minval(model%probs(:, 1:2, :)) > 0.0_dp, 'two class probabilities positive')
      call assert_true(all(model%class_share_se > 0.0_dp), 'two class share standard errors')
      call polca_reorder_probs(model%probs, [2, 1], reordered)
      call assert_close(reordered(1, 1, 1), model%probs(2, 1, 1), 1.0e-14_dp, 'reorder class')
   end subroutine test_two_class

   subroutine test_regression()
      integer, parameter :: n = 1200
      integer :: y(n, 3), k(3), i
      real(dp) :: probs(2, 2, 3), x(n, 2), beta(2, 1), prior(n, 2)
      real(dp) :: u(n), pmat(n, 2), start(2, 2, 3)
      integer :: cls(n)
      type(polca_model) :: model

      k = 2
      probs(1, :, 1) = [0.90_dp, 0.10_dp]
      probs(1, :, 2) = [0.85_dp, 0.15_dp]
      probs(1, :, 3) = [0.80_dp, 0.20_dp]
      probs(2, :, 1) = [0.10_dp, 0.90_dp]
      probs(2, :, 2) = [0.20_dp, 0.80_dp]
      probs(2, :, 3) = [0.15_dp, 0.85_dp]
      beta(:, 1) = [-0.25_dp, 1.1_dp]
      do i = 1, n
         x(i, 1) = 1.0_dp
         x(i, 2) = -1.5_dp + 3.0_dp * real(i - 1, dp) / real(n - 1, dp)
      end do
      call polca_update_prior(beta, x, prior)
      call deterministic_uniform(n, 17, u)
      call polca_rmulti(prior, u, cls)
      do i = 1, 3
         pmat(:, 1) = probs(cls, 1, i)
         pmat(:, 2) = probs(cls, 2, i)
         call deterministic_uniform(n, 40 + i, u)
         call polca_rmulti(pmat, u, y(:, i))
      end do
      start = probs
      call polca_fit(y, k, 2, model, x=x, probs_start=start, maxiter=700, tol=1.0e-10_dp, calc_se=.true.)
      call assert_true(model%has_covariates, 'regression covariate flag')
      call assert_true(model%converged, 'regression convergence')
      call assert_true(size(model%coeff, 1) == 2, 'regression coefficient rows')
      call assert_true(model%coeff(2, 1) > 0.5_dp, 'regression slope direction')
      call assert_true(all(model%coeff_se > 0.0_dp), 'regression coefficient standard errors')
   end subroutine test_regression

   subroutine test_simulation()
      real(dp) :: probs(2, 2, 2), share(2)
      type(polca_simulation) :: sim

      probs(1, :, 1) = [0.8_dp, 0.2_dp]
      probs(1, :, 2) = [0.7_dp, 0.3_dp]
      probs(2, :, 1) = [0.2_dp, 0.8_dp]
      probs(2, :, 2) = [0.3_dp, 0.7_dp]
      share = [0.35_dp, 0.65_dp]
      call polca_simdata(400, probs, sim, class_share=share, missing_fraction=0.1_dp, seed=77)
      call assert_true(size(sim%y, 1) == 400, 'simulation row count')
      call assert_true(all(sim%true_class >= 1 .and. sim%true_class <= 2), 'simulation class range')
      call assert_true(count(sim%y == 0) > 0, 'simulation missing data')
      call assert_close(sum(sim%class_share), 1.0_dp, 1.0e-12_dp, 'simulation realized shares')
   end subroutine test_simulation

   pure real(dp) function mixture_cell(probs, share, cell) result(p)
      real(dp), intent(in) :: probs(:, :, :) !! Class-by-category-by-item probabilities for a test mixture.
      real(dp), intent(in) :: share(:) !! Class shares for the test mixture.
      integer, intent(in) :: cell(:) !! Complete response pattern to evaluate.
      real(dp) :: pr
      integer :: j, r

      p = 0.0_dp
      do r = 1, size(share)
         pr = share(r)
         do j = 1, size(cell)
            pr = pr * probs(r, cell(j), j)
         end do
         p = p + pr
      end do
   end function mixture_cell

   pure subroutine deterministic_uniform(n, offset, u)
      use, intrinsic :: iso_fortran_env, only : int64
      integer, intent(in) :: n !! Number of deterministic pseudo-uniform values to generate.
      integer, intent(in) :: offset !! Integer offset selecting an independent deterministic stream.
      real(dp), intent(out) :: u(:) !! Deterministic values in the open unit interval.
      integer(int64) :: state
      integer :: i

      state = modulo(123457_int64 + 104729_int64 * int(offset, int64), 2147483646_int64) + 1_int64
      do i = 1, n
         state = modulo(48271_int64 * state, 2147483647_int64)
         u(i) = real(state, dp) / 2147483647.0_dp
      end do
   end subroutine deterministic_uniform

   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual !! Computed scalar value.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tol !! Maximum allowed absolute error.
      character(*), intent(in) :: label !! Short assertion description printed on failure.

      if (abs(actual - expected) > tol) then
         write(*, '(a,2(1x,es24.16))') trim(label), actual, expected
         error stop 'assert_close failed'
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Condition that must be true for the test to pass.
      character(*), intent(in) :: label !! Short assertion description printed on failure.

      if (.not. condition) then
         write(*, '(a)') trim(label)
         error stop 'assert_true failed'
      end if
   end subroutine assert_true

end program test_polca
