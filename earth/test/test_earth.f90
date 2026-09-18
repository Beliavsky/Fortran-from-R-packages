program test_earth
   use earth_api, only : dp, earth_model, earth_fit, earth_predict, earth_model_matrix, &
                         earth_residuals, earth_coefficients, earth_hatvalues, &
                         earth_deviance, earth_extract_aic, earth_evimp, expand_bpairs, &
                         contr_earth_response, effective_nbr_of_params
   implicit none

   call test_contrasts
   call test_binomial_pairs
   call test_weighted_multivariate_linear_fit
   call test_piecewise_fit_and_importance
   call test_upstream_trees_regression
   call test_linear_interaction
   print '(a)', 'All earth tests passed.'

contains

   subroutine test_contrasts
      real(dp), allocatable :: contrasts(:, :)
      real(dp) :: expected(3, 3)

      expected = 0.0_dp
      expected(1, 1) = 1.0_dp
      expected(2, 2) = 1.0_dp
      expected(3, 3) = 1.0_dp
      call contr_earth_response(3, contrasts)
      call check(all(shape(contrasts) == [3, 3]), 'contrast matrix shape')
      call check(maxval(abs(contrasts - expected)) < 1.0e-14_dp, 'identity response contrasts')
   end subroutine test_contrasts

   subroutine test_binomial_pairs
      real(dp) :: x(3, 2)
      real(dp) :: pairs(3, 2)
      real(dp), allocatable :: xlong(:, :)
      logical, allocatable :: ylong(:)
      integer, allocatable :: index(:)
      logical :: expected_y(7)
      integer :: ierr

      x = reshape([1.0_dp, 2.0_dp, 3.0_dp, 10.0_dp, 20.0_dp, 30.0_dp], [3, 2])
      pairs = reshape([2.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp], [3, 2])
      expected_y = [.false., .true., .true., .false., .false., .false., .true.]
      call expand_bpairs(x, pairs, xlong, ylong, index, ierr)
      call check(ierr == 0, 'expand_bpairs success')
      call check(all(shape(xlong) == [7, 2]), 'expanded binomial predictor shape')
      call check(all(index == [1, 4, 5]), 'binomial-pair first-row indices')
      call check(all(ylong .eqv. expected_y), 'binomial-pair false/true expansion order')
      call check(maxval(abs(xlong(1:3, :) - spread(x(1, :), 1, 3))) < 1.0e-14_dp, &
                 'first binomial predictor row replication')
      call check(maxval(abs(xlong(4, :) - x(2, :))) < 1.0e-14_dp, &
                 'zero-zero binomial pair retained as one false row')
   end subroutine test_binomial_pairs

   subroutine test_weighted_multivariate_linear_fit
      integer, parameter :: n = 41
      real(dp) :: x(n, 1)
      real(dp) :: y(n, 2)
      real(dp) :: weights(n)
      logical :: linpreds(1)
      type(earth_model) :: model
      real(dp), allocatable :: prediction(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: leverages(:)
      real(dp), allocatable :: bx(:, :)
      real(dp) :: effective_df
      real(dp) :: criterion
      integer :: i
      integer :: ierr
      character(len=:), allocatable :: message

      do i = 1, n
         x(i, 1) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
         y(i, 1) = 2.0_dp + 3.0_dp * x(i, 1)
         y(i, 2) = -1.0_dp + 0.5_dp * x(i, 1)
         weights(i) = 1.0_dp + 0.05_dp * real(mod(i, 7), dp)
      end do
      linpreds = .true.
      call earth_fit(x, y, model, degree=1, linpreds=linpreds, weights=weights, &
                     prune=.false., ierr=ierr, message=message)
      call check(ierr == 0, 'weighted multiple-response linear earth_fit')
      call check(model%n_selected == 2, 'linear fit selects intercept and one linear term')
      call earth_predict(model, x, prediction, ierr)
      call check(ierr == 0, 'linear prediction success')
      call check(maxval(abs(prediction - y)) < 2.0e-12_dp, 'weighted multiple-response exact fit')

      call earth_residuals(model, residuals, ierr)
      call check(ierr == 0, 'residual extraction success')
      call check(maxval(abs(residuals)) < 2.0e-12_dp, 'stored residuals for exact fit')
      call earth_coefficients(model, coefficients, ierr)
      call check(ierr == 0, 'coefficient extraction success')
      call check(all(shape(coefficients) == [2, 2]), 'multiple-response coefficient shape')
      call earth_hatvalues(model, leverages, ierr)
      call check(ierr == 0, 'leverage extraction success')
      call check(all(leverages >= -1.0e-12_dp .and. leverages <= 1.0_dp + 1.0e-12_dp), &
                 'leverage range')
      call check(abs(sum(leverages) - 2.0_dp) < 1.0e-10_dp, 'leverage sum equals selected rank')

      call earth_model_matrix(model, x, bx, ierr)
      call check(ierr == 0, 'selected model matrix construction')
      call check(maxval(abs(matmul(bx, coefficients) - prediction)) < 1.0e-12_dp, &
                 'model matrix and coefficients reproduce predictions')
      call check(abs(earth_deviance(model) - model%rss) < 1.0e-14_dp, 'earth deviance is RSS')
      call earth_extract_aic(model, effective_df, criterion, ierr)
      call check(ierr == 0, 'extractAIC-compatible criterion success')
      call check(abs(effective_df - effective_nbr_of_params(model%n_selected, model%penalty)) < 1.0e-14_dp, &
                 'effective parameter count')
      call check(abs(criterion - model%gcv) < 1.0e-14_dp, 'extractAIC-compatible criterion is GCV')
   end subroutine test_weighted_multivariate_linear_fit

   subroutine test_piecewise_fit_and_importance
      integer, parameter :: n = 81
      real(dp) :: x(n, 2)
      real(dp) :: y(n, 1)
      type(earth_model) :: model
      real(dp), allocatable :: prediction(:, :)
      real(dp), allocatable :: nsubsets(:)
      real(dp), allocatable :: gcv_importance(:)
      real(dp), allocatable :: rss_importance(:)
      logical, allocatable :: used(:)
      integer :: i
      integer :: ierr
      character(len=:), allocatable :: message

      do i = 1, n
         x(i, 1) = real(i - 1, dp) / real(n - 1, dp)
         x(i, 2) = sin(0.73_dp * real(i, dp))
         y(i, 1) = 1.0_dp + 3.0_dp * max(0.0_dp, x(i, 1) - 0.35_dp) &
                   - 2.0_dp * max(0.0_dp, 0.725_dp - x(i, 1))
      end do
      call earth_fit(x, y, model, degree=1, nk=21, penalty=2.0_dp, thresh=1.0e-10_dp, &
                     minspan=1, endspan=1, ierr=ierr, message=message)
      call check(ierr == 0, 'piecewise-linear earth_fit')
      call earth_predict(model, x, prediction, ierr)
      call check(ierr == 0, 'piecewise prediction success')
      call check(maxval(abs(prediction - y)) < 2.0e-11_dp, 'piecewise hinge regression exact fit')
      call check(model%rsq > 1.0_dp - 1.0e-12_dp, 'piecewise model R-squared')

      call earth_evimp(model, nsubsets, gcv_importance, rss_importance, used, ierr)
      call check(ierr == 0, 'earth variable importance success')
      call check(used(1) .and. .not. used(2), 'variable importance identifies the signal predictor')
      call check(abs(gcv_importance(1) - 100.0_dp) < 1.0e-10_dp, 'GCV importance normalization')
      call check(abs(rss_importance(1) - 100.0_dp) < 1.0e-10_dp, 'RSS importance normalization')
      call check(nsubsets(1) > 0.0_dp .and. abs(nsubsets(2)) < 1.0e-14_dp, 'subset-count importance')
   end subroutine test_piecewise_fit_and_importance

   subroutine test_upstream_trees_regression
      real(dp), parameter :: girth(31) = [ &
         8.3_dp, 8.6_dp, 8.8_dp, 10.5_dp, 10.7_dp, 10.8_dp, 11.0_dp, 11.0_dp, &
         11.1_dp, 11.2_dp, 11.3_dp, 11.4_dp, 11.4_dp, 11.7_dp, 12.0_dp, 12.9_dp, &
         12.9_dp, 13.3_dp, 13.7_dp, 13.8_dp, 14.0_dp, 14.2_dp, 14.5_dp, 16.0_dp, &
         16.3_dp, 17.3_dp, 17.5_dp, 17.9_dp, 18.0_dp, 18.0_dp, 20.6_dp ]
      real(dp), parameter :: height(31) = [ &
         70.0_dp, 65.0_dp, 63.0_dp, 72.0_dp, 81.0_dp, 83.0_dp, 66.0_dp, 75.0_dp, &
         80.0_dp, 75.0_dp, 79.0_dp, 76.0_dp, 76.0_dp, 69.0_dp, 75.0_dp, 74.0_dp, &
         85.0_dp, 86.0_dp, 71.0_dp, 64.0_dp, 78.0_dp, 80.0_dp, 74.0_dp, 72.0_dp, &
         77.0_dp, 81.0_dp, 82.0_dp, 80.0_dp, 80.0_dp, 80.0_dp, 87.0_dp ]
      real(dp), parameter :: volume(31) = [ &
         10.3_dp, 10.3_dp, 10.2_dp, 16.4_dp, 18.8_dp, 19.7_dp, 15.6_dp, 18.2_dp, &
         22.6_dp, 19.9_dp, 24.2_dp, 21.0_dp, 21.4_dp, 21.3_dp, 19.1_dp, 22.2_dp, &
         33.8_dp, 27.4_dp, 25.7_dp, 24.9_dp, 34.5_dp, 31.7_dp, 36.3_dp, 38.3_dp, &
         42.6_dp, 55.4_dp, 55.7_dp, 58.3_dp, 51.5_dp, 51.0_dp, 77.0_dp ]
      real(dp) :: x(31, 2)
      real(dp) :: y(31, 1)
      type(earth_model) :: model
      integer :: k
      integer :: term
      integer :: ierr
      logical :: found_intercept
      logical :: found_girth_right
      logical :: found_girth_left
      logical :: found_height_right
      character(len=:), allocatable :: message

      x(:, 1) = girth
      x(:, 2) = height
      y(:, 1) = volume
      call earth_fit(x, y, model, ierr=ierr, message=message)
      call check(ierr == 0, 'upstream trees earth regression fit')
      call check(model%n_selected == 4, 'upstream trees selected-term count')
      call check(abs(model%rss - 209.1138545666_dp) < 1.0e-8_dp, 'upstream trees RSS parity')
      call check(abs(model%rsq - 0.974202850859_dp) < 1.0e-10_dp, 'upstream trees R-squared parity')
      call check(abs(model%gcv - 11.2543914784_dp) < 1.0e-8_dp, 'upstream trees GCV parity')

      found_intercept = .false.
      found_girth_right = .false.
      found_girth_left = .false.
      found_height_right = .false.
      do k = 1, model%n_selected
         term = model%selected_terms(k)
         if (all(model%dirs(term, :) == 0)) then
            found_intercept = abs(model%coefficients(k, 1) - 29.0599535107_dp) < 1.0e-8_dp
         else if (all(model%dirs(term, :) == [1, 0])) then
            found_girth_right = abs(model%cuts(term, 1) - 14.2_dp) < 1.0e-12_dp .and. &
                                abs(model%coefficients(k, 1) - 6.2295143423_dp) < 1.0e-8_dp
         else if (all(model%dirs(term, :) == [-1, 0])) then
            found_girth_left = abs(model%cuts(term, 1) - 14.2_dp) < 1.0e-12_dp .and. &
                               abs(model%coefficients(k, 1) + 3.4198061519_dp) < 1.0e-8_dp
         else if (all(model%dirs(term, :) == [0, 1])) then
            found_height_right = abs(model%cuts(term, 2) - 75.0_dp) < 1.0e-12_dp .and. &
                                 abs(model%coefficients(k, 1) - 0.5813643840_dp) < 1.0e-8_dp
         end if
      end do
      call check(found_intercept, 'upstream trees intercept coefficient parity')
      call check(found_girth_right, 'upstream trees right Girth hinge parity')
      call check(found_girth_left, 'upstream trees left Girth hinge parity')
      call check(found_height_right, 'upstream trees Height hinge parity')
   end subroutine test_upstream_trees_regression

   subroutine test_linear_interaction
      integer, parameter :: n1 = 10
      integer, parameter :: n2 = 10
      integer, parameter :: n = n1 * n2
      real(dp) :: x(n, 2)
      real(dp) :: y(n, 1)
      logical :: linpreds(2)
      type(earth_model) :: model
      real(dp), allocatable :: prediction(:, :)
      integer :: i
      integer :: j
      integer :: k
      integer :: ierr
      character(len=:), allocatable :: message

      k = 0
      do j = 1, n2
         do i = 1, n1
            k = k + 1
            x(k, 1) = real(i - 1, dp) / real(n1 - 1, dp)
            x(k, 2) = real(j - 1, dp) / real(n2 - 1, dp)
            y(k, 1) = 1.0_dp + 4.0_dp * x(k, 1) * x(k, 2)
         end do
      end do
      linpreds = .true.
      call earth_fit(x, y, model, degree=2, nk=10, penalty=3.0_dp, thresh=1.0e-10_dp, &
                     linpreds=linpreds, prune=.false., ierr=ierr, message=message)
      call check(ierr == 0, 'degree-two linear interaction fit')
      call earth_predict(model, x, prediction, ierr)
      call check(ierr == 0, 'interaction prediction success')
      call check(maxval(abs(prediction - y)) < 3.0e-12_dp, 'degree-two interaction exact fit')
      call check(any(all(abs(model%dirs) == 2, dim=2)), 'interaction basis term is retained')
   end subroutine test_linear_interaction

   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion condition that must be true for the deterministic regression test to pass.
      character(len=*), intent(in) :: label !! Short description printed when the associated assertion fails.

      if (.not. condition) then
         write(*, '(a)') 'FAILED: ' // trim(label)
         error stop 1
      end if
   end subroutine check

end program test_earth
