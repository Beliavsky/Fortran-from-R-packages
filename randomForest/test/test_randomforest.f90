program test_randomforest
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
   use randomforest
   implicit none

   call test_classification_numeric()
   call test_classification_categorical()
   call test_classification_stratified()
   call test_regression_numeric()
   call test_regression_categorical()
   call test_utilities()
   call test_unsupervised()
   call test_imputation()
   call test_tuning()
   call test_rfcv()
   call test_mds()
   print '(a)', 'all randomForest tests passed'

contains

   subroutine test_classification_numeric()
      integer, parameter :: n = 20, p = 2
      real(dp) :: x(n, p), prob(n, 2)
      integer :: y(n), pred(n), i, status
      type(rf_options) :: opt
      type(rf_classification_forest) :: forest
      character(len=256) :: message

      do i = 1, n / 2
         x(i, 1) = real(i - 11, dp)
         x(i, 2) = real(mod(i, 3), dp)
         y(i) = 1
      end do
      do i = n / 2 + 1, n
         x(i, 1) = real(i - 10, dp)
         x(i, 2) = real(mod(i, 3), dp)
         y(i) = 2
      end do
      opt%ntree = 31
      opt%mtry = 2
      opt%nodesize = 1
      opt%importance = .true.
      opt%proximity = .true.
      opt%keep_inbag = .true.
      opt%seed = 17
      call fit_classification(x, y, forest, options=opt, status=status, message=message)
      call assert_true(status == 0, 'numeric classification fit: '//trim(message))
      call predict_classification(forest, x, pred, probabilities=prob)
      call assert_true(all(pred == y), 'numeric classification training predictions')
      call assert_true(forest%importance_gini(1) > forest%importance_gini(2), 'numeric classification Gini importance')
      call assert_close(maxval(abs(forest%proximity - transpose(forest%proximity))), 0.0_dp, 1.0e-12_dp, &
         'classification proximity symmetry')
      call assert_close(minval([(forest%proximity(i, i), i = 1, n)]), 1.0_dp, 1.0e-12_dp, &
         'classification proximity diagonal')
      call assert_true(all(abs(sum(prob, dim=2) - 1.0_dp) < 1.0e-12_dp), 'classification probabilities sum to one')
   end subroutine test_classification_numeric

   subroutine test_classification_categorical()
      integer, parameter :: n = 16
      real(dp) :: x(n, 1)
      integer :: y(n), pred(n), cats(1), i, status
      type(rf_options) :: opt
      type(rf_classification_forest) :: forest
      character(len=256) :: message

      do i = 1, n
         x(i, 1) = real(1 + mod(i - 1, 4), dp)
         if (nint(x(i, 1)) <= 2) then
            y(i) = 1
         else
            y(i) = 2
         end if
      end do
      cats = 4
      opt%ntree = 9
      opt%mtry = 1
      opt%nodesize = 1
      opt%replace = .false.
      opt%sample_size = n
      opt%seed = 31
      call fit_classification(x, y, forest, ncat=cats, options=opt, status=status, message=message)
      call assert_true(status == 0, 'categorical classification fit: '//trim(message))
      call predict_classification(forest, x, pred)
      call assert_true(all(pred == y), 'categorical classification predictions')
      call assert_true(forest%trees(1)%split_var(1) == 1, 'categorical root split variable')
   end subroutine test_classification_categorical

   subroutine test_classification_stratified()
      integer, parameter :: n = 12
      real(dp) :: x(n, 1)
      integer :: y(n), strata(n), strata_size(2), i, status
      type(rf_options) :: opt
      type(rf_classification_forest) :: forest
      character(len=256) :: message

      do i = 1, n
         x(i, 1) = real(i, dp)
         y(i) = merge(1, 2, i <= 9)
      end do
      strata = y
      strata_size = [3, 3]
      opt%ntree = 5
      opt%mtry = 1
      opt%keep_inbag = .true.
      opt%seed = 101
      call fit_classification(x, y, forest, options=opt, strata=strata, strata_sample_size=strata_size, &
         status=status, message=message)
      call assert_true(status == 0, 'stratified classification fit: '//trim(message))
      do i = 1, opt%ntree
         call assert_true(sum(forest%inbag(:, i)) == 6, 'stratified inbag sample size')
         call assert_true(sum(forest%inbag(1:9, i)) == 3, 'stratified class one size')
         call assert_true(sum(forest%inbag(10:12, i)) == 3, 'stratified class two size')
      end do
   end subroutine test_classification_stratified

   subroutine test_regression_numeric()
      integer, parameter :: n = 15
      real(dp) :: x(n, 1), y(n), pred(n), grid(3), pd(3)
      integer :: i, status
      type(rf_options) :: opt
      type(rf_regression_forest) :: forest
      character(len=256) :: message

      do i = 1, n
         x(i, 1) = real(i - 8, dp)
         y(i) = 2.0_dp + 3.0_dp * x(i, 1)
      end do
      opt%ntree = 5
      opt%mtry = 1
      opt%nodesize = 1
      opt%replace = .false.
      opt%sample_size = n
      opt%seed = 43
      call fit_regression(x, y, forest, options=opt, status=status, message=message)
      call assert_true(status == 0, 'numeric regression fit: '//trim(message))
      call predict_regression(forest, x, pred)
      call assert_close(maxval(abs(pred - y)), 0.0_dp, 1.0e-10_dp, 'numeric regression exact full-sample tree')
      grid = [-4.0_dp, 0.0_dp, 4.0_dp]
      call partial_dependence_regression(forest, x, 1, grid, pd)
      call assert_true(pd(1) < pd(2) .and. pd(2) < pd(3), 'regression partial dependence monotone')
   end subroutine test_regression_numeric

   subroutine test_regression_categorical()
      integer, parameter :: n = 12
      real(dp) :: x(n, 1), y(n), pred(n)
      integer :: cats(1), i, status
      type(rf_options) :: opt
      type(rf_regression_forest) :: forest
      character(len=256) :: message

      do i = 1, n
         x(i, 1) = real(1 + mod(i - 1, 4), dp)
         if (nint(x(i, 1)) <= 2) then
            y(i) = 10.0_dp
         else
            y(i) = 20.0_dp
         end if
      end do
      cats = 4
      opt%ntree = 3
      opt%mtry = 1
      opt%nodesize = 1
      opt%replace = .false.
      opt%sample_size = n
      opt%seed = 59
      call fit_regression(x, y, forest, ncat=cats, options=opt, status=status, message=message)
      call assert_true(status == 0, 'categorical regression fit: '//trim(message))
      call predict_regression(forest, x, pred)
      call assert_close(maxval(abs(pred - y)), 0.0_dp, 1.0e-10_dp, 'categorical regression predictions')
   end subroutine test_regression_categorical

   subroutine test_utilities()
      real(dp) :: votes(3, 2), margins(3), prox(4, 4), score(4)
      real(dp) :: a(4, 2), centers(2, 2)
      integer :: observed(3), cls(4)

      votes = reshape([8.0_dp, 3.0_dp, 5.0_dp, 2.0_dp, 7.0_dp, 5.0_dp], shape(votes))
      observed = [1, 2, 1]
      call classification_margin(votes, observed, margins)
      call assert_close(margins(1), 0.6_dp, 1.0e-12_dp, 'margin first case')

      prox = reshape([1.0_dp, 0.8_dp, 0.1_dp, 0.1_dp, &
                      0.8_dp, 1.0_dp, 0.1_dp, 0.1_dp, &
                      0.1_dp, 0.1_dp, 1.0_dp, 0.7_dp, &
                      0.1_dp, 0.1_dp, 0.7_dp, 1.0_dp], shape(prox))
      cls = [1, 1, 2, 2]
      call outlier_scores(prox, cls, score)
      call assert_true(all(abs(score) < 1.0e-12_dp), 'symmetric outlier scores')

      a = reshape([1.0_dp, 2.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 4.0_dp, &
                   1.0_dp, 1.0_dp, 9.0_dp, 9.0_dp], shape(a))
      call roughfix_numeric(a)
      call assert_close(a(3, 1), 2.0_dp, 1.0e-12_dp, 'roughfix median')
      call class_centers(a, cls, prox, centers, n_neighbors=2)
      call assert_true(all(ieee_is_finite(centers)), 'class centers finite')
   end subroutine test_utilities

   subroutine test_unsupervised()
      integer, parameter :: n = 10
      real(dp) :: x(n, 2)
      integer :: i, status
      type(rf_options) :: opt
      type(rf_classification_forest) :: forest
      character(len=256) :: message

      do i = 1, n
         x(i, 1) = real(i, dp)
         x(i, 2) = real(mod(i, 3), dp)
      end do
      opt%ntree = 7
      opt%mtry = 1
      opt%proximity = .true.
      opt%keep_inbag = .true.
      opt%seed = 71
      call fit_unsupervised(x, forest, options=opt, status=status, message=message)
      call assert_true(status == 0, 'unsupervised fit: '//trim(message))
      call assert_true(size(forest%trees) == 7, 'unsupervised tree count')
      call assert_close(maxval(abs(forest%proximity - transpose(forest%proximity))), 0.0_dp, 1.0e-12_dp, &
         'unsupervised proximity symmetry')
   end subroutine test_unsupervised


   subroutine test_imputation()
      integer, parameter :: n = 12
      real(dp) :: x(n, 2), imputed(n, 2)
      integer :: y(n), cats(2), i, status
      type(rf_options) :: opt
      character(len=256) :: message

      do i = 1, n
         x(i, 1) = real(i, dp)
         x(i, 2) = real(1 + mod(i - 1, 3), dp)
         y(i) = merge(1, 2, i <= n / 2)
      end do
      x(3, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
      x(9, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
      cats = [1, 3]
      opt%seed = 211
      opt%mtry = 2
      call rf_impute_classification(x, y, imputed, ncat=cats, iterations=2, ntree=15, options=opt, &
         status=status, message=message)
      call assert_true(status == 0, 'classification proximity imputation: '//trim(message))
      call assert_true(all(ieee_is_finite(imputed)), 'imputation returns finite values')
      call assert_true(nint(imputed(9, 2)) >= 1 .and. nint(imputed(9, 2)) <= 3, 'categorical imputation is valid')
      call assert_close(imputed(1, 1), x(1, 1), 0.0_dp, 'imputation leaves observed values unchanged')
   end subroutine test_imputation

   subroutine test_tuning()
      integer, parameter :: n = 24, p = 4
      real(dp) :: x(n, p), yr(n)
      integer :: yc(n), i, status, best_mtry
      integer, allocatable :: mtry_values(:)
      real(dp), allocatable :: errors(:)
      type(rf_options) :: opt
      character(len=256) :: message

      do i = 1, n
         x(i, 1) = real(i - 12, dp)
         x(i, 2) = real(mod(i, 5), dp)
         x(i, 3) = real(mod(2 * i, 7), dp)
         x(i, 4) = real(mod(3 * i, 11), dp)
         yc(i) = merge(1, 2, x(i, 1) < 0.0_dp)
         yr(i) = 1.0_dp + 2.0_dp * x(i, 1) + 0.1_dp * x(i, 2)
      end do
      opt%seed = 307
      call tune_classification_mtry(x, yc, mtry_values, errors, best_mtry, base_options=opt, &
         mtry_start=2, ntree_try=15, improve=0.0_dp, status=status, message=message)
      call assert_true(status == 0, 'classification mtry tuning: '//trim(message))
      call assert_true(size(mtry_values) >= 2 .and. all(mtry_values >= 1) .and. all(mtry_values <= p), &
         'classification tuning mtry range')
      call assert_true(best_mtry >= 1 .and. best_mtry <= p, 'classification tuning best mtry')
      call assert_true(all(ieee_is_finite(errors)), 'classification tuning finite errors')

      call tune_regression_mtry(x, yr, mtry_values, errors, best_mtry, base_options=opt, &
         mtry_start=1, ntree_try=15, improve=0.0_dp, status=status, message=message)
      call assert_true(status == 0, 'regression mtry tuning: '//trim(message))
      call assert_true(best_mtry >= 1 .and. best_mtry <= p, 'regression tuning best mtry')
      call assert_true(all(ieee_is_finite(errors)), 'regression tuning finite errors')
   end subroutine test_tuning

   subroutine test_rfcv()
      integer, parameter :: n = 30, p = 4
      real(dp) :: x(n, p), yr(n)
      integer :: yc(n), i, status
      integer, allocatable :: nvar(:), cpred(:,:)
      real(dp), allocatable :: errors(:), rpred(:,:)
      type(rf_options) :: opt
      character(len=256) :: message

      do i = 1, n
         x(i, 1) = real(i - 15, dp)
         x(i, 2) = real(mod(i, 4), dp)
         x(i, 3) = real(mod(i, 7), dp)
         x(i, 4) = real(mod(i, 9), dp)
         yc(i) = merge(1, 2, x(i, 1) < 0.0_dp)
         yr(i) = 4.0_dp - 1.5_dp * x(i, 1) + 0.2_dp * x(i, 2)
      end do
      opt%ntree = 13
      opt%seed = 401
      call rfcv_classification(x, yc, nvar, errors, cpred, cv_fold=3, step=0.5_dp, &
         base_options=opt, status=status, message=message)
      call assert_true(status == 0, 'classification rfcv: '//trim(message))
      call assert_true(nvar(1) == p .and. nvar(size(nvar)) == 1, 'classification rfcv feature path')
      call assert_true(all(errors >= 0.0_dp .and. errors <= 1.0_dp), 'classification rfcv error range')
      call assert_true(all(cpred >= 1 .and. cpred <= 2), 'classification rfcv predictions')

      call rfcv_regression(x, yr, nvar, errors, rpred, cv_fold=3, step=0.5_dp, &
         base_options=opt, status=status, message=message)
      call assert_true(status == 0, 'regression rfcv: '//trim(message))
      call assert_true(nvar(1) == p .and. nvar(size(nvar)) == 1, 'regression rfcv feature path')
      call assert_true(all(errors >= 0.0_dp .and. ieee_is_finite(errors)), 'regression rfcv finite errors')
      call assert_true(all(ieee_is_finite(rpred)), 'regression rfcv predictions finite')
   end subroutine test_rfcv

   subroutine test_mds()
      real(dp) :: proximity(3, 3)
      real(dp), allocatable :: coords(:,:), evals(:)
      integer :: status
      character(len=256) :: message

      proximity = reshape([1.0_dp, 0.5_dp, 0.0_dp, &
                           0.5_dp, 1.0_dp, 0.5_dp, &
                           0.0_dp, 0.5_dp, 1.0_dp], shape(proximity))
      call mds_coordinates(proximity, 1, coords, evals, status=status, message=message)
      call assert_true(status == 0, 'classical MDS: '//trim(message))
      call assert_true(size(coords, 1) == 3 .and. size(coords, 2) == 1, 'MDS coordinate shape')
      call assert_close(abs(coords(1, 1) - coords(3, 1)), 1.0_dp, 1.0e-10_dp, 'MDS end-point distance')
      call assert_close(abs(coords(1, 1) - coords(2, 1)), 0.5_dp, 1.0e-10_dp, 'MDS adjacent distance')
      call assert_true(evals(1) > 0.0_dp, 'MDS leading eigenvalue positive')
   end subroutine test_mds

   subroutine assert_true(condition, name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not. condition) then
         print '(a)', 'FAIL: '//trim(name)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, name)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: name
      if (abs(actual - expected) > tolerance) then
         print '(a,2es24.14)', 'FAIL: '//trim(name)//' actual/expected: ', actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_randomforest
