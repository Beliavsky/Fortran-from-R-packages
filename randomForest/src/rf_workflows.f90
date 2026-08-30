! SPDX-License-Identifier: GPL-2.0-or-later
module rf_workflows
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp, i64
   use r_linalg, only : symmetric_eigen
   use r_quantiles, only : r_median
   use rf_rng, only : rf_rng_state, shuffle_int
   use rf_types, only : rf_options, rf_classification_forest, rf_regression_forest
   use rf_classification, only : fit_classification, predict_classification
   use rf_regression, only : fit_regression, predict_regression
   implicit none
   private

   public :: rf_impute_classification, rf_impute_regression
   public :: tune_classification_mtry, tune_regression_mtry
   public :: rfcv_classification, rfcv_regression
   public :: mds_coordinates

contains

   subroutine rf_impute_classification(x, y, imputed, ncat, iterations, ntree, options, status, message)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      real(dp), intent(out) :: imputed(:,:)
      integer, intent(in), optional :: ncat(:), iterations, ntree
      type(rf_options), intent(in), optional :: options
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(rf_options) :: opt
      type(rf_classification_forest) :: forest
      integer, allocatable :: cats(:)
      logical, allocatable :: missing(:,:)
      integer :: iter, niter, ntrees, rc
      character(len=256) :: msg

      call set_ok(status, message)
      if (size(imputed, 1) /= size(x, 1) .or. size(imputed, 2) /= size(x, 2) .or. size(y) /= size(x, 1)) then
         call set_error(1, 'rf_impute_classification: incompatible shapes', status, message)
         return
      end if
      call prepare_ncat(x, ncat, cats, rc)
      if (rc /= 0) then
         call set_error(2, 'rf_impute_classification: invalid categorical coding', status, message)
         return
      end if
      missing = ieee_is_nan(x)
      if (.not. any(missing)) then
         call set_error(3, 'rf_impute_classification: x contains no missing values', status, message)
         return
      end if
      call rough_initial_impute(x, cats, missing, imputed, rc)
      if (rc /= 0) then
         call set_error(4, 'rf_impute_classification: a predictor has no observed values', status, message)
         return
      end if

      niter = 5
      if (present(iterations)) niter = max(1, iterations)
      ntrees = 300
      if (present(ntree)) ntrees = max(1, ntree)
      opt = rf_options()
      if (present(options)) opt = options
      opt%ntree = ntrees
      opt%proximity = .true.
      opt%oob_proximity = .false.
      opt%keep_inbag = .false.
      do iter = 1, niter
         opt%seed = opt%seed + int(iter - 1, i64)
         call fit_classification(imputed, y, forest, ncat=cats, options=opt, status=rc, message=msg)
         if (rc /= 0) then
            call set_error(10 + rc, 'rf_impute_classification: forest fit failed: '//trim(msg), status, message)
            return
         end if
         call proximity_impute(x, missing, cats, forest%proximity, imputed)
      end do
   end subroutine rf_impute_classification

   subroutine rf_impute_regression(x, y, imputed, ncat, iterations, ntree, options, status, message)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), intent(out) :: imputed(:,:)
      integer, intent(in), optional :: ncat(:), iterations, ntree
      type(rf_options), intent(in), optional :: options
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(rf_options) :: opt
      type(rf_regression_forest) :: forest
      integer, allocatable :: cats(:)
      logical, allocatable :: missing(:,:)
      integer :: iter, niter, ntrees, rc
      character(len=256) :: msg

      call set_ok(status, message)
      if (size(imputed, 1) /= size(x, 1) .or. size(imputed, 2) /= size(x, 2) .or. size(y) /= size(x, 1)) then
         call set_error(1, 'rf_impute_regression: incompatible shapes', status, message)
         return
      end if
      call prepare_ncat(x, ncat, cats, rc)
      if (rc /= 0) then
         call set_error(2, 'rf_impute_regression: invalid categorical coding', status, message)
         return
      end if
      missing = ieee_is_nan(x)
      if (.not. any(missing)) then
         call set_error(3, 'rf_impute_regression: x contains no missing values', status, message)
         return
      end if
      call rough_initial_impute(x, cats, missing, imputed, rc)
      if (rc /= 0) then
         call set_error(4, 'rf_impute_regression: a predictor has no observed values', status, message)
         return
      end if

      niter = 5
      if (present(iterations)) niter = max(1, iterations)
      ntrees = 300
      if (present(ntree)) ntrees = max(1, ntree)
      opt = rf_options()
      if (present(options)) opt = options
      opt%ntree = ntrees
      opt%proximity = .true.
      opt%oob_proximity = .false.
      opt%keep_inbag = .false.
      do iter = 1, niter
         opt%seed = opt%seed + int(iter - 1, i64)
         call fit_regression(imputed, y, forest, ncat=cats, options=opt, status=rc, message=msg)
         if (rc /= 0) then
            call set_error(10 + rc, 'rf_impute_regression: forest fit failed: '//trim(msg), status, message)
            return
         end if
         call proximity_impute(x, missing, cats, forest%proximity, imputed)
      end do
   end subroutine rf_impute_regression

   subroutine rough_initial_impute(x, ncat, missing, imputed, status)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: ncat(:)
      logical, intent(in) :: missing(:,:)
      real(dp), intent(out) :: imputed(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: observed(:)
      integer :: j, c, nobs, best_count, best_cat, current_count

      imputed = x
      status = 0
      do j = 1, size(x, 2)
         if (.not. any(missing(:, j))) cycle
         nobs = count(.not. missing(:, j))
         if (nobs == 0) then
            status = 1
            return
         end if
         if (ncat(j) <= 1) then
            allocate(observed(nobs))
            observed = pack(x(:, j), .not. missing(:, j))
            where (missing(:, j)) imputed(:, j) = r_median(observed)
            deallocate(observed)
         else
            best_count = -1
            best_cat = 1
            do c = 1, ncat(j)
               current_count = count((.not. missing(:, j)) .and. nint(x(:, j)) == c)
               if (current_count > best_count) then
                  best_count = current_count
                  best_cat = c
               end if
            end do
            where (missing(:, j)) imputed(:, j) = real(best_cat, dp)
         end if
      end do
   end subroutine rough_initial_impute

   subroutine proximity_impute(original, missing, ncat, proximity, imputed)
      real(dp), intent(in) :: original(:,:), proximity(:,:)
      logical, intent(in) :: missing(:,:)
      integer, intent(in) :: ncat(:)
      real(dp), intent(inout) :: imputed(:,:)
      integer :: i, j, k, c, best_cat, count_cat
      real(dp) :: weight_sum, weighted_sum, best_mean, current_mean

      do j = 1, size(original, 2)
         if (.not. any(missing(:, j))) cycle
         do i = 1, size(original, 1)
            if (.not. missing(i, j)) cycle
            if (ncat(j) <= 1) then
               weight_sum = 0.0_dp
               weighted_sum = 0.0_dp
               do k = 1, size(original, 1)
                  if (missing(k, j)) cycle
                  weight_sum = weight_sum + proximity(i, k)
                  weighted_sum = weighted_sum + proximity(i, k) * imputed(k, j)
               end do
               imputed(i, j) = weighted_sum / (1.0e-8_dp + weight_sum)
            else
               best_cat = max(1, nint(imputed(i, j)))
               best_mean = -huge(1.0_dp)
               do c = 1, ncat(j)
                  weight_sum = 0.0_dp
                  count_cat = 0
                  do k = 1, size(original, 1)
                     if (missing(k, j)) cycle
                     if (nint(original(k, j)) == c) then
                        weight_sum = weight_sum + proximity(i, k)
                        count_cat = count_cat + 1
                     end if
                  end do
                  if (count_cat > 0) then
                     current_mean = weight_sum / real(count_cat, dp)
                     if (current_mean > best_mean) then
                        best_mean = current_mean
                        best_cat = c
                     end if
                  end if
               end do
               imputed(i, j) = real(best_cat, dp)
            end if
         end do
      end do
   end subroutine proximity_impute

   subroutine tune_classification_mtry(x, y, mtry_values, oob_error, best_mtry, ncat, base_options, &
      mtry_start, ntree_try, step_factor, improve, status, message)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      integer, allocatable, intent(out) :: mtry_values(:)
      real(dp), allocatable, intent(out) :: oob_error(:)
      integer, intent(out) :: best_mtry
      integer, intent(in), optional :: ncat(:), mtry_start, ntree_try
      type(rf_options), intent(in), optional :: base_options
      real(dp), intent(in), optional :: step_factor, improve
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(rf_options) :: opt
      integer, allocatable :: mtmp(:)
      real(dp), allocatable :: etmp(:)
      integer :: p, start, ntry, nstored, direction, current, old, rc
      real(dp) :: factor, threshold, error_old, error_current, gain
      type(rf_classification_forest) :: forest
      character(len=256) :: msg

      call set_ok(status, message)
      p = size(x, 2)
      if (p <= 0 .or. size(y) /= size(x, 1)) then
         call set_error(1, 'tune_classification_mtry: incompatible shapes', status, message)
         allocate(mtry_values(0), oob_error(0))
         best_mtry = 0
         return
      end if
      start = max(1, int(sqrt(real(p, dp))))
      if (present(mtry_start)) start = max(1, min(p, mtry_start))
      ntry = 50
      if (present(ntree_try)) ntry = max(1, ntree_try)
      factor = 2.0_dp
      if (present(step_factor)) factor = step_factor
      threshold = 0.05_dp
      if (present(improve)) threshold = improve
      if (factor <= 1.0_dp .or. threshold < 0.0_dp) then
         call set_error(2, 'tune_classification_mtry: step_factor must exceed one and improve must be nonnegative', &
            status, message)
         allocate(mtry_values(0), oob_error(0))
         best_mtry = 0
         return
      end if
      opt = rf_options()
      if (present(base_options)) opt = base_options
      opt%ntree = ntry
      opt%mtry = start
      opt%importance = .false.
      opt%proximity = .false.
      call fit_classification(x, y, forest, ncat=ncat, options=opt, status=rc, message=msg)
      if (rc /= 0) then
         call set_error(10 + rc, 'tune_classification_mtry: initial fit failed: '//trim(msg), status, message)
         allocate(mtry_values(0), oob_error(0))
         best_mtry = 0
         return
      end if
      error_old = forest%error_curve(size(forest%trees), 1)
      allocate(mtmp(2 * p + 1), etmp(2 * p + 1))
      nstored = 1
      mtmp(1) = start
      etmp(1) = error_old

      do direction = -1, 1, 2
         gain = 1.1_dp * threshold
         current = start
         do while (gain >= threshold)
            old = current
            if (direction < 0) then
               current = max(1, ceiling(real(current, dp) / factor))
            else
               current = min(p, floor(real(current, dp) * factor))
            end if
            if (current == old) exit
            opt%mtry = current
            opt%seed = opt%seed + int(nstored, i64)
            call fit_classification(x, y, forest, ncat=ncat, options=opt, status=rc, message=msg)
            if (rc /= 0) then
               call set_error(20 + rc, 'tune_classification_mtry: trial fit failed: '//trim(msg), status, message)
               allocate(mtry_values(0), oob_error(0))
               best_mtry = 0
               return
            end if
            error_current = forest%error_curve(size(forest%trees), 1)
            call append_unique_trial(current, error_current, mtmp, etmp, nstored)
            if (error_old > 0.0_dp) then
               gain = 1.0_dp - error_current / error_old
            else
               gain = -huge(1.0_dp)
            end if
            if (gain > threshold) error_old = error_current
         end do
      end do
      call finish_tuning(mtmp, etmp, nstored, mtry_values, oob_error, best_mtry)
   end subroutine tune_classification_mtry

   subroutine tune_regression_mtry(x, y, mtry_values, oob_error, best_mtry, ncat, base_options, &
      mtry_start, ntree_try, step_factor, improve, status, message)
      real(dp), intent(in) :: x(:,:), y(:)
      integer, allocatable, intent(out) :: mtry_values(:)
      real(dp), allocatable, intent(out) :: oob_error(:)
      integer, intent(out) :: best_mtry
      integer, intent(in), optional :: ncat(:), mtry_start, ntree_try
      type(rf_options), intent(in), optional :: base_options
      real(dp), intent(in), optional :: step_factor, improve
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(rf_options) :: opt
      integer, allocatable :: mtmp(:)
      real(dp), allocatable :: etmp(:)
      integer :: p, start, ntry, nstored, direction, current, old, rc
      real(dp) :: factor, threshold, error_old, error_current, gain
      type(rf_regression_forest) :: forest
      character(len=256) :: msg

      call set_ok(status, message)
      p = size(x, 2)
      if (p <= 0 .or. size(y) /= size(x, 1)) then
         call set_error(1, 'tune_regression_mtry: incompatible shapes', status, message)
         allocate(mtry_values(0), oob_error(0))
         best_mtry = 0
         return
      end if
      start = max(1, p / 3)
      if (present(mtry_start)) start = max(1, min(p, mtry_start))
      ntry = 50
      if (present(ntree_try)) ntry = max(1, ntree_try)
      factor = 2.0_dp
      if (present(step_factor)) factor = step_factor
      threshold = 0.05_dp
      if (present(improve)) threshold = improve
      if (factor <= 1.0_dp .or. threshold < 0.0_dp) then
         call set_error(2, 'tune_regression_mtry: step_factor must exceed one and improve must be nonnegative', status, message)
         allocate(mtry_values(0), oob_error(0))
         best_mtry = 0
         return
      end if
      opt = rf_options()
      if (present(base_options)) opt = base_options
      opt%ntree = ntry
      opt%mtry = start
      opt%importance = .false.
      opt%proximity = .false.
      call fit_regression(x, y, forest, ncat=ncat, options=opt, status=rc, message=msg)
      if (rc /= 0) then
         call set_error(10 + rc, 'tune_regression_mtry: initial fit failed: '//trim(msg), status, message)
         allocate(mtry_values(0), oob_error(0))
         best_mtry = 0
         return
      end if
      error_old = forest%mse_curve(size(forest%trees))
      allocate(mtmp(2 * p + 1), etmp(2 * p + 1))
      nstored = 1
      mtmp(1) = start
      etmp(1) = error_old

      do direction = -1, 1, 2
         gain = 1.1_dp * threshold
         current = start
         do while (gain >= threshold)
            old = current
            if (direction < 0) then
               current = max(1, ceiling(real(current, dp) / factor))
            else
               current = min(p, floor(real(current, dp) * factor))
            end if
            if (current == old) exit
            opt%mtry = current
            opt%seed = opt%seed + int(nstored, i64)
            call fit_regression(x, y, forest, ncat=ncat, options=opt, status=rc, message=msg)
            if (rc /= 0) then
               call set_error(20 + rc, 'tune_regression_mtry: trial fit failed: '//trim(msg), status, message)
               allocate(mtry_values(0), oob_error(0))
               best_mtry = 0
               return
            end if
            error_current = forest%mse_curve(size(forest%trees))
            call append_unique_trial(current, error_current, mtmp, etmp, nstored)
            if (error_old > 0.0_dp) then
               gain = 1.0_dp - error_current / error_old
            else
               gain = -huge(1.0_dp)
            end if
            if (gain > threshold) error_old = error_current
         end do
      end do
      call finish_tuning(mtmp, etmp, nstored, mtry_values, oob_error, best_mtry)
   end subroutine tune_regression_mtry

   subroutine append_unique_trial(mtry, err, mvalues, errors, nstored)
      integer, intent(in) :: mtry
      real(dp), intent(in) :: err
      integer, intent(inout) :: mvalues(:), nstored
      real(dp), intent(inout) :: errors(:)
      integer :: i

      do i = 1, nstored
         if (mvalues(i) == mtry) then
            errors(i) = err
            return
         end if
      end do
      nstored = nstored + 1
      mvalues(nstored) = mtry
      errors(nstored) = err
   end subroutine append_unique_trial

   subroutine finish_tuning(mtmp, etmp, nstored, mtry_values, oob_error, best_mtry)
      integer, intent(in) :: mtmp(:), nstored
      real(dp), intent(in) :: etmp(:)
      integer, allocatable, intent(out) :: mtry_values(:)
      real(dp), allocatable, intent(out) :: oob_error(:)
      integer, intent(out) :: best_mtry
      integer :: i, j, mi
      real(dp) :: ei

      allocate(mtry_values(nstored), oob_error(nstored))
      mtry_values = mtmp(1:nstored)
      oob_error = etmp(1:nstored)
      do i = 2, nstored
         mi = mtry_values(i)
         ei = oob_error(i)
         j = i - 1
         do while (j >= 1)
            if (mtry_values(j) <= mi) exit
            mtry_values(j + 1) = mtry_values(j)
            oob_error(j + 1) = oob_error(j)
            j = j - 1
         end do
         mtry_values(j + 1) = mi
         oob_error(j + 1) = ei
      end do
      best_mtry = mtry_values(minloc(oob_error, dim=1))
   end subroutine finish_tuning

   subroutine rfcv_classification(x, y, nvar_sizes, error_cv, predicted, ncat, cv_fold, step, recursive, &
      base_options, status, message)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      integer, allocatable, intent(out) :: nvar_sizes(:), predicted(:,:)
      real(dp), allocatable, intent(out) :: error_cv(:)
      integer, intent(in), optional :: ncat(:), cv_fold
      real(dp), intent(in), optional :: step
      logical, intent(in), optional :: recursive
      type(rf_options), intent(in), optional :: base_options
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      integer, allocatable :: cats(:), fold(:), train_idx(:), test_idx(:), ranking(:), subset(:), subcats(:), pred(:)
      real(dp), allocatable :: xtrain(:,:), xtest(:,:)
      type(rf_options) :: opt
      type(rf_classification_forest) :: forest
      type(rf_rng_state) :: rng
      integer :: n, p, folds, f, j, ntrain, ntest, rc, k, c, pos, nclass
      real(dp) :: reduction
      logical :: rec
      character(len=256) :: msg

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (n <= 1 .or. p <= 0 .or. size(y) /= n) then
         call set_error(1, 'rfcv_classification: incompatible shapes', status, message)
         allocate(nvar_sizes(0), error_cv(0), predicted(0, 0))
         return
      end if
      call prepare_ncat(x, ncat, cats, rc)
      if (rc /= 0) then
         call set_error(2, 'rfcv_classification: invalid categorical coding', status, message)
         allocate(nvar_sizes(0), error_cv(0), predicted(0, 0))
         return
      end if
      folds = 5
      if (present(cv_fold)) folds = max(2, min(n, cv_fold))
      reduction = 0.5_dp
      if (present(step)) reduction = step
      if (reduction <= 0.0_dp .or. reduction >= 1.0_dp) then
         call set_error(3, 'rfcv_classification: step must be in (0,1)', status, message)
         allocate(nvar_sizes(0), error_cv(0), predicted(0, 0))
         return
      end if
      rec = .false.
      if (present(recursive)) rec = recursive
      call make_feature_sizes(p, reduction, nvar_sizes)
      allocate(predicted(n, size(nvar_sizes)), fold(n))
      predicted = spread(y, dim=2, ncopies=size(nvar_sizes))
      opt = rf_options()
      if (present(base_options)) opt = base_options
      call rng%seed(opt%seed)
      nclass = maxval(y)
      fold = 0
      do c = 1, nclass
         k = count(y == c)
         allocate(subset(k))
         pos = 0
         do j = 1, n
            if (y(j) == c) then
               pos = pos + 1
               subset(pos) = j
            end if
         end do
         call shuffle_int(rng, subset)
         do j = 1, k
            fold(subset(j)) = 1 + mod(j - 1, folds)
         end do
         deallocate(subset)
      end do

      do f = 1, folds
         ntrain = count(fold /= f)
         ntest = count(fold == f)
         allocate(train_idx(ntrain), test_idx(ntest))
         train_idx = pack([(j, j = 1, n)], fold /= f)
         test_idx = pack([(j, j = 1, n)], fold == f)
         allocate(ranking(p))
         ranking = [(j, j = 1, p)]
         do k = 1, size(nvar_sizes)
            allocate(subset(nvar_sizes(k)), subcats(nvar_sizes(k)))
            subset = ranking(1:nvar_sizes(k))
            subcats = cats(subset)
            allocate(xtrain(ntrain, nvar_sizes(k)), xtest(ntest, nvar_sizes(k)), pred(ntest))
            xtrain = x(train_idx, subset)
            xtest = x(test_idx, subset)
            opt%mtry = max(1, int(sqrt(real(nvar_sizes(k), dp))))
            opt%importance = (k == 1) .or. rec
            opt%proximity = .false.
            opt%seed = opt%seed + int(f * 1000 + k, i64)
            call fit_classification(xtrain, y(train_idx), forest, ncat=subcats, options=opt, status=rc, message=msg)
            if (rc /= 0) then
               call set_error(10 + rc, 'rfcv_classification: forest fit failed: '//trim(msg), status, message)
               deallocate(train_idx, test_idx, ranking, subset, subcats, xtrain, xtest, pred)
               return
            end if
            call predict_classification(forest, xtest, pred)
            predicted(test_idx, k) = pred
            if (k == 1 .or. rec) then
               call rank_class_importance(forest, subset, ranking)
            end if
            deallocate(subset, subcats, xtrain, xtest, pred)
         end do
         deallocate(train_idx, test_idx, ranking)
      end do
      allocate(error_cv(size(nvar_sizes)))
      do k = 1, size(nvar_sizes)
         error_cv(k) = real(count(predicted(:, k) /= y), dp) / real(n, dp)
      end do
   end subroutine rfcv_classification

   subroutine rfcv_regression(x, y, nvar_sizes, error_cv, predicted, ncat, cv_fold, step, recursive, &
      base_options, status, message)
      real(dp), intent(in) :: x(:,:), y(:)
      integer, allocatable, intent(out) :: nvar_sizes(:)
      real(dp), allocatable, intent(out) :: error_cv(:), predicted(:,:)
      integer, intent(in), optional :: ncat(:), cv_fold
      real(dp), intent(in), optional :: step
      logical, intent(in), optional :: recursive
      type(rf_options), intent(in), optional :: base_options
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      integer, allocatable :: cats(:), fold(:), train_idx(:), test_idx(:), ranking(:), subset(:), subcats(:), order(:)
      real(dp), allocatable :: xtrain(:,:), xtest(:,:), pred(:)
      type(rf_options) :: opt
      type(rf_regression_forest) :: forest
      type(rf_rng_state) :: rng
      integer :: n, p, folds, f, j, ntrain, ntest, rc, k, group, ngroup, pos
      real(dp) :: reduction
      logical :: rec
      character(len=256) :: msg

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (n <= 1 .or. p <= 0 .or. size(y) /= n) then
         call set_error(1, 'rfcv_regression: incompatible shapes', status, message)
         allocate(nvar_sizes(0), error_cv(0), predicted(0, 0))
         return
      end if
      call prepare_ncat(x, ncat, cats, rc)
      if (rc /= 0) then
         call set_error(2, 'rfcv_regression: invalid categorical coding', status, message)
         allocate(nvar_sizes(0), error_cv(0), predicted(0, 0))
         return
      end if
      folds = 5
      if (present(cv_fold)) folds = max(2, min(n, cv_fold))
      reduction = 0.5_dp
      if (present(step)) reduction = step
      if (reduction <= 0.0_dp .or. reduction >= 1.0_dp) then
         call set_error(3, 'rfcv_regression: step must be in (0,1)', status, message)
         allocate(nvar_sizes(0), error_cv(0), predicted(0, 0))
         return
      end if
      rec = .false.
      if (present(recursive)) rec = recursive
      call make_feature_sizes(p, reduction, nvar_sizes)
      allocate(predicted(n, size(nvar_sizes)), fold(n), order(n))
      predicted = spread(y, dim=2, ncopies=size(nvar_sizes))
      order = [(j, j = 1, n)]
      call order_by_response(y, order)
      opt = rf_options()
      if (present(base_options)) opt = base_options
      call rng%seed(opt%seed)
      fold = 0
      do group = 1, min(5, n)
         ngroup = count([(mod(j - 1, min(5, n)) + 1 == group, j = 1, n)])
         allocate(subset(ngroup))
         pos = 0
         do j = 1, n
            if (mod(j - 1, min(5, n)) + 1 == group) then
               pos = pos + 1
               subset(pos) = order(j)
            end if
         end do
         call shuffle_int(rng, subset)
         do j = 1, ngroup
            fold(subset(j)) = 1 + mod(j - 1, folds)
         end do
         deallocate(subset)
      end do

      do f = 1, folds
         ntrain = count(fold /= f)
         ntest = count(fold == f)
         allocate(train_idx(ntrain), test_idx(ntest))
         train_idx = pack([(j, j = 1, n)], fold /= f)
         test_idx = pack([(j, j = 1, n)], fold == f)
         allocate(ranking(p))
         ranking = [(j, j = 1, p)]
         do k = 1, size(nvar_sizes)
            allocate(subset(nvar_sizes(k)), subcats(nvar_sizes(k)))
            subset = ranking(1:nvar_sizes(k))
            subcats = cats(subset)
            allocate(xtrain(ntrain, nvar_sizes(k)), xtest(ntest, nvar_sizes(k)), pred(ntest))
            xtrain = x(train_idx, subset)
            xtest = x(test_idx, subset)
            opt%mtry = max(1, int(sqrt(real(nvar_sizes(k), dp))))
            opt%importance = (k == 1) .or. rec
            opt%proximity = .false.
            opt%seed = opt%seed + int(f * 1000 + k, i64)
            call fit_regression(xtrain, y(train_idx), forest, ncat=subcats, options=opt, status=rc, message=msg)
            if (rc /= 0) then
               call set_error(10 + rc, 'rfcv_regression: forest fit failed: '//trim(msg), status, message)
               deallocate(train_idx, test_idx, ranking, subset, subcats, xtrain, xtest, pred)
               return
            end if
            call predict_regression(forest, xtest, pred)
            predicted(test_idx, k) = pred
            if (k == 1 .or. rec) then
               call rank_reg_importance(forest, subset, ranking)
            end if
            deallocate(subset, subcats, xtrain, xtest, pred)
         end do
         deallocate(train_idx, test_idx, ranking)
      end do
      allocate(error_cv(size(nvar_sizes)))
      do k = 1, size(nvar_sizes)
         error_cv(k) = sum((predicted(:, k) - y) ** 2) / real(n, dp)
      end do
   end subroutine rfcv_regression

   subroutine make_feature_sizes(p, step, sizes)
      integer, intent(in) :: p
      real(dp), intent(in) :: step
      integer, allocatable, intent(out) :: sizes(:)
      integer, allocatable :: temp(:)
      integer :: count_sizes, level, candidate

      allocate(temp(p + 2))
      count_sizes = 1
      temp(1) = p
      level = 1
      do
         candidate = max(1, nint(real(p, dp) * step ** level))
         if (candidate /= temp(count_sizes)) then
            count_sizes = count_sizes + 1
            temp(count_sizes) = candidate
         end if
         if (candidate == 1) exit
         level = level + 1
      end do
      allocate(sizes(count_sizes))
      sizes = temp(1:count_sizes)
   end subroutine make_feature_sizes

   subroutine rank_class_importance(forest, current_subset, ranking)
      type(rf_classification_forest), intent(in) :: forest
      integer, intent(in) :: current_subset(:)
      integer, intent(inout) :: ranking(:)
      real(dp), allocatable :: importance(:)
      integer, allocatable :: local_order(:)
      integer :: i

      if (.not. allocated(forest%importance_accuracy)) return
      allocate(importance(size(current_subset)), local_order(size(current_subset)))
      importance = forest%importance_accuracy(:, forest%nclass + 1)
      local_order = [(i, i = 1, size(current_subset))]
      call order_descending(importance, local_order)
      ranking(1:size(current_subset)) = current_subset(local_order)
   end subroutine rank_class_importance

   subroutine rank_reg_importance(forest, current_subset, ranking)
      type(rf_regression_forest), intent(in) :: forest
      integer, intent(in) :: current_subset(:)
      integer, intent(inout) :: ranking(:)
      integer, allocatable :: local_order(:)
      integer :: i

      if (.not. allocated(forest%importance_accuracy)) return
      allocate(local_order(size(current_subset)))
      local_order = [(i, i = 1, size(current_subset))]
      call order_descending(forest%importance_accuracy, local_order)
      ranking(1:size(current_subset)) = current_subset(local_order)
   end subroutine rank_reg_importance

   subroutine order_descending(values, order)
      real(dp), intent(in) :: values(:)
      integer, intent(inout) :: order(:)
      integer :: i, j, key

      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (values(order(j)) >= values(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine order_descending

   subroutine order_by_response(y, order)
      real(dp), intent(in) :: y(:)
      integer, intent(inout) :: order(:)
      integer :: i, j, key

      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (y(order(j)) <= y(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine order_by_response

   subroutine mds_coordinates(proximity, k, coordinates, eigenvalues, status, message)
      real(dp), intent(in) :: proximity(:,:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: coordinates(:,:), eigenvalues(:)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      real(dp), allocatable :: d2(:,:), b(:,:), row_mean(:), vectors(:,:)
      real(dp) :: grand_mean
      integer :: n, ndim, i, j, info

      call set_ok(status, message)
      n = size(proximity, 1)
      if (size(proximity, 2) /= n .or. n == 0 .or. k < 1) then
         call set_error(1, 'mds_coordinates: proximity must be nonempty square and k positive', status, message)
         allocate(coordinates(0, 0), eigenvalues(0))
         return
      end if
      if (maxval(abs(proximity - transpose(proximity))) > 100.0_dp * epsilon(1.0_dp)) then
         call set_error(2, 'mds_coordinates: proximity must be symmetric', status, message)
         allocate(coordinates(0, 0), eigenvalues(0))
         return
      end if
      allocate(d2(n, n), b(n, n), row_mean(n))
      d2 = (1.0_dp - proximity) ** 2
      row_mean = sum(d2, dim=2) / real(n, dp)
      grand_mean = sum(d2) / real(n * n, dp)
      do j = 1, n
         do i = 1, n
            b(i, j) = -0.5_dp * (d2(i, j) - row_mean(i) - row_mean(j) + grand_mean)
         end do
      end do
      call symmetric_eigen(b, eigenvalues, vectors, info, descending=.true.)
      if (info /= 0) then
         call set_error(3, 'mds_coordinates: symmetric eigensolver failed', status, message)
         allocate(coordinates(0, 0))
         if (.not. allocated(eigenvalues)) allocate(eigenvalues(0))
         return
      end if
      ndim = min(k, n)
      allocate(coordinates(n, ndim))
      coordinates = 0.0_dp
      do j = 1, ndim
         if (eigenvalues(j) > 0.0_dp) coordinates(:, j) = vectors(:, j) * sqrt(eigenvalues(j))
      end do
   end subroutine mds_coordinates

   subroutine prepare_ncat(x, requested, ncat, status)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in), optional :: requested(:)
      integer, allocatable, intent(out) :: ncat(:)
      integer, intent(out) :: status
      integer :: j, i

      allocate(ncat(size(x, 2)))
      if (present(requested)) then
         if (size(requested) /= size(x, 2) .or. any(requested < 1) .or. any(requested > 53)) then
            status = 1
            return
         end if
         ncat = requested
      else
         ncat = 1
      end if
      status = 0
      do j = 1, size(x, 2)
         if (ncat(j) <= 1) cycle
         do i = 1, size(x, 1)
            if (ieee_is_nan(x(i, j))) cycle
            if (abs(x(i, j) - real(nint(x(i, j)), dp)) > 100.0_dp * epsilon(1.0_dp)) then
               status = 2
               return
            end if
            if (nint(x(i, j)) < 1 .or. nint(x(i, j)) > ncat(j)) then
               status = 3
               return
            end if
         end do
      end do
   end subroutine prepare_ncat

   subroutine set_ok(status, message)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      if (present(status)) status = 0
      if (present(message)) message = ''
   end subroutine set_ok

   subroutine set_error(code, text, status, message)
      integer, intent(in) :: code
      character(len=*), intent(in) :: text
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      if (present(status)) status = code
      if (present(message)) message = text
   end subroutine set_error

end module rf_workflows
