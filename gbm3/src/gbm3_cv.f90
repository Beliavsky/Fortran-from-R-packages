! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_cv
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gbm3_kinds, only : dp
   use gbm3_constants
   use gbm3_types, only : gbm_options, gbm_model, gbm_cv_result
   use gbm3_math, only : shuffle_int
   use gbm3_core, only : gbm_fit, gbm_predict
   implicit none
   private

   interface gbm_cross_validate
      module procedure gbm_cross_validate_vector
      module procedure gbm_cross_validate_cox
   end interface gbm_cross_validate

   public :: gbm_cross_validate, gbm_cross_validate_vector, gbm_cross_validate_cox
   public :: gbm_cv_best_iteration

contains

   subroutine gbm_cross_validate_vector(x, y, result, n_folds, options, weight, offset, &
                                        var_classes, monotone, id, group, fold_id, stratify, full_model)
      real(dp), intent(in) :: x(:, :), y(:)
      type(gbm_cv_result), intent(out) :: result
      integer, intent(in) :: n_folds
      type(gbm_options), intent(in), optional :: options
      real(dp), intent(in), optional :: weight(:), offset(:)
      integer, intent(in), optional :: var_classes(:), monotone(:), id(:), group(:), fold_id(:)
      logical, intent(in), optional :: stratify
      type(gbm_model), intent(out), optional :: full_model

      type(gbm_options) :: opt, fold_opt
      type(gbm_model) :: full_fit
      type(gbm_model), allocatable :: fold_models(:)
      real(dp), allocatable :: w(:), off(:), xfold(:, :), yfold(:), wfold(:), offfold(:), pred(:)
      integer, allocatable :: vc(:), mono(:), ids(:), grp(:), folds(:), keys(:), idx_train(:), idx_valid(:)
      integer, allocatable :: idfold(:), grpfold(:)
      integer :: n, p, nbase, fold, ntr, nval, t, i
      logical :: do_stratify

      n = size(x, 1)
      p = size(x, 2)
      if (size(y) /= n) error stop "gbm_cross_validate: response length does not match x"
      if (n <= 1 .or. p <= 0) error stop "gbm_cross_validate: x must have at least two rows"
      if (n_folds < 2) error stop "gbm_cross_validate: n_folds must be at least 2"

      opt = gbm_options()
      if (present(options)) opt = options
      if (opt%distribution == GBM_COXPH) error stop "gbm_cross_validate: Cox PH requires survival response"
      if (opt%distribution == GBM_PAIRWISE .and. .not. present(group)) &
         error stop "gbm_cross_validate: pairwise distribution requires group"
      nbase = opt%num_train
      if (nbase <= 0) nbase = n
      if (nbase > n .or. nbase < 2) error stop "gbm_cross_validate: invalid num_train"

      allocate(w(n), off(n), vc(p), mono(p), ids(n), grp(n))
      w = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n) error stop "gbm_cross_validate: weight length mismatch"
         w = weight
      end if
      off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_cross_validate: offset length mismatch"
         off = offset
      end if
      vc = 0
      if (present(var_classes)) then
         if (size(var_classes) /= p) error stop "gbm_cross_validate: var_classes length mismatch"
         vc = var_classes
      end if
      mono = 0
      if (present(monotone)) then
         if (size(monotone) /= p) error stop "gbm_cross_validate: monotone length mismatch"
         mono = monotone
      end if
      do i = 1, n
         ids(i) = i
         grp(i) = i
      end do
      if (present(id)) then
         if (size(id) /= n) error stop "gbm_cross_validate: id length mismatch"
         ids = id
      end if
      if (present(group)) then
         if (size(group) /= n) error stop "gbm_cross_validate: group length mismatch"
         grp = group
      end if

      allocate(keys(nbase), folds(nbase))
      if (opt%distribution == GBM_PAIRWISE) then
         keys = grp(1:nbase)
      else
         keys = ids(1:nbase)
      end if
      do_stratify = .false.
      if (present(stratify)) do_stratify = stratify
      if (do_stratify .and. opt%distribution /= GBM_BERNOULLI) &
         error stop "gbm_cross_validate: stratify is only supported for Bernoulli"
      if (present(fold_id)) then
         if (size(fold_id) /= nbase) error stop "gbm_cross_validate: fold_id length must equal num_train"
         call make_cv_folds(keys, n_folds, folds, requested=fold_id, stratify=do_stratify, y=y(1:nbase))
      else
         call make_cv_folds(keys, n_folds, folds, stratify=do_stratify, y=y(1:nbase))
      end if

      ! gbm3 creates CV groups first, then fits the full model, then the folds.
      call gbm_fit(x, y, full_fit, opt, weight=w, offset=off, var_classes=vc, monotone=mono, id=ids, group=grp)
      if (present(full_model)) full_model = full_fit

      result%n_folds = n_folds
      allocate(result%error(opt%num_trees), result%fitted(nbase), result%fold_id(nbase), fold_models(n_folds))
      result%error = 0.0_dp
      result%fitted = 0.0_dp
      result%fold_id = folds

      do fold = 1, n_folds
         call partition_rows(folds, fold, idx_train, idx_valid)
         ntr = size(idx_train)
         nval = size(idx_valid)
         if (ntr < 1 .or. nval < 1) error stop "gbm_cross_validate: each fold needs training and validation rows"

         allocate(xfold(nbase, p), yfold(nbase), wfold(nbase), offfold(nbase), idfold(nbase), grpfold(nbase))
         xfold(1:ntr, :) = x(idx_train, :)
         xfold(ntr + 1:nbase, :) = x(idx_valid, :)
         yfold(1:ntr) = y(idx_train)
         yfold(ntr + 1:nbase) = y(idx_valid)
         wfold(1:ntr) = w(idx_train)
         wfold(ntr + 1:nbase) = w(idx_valid)
         offfold(1:ntr) = off(idx_train)
         offfold(ntr + 1:nbase) = off(idx_valid)
         idfold(1:ntr) = ids(idx_train)
         idfold(ntr + 1:nbase) = ids(idx_valid)
         grpfold(1:ntr) = grp(idx_train)
         grpfold(ntr + 1:nbase) = grp(idx_valid)

         fold_opt = opt
         fold_opt%num_train = ntr
         call gbm_fit(xfold, yfold, fold_models(fold), fold_opt, weight=wfold, offset=offfold, &
                      var_classes=vc, monotone=mono, id=idfold, group=grpfold)
         do t = 1, opt%num_trees
            result%error(t) = result%error(t) + fold_models(fold)%validation_error(t) * real(nval, dp)
         end do
         deallocate(xfold, yfold, wfold, offfold, idfold, grpfold, idx_train, idx_valid)
      end do
      result%error = result%error / real(nbase, dp)
      result%best_iteration = argmin_finite(result%error)
      if (result%best_iteration == 0) error stop "gbm_cross_validate: no finite CV error"

      do fold = 1, n_folds
         idx_valid = pack([(i, i=1,nbase)], folds == fold)
         pred = gbm_predict(fold_models(fold), x(idx_valid, :), n_trees=result%best_iteration)
         result%fitted(idx_valid) = pred
         deallocate(idx_valid, pred)
      end do
   end subroutine gbm_cross_validate_vector

   subroutine gbm_cross_validate_cox(x, surv, result, n_folds, options, weight, offset, &
                                     var_classes, monotone, id, strata, fold_id, full_model)
      real(dp), intent(in) :: x(:, :), surv(:, :)
      type(gbm_cv_result), intent(out) :: result
      integer, intent(in) :: n_folds
      type(gbm_options), intent(in), optional :: options
      real(dp), intent(in), optional :: weight(:), offset(:)
      integer, intent(in), optional :: var_classes(:), monotone(:), id(:), strata(:), fold_id(:)
      type(gbm_model), intent(out), optional :: full_model

      type(gbm_options) :: opt, fold_opt
      type(gbm_model) :: full_fit
      type(gbm_model), allocatable :: fold_models(:)
      real(dp), allocatable :: w(:), off(:), xfold(:, :), sfold(:, :), wfold(:), offfold(:), pred(:)
      integer, allocatable :: vc(:), mono(:), ids(:), str(:), folds(:), keys(:), idx_train(:), idx_valid(:)
      integer, allocatable :: idfold(:), strfold(:)
      integer :: n, p, nbase, fold, ntr, nval, t, i

      n = size(x, 1)
      p = size(x, 2)
      if (size(surv, 1) /= n) error stop "gbm_cross_validate: survival response row mismatch"
      if (size(surv, 2) /= 2 .and. size(surv, 2) /= 3) &
         error stop "gbm_cross_validate: survival response must have 2 or 3 columns"
      if (n <= 1 .or. p <= 0) error stop "gbm_cross_validate: x must have at least two rows"
      if (n_folds < 2) error stop "gbm_cross_validate: n_folds must be at least 2"

      opt = gbm_options(distribution=GBM_COXPH)
      if (present(options)) opt = options
      if (opt%distribution /= GBM_COXPH) error stop "gbm_cross_validate: survival response requires GBM_COXPH"
      nbase = opt%num_train
      if (nbase <= 0) nbase = n
      if (nbase > n .or. nbase < 2) error stop "gbm_cross_validate: invalid num_train"

      allocate(w(n), off(n), vc(p), mono(p), ids(n), str(n))
      w = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n) error stop "gbm_cross_validate: weight length mismatch"
         w = weight
      end if
      off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_cross_validate: offset length mismatch"
         off = offset
      end if
      vc = 0
      if (present(var_classes)) then
         if (size(var_classes) /= p) error stop "gbm_cross_validate: var_classes length mismatch"
         vc = var_classes
      end if
      mono = 0
      if (present(monotone)) then
         if (size(monotone) /= p) error stop "gbm_cross_validate: monotone length mismatch"
         mono = monotone
      end if
      do i = 1, n
         ids(i) = i
      end do
      if (present(id)) then
         if (size(id) /= n) error stop "gbm_cross_validate: id length mismatch"
         ids = id
      end if
      str = 1
      if (present(strata)) then
         if (size(strata) /= n) error stop "gbm_cross_validate: strata length mismatch"
         str = strata
      end if

      allocate(keys(nbase), folds(nbase))
      keys = ids(1:nbase)
      if (present(fold_id)) then
         if (size(fold_id) /= nbase) error stop "gbm_cross_validate: fold_id length must equal num_train"
         call make_cv_folds(keys, n_folds, folds, requested=fold_id)
      else
         call make_cv_folds(keys, n_folds, folds)
      end if

      call gbm_fit(x, surv, full_fit, opt, weight=w, offset=off, var_classes=vc, monotone=mono, id=ids, strata=str)
      if (present(full_model)) full_model = full_fit

      result%n_folds = n_folds
      allocate(result%error(opt%num_trees), result%fitted(nbase), result%fold_id(nbase), fold_models(n_folds))
      result%error = 0.0_dp
      result%fitted = 0.0_dp
      result%fold_id = folds

      do fold = 1, n_folds
         call partition_rows(folds, fold, idx_train, idx_valid)
         ntr = size(idx_train)
         nval = size(idx_valid)
         allocate(xfold(nbase, p), sfold(nbase, size(surv, 2)), wfold(nbase), offfold(nbase), &
                  idfold(nbase), strfold(nbase))
         xfold(1:ntr, :) = x(idx_train, :)
         xfold(ntr + 1:nbase, :) = x(idx_valid, :)
         sfold(1:ntr, :) = surv(idx_train, :)
         sfold(ntr + 1:nbase, :) = surv(idx_valid, :)
         wfold(1:ntr) = w(idx_train)
         wfold(ntr + 1:nbase) = w(idx_valid)
         offfold(1:ntr) = off(idx_train)
         offfold(ntr + 1:nbase) = off(idx_valid)
         idfold(1:ntr) = ids(idx_train)
         idfold(ntr + 1:nbase) = ids(idx_valid)
         strfold(1:ntr) = str(idx_train)
         strfold(ntr + 1:nbase) = str(idx_valid)

         fold_opt = opt
         fold_opt%num_train = ntr
         call gbm_fit(xfold, sfold, fold_models(fold), fold_opt, weight=wfold, offset=offfold, &
                      var_classes=vc, monotone=mono, id=idfold, strata=strfold)
         do t = 1, opt%num_trees
            result%error(t) = result%error(t) + fold_models(fold)%validation_error(t) * real(nval, dp)
         end do
         deallocate(xfold, sfold, wfold, offfold, idfold, strfold, idx_train, idx_valid)
      end do
      result%error = result%error / real(nbase, dp)
      result%best_iteration = argmin_finite(result%error)
      if (result%best_iteration == 0) error stop "gbm_cross_validate: no finite CV error"

      do fold = 1, n_folds
         idx_valid = pack([(i, i=1,nbase)], folds == fold)
         pred = gbm_predict(fold_models(fold), x(idx_valid, :), n_trees=result%best_iteration)
         result%fitted(idx_valid) = pred
         deallocate(idx_valid, pred)
      end do
   end subroutine gbm_cross_validate_cox

   integer function gbm_cv_best_iteration(result) result(iteration)
      type(gbm_cv_result), intent(in) :: result
      if (result%best_iteration > 0) then
         iteration = result%best_iteration
      else if (allocated(result%error)) then
         iteration = argmin_finite(result%error)
      else
         iteration = 0
      end if
   end function gbm_cv_best_iteration

   subroutine make_cv_folds(keys, n_folds, folds, requested, stratify, y)
      integer, intent(in) :: keys(:), n_folds
      integer, intent(out) :: folds(size(keys))
      integer, intent(in), optional :: requested(:)
      logical, intent(in), optional :: stratify
      real(dp), intent(in), optional :: y(:)

      integer, allocatable :: unique_keys(:), unit_of_row(:), unit_fold(:), order(:), class0(:), class1(:)
      integer :: nunit, i, j, pos, n0, n1, cls
      logical :: found, do_stratify

      if (size(keys) < n_folds) error stop "gbm_cross_validate: more folds than training rows"
      if (present(requested)) then
         if (size(requested) /= size(keys)) error stop "make_cv_folds: requested fold length mismatch"
         if (any(requested < 1) .or. any(requested > n_folds)) &
            error stop "gbm_cross_validate: fold_id values must be in 1..n_folds"
         folds = requested
         do i = 1, size(keys)
            do j = i + 1, size(keys)
               if (keys(j) == keys(i) .and. folds(j) /= folds(i)) &
                  error stop "gbm_cross_validate: rows with the same id/group must share a fold"
            end do
         end do
         do i = 1, n_folds
            if (count(folds == i) == 0) error stop "gbm_cross_validate: fold_id contains an empty fold"
         end do
         return
      end if

      allocate(unique_keys(size(keys)), unit_of_row(size(keys)))
      nunit = 0
      do i = 1, size(keys)
         found = .false.
         do j = 1, nunit
            if (keys(i) == unique_keys(j)) then
               unit_of_row(i) = j
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            nunit = nunit + 1
            unique_keys(nunit) = keys(i)
            unit_of_row(i) = nunit
         end if
      end do
      if (nunit < n_folds) error stop "gbm_cross_validate: fewer independent ids/groups than folds"
      allocate(unit_fold(nunit))

      do_stratify = .false.
      if (present(stratify)) do_stratify = stratify
      if (do_stratify) then
         if (.not. present(y)) error stop "make_cv_folds: stratification requires response"
         if (size(y) /= size(keys)) error stop "make_cv_folds: stratification response length mismatch"
         allocate(class0(nunit), class1(nunit))
         n0 = 0
         n1 = 0
         do j = 1, nunit
            cls = -1
            do i = 1, size(keys)
               if (unit_of_row(i) /= j) cycle
               if (abs(y(i)) <= 1.0e-12_dp) then
                  if (cls == 1) error stop "gbm_cross_validate: one id spans both Bernoulli classes"
                  cls = 0
               else if (abs(y(i) - 1.0_dp) <= 1.0e-12_dp) then
                  if (cls == 0) error stop "gbm_cross_validate: one id spans both Bernoulli classes"
                  cls = 1
               else
                  error stop "gbm_cross_validate: Bernoulli stratification requires 0/1 responses"
               end if
            end do
            if (cls == 0) then
               n0 = n0 + 1
               class0(n0) = j
            else
               n1 = n1 + 1
               class1(n1) = j
            end if
         end do
         if (min(n0, n1) < n_folds) &
            error stop "gbm_cross_validate: smallest Bernoulli class has fewer ids than folds"
         call shuffle_int(class0(1:n0))
         call shuffle_int(class1(1:n1))
         do pos = 1, n0
            unit_fold(class0(pos)) = modulo(pos - 1, n_folds) + 1
         end do
         do pos = 1, n1
            unit_fold(class1(pos)) = modulo(pos - 1, n_folds) + 1
         end do
      else
         allocate(order(nunit))
         order = [(i, i=1,nunit)]
         call shuffle_int(order)
         do pos = 1, nunit
            unit_fold(order(pos)) = modulo(pos - 1, n_folds) + 1
         end do
      end if

      do i = 1, size(keys)
         folds(i) = unit_fold(unit_of_row(i))
      end do
   end subroutine make_cv_folds

   subroutine partition_rows(folds, target, train_rows, valid_rows)
      integer, intent(in) :: folds(:), target
      integer, allocatable, intent(out) :: train_rows(:), valid_rows(:)
      integer :: i
      train_rows = pack([(i, i=1,size(folds))], folds /= target)
      valid_rows = pack([(i, i=1,size(folds))], folds == target)
   end subroutine partition_rows

   integer function argmin_finite(x) result(idx)
      real(dp), intent(in) :: x(:)
      real(dp) :: best
      integer :: i
      idx = 0
      best = huge(1.0_dp)
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) cycle
         if (idx == 0 .or. x(i) < best) then
            idx = i
            best = x(i)
         end if
      end do
   end function argmin_finite

end module gbm3_cv
