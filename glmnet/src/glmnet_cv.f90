! SPDX-License-Identifier: GPL-2.0-only
module glmnet_cv
   use glmnet_kinds, only : dp, glmnet_eps
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument
   use glmnet_types, only : glmnet_control_type, glmnet_cv_result, &
      glmnet_path_result, glmnet_assessment_result, glmnet_family_gaussian, &
      glmnet_family_binomial, glmnet_family_poisson
   use glmnet_utils, only : deterministic_fold_ids
   use glmnet_gaussian, only : fit_gaussian_path
   use glmnet_glm, only : fit_binomial_path, fit_poisson_path
   use glmnet_multinomial, only : fit_multinomial_path
   use glmnet_cox, only : fit_cox_path, concordance_index
   use glmnet_predict, only : predict_glmnet
   use glmnet_assess, only : assess_glmnet, assess_multinomial
   implicit none
   private
   public :: cv_glmnet, cv_multinomial, cv_cox, build_predmat
contains
   subroutine cv_glmnet(x, y, family_code, result, control, nfolds, fold_id_in, &
      measure, weights, offset, seed)
      real(dp), intent(in) :: x(:,:), y(:)
      integer, intent(in) :: family_code
      type(glmnet_cv_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      integer, intent(in), optional :: nfolds, fold_id_in(:), seed
      character(len=*), intent(in), optional :: measure
      real(dp), intent(in), optional :: weights(:), offset(:)
      type(glmnet_control_type) :: ctl
      type(glmnet_path_result) :: fold_fit
      type(glmnet_assessment_result) :: assessment
      real(dp), allocatable :: x_train(:,:), x_test(:,:), y_train(:), y_test(:)
      real(dp), allocatable :: w_all(:), w_train(:), w_test(:), o_all(:), o_train(:), o_test(:)
      real(dp), allocatable :: pred(:,:,:)
      integer, allocatable :: train_index(:), test_index(:), fold_id(:)
      integer :: n, kfold, f, status, i, l, ntrain, ntest
      character(len=16) :: metric

      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      n = size(x, 1)
      result%status = glmnet_success
      if (size(y) /= n .or. n < 3) then
         result%status = glmnet_invalid_argument
         return
      end if
      kfold = min(10, n)
      if (present(nfolds)) kfold = nfolds
      if (kfold < 2 .or. kfold > n) then
         result%status = glmnet_invalid_argument
         return
      end if
      if (present(fold_id_in)) then
         if (size(fold_id_in) /= n .or. any(fold_id_in < 1) .or. any(fold_id_in > kfold)) then
            result%status = glmnet_invalid_argument
            return
         end if
         allocate(fold_id(n)); fold_id = fold_id_in
      else
         if (present(seed)) then
            call deterministic_fold_ids(n, kfold, fold_id, seed)
         else
            call deterministic_fold_ids(n, kfold, fold_id)
         end if
      end if
      allocate(w_all(n), o_all(n))
      if (present(weights)) then
         if (size(weights) /= n) then
            result%status = glmnet_invalid_argument
            return
         end if
         w_all = weights
      else
         w_all = 1.0_dp
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = glmnet_invalid_argument
            return
         end if
         o_all = offset
      else
         o_all = 0.0_dp
      end if
      select case (family_code)
      case (glmnet_family_gaussian)
         call fit_gaussian_path(x, y, result%fit, ctl, weights_in=w_all, offset_in=o_all)
         metric = 'mse'
      case (glmnet_family_binomial)
         call fit_binomial_path(x, y, result%fit, ctl, weights_in=w_all, offset_in=o_all)
         metric = 'deviance'
      case (glmnet_family_poisson)
         call fit_poisson_path(x, y, result%fit, ctl, weights_in=w_all, offset_in=o_all)
         metric = 'deviance'
      case default
         result%status = glmnet_invalid_argument
         return
      end select
      if (present(measure)) metric = adjustl(measure)
      if (result%fit%nlambda < 1) then
         result%status = result%fit%status
         return
      end if
      result%nfolds = kfold
      result%measure = metric
      result%fold_id = fold_id
      allocate(result%fold_values(kfold, result%fit%nlambda), &
         result%predictions(n, 1, result%fit%nlambda))
      result%predictions = 0.0_dp
      do f = 1, kfold
         ntest = count(fold_id == f)
         ntrain = n - ntest
         allocate(train_index(ntrain), test_index(ntest))
         ntrain = 0; ntest = 0
         do i = 1, n
            if (fold_id(i) == f) then
               ntest = ntest + 1; test_index(ntest) = i
            else
               ntrain = ntrain + 1; train_index(ntrain) = i
            end if
         end do
         allocate(x_train(ntrain, size(x, 2)), x_test(ntest, size(x, 2)), &
            y_train(ntrain), y_test(ntest), w_train(ntrain), w_test(ntest), &
            o_train(ntrain), o_test(ntest))
         x_train = x(train_index, :); x_test = x(test_index, :)
         y_train = y(train_index); y_test = y(test_index)
         w_train = w_all(train_index); w_test = w_all(test_index)
         o_train = o_all(train_index); o_test = o_all(test_index)
         select case (family_code)
         case (glmnet_family_gaussian)
            call fit_gaussian_path(x_train, y_train, fold_fit, ctl, weights_in=w_train, &
               offset_in=o_train, lambda_in=result%fit%lambda)
         case (glmnet_family_binomial)
            call fit_binomial_path(x_train, y_train, fold_fit, ctl, weights_in=w_train, &
               offset_in=o_train, lambda_in=result%fit%lambda)
         case (glmnet_family_poisson)
            call fit_poisson_path(x_train, y_train, fold_fit, ctl, weights_in=w_train, &
               offset_in=o_train, lambda_in=result%fit%lambda)
         end select
         call assess_glmnet(fold_fit, x_test, y_test, assessment, metric, w_test, o_test)
         if (assessment%status /= glmnet_success) then
            result%status = assessment%status
            return
         end if
         result%fold_values(f, :) = assessment%value
         call predict_with_vector_offset(fold_fit, x_test, o_test, pred, status)
         if (status /= glmnet_success) then
            result%status = status
            return
         end if
         do l = 1, result%fit%nlambda
            result%predictions(test_index, 1, l) = pred(:, 1, l)
         end do
         deallocate(train_index, test_index, x_train, x_test, y_train, y_test, &
            w_train, w_test, o_train, o_test, pred)
      end do
      call finalize_cv(result, higher_is_better=(trim(metric) == 'auc'))
   end subroutine cv_glmnet

   subroutine cv_multinomial(x, class_id, result, control, nfolds, fold_id_in, &
      measure, weights, offset, seed, nclass)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: class_id(:)
      type(glmnet_cv_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      integer, intent(in), optional :: nfolds, fold_id_in(:), seed, nclass
      character(len=*), intent(in), optional :: measure
      real(dp), intent(in), optional :: weights(:), offset(:,:)
      type(glmnet_control_type) :: ctl
      type(glmnet_path_result) :: fold_fit
      type(glmnet_assessment_result) :: assessment
      real(dp), allocatable :: x_train(:,:), x_test(:,:), w_all(:), w_train(:), w_test(:)
      real(dp), allocatable :: o_all(:,:), o_train(:,:), o_test(:,:), pred(:,:,:)
      integer, allocatable :: y_train(:), y_test(:), train_index(:), test_index(:), fold_id(:)
      integer :: n, k, kfold, f, status, i, l, ntrain, ntest
      character(len=16) :: metric
      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      n = size(x, 1)
      if (present(nclass)) then
         k = nclass
      else if (size(class_id) > 0) then
         k = maxval(class_id)
      else
         k = 0
      end if
      result%status = glmnet_success
      if (size(class_id) /= n .or. k < 2) then
         result%status = glmnet_invalid_argument
         return
      end if
      kfold = min(10, n)
      if (present(nfolds)) kfold = nfolds
      if (kfold < 2 .or. kfold > n) then
         result%status = glmnet_invalid_argument
         return
      end if
      if (present(fold_id_in)) then
         if (size(fold_id_in) /= n) then
            result%status = glmnet_invalid_argument
            return
         end if
         allocate(fold_id(n)); fold_id = fold_id_in
      else
         if (present(seed)) then
            call deterministic_fold_ids(n, kfold, fold_id, seed)
         else
            call deterministic_fold_ids(n, kfold, fold_id)
         end if
      end if
      allocate(w_all(n), o_all(n, k))
      if (present(weights)) then
         if (size(weights) /= n) then
            result%status = glmnet_invalid_argument
            return
         end if
         w_all = weights
      else
         w_all = 1.0_dp
      end if
      o_all = 0.0_dp
      if (present(offset)) then
         if (size(offset, 1) /= n .or. size(offset, 2) /= k) then
            result%status = glmnet_invalid_argument
            return
         end if
         o_all = offset
      end if
      call fit_multinomial_path(x, class_id, result%fit, ctl, weights_in=w_all, &
         offset_in=o_all, nclass=k)
      metric = 'deviance'
      if (present(measure)) metric = adjustl(measure)
      result%nfolds = kfold
      result%measure = metric
      result%fold_id = fold_id
      allocate(result%fold_values(kfold, result%fit%nlambda), &
         result%predictions(n, k, result%fit%nlambda))
      result%predictions = 0.0_dp
      do f = 1, kfold
         ntest = count(fold_id == f); ntrain = n - ntest
         allocate(train_index(ntrain), test_index(ntest))
         ntrain = 0; ntest = 0
         do i = 1, n
            if (fold_id(i) == f) then
               ntest = ntest + 1; test_index(ntest) = i
            else
               ntrain = ntrain + 1; train_index(ntrain) = i
            end if
         end do
         allocate(x_train(ntrain, size(x, 2)), x_test(ntest, size(x, 2)), &
            y_train(ntrain), y_test(ntest), w_train(ntrain), w_test(ntest), &
            o_train(ntrain, k), o_test(ntest, k))
         x_train = x(train_index, :); x_test = x(test_index, :)
         y_train = class_id(train_index); y_test = class_id(test_index)
         w_train = w_all(train_index); w_test = w_all(test_index)
         o_train = o_all(train_index, :); o_test = o_all(test_index, :)
         call fit_multinomial_path(x_train, y_train, fold_fit, ctl, weights_in=w_train, &
            offset_in=o_train, lambda_in=result%fit%lambda, nclass=k)
         call assess_multinomial(fold_fit, x_test, y_test, assessment, metric, w_test, o_test)
         if (assessment%status /= glmnet_success) then
            result%status = assessment%status
            return
         end if
         result%fold_values(f, :) = assessment%value
         call predict_glmnet(fold_fit, x_test, pred, status, prediction_type='response', offset=o_test)
         if (status /= glmnet_success) then
            result%status = status
            return
         end if
         do l = 1, result%fit%nlambda
            result%predictions(test_index, :, l) = pred(:, :, l)
         end do
         deallocate(train_index, test_index, x_train, x_test, y_train, y_test, &
            w_train, w_test, o_train, o_test, pred)
      end do
      call finalize_cv(result, higher_is_better=.false.)
   end subroutine cv_multinomial

   subroutine cv_cox(x, start_time, stop_time, event, result, control, nfolds, &
      fold_id_in, weights, offset, strata, seed, efron)
      real(dp), intent(in) :: x(:,:), start_time(:), stop_time(:)
      integer, intent(in) :: event(:)
      type(glmnet_cv_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      integer, intent(in), optional :: nfolds, fold_id_in(:), strata(:), seed
      real(dp), intent(in), optional :: weights(:), offset(:)
      logical, intent(in), optional :: efron
      type(glmnet_control_type) :: ctl
      type(glmnet_path_result) :: fold_fit
      real(dp), allocatable :: w_all(:), o_all(:), x_train(:,:), x_test(:,:), pred(:,:,:)
      real(dp), allocatable :: st_train(:), sp_train(:), st_test(:), sp_test(:)
      real(dp), allocatable :: w_train(:), w_test(:), o_train(:), o_test(:)
      integer, allocatable :: e_train(:), e_test(:), s_all(:), s_train(:), s_test(:)
      integer, allocatable :: train_index(:), test_index(:), fold_id(:)
      integer :: n, kfold, f, status, i, l, ntrain, ntest
      logical :: ef
      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      n = size(x, 1)
      result%status = glmnet_success
      if (size(start_time) /= n .or. size(stop_time) /= n .or. size(event) /= n) then
         result%status = glmnet_invalid_argument
         return
      end if
      kfold = min(10, n)
      if (present(nfolds)) kfold = nfolds
      if (kfold < 2 .or. kfold > n) then
         result%status = glmnet_invalid_argument
         return
      end if
      if (present(fold_id_in)) then
         if (size(fold_id_in) /= n) then
            result%status = glmnet_invalid_argument
            return
         end if
         allocate(fold_id(n)); fold_id = fold_id_in
      else
         if (present(seed)) then
            call deterministic_fold_ids(n, kfold, fold_id, seed)
         else
            call deterministic_fold_ids(n, kfold, fold_id)
         end if
      end if
      allocate(w_all(n), o_all(n), s_all(n))
      w_all = 1.0_dp; o_all = 0.0_dp; s_all = 1
      if (present(weights)) w_all = weights
      if (present(offset)) o_all = offset
      if (present(strata)) s_all = strata
      ef = .false.; if (present(efron)) ef = efron
      call fit_cox_path(x, start_time, stop_time, event, result%fit, ctl, &
         weights_in=w_all, offset_in=o_all, strata_in=s_all, efron=ef)
      result%nfolds = kfold
      result%measure = 'C'
      result%fold_id = fold_id
      allocate(result%fold_values(kfold, result%fit%nlambda), &
         result%predictions(n, 1, result%fit%nlambda))
      result%predictions = 0.0_dp
      do f = 1, kfold
         ntest = count(fold_id == f); ntrain = n - ntest
         allocate(train_index(ntrain), test_index(ntest))
         ntrain = 0; ntest = 0
         do i = 1, n
            if (fold_id(i) == f) then
               ntest = ntest + 1; test_index(ntest) = i
            else
               ntrain = ntrain + 1; train_index(ntrain) = i
            end if
         end do
         allocate(x_train(ntrain, size(x, 2)), x_test(ntest, size(x, 2)), &
            st_train(ntrain), sp_train(ntrain), st_test(ntest), sp_test(ntest), &
            e_train(ntrain), e_test(ntest), w_train(ntrain), w_test(ntest), &
            o_train(ntrain), o_test(ntest), s_train(ntrain), s_test(ntest))
         x_train = x(train_index, :); x_test = x(test_index, :)
         st_train = start_time(train_index); sp_train = stop_time(train_index)
         st_test = start_time(test_index); sp_test = stop_time(test_index)
         e_train = event(train_index); e_test = event(test_index)
         w_train = w_all(train_index); w_test = w_all(test_index)
         o_train = o_all(train_index); o_test = o_all(test_index)
         s_train = s_all(train_index); s_test = s_all(test_index)
         call fit_cox_path(x_train, st_train, sp_train, e_train, fold_fit, ctl, &
            weights_in=w_train, offset_in=o_train, strata_in=s_train, &
            lambda_in=result%fit%lambda, efron=ef)
         call predict_with_vector_offset(fold_fit, x_test, o_test, pred, status, 'link')
         if (status /= glmnet_success) then
            result%status = status
            return
         end if
         do l = 1, result%fit%nlambda
            result%fold_values(f, l) = concordance_index(sp_test, e_test, pred(:, 1, l), w_test)
            result%predictions(test_index, 1, l) = pred(:, 1, l)
         end do
         deallocate(train_index, test_index, x_train, x_test, st_train, sp_train, &
            st_test, sp_test, e_train, e_test, w_train, w_test, o_train, o_test, &
            s_train, s_test, pred)
      end do
      call finalize_cv(result, higher_is_better=.true.)
   end subroutine cv_cox

   subroutine predict_with_vector_offset(fit, x, offset, prediction, status, ptype)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:), offset(:)
      real(dp), allocatable, intent(out) :: prediction(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(in), optional :: ptype
      real(dp), allocatable :: om(:,:)
      character(len=16) :: kind
      allocate(om(size(x, 1), fit%nout))
      om = 0.0_dp
      om(:, 1) = offset
      kind = 'response'; if (present(ptype)) kind = ptype
      call predict_glmnet(fit, x, prediction, status, prediction_type=kind, offset=om)
   end subroutine predict_with_vector_offset

   subroutine finalize_cv(result, higher_is_better)
      type(glmnet_cv_result), intent(inout) :: result
      logical, intent(in) :: higher_is_better
      integer :: l, k, best
      real(dp) :: threshold
      k = result%nfolds
      allocate(result%cv_mean(result%fit%nlambda), result%cv_sd(result%fit%nlambda))
      do l = 1, result%fit%nlambda
         result%cv_mean(l) = sum(result%fold_values(:, l)) / real(k, dp)
         if (k > 1) then
            result%cv_sd(l) = sqrt(sum((result%fold_values(:, l) - result%cv_mean(l)) ** 2) / &
               real(k - 1, dp)) / sqrt(real(k, dp))
         else
            result%cv_sd(l) = 0.0_dp
         end if
      end do
      if (higher_is_better) then
         best = maxloc(result%cv_mean, dim=1)
         threshold = result%cv_mean(best) - result%cv_sd(best)
         result%index_1se = best
         do l = 1, best
            if (result%cv_mean(l) >= threshold) then
               result%index_1se = l
               exit
            end if
         end do
      else
         best = minloc(result%cv_mean, dim=1)
         threshold = result%cv_mean(best) + result%cv_sd(best)
         result%index_1se = best
         do l = 1, best
            if (result%cv_mean(l) <= threshold) then
               result%index_1se = l
               exit
            end if
         end do
      end if
      result%index_min = best
      result%lambda_min = result%fit%lambda(best)
      result%lambda_1se = result%fit%lambda(result%index_1se)
   end subroutine finalize_cv

   subroutine build_predmat(result, predictions, status)
      type(glmnet_cv_result), intent(in) :: result
      real(dp), allocatable, intent(out) :: predictions(:,:,:)
      integer, intent(out) :: status
      if (.not. allocated(result%predictions)) then
         status = glmnet_invalid_argument
         allocate(predictions(0, 0, 0))
         return
      end if
      allocate(predictions(size(result%predictions, 1), size(result%predictions, 2), &
         size(result%predictions, 3)))
      predictions = result%predictions
      status = glmnet_success
   end subroutine build_predmat
end module glmnet_cv
