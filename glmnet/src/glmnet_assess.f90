! SPDX-License-Identifier: GPL-2.0-only
module glmnet_assess
   use glmnet_kinds, only : dp, glmnet_eps
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument
   use glmnet_types, only : glmnet_path_result, glmnet_assessment_result, &
      glmnet_roc_result, glmnet_family_gaussian, glmnet_family_mgaussian, &
      glmnet_family_binomial, glmnet_family_poisson, glmnet_family_multinomial, &
      glmnet_family_cox
   use glmnet_predict, only : predict_glmnet
   use glmnet_utils, only : normalize_weights, safe_log, sort_indices_real
   use glmnet_cox, only : concordance_index
   implicit none
   private
   public :: assess_glmnet, assess_multinomial, assess_cox
   public :: auc, roc_glmnet, confusion_glmnet, glmnet_measures
contains
   subroutine assess_glmnet(fit, x, y, assessment, measure, weights, offset)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:), y(:)
      type(glmnet_assessment_result), intent(out) :: assessment
      character(len=*), intent(in), optional :: measure
      real(dp), intent(in), optional :: weights(:), offset(:)
      real(dp), allocatable :: prediction(:,:,:), w(:), offset_matrix(:,:)
      character(len=16) :: metric
      integer :: status, l, n
      n = size(x, 1)
      assessment%status = glmnet_success
      if (size(y) /= n) then
         assessment%status = glmnet_invalid_argument
         return
      end if
      metric = default_measure(fit%family_code)
      if (present(measure)) metric = adjustl(measure)
      assessment%measure = metric
      call normalize_weights(n, weights, w, status)
      if (status /= glmnet_success) then
         assessment%status = status
         return
      end if
      allocate(offset_matrix(n, fit%nout))
      offset_matrix = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) then
            assessment%status = glmnet_invalid_argument
            return
         end if
         offset_matrix(:, 1) = offset
      end if
      call predict_glmnet(fit, x, prediction, status, prediction_type='response', &
         offset=offset_matrix)
      if (status /= glmnet_success) then
         assessment%status = status
         return
      end if
      allocate(assessment%value(fit%nlambda))
      do l = 1, fit%nlambda
         select case (trim(metric))
         case ('mse')
            assessment%value(l) = sum(w * (y - prediction(:, 1, l)) ** 2)
         case ('mae')
            assessment%value(l) = sum(w * abs(y - prediction(:, 1, l)))
         case ('deviance')
            if (fit%family_code == glmnet_family_binomial) then
               assessment%value(l) = binomial_deviance(y, prediction(:, 1, l), w)
            else if (fit%family_code == glmnet_family_poisson) then
               assessment%value(l) = poisson_deviance(y, prediction(:, 1, l), w)
            else
               assessment%value(l) = sum(w * (y - prediction(:, 1, l)) ** 2)
            end if
         case ('class')
            assessment%value(l) = sum(w * merge(1.0_dp, 0.0_dp, &
               (prediction(:, 1, l) >= 0.5_dp) .neqv. (y >= 0.5_dp)))
         case ('auc')
            assessment%value(l) = auc(y, prediction(:, 1, l), w)
         case default
            assessment%value(l) = sum(w * (y - prediction(:, 1, l)) ** 2)
         end select
      end do
   end subroutine assess_glmnet

   subroutine assess_multinomial(fit, x, class_id, assessment, measure, weights, offset)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: class_id(:)
      type(glmnet_assessment_result), intent(out) :: assessment
      character(len=*), intent(in), optional :: measure
      real(dp), intent(in), optional :: weights(:), offset(:,:)
      real(dp), allocatable :: prediction(:,:,:), w(:), offset_matrix(:,:)
      character(len=16) :: metric
      integer :: status, l, i, predicted
      assessment%status = glmnet_success
      if (fit%family_code /= glmnet_family_multinomial .or. size(class_id) /= size(x, 1)) then
         assessment%status = glmnet_invalid_argument
         return
      end if
      metric = 'deviance'
      if (present(measure)) metric = adjustl(measure)
      assessment%measure = metric
      call normalize_weights(size(x, 1), weights, w, status)
      if (status /= glmnet_success) then
         assessment%status = status
         return
      end if
      allocate(offset_matrix(size(x, 1), fit%nout))
      offset_matrix = 0.0_dp
      if (present(offset)) then
         if (size(offset, 1) /= size(x, 1) .or. size(offset, 2) /= fit%nout) then
            assessment%status = glmnet_invalid_argument
            return
         end if
         offset_matrix = offset
      end if
      call predict_glmnet(fit, x, prediction, status, prediction_type='response', &
         offset=offset_matrix)
      if (status /= glmnet_success) then
         assessment%status = status
         return
      end if
      allocate(assessment%value(fit%nlambda))
      assessment%value = 0.0_dp
      do l = 1, fit%nlambda
         do i = 1, size(x, 1)
            if (trim(metric) == 'class') then
               predicted = maxloc(prediction(i, :, l), dim=1)
               if (predicted /= class_id(i)) assessment%value(l) = assessment%value(l) + w(i)
            else
               assessment%value(l) = assessment%value(l) - 2.0_dp * w(i) * &
                  safe_log(prediction(i, class_id(i), l))
            end if
         end do
      end do
   end subroutine assess_multinomial

   subroutine assess_cox(fit, x, time, event, assessment, weights, offset)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:), time(:)
      integer, intent(in) :: event(:)
      type(glmnet_assessment_result), intent(out) :: assessment
      real(dp), intent(in), optional :: weights(:), offset(:)
      real(dp), allocatable :: prediction(:,:,:), w(:), offset_matrix(:,:)
      integer :: status, l, n
      n = size(x, 1)
      assessment%status = glmnet_success
      assessment%measure = 'C'
      if (size(time) /= n .or. size(event) /= n) then
         assessment%status = glmnet_invalid_argument
         return
      end if
      call normalize_weights(n, weights, w, status)
      if (status /= glmnet_success) then
         assessment%status = status
         return
      end if
      allocate(offset_matrix(n, 1))
      offset_matrix = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) then
            assessment%status = glmnet_invalid_argument
            return
         end if
         offset_matrix(:, 1) = offset
      end if
      call predict_glmnet(fit, x, prediction, status, prediction_type='link', &
         offset=offset_matrix)
      if (status /= glmnet_success) then
         assessment%status = status
         return
      end if
      allocate(assessment%value(fit%nlambda))
      do l = 1, fit%nlambda
         assessment%value(l) = concordance_index(time, event, prediction(:, 1, l), w)
      end do
   end subroutine assess_cox

   function auc(y, score, weights) result(value)
      real(dp), intent(in) :: y(:), score(:)
      real(dp), intent(in), optional :: weights(:)
      real(dp) :: value
      real(dp), allocatable :: w(:)
      real(dp) :: positive_weight, negative_weight, concordant
      integer :: i, j, n
      n = size(y)
      allocate(w(n))
      if (present(weights)) then
         if (size(weights) == n) then
            w = weights
         else
            w = 1.0_dp
         end if
      else
         w = 1.0_dp
      end if
      positive_weight = sum(w, mask=y >= 0.5_dp)
      negative_weight = sum(w, mask=y < 0.5_dp)
      if (positive_weight <= glmnet_eps .or. negative_weight <= glmnet_eps) then
         value = 0.5_dp
         return
      end if
      concordant = 0.0_dp
      do i = 1, n
         if (y(i) < 0.5_dp) cycle
         do j = 1, n
            if (y(j) >= 0.5_dp) cycle
            if (abs(score(i) - score(j)) <= 10.0_dp * glmnet_eps) then
               concordant = concordant + 0.5_dp * w(i) * w(j)
            else if (score(i) > score(j)) then
               concordant = concordant + w(i) * w(j)
            end if
         end do
      end do
      value = concordant / (positive_weight * negative_weight)
   end function auc

   subroutine roc_glmnet(y, score, result, weights)
      real(dp), intent(in) :: y(:), score(:)
      type(glmnet_roc_result), intent(out) :: result
      real(dp), intent(in), optional :: weights(:)
      real(dp), allocatable :: w(:)
      integer, allocatable :: order(:)
      real(dp) :: pos, neg, tp, fp
      integer :: i, n
      result%status = glmnet_success
      n = size(y)
      if (size(score) /= n .or. n < 1) then
         result%status = glmnet_invalid_argument
         return
      end if
      allocate(w(n))
      if (present(weights)) then
         if (size(weights) /= n) then
            result%status = glmnet_invalid_argument
            return
         end if
         w = weights
      else
         w = 1.0_dp
      end if
      call sort_indices_real(score, order, descending=.true.)
      allocate(result%threshold(n + 1), result%false_positive_rate(n + 1), &
         result%true_positive_rate(n + 1))
      pos = sum(w, mask=y >= 0.5_dp)
      neg = sum(w, mask=y < 0.5_dp)
      tp = 0.0_dp
      fp = 0.0_dp
      result%threshold(1) = huge(1.0_dp)
      result%true_positive_rate(1) = 0.0_dp
      result%false_positive_rate(1) = 0.0_dp
      do i = 1, n
         if (y(order(i)) >= 0.5_dp) then
            tp = tp + w(order(i))
         else
            fp = fp + w(order(i))
         end if
         result%threshold(i + 1) = score(order(i))
         result%true_positive_rate(i + 1) = tp / max(pos, glmnet_eps)
         result%false_positive_rate(i + 1) = fp / max(neg, glmnet_eps)
      end do
      result%auc = auc(y, score, w)
   end subroutine roc_glmnet

   subroutine confusion_glmnet(y, probability, threshold, table, status, weights)
      real(dp), intent(in) :: y(:), probability(:), threshold
      real(dp), intent(out) :: table(2, 2)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: weights(:)
      real(dp), allocatable :: w(:)
      integer :: i, actual, predicted
      status = glmnet_success
      table = 0.0_dp
      if (size(y) /= size(probability)) then
         status = glmnet_invalid_argument
         return
      end if
      allocate(w(size(y)))
      if (present(weights)) then
         if (size(weights) /= size(y)) then
            status = glmnet_invalid_argument
            return
         end if
         w = weights
      else
         w = 1.0_dp
      end if
      do i = 1, size(y)
         actual = merge(2, 1, y(i) >= 0.5_dp)
         predicted = merge(2, 1, probability(i) >= threshold)
         table(actual, predicted) = table(actual, predicted) + w(i)
      end do
   end subroutine confusion_glmnet

   subroutine glmnet_measures(family_code, names)
      integer, intent(in) :: family_code
      character(len=16), allocatable, intent(out) :: names(:)
      select case (family_code)
      case (glmnet_family_gaussian, glmnet_family_mgaussian)
         allocate(names(2)); names = ['mse             ', 'mae             ']
      case (glmnet_family_binomial)
         allocate(names(4)); names = ['deviance        ', 'class           ', &
            'auc             ', 'mse             ']
      case (glmnet_family_poisson)
         allocate(names(2)); names = ['deviance        ', 'mse             ']
      case (glmnet_family_multinomial)
         allocate(names(2)); names = ['deviance        ', 'class           ']
      case (glmnet_family_cox)
         allocate(names(2)); names = ['deviance        ', 'C               ']
      case default
         allocate(names(1)); names = ['mse             ']
      end select
   end subroutine glmnet_measures

   pure function default_measure(family_code) result(metric)
      integer, intent(in) :: family_code
      character(len=16) :: metric
      select case (family_code)
      case (glmnet_family_binomial, glmnet_family_poisson, glmnet_family_multinomial)
         metric = 'deviance'
      case (glmnet_family_cox)
         metric = 'C'
      case default
         metric = 'mse'
      end select
   end function default_measure

   pure function binomial_deviance(y, probability, weights) result(value)
      real(dp), intent(in) :: y(:), probability(:), weights(:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, size(y)
         if (y(i) > 0.0_dp) value = value - 2.0_dp * weights(i) * y(i) * &
            safe_log(probability(i))
         if (y(i) < 1.0_dp) value = value - 2.0_dp * weights(i) * &
            (1.0_dp - y(i)) * safe_log(1.0_dp - probability(i))
      end do
   end function binomial_deviance

   pure function poisson_deviance(y, mu, weights) result(value)
      real(dp), intent(in) :: y(:), mu(:), weights(:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, size(y)
         if (y(i) > 0.0_dp) then
            value = value + 2.0_dp * weights(i) * &
               (y(i) * safe_log(y(i) / max(mu(i), glmnet_eps)) - (y(i) - mu(i)))
         else
            value = value + 2.0_dp * weights(i) * mu(i)
         end if
      end do
   end function poisson_deviance
end module glmnet_assess
