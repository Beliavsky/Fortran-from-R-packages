! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_da
  use rrcov_kinds, only : dp
  use rrcov_types, only : lda_model, qda_model, covariance_result, &
    rrcov_success, rrcov_invalid_argument, rrcov_dimension_error
  use rrcov_stats, only : mean_vector, covariance_matrix
  use rrcov_linalg, only : symmetric_inverse, make_positive_definite, log_determinant
  use rrcov_robust, only : robust_covariance
  implicit none
  private
  public :: lda_classic_fit, lda_cov_fit, linda_fit, lda_pp_fit, lda_predict
  public :: qda_classic_fit, qda_cov_fit, qda_predict, confusion_matrix
contains
  subroutine lda_classic_fit(x, grouping, model, priors)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(lda_model), intent(out) :: model
    real(dp), intent(in), optional :: priors(:)
    real(dp), allocatable :: pooled(:, :), group_data(:, :), residual(:), group_cov(:, :)
    integer, allocatable :: rows(:)
    integer :: n, p, g, i, count_rows, status
    call initialize_lda(x, grouping, model, priors, status)
    if (status /= rrcov_success) return
    n = size(x, 1)
    p = size(x, 2)
    allocate(pooled(p, p), residual(p))
    pooled = 0.0_dp
    do g = 1, model%n_groups
      count_rows = count(grouping == model%labels(g))
      allocate(rows(count_rows), group_data(count_rows, p))
      rows = pack([(i, i=1, n)], grouping == model%labels(g))
      group_data = x(rows, :)
      model%means(g, :) = mean_vector(group_data)
      group_cov = covariance_matrix(group_data, unbiased=.false., status=status)
      pooled = pooled + real(count_rows, dp) * group_cov
      deallocate(rows, group_data, group_cov)
    end do
    pooled = pooled / real(max(1, n - model%n_groups), dp)
    model%covariance = make_positive_definite(pooled, 1.0e-9_dp)
    model%inverse = symmetric_inverse(model%covariance, status)
    call finalize_lda_coefficients(model)
    model%status = status
    model%method = "Classical linear discriminant analysis"
  end subroutine lda_classic_fit

  subroutine lda_cov_fit(x, grouping, model, covariance_method, priors, alpha, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(lda_model), intent(out) :: model
    character(len=*), intent(in), optional :: covariance_method
    real(dp), intent(in), optional :: priors(:), alpha
    integer, intent(in), optional :: nsamp, seed
    type(covariance_result) :: group_estimate, pooled_estimate
    real(dp), allocatable :: group_data(:, :), residuals(:, :)
    integer, allocatable :: rows(:)
    character(len=24) :: method
    integer :: n, p, g, i, count_rows, status
    method = "mcd"
    if (present(covariance_method)) method = covariance_method
    call initialize_lda(x, grouping, model, priors, status)
    if (status /= rrcov_success) return
    n = size(x, 1)
    p = size(x, 2)
    allocate(residuals(n, p))
    do g = 1, model%n_groups
      count_rows = count(grouping == model%labels(g))
      allocate(rows(count_rows), group_data(count_rows, p))
      rows = pack([(i, i=1, n)], grouping == model%labels(g))
      group_data = x(rows, :)
      call robust_covariance(group_data, method, group_estimate, alpha=alpha, nsamp=nsamp, seed=seed)
      model%means(g, :) = group_estimate%center
      do i = 1, count_rows
        residuals(rows(i), :) = group_data(i, :) - group_estimate%center
      end do
      deallocate(rows, group_data)
    end do
    call robust_covariance(residuals, method, pooled_estimate, alpha=alpha, nsamp=nsamp, seed=seed)
    model%covariance = make_positive_definite(pooled_estimate%covariance, 1.0e-9_dp)
    model%inverse = symmetric_inverse(model%covariance, status)
    call finalize_lda_coefficients(model)
    model%status = pooled_estimate%status
    model%method = "Robust covariance LDA (" // trim(method) // ")"
  end subroutine lda_cov_fit

  subroutine linda_fit(x, grouping, model, priors, alpha, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(lda_model), intent(out) :: model
    real(dp), intent(in), optional :: priors(:), alpha
    integer, intent(in), optional :: nsamp, seed
    call lda_cov_fit(x, grouping, model, "mcd", priors=priors, alpha=alpha, nsamp=nsamp, seed=seed)
    model%method = "Linda robust linear discriminant analysis"
  end subroutine linda_fit

  subroutine lda_pp_fit(x, grouping, model, priors, robust_method, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(lda_model), intent(out) :: model
    real(dp), intent(in), optional :: priors(:)
    character(len=*), intent(in), optional :: robust_method
    integer, intent(in), optional :: nsamp, seed
    character(len=24) :: method
    method = "ogk"
    if (present(robust_method)) method = robust_method
    call lda_cov_fit(x, grouping, model, method, priors=priors, nsamp=nsamp, seed=seed)
    model%method = "Projection-pursuit robust LDA"
  end subroutine lda_pp_fit

  subroutine lda_predict(model, x, predicted, posterior, scores, status)
    type(lda_model), intent(in) :: model
    real(dp), intent(in) :: x(:, :)
    integer, allocatable, intent(out) :: predicted(:)
    real(dp), allocatable, intent(out), optional :: posterior(:, :), scores(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: work(:, :), probabilities(:, :)
    real(dp) :: maximum, denominator
    integer :: n, g, i, best
    n = size(x, 1)
    g = model%n_groups
    allocate(predicted(n), work(n, g))
    if (size(x, 2) /= model%n_features .or. model%status /= rrcov_success) then
      predicted = 0
      if (present(status)) status = rrcov_dimension_error
      if (present(scores)) allocate(scores(0, 0))
      if (present(posterior)) allocate(posterior(0, 0))
      return
    end if
    work = matmul(x, model%coefficients)
    do i = 1, n
      work(i, :) = work(i, :) + model%constants
      best = maxloc(work(i, :), dim=1)
      predicted(i) = model%labels(best)
    end do
    if (present(scores)) scores = work
    if (present(posterior)) then
      allocate(probabilities(n, g))
      do i = 1, n
        maximum = maxval(work(i, :))
        probabilities(i, :) = exp(work(i, :) - maximum)
        denominator = sum(probabilities(i, :))
        probabilities(i, :) = probabilities(i, :) / max(denominator, tiny(1.0_dp))
      end do
      posterior = probabilities
    end if
    if (present(status)) status = rrcov_success
  end subroutine lda_predict

  subroutine qda_classic_fit(x, grouping, model, priors)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(qda_model), intent(out) :: model
    real(dp), intent(in), optional :: priors(:)
    real(dp), allocatable :: group_data(:, :)
    integer, allocatable :: rows(:)
    integer :: n, p, g, i, count_rows, status
    call initialize_qda(x, grouping, model, priors, status)
    if (status /= rrcov_success) return
    n = size(x, 1)
    p = size(x, 2)
    do g = 1, model%n_groups
      count_rows = count(grouping == model%labels(g))
      allocate(rows(count_rows), group_data(count_rows, p))
      rows = pack([(i, i=1, n)], grouping == model%labels(g))
      group_data = x(rows, :)
      model%means(g, :) = mean_vector(group_data)
      model%covariance(:, :, g) = covariance_matrix(group_data, unbiased=.true., status=status)
      model%covariance(:, :, g) = make_positive_definite(model%covariance(:, :, g), 1.0e-9_dp)
      model%inverse(:, :, g) = symmetric_inverse(model%covariance(:, :, g), status)
      model%log_determinant(g) = log_determinant(model%covariance(:, :, g), status)
      deallocate(rows, group_data)
    end do
    model%status = rrcov_success
    model%method = "Classical quadratic discriminant analysis"
  end subroutine qda_classic_fit

  subroutine qda_cov_fit(x, grouping, model, covariance_method, priors, alpha, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(qda_model), intent(out) :: model
    character(len=*), intent(in), optional :: covariance_method
    real(dp), intent(in), optional :: priors(:), alpha
    integer, intent(in), optional :: nsamp, seed
    type(covariance_result) :: estimate
    real(dp), allocatable :: group_data(:, :)
    integer, allocatable :: rows(:)
    character(len=24) :: method
    integer :: n, p, g, i, count_rows, status
    method = "mcd"
    if (present(covariance_method)) method = covariance_method
    call initialize_qda(x, grouping, model, priors, status)
    if (status /= rrcov_success) return
    n = size(x, 1)
    p = size(x, 2)
    do g = 1, model%n_groups
      count_rows = count(grouping == model%labels(g))
      allocate(rows(count_rows), group_data(count_rows, p))
      rows = pack([(i, i=1, n)], grouping == model%labels(g))
      group_data = x(rows, :)
      call robust_covariance(group_data, method, estimate, alpha=alpha, nsamp=nsamp, seed=seed)
      model%means(g, :) = estimate%center
      model%covariance(:, :, g) = make_positive_definite(estimate%covariance, 1.0e-9_dp)
      model%inverse(:, :, g) = symmetric_inverse(model%covariance(:, :, g), status)
      model%log_determinant(g) = log_determinant(model%covariance(:, :, g), status)
      deallocate(rows, group_data)
    end do
    model%status = rrcov_success
    model%method = "Robust covariance QDA (" // trim(method) // ")"
  end subroutine qda_cov_fit

  subroutine qda_predict(model, x, predicted, posterior, scores, status)
    type(qda_model), intent(in) :: model
    real(dp), intent(in) :: x(:, :)
    integer, allocatable, intent(out) :: predicted(:)
    real(dp), allocatable, intent(out), optional :: posterior(:, :), scores(:, :)
    integer, intent(out), optional :: status
    real(dp), allocatable :: work(:, :), probabilities(:, :), delta(:)
    real(dp) :: maximum, denominator
    integer :: n, p, g, i, best
    n = size(x, 1)
    p = size(x, 2)
    allocate(predicted(n), work(n, model%n_groups), delta(p))
    if (p /= model%n_features .or. model%status /= rrcov_success) then
      predicted = 0
      if (present(status)) status = rrcov_dimension_error
      if (present(scores)) allocate(scores(0, 0))
      if (present(posterior)) allocate(posterior(0, 0))
      return
    end if
    do i = 1, n
      do g = 1, model%n_groups
        delta = x(i, :) - model%means(g, :)
        work(i, g) = -0.5_dp * dot_product(delta, matmul(model%inverse(:, :, g), delta)) - &
          0.5_dp * model%log_determinant(g) + log(model%priors(g))
      end do
      best = maxloc(work(i, :), dim=1)
      predicted(i) = model%labels(best)
    end do
    if (present(scores)) scores = work
    if (present(posterior)) then
      allocate(probabilities(n, model%n_groups))
      do i = 1, n
        maximum = maxval(work(i, :))
        probabilities(i, :) = exp(work(i, :) - maximum)
        denominator = sum(probabilities(i, :))
        probabilities(i, :) = probabilities(i, :) / max(denominator, tiny(1.0_dp))
      end do
      posterior = probabilities
    end if
    if (present(status)) status = rrcov_success
  end subroutine qda_predict

  subroutine confusion_matrix(actual, predicted, labels, table, error_rate, status)
    integer, intent(in) :: actual(:), predicted(:)
    integer, allocatable, intent(out) :: labels(:), table(:, :)
    real(dp), intent(out) :: error_rate
    integer, intent(out) :: status
    integer :: i, a, b
    if (size(actual) /= size(predicted) .or. size(actual) == 0) then
      allocate(labels(0), table(0, 0))
      error_rate = 0.0_dp
      status = rrcov_dimension_error
      return
    end if
    call unique_labels([actual, predicted], labels)
    allocate(table(size(labels), size(labels)))
    table = 0
    do i = 1, size(actual)
      a = find_label(labels, actual(i))
      b = find_label(labels, predicted(i))
      table(a, b) = table(a, b) + 1
    end do
    error_rate = real(count(actual /= predicted), dp) / real(size(actual), dp)
    status = rrcov_success
  end subroutine confusion_matrix

  subroutine initialize_lda(x, grouping, model, priors, status)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(lda_model), intent(out) :: model
    real(dp), intent(in), optional :: priors(:)
    integer, intent(out) :: status
    integer :: g, n, p
    n = size(x, 1)
    p = size(x, 2)
    if (size(grouping) /= n .or. n < 2 .or. p < 1) then
      model%status = rrcov_dimension_error
      status = model%status
      return
    end if
    call unique_labels(grouping, model%labels)
    model%n_groups = size(model%labels)
    model%n_features = p
    if (model%n_groups < 2) then
      model%status = rrcov_invalid_argument
      status = model%status
      return
    end if
    allocate(model%means(model%n_groups, p), model%covariance(p, p), model%inverse(p, p))
    allocate(model%priors(model%n_groups), model%coefficients(p, model%n_groups), model%constants(model%n_groups))
    if (present(priors)) then
      if (size(priors) == model%n_groups .and. all(priors > 0.0_dp)) then
        model%priors = priors / sum(priors)
      else
        model%status = rrcov_invalid_argument
        status = model%status
        return
      end if
    else
      do g = 1, model%n_groups
        model%priors(g) = real(count(grouping == model%labels(g)), dp) / real(n, dp)
      end do
    end if
    model%means = 0.0_dp
    model%covariance = 0.0_dp
    model%inverse = 0.0_dp
    model%coefficients = 0.0_dp
    model%constants = 0.0_dp
    model%status = rrcov_success
    status = rrcov_success
  end subroutine initialize_lda

  subroutine finalize_lda_coefficients(model)
    type(lda_model), intent(inout) :: model
    integer :: g
    do g = 1, model%n_groups
      model%coefficients(:, g) = matmul(model%inverse, model%means(g, :))
      model%constants(g) = -0.5_dp * dot_product(model%means(g, :), model%coefficients(:, g)) + &
        log(model%priors(g))
    end do
  end subroutine finalize_lda_coefficients

  subroutine initialize_qda(x, grouping, model, priors, status)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: grouping(:)
    type(qda_model), intent(out) :: model
    real(dp), intent(in), optional :: priors(:)
    integer, intent(out) :: status
    integer :: g, n, p
    n = size(x, 1)
    p = size(x, 2)
    if (size(grouping) /= n .or. n < 2 .or. p < 1) then
      model%status = rrcov_dimension_error
      status = model%status
      return
    end if
    call unique_labels(grouping, model%labels)
    model%n_groups = size(model%labels)
    model%n_features = p
    if (model%n_groups < 2) then
      model%status = rrcov_invalid_argument
      status = model%status
      return
    end if
    allocate(model%means(model%n_groups, p), model%covariance(p, p, model%n_groups))
    allocate(model%inverse(p, p, model%n_groups), model%log_determinant(model%n_groups))
    allocate(model%priors(model%n_groups))
    if (present(priors)) then
      if (size(priors) == model%n_groups .and. all(priors > 0.0_dp)) then
        model%priors = priors / sum(priors)
      else
        model%status = rrcov_invalid_argument
        status = model%status
        return
      end if
    else
      do g = 1, model%n_groups
        model%priors(g) = real(count(grouping == model%labels(g)), dp) / real(n, dp)
      end do
    end if
    model%means = 0.0_dp
    model%covariance = 0.0_dp
    model%inverse = 0.0_dp
    model%log_determinant = 0.0_dp
    model%status = rrcov_success
    status = rrcov_success
  end subroutine initialize_qda

  subroutine unique_labels(grouping, labels)
    integer, intent(in) :: grouping(:)
    integer, allocatable, intent(out) :: labels(:)
    integer, allocatable :: work(:)
    integer :: i, count_labels
    allocate(work(size(grouping)))
    count_labels = 0
    do i = 1, size(grouping)
      if (count_labels == 0 .or. .not. any(work(1:count_labels) == grouping(i))) then
        count_labels = count_labels + 1
        work(count_labels) = grouping(i)
      end if
    end do
    allocate(labels(count_labels))
    labels = work(1:count_labels)
  end subroutine unique_labels

  pure function find_label(labels, value) result(index)
    integer, intent(in) :: labels(:), value
    integer :: index, i
    index = 0
    do i = 1, size(labels)
      if (labels(i) == value) then
        index = i
        return
      end if
    end do
  end function find_label
end module rrcov_da
