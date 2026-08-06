! SPDX-License-Identifier: GPL-3.0-only
module mass_ordinal
  use rrcov_kinds, only : dp
  use rrcov_stats, only : normal_cdf
  use mass_types, only : ordinal_model, mass_success, mass_invalid_argument, mass_dimension_error
  use mass_math, only : pi_dp, logistic_cdf, bfgs_minimize, numerical_hessian, covariance_from_hessian
  implicit none
  private
  public :: polr_fit, polr_predict
contains

  subroutine polr_fit(x, y, model, method, weights, offset, start, maxit, tolerance)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: y(:)
    type(ordinal_model), intent(out) :: model
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: weights(:), offset(:), start(:), tolerance
    integer, intent(in), optional :: maxit
    character(len=16) :: link
    real(dp), allocatable :: w(:), off(:), par(:), hessian(:, :), covariance(:, :)
    real(dp) :: fval, tol
    integer :: n, p, categories, q, i, status, iterations, mit, hstatus

    n = size(y)
    p = size(x, 2)
    if (size(x, 1) /= n .or. n < 1 .or. minval(y) < 1) then
      model%status = mass_dimension_error
      return
    end if
    categories = maxval(y)
    q = categories - 1
    if (categories < 3) then
      model%status = mass_invalid_argument
      return
    end if
    allocate(w(n), off(n))
    w = 1.0_dp
    off = 0.0_dp
    if (present(weights)) then
      if (size(weights) /= n .or. any(weights < 0.0_dp)) then
        model%status = mass_invalid_argument
        return
      end if
      w = weights
    end if
    if (present(offset)) then
      if (size(offset) /= n) then
        model%status = mass_invalid_argument
        return
      end if
      off = offset
    end if
    link = "logistic"
    if (present(method)) link = trim(method)
    mit = 500
    if (present(maxit)) mit = maxit
    tol = 1.0e-7_dp
    if (present(tolerance)) tol = tolerance

    allocate(par(p + q))
    par = 0.0_dp
    if (present(start)) then
      if (size(start) /= p + q) then
        model%status = mass_invalid_argument
        return
      end if
      par(1:p) = start(1:p)
      par(p + 1) = start(p + 1)
      do i = 2, q
        par(p + i) = log(max(start(p + i) - start(p + i - 1), 1.0e-6_dp))
      end do
    else
      par(p + 1) = link_quantile(1.0_dp / real(categories, dp), link)
      do i = 2, q
        par(p + i) = log(max(link_quantile(real(i, dp) / real(categories, dp), link) - &
          link_quantile(real(i - 1, dp) / real(categories, dp), link), 1.0e-3_dp))
      end do
    end if

    call bfgs_minimize(objective, par, fval, status, iterations, mit, tol)
    allocate(model%coefficients(p), model%zeta(q))
    model%coefficients = par(1:p)
    model%zeta(1) = par(p + 1)
    do i = 2, q
      model%zeta(i) = model%zeta(i - 1) + exp(par(p + i))
    end do
    call numerical_hessian(objective, par, hessian)
    covariance = covariance_from_hessian(hessian, hstatus)
    model%covariance = covariance
    model%n_categories = categories
    model%iterations = iterations
    model%status = status
    model%log_likelihood = -fval
    model%aic = 2.0_dp * fval + 2.0_dp * real(p + q, dp)
    model%link = link
  contains
    function objective(parameters) result(value)
      real(dp), intent(in) :: parameters(:)
      real(dp) :: value, eta, lower, upper, probability
      real(dp), allocatable :: zeta(:)
      integer :: observation, cut
      allocate(zeta(q))
      zeta(1) = parameters(p + 1)
      do cut = 2, q
        zeta(cut) = zeta(cut - 1) + exp(parameters(p + cut))
      end do
      value = 0.0_dp
      do observation = 1, n
        eta = off(observation)
        if (p > 0) eta = eta + dot_product(x(observation, :), parameters(1:p))
        if (y(observation) == 1) then
          lower = 0.0_dp
        else
          lower = link_cdf(max(-100.0_dp, zeta(y(observation) - 1) - eta), link)
        end if
        if (y(observation) == categories) then
          upper = 1.0_dp
        else
          upper = link_cdf(min(100.0_dp, zeta(y(observation)) - eta), link)
        end if
        probability = upper - lower
        if (probability <= 0.0_dp) then
          value = huge(1.0_dp) / 100.0_dp
          return
        end if
        value = value - w(observation) * log(probability)
      end do
    end function objective
  end subroutine polr_fit

  subroutine polr_predict(model, x, predicted, probabilities, offset, status)
    type(ordinal_model), intent(in) :: model
    real(dp), intent(in) :: x(:, :)
    integer, allocatable, intent(out) :: predicted(:)
    real(dp), allocatable, intent(out) :: probabilities(:, :)
    real(dp), intent(in), optional :: offset(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: off(:), cumulative(:)
    real(dp) :: eta
    integer :: n, categories, i, j

    n = size(x, 1)
    categories = model%n_categories
    if (size(x, 2) /= size(model%coefficients) .or. categories < 2) then
      allocate(predicted(0), probabilities(0, 0))
      if (present(status)) status = mass_dimension_error
      return
    end if
    allocate(off(n))
    off = 0.0_dp
    if (present(offset)) then
      if (size(offset) /= n) then
        allocate(predicted(0), probabilities(0, 0))
        if (present(status)) status = mass_dimension_error
        return
      end if
      off = offset
    end if
    allocate(predicted(n), probabilities(n, categories), cumulative(categories - 1))
    do i = 1, n
      eta = off(i) + dot_product(x(i, :), model%coefficients)
      do j = 1, categories - 1
        cumulative(j) = link_cdf(model%zeta(j) - eta, model%link)
      end do
      probabilities(i, 1) = cumulative(1)
      do j = 2, categories - 1
        probabilities(i, j) = cumulative(j) - cumulative(j - 1)
      end do
      probabilities(i, categories) = 1.0_dp - cumulative(categories - 1)
      predicted(i) = maxloc(probabilities(i, :), dim=1)
    end do
    if (present(status)) status = mass_success
  end subroutine polr_predict

  pure elemental function link_cdf(x, link) result(value)
    real(dp), intent(in) :: x
    character(len=*), intent(in) :: link
    real(dp) :: value
    select case (trim(link))
    case ("probit")
      value = normal_cdf(x)
    case ("loglog")
      value = exp(-exp(-x))
    case ("cloglog")
      value = 1.0_dp - exp(-exp(x))
    case ("cauchit")
      value = 0.5_dp + atan(x) / pi_dp
    case default
      value = logistic_cdf(x)
    end select
  end function link_cdf

  function link_quantile(p, link) result(value)
    use mass_math, only : normal_quantile
    real(dp), intent(in) :: p
    character(len=*), intent(in) :: link
    real(dp) :: value
    select case (trim(link))
    case ("probit")
      value = normal_quantile(p)
    case ("loglog")
      value = -log(-log(p))
    case ("cloglog")
      value = log(-log(1.0_dp - p))
    case ("cauchit")
      value = tan(pi_dp * (p - 0.5_dp))
    case default
      value = log(p / (1.0_dp - p))
    end select
  end function link_quantile

end module mass_ordinal
