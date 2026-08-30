! SPDX-License-Identifier: GPL-2.0-only
module multcomp_glht
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
  use multcomp_kinds, only : dp
  use multcomp_types, only : parm_type, glht_type, mtest_type, &
    confidence_interval_type, global_test_type, alt_two_sided, alt_less, &
    alt_greater, set_collection
  use multcomp_math, only : covariance_to_correlation_mc, p_adjust, &
    normal_pvalue, student_pvalue, generalized_inverse, chi_square_upper, &
    f_upper, stable_order
  use multcomp_maxsets, only : compute_maxsets
  use mvtnorm, only : probability_control, probability_result, quantile_result, &
    pmvt, qmvt
  implicit none
  private

  public :: glht_fit
  public :: glht_identity
  public :: glht_coefficients
  public :: glht_test
  public :: glht_confint
  public :: glht_global_test
  public :: glht_critical_value
  public :: alternative_code

contains

  subroutine glht_fit(parameters, linfct, result, rhs, alternative)
    type(parm_type), intent(in) :: parameters !! Coefficients, covariance, and degrees of freedom for the fitted model.
    real(dp), intent(in) :: linfct(:, :) !! Linear-hypothesis matrix K with one hypothesis per row.
    type(glht_type), intent(out) :: result !! General linear-hypothesis object containing transformed estimates and covariance.
    real(dp), intent(in), optional :: rhs(:) !! Null-hypothesis right-hand sides; defaults to zero for every row of K.
    character(len=*), intent(in), optional :: alternative !! Alternative: two.sided, less, or greater; default two.sided.

    real(dp), allocatable :: sd(:)
    logical :: corr_ok
    integer :: i
    integer :: m
    integer :: p

    result%ok = .false.
    result%message = ''
    if (.not. parameters%ok) then
      result%message = 'parameter object is not valid'
      return
    end if

    m = size(linfct, 1)
    p = size(linfct, 2)
    if (p /= size(parameters%coef) .or. m < 1) then
      result%message = 'linear-function matrix has incompatible dimensions'
      return
    end if

    allocate(result%linfct(m, p), result%rhs(m), result%estimate(m))
    allocate(result%covariance(m, m), result%standard_error(m), result%statistic(m))
    result%linfct = linfct
    result%rhs = 0.0_dp
    if (present(rhs)) then
      if (size(rhs) /= m) then
        result%message = 'rhs length must equal the number of linear hypotheses'
        return
      end if
      result%rhs = rhs
    end if

    result%alternative = alt_two_sided
    if (present(alternative)) then
      result%alternative = alternative_code(alternative)
      if (result%alternative == 0) then
        result%message = 'alternative must be two.sided, less, or greater'
        return
      end if
    end if

    result%estimate = matmul(linfct, parameters%coef)
    result%covariance = matmul(linfct, matmul(parameters%vcov, transpose(linfct)))
    result%covariance = 0.5_dp * (result%covariance + transpose(result%covariance))
    do i = 1, m
      if (result%covariance(i, i) <= 0.0_dp) then
        result%message = 'a linear hypothesis has nonpositive variance'
        return
      end if
      result%standard_error(i) = sqrt(result%covariance(i, i))
    end do
    result%statistic = (result%estimate - result%rhs) / result%standard_error
    call covariance_to_correlation_mc(result%covariance, result%correlation, sd, corr_ok)
    if (.not. corr_ok) then
      result%message = 'failed to standardize hypothesis covariance matrix'
      return
    end if
    result%df = parameters%df
    result%ok = .true.
  end subroutine glht_fit

  subroutine glht_identity(parameters, result, alternative)
    type(parm_type), intent(in) :: parameters !! Parameter object whose coefficients are each tested against zero.
    type(glht_type), intent(out) :: result !! General linear-hypothesis object using an identity K matrix.
    character(len=*), intent(in), optional :: alternative !! Alternative: two.sided, less, or greater; default two.sided.

    real(dp), allocatable :: identity(:, :)
    integer :: i
    integer :: n

    n = size(parameters%coef)
    allocate(identity(n, n))
    identity = 0.0_dp
    do i = 1, n
      identity(i, i) = 1.0_dp
    end do
    if (present(alternative)) then
      call glht_fit(parameters, identity, result, alternative=alternative)
    else
      call glht_fit(parameters, identity, result)
    end if
  end subroutine glht_identity

  subroutine glht_coefficients(parameters, indices, result, alternative)
    type(parm_type), intent(in) :: parameters !! Parameter object whose selected coefficients are to be tested against zero.
    integer, intent(in) :: indices(:) !! One-based coefficient indices to retain, in the requested reporting order.
    type(glht_type), intent(out) :: result !! GLHT object for the selected identity rows, analogous to multcomp::cftest.
    character(len=*), intent(in), optional :: alternative !! Alternative: two.sided, less, or greater; default two.sided.

    real(dp), allocatable :: linfct(:, :)
    integer :: i

    result%ok = .false.
    result%message = ''
    if (.not. parameters%ok) then
      result%message = 'parameter object is not valid'
      return
    end if
    if (size(indices) < 1 .or. any(indices < 1) .or. any(indices > size(parameters%coef))) then
      result%message = 'coefficient indices must lie within the parameter vector'
      return
    end if

    allocate(linfct(size(indices), size(parameters%coef)))
    linfct = 0.0_dp
    do i = 1, size(indices)
      linfct(i, indices(i)) = 1.0_dp
    end do
    if (present(alternative)) then
      call glht_fit(parameters, linfct, result, alternative=alternative)
    else
      call glht_fit(parameters, linfct, result)
    end if
  end subroutine glht_coefficients

  subroutine glht_test(object, method, result, control)
    type(glht_type), intent(in) :: object !! Valid GLHT object to test simultaneously or marginally.
    character(len=*), intent(in) :: method !! univariate, single-step, free, Shaffer, Westfall, or an R p.adjust method.
    type(mtest_type), intent(out) :: result !! Test estimates, statistics, adjusted p-values, and integration error.
    type(probability_control), intent(in), optional :: control !! mvtnorm integration controls for max-type adjustments.

    character(len=:), allocatable :: key
    real(dp), allocatable :: raw(:)
    logical :: adjust_ok
    integer :: i

    result%ok = .false.
    result%message = ''
    if (.not. object%ok) then
      result%message = 'GLHT object is not valid'
      return
    end if

    allocate(result%estimate(size(object%estimate)))
    allocate(result%standard_error(size(object%standard_error)))
    allocate(result%statistic(size(object%statistic)))
    result%estimate = object%estimate
    result%standard_error = object%standard_error
    result%statistic = object%statistic
    result%method = trim(method)

    key = lowercase(trim(adjustl(method)))
    select case (key)
    case ('univariate', 'none')
      allocate(result%pvalue(size(object%statistic)))
      do i = 1, size(object%statistic)
        result%pvalue(i) = marginal_pvalue(object, i)
      end do
    case ('single-step', 'adjusted')
      call single_step_pvalues(object, result%pvalue, result%error, control)
    case ('free')
      call free_pvalues(object, result%pvalue, result%error, control)
    case ('shaffer', 'westfall')
      call closed_pvalues(object, key, result%pvalue, result%error, result%ok, &
        result%message, control)
      if (.not. result%ok) return
    case default
      allocate(raw(size(object%statistic)))
      do i = 1, size(object%statistic)
        raw(i) = marginal_pvalue(object, i)
      end do
      call p_adjust(raw, key, result%pvalue, adjust_ok)
      if (.not. adjust_ok) then
        result%message = 'unknown p-value adjustment method'
        return
      end if
    end select

    result%ok = .true.
  end subroutine glht_test

  subroutine glht_confint(object, level, result, adjusted, control, critical)
    type(glht_type), intent(in) :: object !! Valid GLHT object for which simultaneous confidence intervals are required.
    real(dp), intent(in) :: level !! Family-wise confidence level strictly between zero and one.
    type(confidence_interval_type), intent(out) :: result !! Confidence limits, estimates, critical value, and integration error.
    logical, intent(in), optional :: adjusted !! If false, use a univariate critical value; default true.
    type(probability_control), intent(in), optional :: control !! mvtnorm integration controls for adjusted critical values.
    real(dp), intent(in), optional :: critical !! User-supplied scalar critical value replacing numerical quantile calculation.

    logical :: use_adjusted
    real(dp) :: crit
    real(dp) :: error
    logical :: critical_ok
    integer :: n

    result%ok = .false.
    result%message = ''
    if (.not. object%ok) then
      result%message = 'GLHT object is not valid'
      return
    end if
    if (level <= 0.0_dp .or. level >= 1.0_dp) then
      result%message = 'confidence level must lie strictly between zero and one'
      return
    end if

    use_adjusted = .true.
    if (present(adjusted)) use_adjusted = adjusted
    if (present(critical)) then
      crit = critical
      error = 0.0_dp
      critical_ok = .true.
    else
      call glht_critical_value(object, level, use_adjusted, crit, error, critical_ok, control)
    end if
    if (.not. critical_ok) then
      result%message = 'failed to compute the requested critical value'
      return
    end if

    n = size(object%estimate)
    allocate(result%estimate(n), result%lower(n), result%upper(n))
    result%estimate = object%estimate
    select case (object%alternative)
    case (alt_less)
      result%lower = ieee_value(0.0_dp, ieee_negative_inf)
      result%upper = object%estimate + crit * object%standard_error
    case (alt_greater)
      result%lower = object%estimate + crit * object%standard_error
      result%upper = ieee_value(0.0_dp, ieee_positive_inf)
    case default
      result%lower = object%estimate - crit * object%standard_error
      result%upper = object%estimate + crit * object%standard_error
    end select
    result%level = level
    result%critical = crit
    result%error = error
    result%adjusted = use_adjusted
    result%ok = .true.
  end subroutine glht_confint

  subroutine glht_critical_value(object, level, adjusted, critical, error, ok, control)
    type(glht_type), intent(in) :: object !! Valid GLHT object defining correlation and alternative direction.
    real(dp), intent(in) :: level !! Probability content used to define the max-type critical value.
    logical, intent(in) :: adjusted !! If true use all hypotheses; otherwise use a one-dimensional reference distribution.
    real(dp), intent(out) :: critical !! Computed scalar critical value.
    real(dp), intent(out) :: error !! Quantile probability error reported by mvtnorm.
    logical, intent(out) :: ok !! True when the quantile solver converged.
    type(probability_control), intent(in), optional :: control !! mvtnorm integration controls for quantile evaluation.

    type(probability_control) :: ctl
    type(quantile_result) :: quant
    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: delta(:)
    character(len=16) :: tail

    ctl = probability_control()
    if (present(control)) ctl = control
    if (adjusted) then
      call integration_correlation(object%correlation, correlation)
      allocate(delta(size(object%estimate)))
    else
      allocate(correlation(1, 1), delta(1))
      correlation(1, 1) = 1.0_dp
    end if
    delta = 0.0_dp

    select case (object%alternative)
    case (alt_less)
      tail = 'lower.tail'
    case (alt_greater)
      tail = 'upper.tail'
    case default
      tail = 'both.tails'
    end select

    quant = qmvt(level, delta, correlation, object%df, tail=tail, control=ctl)
    critical = quant%quantile
    error = quant%error
    ok = quant%converged
  end subroutine glht_critical_value

  subroutine glht_global_test(object, result, request_f)
    type(glht_type), intent(in) :: object !! Valid GLHT object defining the joint null hypothesis.
    type(global_test_type), intent(out) :: result !! Wald chi-square or F test result for all hypotheses jointly.
    logical, intent(in), optional :: request_f !! If true and df is positive, report an F test instead of chi-square.

    real(dp), allocatable :: difference(:)
    real(dp), allocatable :: pinv(:, :)
    logical :: inverse_ok
    logical :: use_f
    integer :: rank

    result%ok = .false.
    result%message = ''
    if (.not. object%ok) then
      result%message = 'GLHT object is not valid'
      return
    end if

    difference = object%estimate - object%rhs
    call generalized_inverse(object%covariance, pinv, rank, inverse_ok)
    if (.not. inverse_ok .or. rank < 1) then
      result%message = 'hypothesis covariance has zero numerical rank'
      return
    end if

    result%ssh = dot_product(difference, matmul(pinv, difference))
    result%rank = rank
    use_f = .false.
    if (present(request_f)) use_f = request_f .and. object%df > 0.0_dp
    result%f_test = use_f
    if (use_f) then
      result%statistic = result%ssh / real(rank, dp)
      result%denominator_df = object%df
      result%pvalue = f_upper(result%statistic, real(rank, dp), object%df)
    else
      result%statistic = result%ssh
      result%pvalue = chi_square_upper(result%ssh, real(rank, dp))
    end if
    result%ok = .true.
  end subroutine glht_global_test

  integer function alternative_code(name) result(code)
    character(len=*), intent(in) :: name !! Alternative text: two.sided, less, greater, or common punctuation variants.

    character(len=:), allocatable :: key

    key = lowercase(trim(adjustl(name)))
    select case (key)
    case ('two.sided', 'two-sided', 'two_sided', '==', '=')
      code = alt_two_sided
    case ('less', '<', '<=')
      code = alt_less
    case ('greater', '>', '>=')
      code = alt_greater
    case default
      code = 0
    end select
  end function alternative_code

  real(dp) function marginal_pvalue(object, index) result(p)
    type(glht_type), intent(in) :: object !! GLHT object supplying statistic, df, and alternative direction.
    integer, intent(in) :: index !! One-based hypothesis index whose marginal p-value is requested.

    if (object%df > 0.0_dp) then
      p = student_pvalue(object%statistic(index), object%df, object%alternative)
    else
      p = normal_pvalue(object%statistic(index), object%alternative)
    end if
  end function marginal_pvalue

  subroutine single_step_pvalues(object, pvalue, error, control)
    type(glht_type), intent(in) :: object !! GLHT object whose max-type adjusted p-values are required.
    real(dp), allocatable, intent(out) :: pvalue(:) !! Single-step family-wise adjusted p-values.
    real(dp), intent(out) :: error !! Largest mvtnorm integration error over hypotheses.
    type(probability_control), intent(in), optional :: control !! mvtnorm integration controls.

    type(probability_control) :: ctl
    type(probability_result) :: probability
    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: delta(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    real(dp) :: q
    integer :: i
    integer :: n

    ctl = probability_control()
    if (present(control)) ctl = control
    n = size(object%statistic)
    call integration_correlation(object%correlation, correlation)
    allocate(pvalue(n), lower(n), upper(n), delta(n))
    delta = 0.0_dp
    error = 0.0_dp
    do i = 1, n
      q = object%statistic(i)
      select case (object%alternative)
      case (alt_less)
        lower = q
        upper = 1.0e100_dp
      case (alt_greater)
        lower = -1.0e100_dp
        upper = q
      case default
        lower = -abs(q)
        upper = abs(q)
      end select
      probability = pmvt(lower, upper, delta, correlation, object%df, ctl)
      pvalue(i) = min(1.0_dp, max(0.0_dp, 1.0_dp - probability%value))
      error = max(error, probability%error)
    end do
  end subroutine single_step_pvalues

  subroutine free_pvalues(object, pvalue, error, control)
    type(glht_type), intent(in) :: object !! GLHT object for the free step-down max-test adjustment.
    real(dp), allocatable, intent(out) :: pvalue(:) !! Free step-down adjusted p-values in original order.
    real(dp), intent(out) :: error !! Largest mvtnorm integration error encountered.
    type(probability_control), intent(in), optional :: control !! mvtnorm integration controls.

    type(glht_type) :: subset
    real(dp), allocatable :: subp(:)
    integer, allocatable :: remaining(:)
    integer, allocatable :: next_remaining(:)
    real(dp) :: suberror
    real(dp) :: step_p
    integer :: i
    integer :: j
    integer :: n
    integer :: remove_at

    n = size(object%estimate)
    allocate(pvalue(n), remaining(n))
    pvalue = 0.0_dp
    remaining = [(i, i = 1, n)]
    error = 0.0_dp

    do while (size(remaining) > 0)
      call subset_glht(object, remaining, subset)
      call single_step_pvalues(subset, subp, suberror, control)
      error = max(error, suberror)
      remove_at = minloc(subp, dim=1)
      step_p = subp(remove_at)
      do i = 1, size(remaining)
        pvalue(remaining(i)) = max(pvalue(remaining(i)), step_p)
      end do
      if (size(remaining) == 1) exit
      allocate(next_remaining(size(remaining) - 1))
      j = 0
      do i = 1, size(remaining)
        if (i == remove_at) cycle
        j = j + 1
        next_remaining(j) = remaining(i)
      end do
      call move_alloc(next_remaining, remaining)
      if (allocated(subp)) deallocate(subp)
    end do
    pvalue = min(1.0_dp, max(0.0_dp, pvalue))
  end subroutine free_pvalues

  subroutine closed_pvalues(object, method, pvalue, error, ok, message, control)
    type(glht_type), intent(in) :: object !! GLHT object for Shaffer or Westfall closed testing.
    character(len=*), intent(in) :: method !! Closed-test method name, either shaffer or westfall.
    real(dp), allocatable, intent(out) :: pvalue(:) !! Closed-testing adjusted p-values in original hypothesis order.
    real(dp), intent(out) :: error !! Largest mvtnorm integration error encountered.
    logical, intent(out) :: ok !! True when maximal sets and all local tests are computed successfully.
    character(len=*), intent(out) :: message !! Failure explanation; blank on success.
    type(probability_control), intent(in), optional :: control !! mvtnorm integration controls for Westfall local tests.

    type(set_collection), allocatable :: maximal(:)
    type(glht_type) :: subset
    real(dp), allocatable :: local_p(:)
    real(dp), allocatable :: ordered_stat(:)
    real(dp), allocatable :: ordered_rhs(:)
    real(dp), allocatable :: ordered_k(:, :)
    real(dp), allocatable :: ordered_p(:)
    integer, allocatable :: order(:)
    integer :: i
    integer :: j
    integer :: m
    integer :: local_index
    real(dp) :: local_error
    real(dp) :: set_value
    real(dp) :: value
    logical :: max_ok
    character(len=256) :: max_message

    m = size(object%estimate)
    allocate(ordered_stat(m))
    select case (object%alternative)
    case (alt_less)
      ordered_stat = object%statistic
    case (alt_greater)
      ordered_stat = -object%statistic
    case default
      ordered_stat = -abs(object%statistic)
    end select
    do i = 1, m
      ordered_stat(i) = round_significant(ordered_stat(i), 10)
    end do
    call stable_order(ordered_stat, order)
    ordered_k = object%linfct(order, :)
    ordered_rhs = object%rhs(order)
    call compute_maxsets(ordered_k, maximal, max_ok, max_message)
    if (.not. max_ok) then
      ok = .false.
      message = trim(max_message)
      allocate(pvalue(m))
      pvalue = 1.0_dp
      error = 0.0_dp
      return
    end if

    allocate(ordered_p(m))
    ordered_p = 0.0_dp
    error = 0.0_dp
    do i = 1, m
      value = 0.0_dp
      do j = 1, size(maximal(i)%set)
        call subset_ordered_glht(object, order, maximal(i)%set(j)%value, subset)
        if (lowercase(trim(method)) == 'westfall') then
          call single_step_pvalues(subset, local_p, local_error, control)
          error = max(error, local_error)
        else
          allocate(local_p(size(subset%statistic)))
          local_p = [(marginal_pvalue(subset, local_index), &
            local_index = 1, size(subset%statistic))]
          local_p = min(1.0_dp, real(size(local_p), dp) * local_p)
        end if
        set_value = minval(local_p)
        value = max(value, set_value)
        deallocate(local_p)
      end do
      ordered_p(i) = value
      if (i > 1) ordered_p(i) = max(ordered_p(i - 1), ordered_p(i))
    end do

    allocate(pvalue(m))
    do i = 1, m
      pvalue(order(i)) = ordered_p(i)
    end do
    pvalue = min(1.0_dp, max(0.0_dp, pvalue))
    ok = .true.
    message = ''
  end subroutine closed_pvalues

  subroutine subset_glht(object, indices, subset)
    type(glht_type), intent(in) :: object !! Parent GLHT object to subset by hypothesis index.
    integer, intent(in) :: indices(:) !! One-based original hypothesis indices to retain.
    type(glht_type), intent(out) :: subset !! GLHT object restricted to the requested hypotheses.

    allocate(subset%linfct(size(indices), size(object%linfct, 2)))
    allocate(subset%rhs(size(indices)), subset%estimate(size(indices)))
    allocate(subset%covariance(size(indices), size(indices)))
    allocate(subset%standard_error(size(indices)), subset%statistic(size(indices)))
    allocate(subset%correlation(size(indices), size(indices)))
    subset%linfct = object%linfct(indices, :)
    subset%rhs = object%rhs(indices)
    subset%estimate = object%estimate(indices)
    subset%covariance = object%covariance(indices, indices)
    subset%standard_error = object%standard_error(indices)
    subset%statistic = object%statistic(indices)
    subset%correlation = object%correlation(indices, indices)
    subset%df = object%df
    subset%alternative = object%alternative
    subset%ok = .true.
  end subroutine subset_glht

  subroutine subset_ordered_glht(object, order, local_indices, subset)
    type(glht_type), intent(in) :: object !! Parent GLHT object in original hypothesis order.
    integer, intent(in) :: order(:) !! Permutation defining the closed-testing ordered hypotheses.
    integer, intent(in) :: local_indices(:) !! One-based positions within order to retain.
    type(glht_type), intent(out) :: subset !! GLHT object for the corresponding original hypotheses.

    integer, allocatable :: original_indices(:)

    original_indices = order(local_indices)
    call subset_glht(object, original_indices, subset)
  end subroutine subset_ordered_glht

  subroutine integration_correlation(input, output)
    real(dp), intent(in) :: input(:, :) !! Correlation matrix, which may be positive semidefinite for dependent contrasts.
    real(dp), allocatable, intent(out) :: output(:, :) !! Positive-definite correlation used by the numerical mvtnorm integrator.

    real(dp), allocatable :: inverse(:, :)
    real(dp), parameter :: nugget = 1.0e-10_dp
    real(dp) :: shrink
    integer :: i
    integer :: rank
    logical :: inverse_ok

    output = input
    call generalized_inverse(input, inverse, rank, inverse_ok)
    if (.not. inverse_ok) return
    if (rank >= size(input, 1)) return

    shrink = 1.0_dp / (1.0_dp + nugget)
    output = shrink * input
    do i = 1, size(input, 1)
      output(i, i) = 1.0_dp
    end do
  end subroutine integration_correlation

  real(dp) function round_significant(x, digits) result(value)
    real(dp), intent(in) :: x !! Number to round to a requested count of significant decimal digits.
    integer, intent(in) :: digits !! Positive number of significant decimal digits to retain.

    real(dp) :: exponent
    real(dp) :: scale

    if (abs(x) <= tiny(1.0_dp) .or. digits <= 0) then
      value = x
      return
    end if
    exponent = floor(log10(abs(x)))
    scale = 10.0_dp ** real(digits - 1, dp) / (10.0_dp ** exponent)
    value = anint(x * scale) / scale
  end function round_significant

  pure function lowercase(text) result(out)
    character(len=*), intent(in) :: text !! Method name to normalize for case-insensitive matching.
    character(len=len(text)) :: out

    integer :: c
    integer :: i

    out = text
    do i = 1, len(text)
      c = iachar(out(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
    end do
  end function lowercase

end module multcomp_glht
