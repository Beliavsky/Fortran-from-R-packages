! SPDX-License-Identifier: MIT
! Copyright (c) 2023 Bernardo Reckziegel
module epo_core
  use epo_kinds, only : dp
  use epo_linalg, only : quadratic_form, solve_spd
  use epo_statistics, only : all_finite_matrix, all_finite_vector, &
    covariance_to_correlation, sample_covariance
  use epo_types, only : epo_invalid_input, epo_normalization_failure, &
    epo_result, epo_singular_matrix, set_epo_error
  implicit none
  private

  public :: epo_optimize
  public :: epo_from_covariance
  public :: simple_epo
  public :: anchored_epo

contains

  function epo_optimize(returns, signal, lambda, method, w, anchor, &
      normalize, endogenous) result(result)
    real(dp), intent(in) :: returns(:,:)
    real(dp), intent(in) :: signal(:)
    real(dp), intent(in) :: lambda
    character(len=*), intent(in) :: method
    real(dp), intent(in) :: w
    real(dp), intent(in), optional :: anchor(:)
    logical, intent(in), optional :: normalize
    logical, intent(in), optional :: endogenous
    type(epo_result) :: result

    real(dp), allocatable :: means(:), covariance(:,:)
    logical :: ok
    integer :: n_assets

    n_assets = size(returns,2)
    if (size(returns,1) < 2 .or. n_assets < 1) then
      call set_epo_error(result, epo_invalid_input, &
        'returns must contain at least two observations and one asset')
      return
    end if

    if (size(signal) /= n_assets) then
      call set_epo_error(result, epo_invalid_input, &
        'signal length must equal the number of assets')
      return
    end if

    allocate(means(n_assets), covariance(n_assets,n_assets))
    call sample_covariance(returns, means, covariance, ok)
    if (.not. ok) then
      call set_epo_error(result, epo_invalid_input, &
        'unable to calculate a finite sample covariance matrix')
      return
    end if

    result = epo_from_covariance(covariance, signal, lambda, method, w, &
      anchor, normalize, endogenous)
  end function epo_optimize

  function epo_from_covariance(covariance, signal, lambda, method, w, &
      anchor, normalize, endogenous) result(result)
    real(dp), intent(in) :: covariance(:,:)
    real(dp), intent(in) :: signal(:)
    real(dp), intent(in) :: lambda
    character(len=*), intent(in) :: method
    real(dp), intent(in) :: w
    real(dp), intent(in), optional :: anchor(:)
    logical, intent(in), optional :: normalize
    logical, intent(in), optional :: endogenous
    type(epo_result) :: result

    character(len=:), allocatable :: selected_method
    real(dp), allocatable :: diagonal_covariance(:,:), rhs(:), solved_signal(:)
    real(dp) :: denominator, numerator, scale, sum_weights, tol
    logical :: do_normalize, endogenous_risk, ok
    integer :: i, n_assets

    n_assets = size(covariance,1)
    selected_method = lowercase(trim(method))
    do_normalize = .true.
    if (present(normalize)) do_normalize = normalize
    endogenous_risk = .true.
    if (present(endogenous)) endogenous_risk = endogenous

    result%method = selected_method
    result%normalized = do_normalize
    result%endogenous = endogenous_risk

    if (n_assets < 1 .or. size(covariance,2) /= n_assets) then
      call set_epo_error(result, epo_invalid_input, &
        'covariance must be a nonempty square matrix')
      return
    end if
    if (size(signal) /= n_assets) then
      call set_epo_error(result, epo_invalid_input, &
        'signal length must equal the covariance dimension')
      return
    end if
    if (.not. all_finite_matrix(covariance) .or. &
        .not. all_finite_vector(signal)) then
      call set_epo_error(result, epo_invalid_input, &
        'covariance and signal must contain only finite values')
      return
    end if
    if (w < 0.0_dp .or. w > 1.0_dp) then
      call set_epo_error(result, epo_invalid_input, &
        'w must lie between zero and one')
      return
    end if
    if (selected_method /= 'simple' .and. selected_method /= 'anchored') then
      call set_epo_error(result, epo_invalid_input, &
        'method must be simple or anchored')
      return
    end if
    if (selected_method == 'simple' .or. .not. endogenous_risk) then
      if (lambda <= 0.0_dp) then
        call set_epo_error(result, epo_invalid_input, &
          'lambda must be positive when it is used')
        return
      end if
    end if
    if (selected_method == 'anchored') then
      if (.not. present(anchor)) then
        call set_epo_error(result, epo_invalid_input, &
          'the anchored method requires an anchor vector')
        return
      end if
      if (size(anchor) /= n_assets) then
        call set_epo_error(result, epo_invalid_input, &
          'anchor length must equal the covariance dimension')
        return
      end if
      if (.not. all_finite_vector(anchor)) then
        call set_epo_error(result, epo_invalid_input, &
          'anchor must contain only finite values')
        return
      end if
    end if

    allocate(result%weights(n_assets))
    allocate(result%covariance(n_assets,n_assets))
    allocate(result%correlation(n_assets,n_assets))
    allocate(result%shrunk_correlation(n_assets,n_assets))
    allocate(result%shrunk_covariance(n_assets,n_assets))
    allocate(diagonal_covariance(n_assets,n_assets))
    allocate(rhs(n_assets), solved_signal(n_assets))

    result%covariance = 0.5_dp * (covariance + transpose(covariance))
    call covariance_to_correlation(result%covariance, result%correlation, ok)
    if (.not. ok) then
      call set_epo_error(result, epo_invalid_input, &
        'all asset variances must be strictly positive')
      return
    end if

    diagonal_covariance = 0.0_dp
    do i = 1, n_assets
      diagonal_covariance(i,i) = result%covariance(i,i)
    end do

    result%shrunk_correlation = (1.0_dp - w) * result%correlation
    do i = 1, n_assets
      result%shrunk_correlation(i,i) = &
        result%shrunk_correlation(i,i) + w
    end do

    result%shrunk_covariance = (1.0_dp - w) * result%covariance + &
      w * diagonal_covariance

    call solve_spd(result%shrunk_covariance, signal, solved_signal, ok)
    if (.not. ok) then
      call set_epo_error(result, epo_singular_matrix, &
        'the shrunk covariance matrix is not positive definite')
      return
    end if

    if (selected_method == 'simple') then
      result%weights = solved_signal / lambda
    else
      if (endogenous_risk) then
        numerator = quadratic_form(anchor, result%shrunk_covariance)
        denominator = dot_product(signal, solved_signal)
        tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(numerator))
        if (numerator < -tol .or. denominator <= tol) then
          call set_epo_error(result, epo_invalid_input, &
            'unable to calculate endogenous risk-aversion scaling')
          return
        end if
        result%gamma = sqrt(max(0.0_dp, numerator) / denominator)
        scale = result%gamma
      else
        scale = 1.0_dp / lambda
      end if

      rhs = (1.0_dp - w) * scale * signal + &
        w * matmul(diagonal_covariance, anchor)
      call solve_spd(result%shrunk_covariance, rhs, result%weights, ok)
      if (.not. ok) then
        call set_epo_error(result, epo_singular_matrix, &
          'the shrunk covariance matrix is not positive definite')
        return
      end if
    end if

    sum_weights = sum(result%weights)
    result%weight_sum_before_normalization = sum_weights
    if (do_normalize) then
      tol = 1000.0_dp * epsilon(1.0_dp) * &
        max(1.0_dp, sum(abs(result%weights)))
      if (abs(sum_weights) <= tol) then
        call set_epo_error(result, epo_normalization_failure, &
          'portfolio weights cannot be normalized because their sum is zero')
        return
      end if
      result%weights = result%weights / sum_weights
    end if

    result%ok = .true.
    result%status = 0
    result%message = 'success'
  end function epo_from_covariance

  function simple_epo(returns, signal, lambda, w, normalize) result(result)
    real(dp), intent(in) :: returns(:,:)
    real(dp), intent(in) :: signal(:)
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: w
    logical, intent(in), optional :: normalize
    type(epo_result) :: result

    result = epo_optimize(returns, signal, lambda, 'simple', w, &
      normalize=normalize)
  end function simple_epo

  function anchored_epo(returns, signal, lambda, w, anchor, normalize, &
      endogenous) result(result)
    real(dp), intent(in) :: returns(:,:)
    real(dp), intent(in) :: signal(:)
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: w
    real(dp), intent(in) :: anchor(:)
    logical, intent(in), optional :: normalize
    logical, intent(in), optional :: endogenous
    type(epo_result) :: result

    result = epo_optimize(returns, signal, lambda, 'anchored', w, &
      anchor=anchor, normalize=normalize, endogenous=endogenous)
  end function anchored_epo

  pure function lowercase(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: code, i

    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end function lowercase

end module epo_core
