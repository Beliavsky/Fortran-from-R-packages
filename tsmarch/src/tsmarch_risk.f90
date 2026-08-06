! SPDX-License-Identifier: GPL-2.0-only
module tsmarch_risk
  use ghyp_kinds, only : dp
  use tsd_math, only : normal_quantile, normal_pdf
  use tsmarch_types, only : risk_result, tsm_success, tsm_invalid_argument
  use tsmarch_linalg, only : type7_quantile
  implicit none
  private
  public :: simulation_risk, gaussian_risk, portfolio_paths
  public :: value_at_risk, expected_shortfall

contains

  function portfolio_paths(simulated, weights) result(paths)
    real(dp), intent(in) :: simulated(:, :, :)
    real(dp), intent(in) :: weights(:)
    real(dp), allocatable :: paths(:, :)
    integer :: h, p
    allocate(paths(size(simulated, 1), size(simulated, 3)))
    do p = 1, size(simulated, 3)
      do h = 1, size(simulated, 1)
        paths(h, p) = dot_product(simulated(h, :, p), weights)
      end do
    end do
  end function portfolio_paths

  function simulation_risk(simulated, probabilities, weights) result(out)
    real(dp), intent(in) :: simulated(:, :, :)
    real(dp), intent(in) :: probabilities(:)
    real(dp), intent(in), optional :: weights(:)
    type(risk_result) :: out
    real(dp), allocatable :: paths(:, :), tail(:)
    real(dp) :: q
    integer :: h, j, p, count_tail, nassets
    nassets = size(simulated, 2)
    if (size(simulated, 1) < 1 .or. size(simulated, 3) < 2 .or. any(probabilities <= 0.0_dp) .or. &
        any(probabilities >= 1.0_dp)) then
      out%status = tsm_invalid_argument
      out%message = 'simulation risk requires paths and probabilities strictly between zero and one'
      return
    end if
    if (present(weights)) then
      if (size(weights) /= nassets) then
        out%status = tsm_invalid_argument
        out%message = 'portfolio weights do not conform to the simulated asset dimension'
        return
      end if
      paths = portfolio_paths(simulated, weights)
    else
      allocate(paths(size(simulated, 1), size(simulated, 3)))
      do p = 1, size(simulated, 3)
        paths(:, p) = simulated(:, 1, p)
      end do
    end if
    allocate(out%probabilities(size(probabilities)))
    out%probabilities = probabilities
    allocate(out%value_at_risk(size(paths, 1), size(probabilities)))
    allocate(out%expected_shortfall(size(paths, 1), size(probabilities)))
    do h = 1, size(paths, 1)
      do j = 1, size(probabilities)
        q = type7_quantile(paths(h, :), probabilities(j))
        out%value_at_risk(h, j) = q
        count_tail = count(paths(h, :) <= q)
        if (count_tail > 0) then
          allocate(tail(count_tail))
          tail = pack(paths(h, :), paths(h, :) <= q)
          out%expected_shortfall(h, j) = sum(tail) / real(count_tail, dp)
          deallocate(tail)
        else
          out%expected_shortfall(h, j) = q
        end if
      end do
    end do
    out%status = tsm_success
    out%message = 'ok'
  end function simulation_risk

  function gaussian_risk(mean, sigma, probabilities) result(out)
    real(dp), intent(in) :: mean(:), sigma(:), probabilities(:)
    type(risk_result) :: out
    real(dp) :: z
    integer :: h, j
    if (size(mean) /= size(sigma) .or. any(sigma < 0.0_dp) .or. any(probabilities <= 0.0_dp) .or. &
        any(probabilities >= 1.0_dp)) then
      out%status = tsm_invalid_argument
      out%message = 'invalid Gaussian risk arguments'
      return
    end if
    allocate(out%probabilities(size(probabilities)))
    out%probabilities = probabilities
    allocate(out%value_at_risk(size(mean), size(probabilities)))
    allocate(out%expected_shortfall(size(mean), size(probabilities)))
    do h = 1, size(mean)
      do j = 1, size(probabilities)
        z = normal_quantile(probabilities(j))
        out%value_at_risk(h, j) = mean(h) + sigma(h) * z
        out%expected_shortfall(h, j) = mean(h) - sigma(h) * normal_pdf(z) / probabilities(j)
      end do
    end do
    out%status = tsm_success
    out%message = 'ok'
  end function gaussian_risk

  function value_at_risk(simulated, probabilities, weights) result(values)
    real(dp), intent(in) :: simulated(:, :, :), probabilities(:)
    real(dp), intent(in), optional :: weights(:)
    real(dp), allocatable :: values(:, :)
    type(risk_result) :: out
    if (present(weights)) then
      out = simulation_risk(simulated, probabilities, weights)
    else
      out = simulation_risk(simulated, probabilities)
    end if
    if (out%status == tsm_success) then
      values = out%value_at_risk
    else
      allocate(values(0, 0))
    end if
  end function value_at_risk

  function expected_shortfall(simulated, probabilities, weights) result(values)
    real(dp), intent(in) :: simulated(:, :, :), probabilities(:)
    real(dp), intent(in), optional :: weights(:)
    real(dp), allocatable :: values(:, :)
    type(risk_result) :: out
    if (present(weights)) then
      out = simulation_risk(simulated, probabilities, weights)
    else
      out = simulation_risk(simulated, probabilities)
    end if
    if (out%status == tsm_success) then
      values = out%expected_shortfall
    else
      allocate(values(0, 0))
    end if
  end function expected_shortfall

end module tsmarch_risk
