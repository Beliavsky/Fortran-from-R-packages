module trading_climate
  use trading_kinds, only : dp
  implicit none
  private

  public :: carbon_footprint
  public :: carbon_intensity
  public :: total_carbon_emissions
  public :: weighted_average_carbon_intensity

contains

  pure real(dp) function carbon_footprint(exposures, emissions, capitalization) result(value)
    real(dp), intent(in) :: exposures(:)
    real(dp), intent(in) :: emissions(:)
    real(dp), intent(in) :: capitalization(:)
    integer :: n

    n = min(size(exposures), min(size(emissions), size(capitalization)))
    if (n == 0 .or. abs(sum(exposures(:n))) <= tiny(1.0_dp)) then
      value = 0.0_dp
      return
    end if
    value = sum(exposures(:n) / capitalization(:n) * emissions(:n)) / &
      sum(exposures(:n))
  end function carbon_footprint

  pure real(dp) function carbon_intensity(exposures, emissions, capitalization, revenue) &
      result(value)
    real(dp), intent(in) :: exposures(:)
    real(dp), intent(in) :: emissions(:)
    real(dp), intent(in) :: capitalization(:)
    real(dp), intent(in) :: revenue(:)
    real(dp), allocatable :: ownership_weights(:)
    real(dp) :: denominator
    integer :: n

    n = min(size(exposures), min(size(emissions), &
      min(size(capitalization), size(revenue))))
    if (n == 0) then
      value = 0.0_dp
      return
    end if

    allocate(ownership_weights(n))
    ownership_weights = exposures(:n) / capitalization(:n)
    denominator = sum(ownership_weights * revenue(:n))
    if (abs(denominator) <= tiny(1.0_dp)) then
      value = 0.0_dp
    else
      value = sum(ownership_weights * emissions(:n)) / denominator
    end if
  end function carbon_intensity

  pure real(dp) function total_carbon_emissions(exposures, emissions, capitalization) &
      result(value)
    real(dp), intent(in) :: exposures(:)
    real(dp), intent(in) :: emissions(:)
    real(dp), intent(in) :: capitalization(:)
    integer :: n

    n = min(size(exposures), min(size(emissions), size(capitalization)))
    if (n == 0) then
      value = 0.0_dp
    else
      value = sum(exposures(:n) / capitalization(:n) * emissions(:n))
    end if
  end function total_carbon_emissions

  pure real(dp) function weighted_average_carbon_intensity(exposures, emissions, revenue) &
      result(value)
    real(dp), intent(in) :: exposures(:)
    real(dp), intent(in) :: emissions(:)
    real(dp), intent(in) :: revenue(:)
    real(dp) :: total_exposure
    integer :: n

    n = min(size(exposures), min(size(emissions), size(revenue)))
    if (n == 0) then
      value = 0.0_dp
      return
    end if

    total_exposure = sum(exposures(:n))
    if (abs(total_exposure) <= tiny(1.0_dp)) then
      value = 0.0_dp
    else
      value = sum((exposures(:n) / total_exposure) * &
        emissions(:n) / revenue(:n))
    end if
  end function weighted_average_carbon_intensity

end module trading_climate
