module trading_csa
  use trading_kinds, only : dp, str_len
  implicit none
  private

  type, public :: csa_t
    character(len=str_len) :: id = ""
    real(dp) :: threshold_counterparty = 0.0_dp
    real(dp) :: threshold_processing_organization = 0.0_dp
    real(dp) :: minimum_transfer_counterparty = 0.0_dp
    real(dp) :: minimum_transfer_processing_organization = 0.0_dp
    real(dp) :: initial_margin_counterparty = 0.0_dp
    real(dp) :: initial_margin_processing_organization = 0.0_dp
    real(dp) :: mpor_days = 0.0_dp
    real(dp) :: remargin_frequency = 0.0_dp
    real(dp) :: rounding = 0.0_dp
    character(len=str_len) :: counterparty = ""
    character(len=str_len), allocatable :: currencies(:)
    character(len=str_len), allocatable :: trade_groups(:)
    character(len=str_len) :: values_type = ""
  contains
    procedure :: apply_threshold
    procedure :: maturity_factor
  end type csa_t

  type, public :: collateral_t
    character(len=str_len) :: id = ""
    real(dp) :: amount = 0.0_dp
    character(len=str_len) :: csa_id = ""
    character(len=str_len) :: collateral_type = ""
  end type collateral_t

contains

  subroutine apply_threshold(self, mtm_vector, collateralized_mtm)
    class(csa_t), intent(in) :: self
    real(dp), intent(in) :: mtm_vector(:)
    real(dp), intent(out) :: collateralized_mtm(:)
    real(dp) :: collateral
    real(dp) :: mtm_difference
    integer :: odd_index
    integer :: n

    n = size(mtm_vector)
    if (size(collateralized_mtm) /= n) then
      error stop "apply_threshold: result has the wrong size"
    end if
    if (n == 0) return

    collateralized_mtm = mtm_vector
    collateral = 0.0_dp
    mtm_difference = 0.0_dp
    collateralized_mtm(1) = mtm_vector(1)

    odd_index = 1
    do while (odd_index <= n)
      if (collateralized_mtm(odd_index) > &
          self%threshold_counterparty + self%minimum_transfer_counterparty) then
        collateral = collateral + collateralized_mtm(odd_index) - &
          self%threshold_counterparty - self%minimum_transfer_counterparty
      else if (collateralized_mtm(odd_index) < &
          -self%threshold_processing_organization - &
          self%minimum_transfer_processing_organization) then
        collateral = collateral + collateralized_mtm(odd_index) + &
          self%threshold_processing_organization + &
          self%minimum_transfer_processing_organization
      else if (odd_index > 3 .and. collateral > 0.0_dp .and. &
          mtm_difference < 0.0_dp) then
        collateral = max(collateral + mtm_difference, 0.0_dp)
      else if (odd_index > 3 .and. collateral < 0.0_dp .and. &
          mtm_difference > 0.0_dp) then
        collateral = min(collateral + mtm_difference, 0.0_dp)
      end if

      if (odd_index + 1 <= n) then
        collateralized_mtm(odd_index + 1) = mtm_vector(odd_index + 1) - collateral
      end if
      if (odd_index + 2 <= n) then
        collateralized_mtm(odd_index + 2) = mtm_vector(odd_index + 2) - collateral
        mtm_difference = mtm_vector(odd_index + 2) - mtm_vector(odd_index)
      end if
      odd_index = odd_index + 2
    end do
  end subroutine apply_threshold

  pure real(dp) function maturity_factor(self, simplified) result(value)
    class(csa_t), intent(in) :: self
    logical, intent(in) :: simplified
    real(dp) :: margin_period

    if (simplified) then
      value = 0.42_dp
    else
      margin_period = 10.0_dp + self%remargin_frequency - 1.0_dp
      value = 1.5_dp * sqrt(max(margin_period, 0.0_dp) / 250.0_dp)
    end if
  end function maturity_factor

end module trading_csa
