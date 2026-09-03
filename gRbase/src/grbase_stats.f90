module grbase_stats
  use r_kinds, only : dp
  use r_linalg, only : solve_spd
  implicit none
  private

  public :: inverse_spd
  public :: concentration_to_partial_correlation
  public :: covariance_to_partial_correlation

contains

  pure subroutine inverse_spd(a, ainv, info)
    real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix to invert.
    real(dp), allocatable, intent(out) :: ainv(:, :) !! Inverse of `a`; allocated to zero size when dimensions are invalid.
    integer, intent(out) :: info !! Zero on success; nonzero on invalid dimensions or failed SPD factorization.
    real(dp), allocatable :: identity(:, :)
    integer :: i
    integer :: n

    if (size(a, 1) /= size(a, 2)) then
      allocate(ainv(0, 0))
      info = -1
      return
    end if
    n = size(a, 1)
    allocate(identity(n, n), ainv(n, n))
    identity = 0.0_dp
    do i = 1, n
      identity(i, i) = 1.0_dp
    end do
    call solve_spd(a, identity, ainv, info)
  end subroutine inverse_spd

  pure subroutine concentration_to_partial_correlation(k, pcor, info)
    real(dp), intent(in) :: k(:, :) !! Symmetric concentration/precision matrix with positive diagonal entries.
    real(dp), allocatable, intent(out) :: pcor(:, :) !! Partial-correlation matrix with unit diagonal.
    integer, intent(out) :: info !! Zero on success; nonzero for nonsquare input or nonpositive diagonal entries.
    integer :: i
    integer :: j
    integer :: n

    if (size(k, 1) /= size(k, 2)) then
      allocate(pcor(0, 0))
      info = -1
      return
    end if
    n = size(k, 1)
    allocate(pcor(n, n))
    if (any([(k(i, i) <= 0.0_dp, i = 1, n)])) then
      pcor = 0.0_dp
      info = 1
      return
    end if
    do j = 1, n
      do i = 1, n
        if (i == j) then
          pcor(i, j) = 1.0_dp
        else
          pcor(i, j) = -k(i, j) / sqrt(k(i, i) * k(j, j))
        end if
      end do
    end do
    info = 0
  end subroutine concentration_to_partial_correlation

  pure subroutine covariance_to_partial_correlation(covariance, pcor, info)
    real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-definite covariance matrix.
    real(dp), allocatable, intent(out) :: pcor(:, :) !! Partial-correlation matrix implied by the inverse covariance.
    integer, intent(out) :: info !! Zero on success; otherwise the SPD inversion or diagonal-validation status.
    real(dp), allocatable :: precision(:, :)

    call inverse_spd(covariance, precision, info)
    if (info /= 0) then
      allocate(pcor(0, 0))
      return
    end if
    call concentration_to_partial_correlation(precision, pcor, info)
  end subroutine covariance_to_partial_correlation

end module grbase_stats
