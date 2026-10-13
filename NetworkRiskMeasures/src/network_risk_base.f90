! SPDX-License-Identifier: GPL-3.0-only
! Derived from NetworkRiskMeasures 0.1.7, Carlos Cinelli and Thiago Cristiano Silva.
module network_risk_base
   use r_kinds, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   integer, parameter, public :: nr_ok = 0, nr_shape = 1, nr_domain = 2
   integer, parameter, public :: nr_iterations = 3, nr_infeasible = 4, nr_random_exhausted = 5, nr_numeric = 6
   public :: dp, risk_matrix, impact_matrix, vulnerability_matrix, criticality, nan_value
contains
   pure function nan_value() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_value

   pure subroutine risk_matrix(exposures, result, status, buffer, binary, exposure_type, returns)
      real(dp), intent(in) :: exposures(:, :) !! Square nonnegative finite matrix, assets by default.
      real(dp), intent(out) :: result(:, :) !! Output square impact matrix, or vulnerability when requested.
      integer, intent(out) :: status !! Zero on success, nr_shape for dimensions, nr_domain for invalid data.
      real(dp), intent(in), optional :: buffer(:) !! Positive capital buffers of length n, required for assets/liabilities.
      logical, intent(in), optional :: binary !! If true, threshold relative exposures at one; default false.
      character(*), intent(in), optional :: exposure_type !! One of assets, liabilities, impact, vulnerability; default assets.
      character(*), intent(in), optional :: returns !! Output orientation: impact (default) or vulnerability.
      character(:), allocatable :: kind, target
      logical :: bin
      integer :: n, i
      result = nan_value()
      status = nr_shape
      n = size(exposures, 1)
      if (n < 1 .or. size(exposures, 2) /= n) return
      if (any(shape(result) /= [n, n])) return
      status = nr_domain
      if (any(.not. ieee_is_finite(exposures))) return
      if (any(exposures < 0.0_dp)) return
      kind = 'assets'
      if (present(exposure_type)) kind = exposure_type
      target = 'impact'
      if (present(returns)) target = returns
      if (target /= 'impact' .and. target /= 'vulnerability') return
      bin = .false.
      if (present(binary)) bin = binary
      if (present(buffer)) then
         status = nr_shape
         if (size(buffer) /= n) return
         status = nr_domain
         if (any(.not. ieee_is_finite(buffer))) return
         if (any(buffer < sqrt(epsilon(1.0_dp)))) return
      end if
      select case (kind)
      case ('impact')
         result = exposures
         if (target == 'vulnerability') result = transpose(exposures)
      case ('vulnerability')
         result = exposures
         if (target == 'impact') result = transpose(exposures)
      case ('assets', 'liabilities')
         if (.not. present(buffer)) return
         result = exposures
         if (kind == 'liabilities') result = transpose(exposures)
         do i = 1, n
            result(i, :) = result(i, :) / buffer(i)
         end do
         if (target == 'impact') result = transpose(result)
      case default
         return
      end select
      if (any(.not. ieee_is_finite(result))) then
         status = nr_numeric
         result = nan_value()
         return
      end if
      if (bin) result = merge(1.0_dp, 0.0_dp, result >= 1.0_dp)
      status = nr_ok
   end subroutine risk_matrix

   pure subroutine impact_matrix(exposures, result, status, buffer, binary, exposure_type)
      real(dp), intent(in) :: exposures(:, :) !! Square finite nonnegative bilateral exposures.
      real(dp), intent(out) :: result(:, :) !! Impact matrix: entry i,j is impact of default i on vertex j.
      integer, intent(out) :: status !! Zero on success; otherwise a public nr_* status code.
      real(dp), intent(in), optional :: buffer(:) !! Capital buffers in vertex order, needed for assets or liabilities.
      logical, intent(in), optional :: binary !! Threshold impact entries at one when true; default false.
      character(*), intent(in), optional :: exposure_type !! Orientation: assets (default), liabilities, impact, vulnerability.
      call risk_matrix(exposures, result, status, buffer, binary, exposure_type, 'impact')
   end subroutine impact_matrix

   pure subroutine vulnerability_matrix(exposures, result, status, buffer, binary, exposure_type)
      real(dp), intent(in) :: exposures(:, :) !! Square finite nonnegative bilateral exposures.
      real(dp), intent(out) :: result(:, :) !! Vulnerability matrix: entry i,j is damage to i from default j.
      integer, intent(out) :: status !! Zero on success; otherwise a public nr_* status code.
      real(dp), intent(in), optional :: buffer(:) !! Capital buffers in vertex order, needed for assets or liabilities.
      logical, intent(in), optional :: binary !! Threshold vulnerability entries at one when true; default false.
      character(*), intent(in), optional :: exposure_type !! Orientation: assets (default), liabilities, impact, vulnerability.
      call risk_matrix(exposures, result, status, buffer, binary, exposure_type, 'vulnerability')
   end subroutine vulnerability_matrix

   pure subroutine criticality(exposures, result, status, buffer, binary, exposure_type)
      real(dp), intent(in) :: exposures(:, :) !! Square finite nonnegative bilateral exposures.
      real(dp), intent(out) :: result(:) !! Length-n row sums of impact entries capped above at one.
      integer, intent(out) :: status !! Zero on success; otherwise a public nr_* status code.
      real(dp), intent(in), optional :: buffer(:) !! Positive capital buffer per vertex for assets or liabilities.
      logical, intent(in), optional :: binary !! Threshold impact entries before summation; default false.
      character(*), intent(in), optional :: exposure_type !! Orientation: assets (default), liabilities, impact, vulnerability.
      real(dp) :: impact(size(exposures, 1), size(exposures, 2))
      result = nan_value()
      status = nr_shape
      if (size(result) /= size(exposures, 1)) return
      call impact_matrix(exposures, impact, status, buffer, binary, exposure_type)
      if (status /= nr_ok) return
      result = sum(min(impact, 1.0_dp), dim=2)
   end subroutine criticality
end module network_risk_base
