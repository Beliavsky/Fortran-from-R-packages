! SPDX-License-Identifier: GPL-3.0-only
! Derived from NetworkRiskMeasures 0.1.7, Carlos Cinelli and Thiago Cristiano Silva.
module network_communicability
   use network_risk_base, only: dp, nan_value, nr_ok, nr_shape, nr_domain, nr_numeric, impact_matrix
   use expm_module, only: expm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: communicability_matrix, impact_susceptibility, impact_fluidity, impact_diffusion
contains
   subroutine communicability_matrix(x, result, status, terms)
      real(dp), intent(in) :: x(:, :) !! Square finite matrix; signed values are accepted for a general matrix exponential.
      real(dp), intent(out) :: result(:, :) !! Communicability exp(x), or Taylor polynomial through the requested degree.
      integer, intent(out) :: status !! Zero on success; nr_shape, nr_domain or nr_numeric on invalid or overflowed data.
      integer, intent(in), optional :: terms !! Taylor degree >=0, or -1 (default) for the shared matrix exponential.
      real(dp) :: term(size(x, 1), size(x, 2))
      real(dp), allocatable :: exponential(:, :)
      integer :: n, degree, i
      result = nan_value()
      status = nr_shape
      n = size(x, 1)
      if (n < 1 .or. size(x, 2) /= n) return
      if (any(shape(result) /= [n, n])) return
      status = nr_domain
      if (any(.not. ieee_is_finite(x))) return
      degree = -1
      if (present(terms)) degree = terms
      if (degree < -1) return
      if (degree == -1) then
         exponential = expm(x)
         result = exponential
      else
         result = 0.0_dp
         do i = 1, n
            result(i, i) = 1.0_dp
         end do
         term = result
         do i = 1, degree
            term = matmul(term, x / real(i, dp))
            result = result + term
            if (any(.not. ieee_is_finite(result))) exit
         end do
      end if
      status = nr_ok
      if (any(.not. ieee_is_finite(result))) status = nr_numeric
   end subroutine communicability_matrix

   subroutine impact_susceptibility(exposures, result, status, buffer, weights, terms, binary, exposure_type)
      real(dp), intent(in) :: exposures(:, :) !! Square finite nonnegative exposures with directed orientation specified below.
      real(dp), intent(out) :: result(:) !! Length-n susceptibility, zero for isolated vertices in the unweighted calculation.
      integer, intent(out) :: status !! Zero on success; otherwise a public nr_* status code.
      real(dp), intent(in), optional :: buffer(:) !! Positive capital buffer per vertex, required for assets or liabilities.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative finite weights multiplying column communicability sums.
      integer, intent(in), optional :: terms !! Taylor degree >=0, or -1 (default) for shared full matrix exponential.
      logical, intent(in), optional :: binary !! Threshold impact entries at one, default true.
      character(*), intent(in), optional :: exposure_type !! Orientation: assets (default), liabilities, impact or vulnerability.
      real(dp) :: v(size(exposures, 1), size(exposures, 2)), comm(size(exposures, 1), size(exposures, 2))
      real(dp) :: incoming(size(exposures, 1))
      integer :: i, n
      logical :: bin
      result = nan_value()
      n = size(exposures, 1)
      status = nr_shape
      if (size(result) /= n) return
      if (present(weights)) then
         if (size(weights) /= n) return
         status = nr_domain
         if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
      end if
      bin = .true.
      if (present(binary)) bin = binary
      call impact_matrix(exposures, v, status, buffer, bin, exposure_type)
      if (status /= nr_ok) return
      v = min(v, 1.0_dp)
      call communicability_matrix(v, comm, status, terms)
      if (status /= nr_ok) return
      do i = 1, n
         comm(i, i) = 0.0_dp
      end do
      result = sum(comm, dim=1)
      if (present(weights)) then
         result = result * weights
      else
         incoming = sum(v, dim=1)
         do i = 1, n
            if (incoming(i) == 0.0_dp) then
               result(i) = 0.0_dp
            else
               result(i) = result(i) / incoming(i)
            end if
         end do
      end if
      if (any(.not. ieee_is_finite(result))) status = nr_numeric
   end subroutine impact_susceptibility

   subroutine impact_fluidity(exposures, result, status, buffer, weights, terms, binary, exposure_type)
      real(dp), intent(in) :: exposures(:, :) !! Square nonnegative finite bilateral exposure matrix.
      real(dp), intent(out) :: result !! Arithmetic mean of all vertex susceptibilities, including isolated zeros.
      integer, intent(out) :: status !! Zero on success; otherwise a public nr_* status code.
      real(dp), intent(in), optional :: buffer(:) !! Positive vertex capital buffers for assets/liabilities inputs.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative finite node weights, used without normalization.
      integer, intent(in), optional :: terms !! Taylor degree >=0, or -1 (default) for full matrix exponential.
      logical, intent(in), optional :: binary !! Threshold impact entries at one before capping, default true.
      character(*), intent(in), optional :: exposure_type !! Input orientation, one of assets, liabilities, impact, vulnerability.
      real(dp) :: susceptibility(size(exposures, 1))
      result = nan_value()
      call impact_susceptibility(exposures, susceptibility, status, buffer, weights, terms, binary, exposure_type)
      if (status /= nr_ok) return
      result = sum(susceptibility / real(size(susceptibility), dp))
      if (.not. ieee_is_finite(result)) status = nr_numeric
   end subroutine impact_fluidity

   subroutine impact_diffusion(exposures, result, status, buffer, weights, terms, binary, exposure_type)
      real(dp), intent(in) :: exposures(:, :) !! Square finite nonnegative bilateral exposures.
      real(dp), intent(out) :: result(:, :) !! n-by-3 values: start, intermediate, and total diffusion in that column order.
      integer, intent(out) :: status !! Zero on success; otherwise a public nr_* status code.
      real(dp), intent(in), optional :: buffer(:) !! Positive capital buffers for assets or liabilities exposures.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative finite destination weights for changed communicability paths.
      integer, intent(in), optional :: terms !! Taylor degree >=0, or -1 (default) for shared full matrix exponential.
      logical, intent(in), optional :: binary !! Threshold impact entries at one before capping; default true.
      character(*), intent(in), optional :: exposure_type !! Orientation: assets (default), liabilities, impact or vulnerability.
      real(dp) :: v(size(exposures, 1), size(exposures, 2)), reduced(size(exposures, 1), size(exposures, 2))
      real(dp) :: comm(size(exposures, 1), size(exposures, 2)), alternate(size(exposures, 1), size(exposures, 2))
      real(dp) :: difference(size(exposures, 1), size(exposures, 2)), row_values(size(exposures, 1)), degree
      integer :: n, q, i
      logical :: bin
      result = nan_value()
      n = size(exposures, 1)
      status = nr_shape
      if (any(shape(result) /= [n, 3])) return
      if (present(weights)) then
         if (size(weights) /= n) return
         status = nr_domain
         if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
      end if
      bin = .true.
      if (present(binary)) bin = binary
      call impact_matrix(exposures, v, status, buffer, bin, exposure_type)
      if (status /= nr_ok) return
      v = min(v, 1.0_dp)
      call communicability_matrix(v, comm, status, terms)
      if (status /= nr_ok) return
      do i = 1, n
         comm(i, i) = 0.0_dp
      end do
      result = 0.0_dp
      do q = 1, n
         degree = sum(v(q, :))
         if (degree == 0.0_dp) cycle
         reduced = v
         reduced(q, :) = 0.0_dp
         call communicability_matrix(reduced, alternate, status, terms)
         if (status /= nr_ok) then
            result = nan_value()
            return
         end if
         do i = 1, n
            alternate(i, i) = 0.0_dp
         end do
         difference = max(0.0_dp, comm - alternate)
         if (present(weights)) then
            row_values = matmul(difference, weights)
         else
            row_values = sum(difference, dim=2) / degree
         end if
         result(q, 1) = row_values(q)
         result(q, 2) = sum(row_values) - row_values(q)
         result(q, 3) = sum(row_values)
      end do
      if (any(.not. ieee_is_finite(result))) status = nr_numeric
   end subroutine impact_diffusion
end module network_communicability
