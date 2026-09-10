! SPDX-License-Identifier: GPL-2.0-or-later
! Cadzow low-rank approximation translated from Rssa 1.1.
module rssa_cadzow
   use rssa_decomposition, only : ssa_1d
   use rssa_kinds, only : dp
   use rssa_reconstruction, only : reconstruct_ssa
   use rssa_types, only : rssa_invalid_input, rssa_success, ssa_result
   implicit none
   private
   public :: cadzow, cadzow_ssa
contains
   function cadzow(series, window, rank, iterations, tolerance, info) result(filtered)
      real(dp), intent(in) :: series(:) !! Input series to project repeatedly onto the rank-constrained Hankel set.
      integer, intent(in), optional :: window !! SSA window length; default `(N+1)/2`.
      integer, intent(in), optional :: rank !! Retained rank; default two.
      integer, intent(in), optional :: iterations !! Maximum Cadzow iterations; default ten.
      real(dp), intent(in), optional :: tolerance !! Relative stopping tolerance; default `1e-8`.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: filtered(:), previous(:), reconstructed(:)
      type(ssa_result) :: fit
      integer :: i, iter, l, nrank, limit, status
      real(dp) :: tol, scale
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      nrank = 2
      if (present(rank)) nrank = rank
      limit = 10
      if (present(iterations)) limit = max(1, iterations)
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      filtered = series
      status = rssa_success
      if (l < 1 .or. l > size(series) .or. nrank < 1) status = rssa_invalid_input
      if (status == rssa_success) then
         do iter = 1, limit
            previous = filtered
            fit = ssa_1d(filtered, l, nrank)
            if (fit%info /= rssa_success) then
               status = fit%info
               exit
            end if
            reconstructed = reconstruct_ssa(fit, [(i, i=1, min(nrank, size(fit%sigma)))])
            filtered = reconstructed
            scale = max(1.0_dp, maxval(abs(previous)))
            if (maxval(abs(filtered - previous)) <= tol * scale) exit
         end do
      end if
      if (present(info)) info = status
   end function cadzow

   function cadzow_ssa(object, rank, iterations, tolerance, info) result(filtered)
      type(ssa_result), intent(in) :: object !! Existing SSA object supplying the original series and window.
      integer, intent(in), optional :: rank !! Retained rank; defaults to the stored component count.
      integer, intent(in), optional :: iterations !! Maximum Cadzow iterations.
      real(dp), intent(in), optional :: tolerance !! Relative stopping tolerance.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: filtered(:)
      integer :: nrank, status
      if (.not. allocated(object%series)) then
         allocate(filtered(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      nrank = 2
      if (allocated(object%sigma)) nrank = max(1, size(object%sigma))
      if (present(rank)) nrank = rank
      if (present(tolerance)) then
         filtered = cadzow(object%series, object%window, nrank, iterations, tolerance, status)
      else
         filtered = cadzow(object%series, object%window, nrank, iterations, info=status)
      end if
      if (present(info)) info = status
   end function cadzow_ssa
end module rssa_cadzow
