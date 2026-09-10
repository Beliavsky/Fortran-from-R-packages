! SPDX-License-Identifier: GPL-2.0-or-later
! Missing-value classification and SSA gap-filling kernels translated from Rssa 1.1.
module rssa_gapfill
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use rssa_decomposition, only : ssa_1d, ssa_2d, ssa_complex, ssa_mssa
   use rssa_forecast, only : apply_lrr, lrr_default
   use rssa_kinds, only : dp
   use rssa_reconstruction, only : reconstruct_2d, reconstruct_complex, reconstruct_mssa, reconstruct_ssa
   use rssa_types, only : cssa_result, gap_summary, mssa_result, rssa_invalid_input, rssa_success
   use rssa_types, only : ssa2d_result, ssa_result
   implicit none
   private
   public :: summarize_gaps, summarize_gaps_complex
   public :: gapfill_ssa, gapfill_mssa_channel, gapfill_complex
   public :: igapfill, igapfill_ssa, igapfill_2d, igapfill_mssa, igapfill_complex
contains
   function summarize_gaps(series, window) result(summary)
      real(dp), intent(in) :: series(:) !! Input series whose NaNs mark missing values.
      integer, intent(in), optional :: window !! SSA window length; default `(N+1)/2`.
      type(gap_summary) :: summary
      integer, allocatable :: missing(:)
      integer :: i, l, p, first, last
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      allocate(missing(count(ieee_is_nan(series))))
      p = 0
      do i = 1, size(series)
         if (ieee_is_nan(series(i))) then
            p = p + 1
            missing(p) = i
         end if
      end do
      summary%n_missing = size(missing)
      summary%missing_indices = missing
      if (size(missing) == 0) return
      first = missing(1)
      last = missing(size(missing))
      summary%n_left = count(missing < l)
      summary%n_right = count(missing > size(series) - l + 1)
      summary%n_internal = summary%n_missing - summary%n_left - summary%n_right
      if (summary%n_missing == last - first + 1) then
         summary%n_dense = summary%n_missing
      else
         summary%n_sparse = summary%n_missing
      end if
   end function summarize_gaps

   function summarize_gaps_complex(series, window) result(summary)
      complex(dp), intent(in) :: series(:) !! Complex series; NaN in either part marks a missing value.
      integer, intent(in), optional :: window !! SSA window length.
      type(gap_summary) :: summary
      real(dp), allocatable :: proxy(:)
      integer :: i
      allocate(proxy(size(series)))
      do i = 1, size(series)
         if (ieee_is_nan(real(series(i), dp)) .or. ieee_is_nan(aimag(series(i)))) then
            proxy(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            proxy(i) = real(series(i), dp)
         end if
      end do
      if (present(window)) then
         summary = summarize_gaps(proxy, window)
      else
         summary = summarize_gaps(proxy)
      end if
   end function summarize_gaps_complex

   function gapfill_ssa(object, incomplete, indices, alpha, info) result(filled)
      type(ssa_result), intent(in) :: object !! SSA template supplying a signal subspace/window.
      real(dp), intent(in) :: incomplete(:) !! Series with NaNs at missing positions.
      integer, intent(in), optional :: indices(:) !! Components used to obtain recurrence coefficients.
      real(dp), intent(in), optional :: alpha !! Blend of forward/backward recurrence for internal gaps; default 0.5.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: filled(:), flrr(:), blrr(:), forward(:), backward(:)
      logical, allocatable :: missing(:)
      integer :: i, left, right, status
      real(dp) :: blend
      filled = incomplete
      missing = ieee_is_nan(incomplete)
      blend = 0.5_dp
      if (present(alpha)) blend = min(1.0_dp, max(0.0_dp, alpha))
      if (present(indices)) then
         flrr = lrr_default(object%u(:, indices), orthonormalize=.false., info=status)
         blrr = lrr_default(object%u(:, indices), reverse=.true., orthonormalize=.false., info=status)
      else
         flrr = lrr_default(object%u, orthonormalize=.false., info=status)
         blrr = lrr_default(object%u, reverse=.true., orthonormalize=.false., info=status)
      end if
      if (status /= rssa_success) then
         if (present(info)) info = status
         return
      end if
      do i = 1, size(incomplete)
         if (.not. missing(i)) cycle
         left = i - 1
         right = i + 1
         if (left >= size(flrr)) then
            forward = apply_lrr(filled(max(1, left - size(flrr) + 1):left), flrr, 1, only_new=.true., info=status)
         else
            allocate(forward(0))
         end if
         if (right + size(blrr) - 1 <= size(filled)) then
            backward = apply_lrr(filled(right:right + size(blrr) - 1), blrr, 1, reverse=.true., only_new=.true., info=status)
         else
            allocate(backward(0))
         end if
         if (size(forward) == 1 .and. size(backward) == 1) then
            filled(i) = (1.0_dp - blend) * forward(1) + blend * backward(1)
         else if (size(forward) == 1) then
            filled(i) = forward(1)
         else if (size(backward) == 1) then
            filled(i) = backward(1)
         end if
         if (allocated(forward)) deallocate(forward)
         if (allocated(backward)) deallocate(backward)
      end do
      if (present(info)) info = rssa_success
   end function gapfill_ssa

   function gapfill_mssa_channel(object, incomplete, channel, indices, alpha, info) result(filled)
      type(mssa_result), intent(in) :: object !! MSSA template supplying the common left subspace.
      real(dp), intent(in) :: incomplete(:) !! One channel with NaNs.
      integer, intent(in) :: channel !! One-based channel index used for validation.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples.
      real(dp), intent(in), optional :: alpha !! Forward/backward blend.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      real(dp), allocatable :: filled(:)
      type(ssa_result) :: proxy
      if (channel < 1 .or. channel > size(object%series, 2)) then
         allocate(filled(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      proxy%window = object%window
      proxy%u = object%u
      proxy%sigma = object%sigma
      filled = gapfill_ssa(proxy, incomplete, indices, alpha, info)
   end function gapfill_mssa_channel

   function gapfill_complex(object, incomplete, indices, info) result(filled)
      type(cssa_result), intent(in) :: object !! CSSA template.
      complex(dp), intent(in) :: incomplete(:) !! Complex series with NaNs in either part.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples; accepted for API parity.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      complex(dp), allocatable :: filled(:)
      integer :: rank, status
      rank = 2
      if (present(indices)) rank = max(1, size(indices))
      filled = igapfill_complex(incomplete, object%window, rank, info=status)
      if (present(info)) info = status
   end function gapfill_complex

   function igapfill(incomplete, window, rank, fill, tolerance, maxiter, info) result(filled)
      real(dp), intent(in) :: incomplete(:) !! Input series with NaNs at values to estimate.
      integer, intent(in), optional :: window !! SSA window length.
      integer, intent(in), optional :: rank !! Signal rank; default two.
      real(dp), intent(in), optional :: fill !! Initial scalar filler; default finite mean.
      real(dp), intent(in), optional :: tolerance !! Maximum missing-cell update tolerance; default `1e-6`.
      integer, intent(in), optional :: maxiter !! Iteration limit; default 100.
      integer, intent(out), optional :: info !! Zero on convergence or Rssa status.
      real(dp), allocatable :: filled(:), previous(:), reconstruction(:)
      logical, allocatable :: missing(:)
      type(ssa_result) :: fit
      real(dp) :: initial, tol, change
      integer :: i, iter, l, nrank, limit, status
      l = (size(incomplete) + 1) / 2
      if (present(window)) l = window
      nrank = 2
      if (present(rank)) nrank = rank
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      limit = 100
      if (present(maxiter)) then
         if (maxiter > 0) limit = maxiter
      end if
      missing = ieee_is_nan(incomplete)
      filled = incomplete
      if (count(.not. missing) == 0 .or. l < 1 .or. l > size(incomplete) .or. nrank < 1) then
         if (present(info)) info = rssa_invalid_input
         return
      end if
      if (present(fill)) then
         initial = fill
      else
         initial = sum(pack(incomplete, .not. missing)) / real(count(.not. missing), dp)
      end if
      where (missing) filled = initial
      status = rssa_success
      do iter = 1, limit
         fit = ssa_1d(filled, l, nrank)
         if (fit%info /= rssa_success) then
            status = fit%info
            exit
         end if
         reconstruction = reconstruct_ssa(fit, [(i, i=1, min(nrank, size(fit%sigma)))])
         previous = filled
         where (missing) filled = reconstruction
         change = 0.0_dp
         if (count(missing) > 0) change = maxval(abs(pack(filled - previous, missing)))
         if (change < tol) exit
      end do
      if (present(info)) info = status
   end function igapfill

   function igapfill_ssa(object, incomplete, indices, fill, tolerance, maxiter, info) result(filled)
      type(ssa_result), intent(in) :: object !! SSA template supplying the window and default rank.
      real(dp), intent(in) :: incomplete(:) !! Series with NaNs.
      integer, intent(in), optional :: indices(:) !! Components whose count defines rank.
      real(dp), intent(in), optional :: fill !! Initial scalar filler.
      real(dp), intent(in), optional :: tolerance !! Update tolerance.
      integer, intent(in), optional :: maxiter !! Iteration limit.
      integer, intent(out), optional :: info !! Zero on convergence or Rssa status.
      real(dp), allocatable :: filled(:)
      integer :: rank, status
      rank = 2
      if (allocated(object%sigma)) rank = max(1, size(object%sigma))
      if (present(indices)) rank = max(1, size(indices))
      if (present(fill)) then
         filled = igapfill(incomplete, object%window, rank, fill, tolerance, maxiter, status)
      else
         filled = igapfill(incomplete, object%window, rank, tolerance=tolerance, maxiter=maxiter, info=status)
      end if
      if (present(info)) info = status
   end function igapfill_ssa

   function igapfill_2d(field, window_shape, rank, fill, tolerance, maxiter, info) result(filled)
      real(dp), intent(in) :: field(:, :) !! Field with NaNs at missing cells.
      integer, intent(in) :: window_shape(2) !! Rectangular SSA window.
      integer, intent(in) :: rank !! Retained rank.
      real(dp), intent(in), optional :: fill !! Initial scalar filler.
      real(dp), intent(in), optional :: tolerance !! Update tolerance; default `1e-6`.
      integer, intent(in), optional :: maxiter !! Iteration limit; default 100.
      integer, intent(out), optional :: info !! Zero on convergence or Rssa status.
      real(dp), allocatable :: filled(:, :), previous(:, :), reconstruction(:, :)
      logical, allocatable :: missing(:, :)
      type(ssa2d_result) :: fit
      real(dp) :: initial, tol, change
      integer :: i, iter, limit, status
      missing = ieee_is_nan(field)
      filled = field
      if (rank < 1 .or. any(window_shape < 1) .or. any(window_shape > shape(field))) then
         if (present(info)) info = rssa_invalid_input
         return
      end if
      initial = 0.0_dp
      if (count(.not. missing) > 0) initial = sum(pack(field, .not. missing)) / real(count(.not. missing), dp)
      if (present(fill)) initial = fill
      where (missing) filled = initial
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      limit = 100
      if (present(maxiter)) then
         if (maxiter > 0) limit = maxiter
      end if
      status = rssa_success
      do iter = 1, limit
         fit = ssa_2d(filled, window_shape, rank)
         if (fit%info /= rssa_success) then
            status = fit%info
            exit
         end if
         reconstruction = reconstruct_2d(fit, [(i, i=1, min(rank, size(fit%sigma)))])
         previous = filled
         where (missing) filled = reconstruction
         change = 0.0_dp
         if (count(missing) > 0) change = maxval(abs(pack(filled - previous, missing)))
         if (change < tol) exit
      end do
      if (present(info)) info = status
   end function igapfill_2d

   function igapfill_mssa(series, window, rank, lengths, fill, tolerance, maxiter, info) result(filled)
      real(dp), intent(in) :: series(:, :) !! MSSA channels with NaNs within active lengths.
      integer, intent(in) :: window !! Common MSSA window.
      integer, intent(in) :: rank !! Retained rank.
      integer, intent(in), optional :: lengths(:) !! Active lengths per channel.
      real(dp), intent(in), optional :: fill !! Initial scalar filler.
      real(dp), intent(in), optional :: tolerance !! Update tolerance.
      integer, intent(in), optional :: maxiter !! Iteration limit.
      integer, intent(out), optional :: info !! Zero on convergence or Rssa status.
      real(dp), allocatable :: filled(:, :), previous(:, :), reconstruction(:, :)
      logical, allocatable :: missing(:, :)
      integer, allocatable :: lens(:)
      type(mssa_result) :: fit
      real(dp) :: initial, tol, change
      integer :: i, iter, limit, status
      allocate(lens(size(series, 2)))
      lens = size(series, 1)
      if (present(lengths)) lens = lengths
      missing = ieee_is_nan(series)
      filled = series
      initial = 0.0_dp
      if (count(.not. missing) > 0) initial = sum(pack(series, .not. missing)) / &
            real(count(.not. missing), dp)
      if (present(fill)) initial = fill
      where (missing) filled = initial
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      limit = 100
      if (present(maxiter)) then
         if (maxiter > 0) limit = maxiter
      end if
      status = rssa_success
      do iter = 1, limit
         fit = ssa_mssa(filled, window, rank, lens)
         if (fit%info /= rssa_success) then
            status = fit%info
            exit
         end if
         reconstruction = reconstruct_mssa(fit, [(i, i=1, min(rank, size(fit%sigma)))])
         previous = filled
         where (missing) filled = reconstruction
         change = 0.0_dp
         if (count(missing) > 0) change = maxval(abs(pack(filled - previous, missing)))
         if (change < tol) exit
      end do
      if (present(info)) info = status
   end function igapfill_mssa

   function igapfill_complex(series, window, rank, fill, tolerance, maxiter, info) result(filled)
      complex(dp), intent(in) :: series(:) !! Complex series; NaN in either component marks missing.
      integer, intent(in) :: window !! CSSA window length.
      integer, intent(in) :: rank !! Retained rank.
      complex(dp), intent(in), optional :: fill !! Initial complex filler.
      real(dp), intent(in), optional :: tolerance !! Update tolerance.
      integer, intent(in), optional :: maxiter !! Iteration limit.
      integer, intent(out), optional :: info !! Zero on convergence or Rssa status.
      complex(dp), allocatable :: filled(:), previous(:), reconstruction(:)
      logical, allocatable :: missing(:)
      type(cssa_result) :: fit
      complex(dp) :: initial
      real(dp) :: tol, change
      integer :: i, iter, limit, status
      missing = ieee_is_nan(real(series, dp)) .or. ieee_is_nan(aimag(series))
      filled = series
      initial = (0.0_dp, 0.0_dp)
      if (count(.not. missing) > 0) initial = sum(pack(series, .not. missing)) / &
            cmplx(real(count(.not. missing), dp), 0.0_dp, dp)
      if (present(fill)) initial = fill
      where (missing) filled = initial
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      limit = 100
      if (present(maxiter)) then
         if (maxiter > 0) limit = maxiter
      end if
      status = rssa_success
      do iter = 1, limit
         fit = ssa_complex(filled, window, rank)
         if (fit%info /= rssa_success) then
            status = fit%info
            exit
         end if
         reconstruction = reconstruct_complex(fit, [(i, i=1, min(rank, size(fit%sigma)))])
         previous = filled
         where (missing) filled = reconstruction
         change = 0.0_dp
         if (count(missing) > 0) change = maxval(abs(pack(filled - previous, missing)))
         if (change < tol) exit
      end do
      if (present(info)) info = status
   end function igapfill_complex
end module rssa_gapfill
