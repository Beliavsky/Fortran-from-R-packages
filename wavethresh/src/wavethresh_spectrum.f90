! SPDX-License-Identifier: GPL-2.0-or-later
! Discrete autocorrelation-wavelet utilities translated from wavethresh 4.7.3.
module wavethresh_spectrum
   use wavethresh_types, only : dp, wt_vector_t, wt_filter_t, wd_t, ewspec_result_t
   use wavethresh_filters, only : filter_select
   use r_linalg, only : inverse_matrix
   implicit none
   private

   public :: psi_j, psi_j_mat, ipndacw, local_spec_wd, local_spec_wst, ewspec

contains

   pure function local_spec_wd(object) result(periodogram)
      type(wd_t), intent(in) :: object !! Stationary wavelet transform whose detail coefficients are converted to raw local
                                       !! energies without linear or nonlinear smoothing.
      type(wd_t) :: periodogram
      integer :: level

      periodogram = object
      if (.not. object%ok) then
         periodogram%ok = .false.
         periodogram%message = "LocalSpec requires a valid stationary wavelet object"
         return
      end if
      if (trim(object%transform_type) /= "station") then
         periodogram%ok = .false.
         periodogram%message = "LocalSpec requires a stationary wavelet transform"
         return
      end if
      if (.not. allocated(object%detail)) then
         periodogram%ok = .false.
         periodogram%message = "LocalSpec requires stored stationary detail coefficients"
         return
      end if
      do level = lbound(periodogram%detail, 1), ubound(periodogram%detail, 1)
         if (.not. allocated(periodogram%detail(level)%values)) cycle
         periodogram%detail(level)%values = periodogram%detail(level)%values**2
      end do
      periodogram%message = "ok"
   end function local_spec_wd

   pure function local_spec_wst(object) result(periodogram)
      type(wd_t), intent(in) :: object !! Stationary wavelet transform passed through the LocalSpec.wst raw-periodogram
                                       !! computational dispatch.
      type(wd_t) :: periodogram

      periodogram = local_spec_wd(object)
   end function local_spec_wst

   function ewspec(object, filter_number, family) result(result)
      type(wd_t), intent(in) :: object !! Stationary wavelet transform used by the DoSWT=FALSE, WPsmooth=FALSE ewspec path.
      real(dp), intent(in), optional :: filter_number !! Real wavelet filter number used for the correction matrix; default is 10.
      character(len=*), intent(in), optional :: family !! Real wavelet family used for the correction matrix; default is
                                                    !! DaubLeAsymm.
      type(ewspec_result_t) :: result
      real(dp), allocatable :: wavelet_periodogram(:,:)
      real(dp), allocatable :: corrected(:,:)
      real(dp), allocatable :: correction_matrix(:,:)
      real(dp) :: fnum
      integer :: info
      integer :: j
      integer :: level
      integer :: n
      character(len=24) :: fam

      result%periodogram = local_spec_wd(object)
      if (.not. result%periodogram%ok) then
         result%message = result%periodogram%message
         return
      end if
      j = result%periodogram%nlevels
      if (j < 1) then
         result%message = "ewspec requires at least one stationary wavelet level"
         return
      end if
      n = result%periodogram%n_original
      if (n < 1) then
         result%message = "ewspec requires a nonempty stationary wavelet object"
         return
      end if
      do level = 0, j - 1
         if (.not. allocated(result%periodogram%detail(level)%values)) then
            result%message = "ewspec requires every stationary detail level"
            return
         end if
         if (size(result%periodogram%detail(level)%values) /= n) then
            result%message = "ewspec requires equal-length stationary detail levels"
            return
         end if
      end do

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      correction_matrix = ipndacw(-j, fnum, trim(fam))
      if (any(shape(correction_matrix) /= [j, j])) then
         result%message = "ewspec could not construct the autocorrelation-wavelet matrix"
         return
      end if
      result%rm = correction_matrix
      call inverse_matrix(correction_matrix, result%irm, info)
      if (info /= 0) then
         result%message = "ewspec autocorrelation-wavelet matrix inversion failed"
         return
      end if

      allocate(wavelet_periodogram(j, n))
      do level = 0, j - 1
         wavelet_periodogram(j - level, :) = result%periodogram%detail(level)%values
      end do
      corrected = matmul(result%irm, wavelet_periodogram)
      result%spectrum = result%periodogram
      do level = 0, j - 1
         result%spectrum%detail(level)%values = corrected(j - level, :)
      end do
      result%ok = .true.
      result%message = "ok"
   end function ewspec

   function psi_j(j, filter_number, family, tol, oplength) result(psi)
      integer, intent(in) :: j !! Finest negative scale requested; values are returned for scales -1 through j.
      real(dp), intent(in), optional :: filter_number !! Wavethresh real-filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh real-filter family; default is DaubLeAsymm.
      real(dp), intent(in), optional :: tol !! Wavelet-response coefficients with magnitude at most this value are discarded.
      integer, intent(in), optional :: oplength !! Maximum total number of autocorrelation-wavelet values; default is 10^7.
      type(wt_vector_t), allocatable :: psi(:)
      type(wt_filter_t) :: filter
      real(dp), allocatable :: scaling(:)
      real(dp), allocatable :: wavelet(:)
      real(dp), allocatable :: packed_wavelet(:)
      real(dp) :: fnum
      real(dp) :: tolerance
      integer :: limit
      integer :: nscale
      integer :: scale
      integer :: factor
      integer :: total_length
      integer :: increment
      character(len=24) :: fam

      allocate(psi(0))
      if (j >= 0) return
      nscale = -j
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      tolerance = 1.0e-100_dp
      if (present(tol)) tolerance = max(0.0_dp, tol)
      limit = 10000000
      if (present(oplength)) limit = max(0, oplength)

      filter = filter_select(fnum, trim(fam))
      if (.not. filter%ok .or. filter%is_complex) return
      if (nscale >= bit_size(factor) - 1) return

      deallocate(psi)
      allocate(psi(nscale))
      scaling = [1.0_dp]
      total_length = 0
      do scale = 1, nscale
         factor = 2**(scale - 1)
         wavelet = linear_convolution(scaling, upsample_filter(filter%high, factor))
         packed_wavelet = pack(wavelet, abs(wavelet) > tolerance)
         if (size(packed_wavelet) == 0) then
            deallocate(psi)
            allocate(psi(0))
            return
         end if
         if (size(packed_wavelet) > ishft(huge(total_length), -1)) then
            deallocate(psi)
            allocate(psi(0))
            return
         end if
         increment = 2 * size(packed_wavelet) - 1
         if (total_length > huge(total_length) - increment) then
            deallocate(psi)
            allocate(psi(0))
            return
         end if
         total_length = total_length + increment
         if (total_length > limit) then
            deallocate(psi)
            allocate(psi(0))
            return
         end if
         psi(scale)%values = autocorrelation_full(packed_wavelet)
         scaling = linear_convolution(scaling, upsample_filter(filter%low, factor))
      end do
   end function psi_j

   function psi_j_mat(j, filter_number, family, tol, oplength) result(matrix)
      integer, intent(in) :: j !! Finest negative scale requested; matrix rows correspond to scales -1 through j.
      real(dp), intent(in), optional :: filter_number !! Wavethresh real-filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh real-filter family; default is DaubLeAsymm.
      real(dp), intent(in), optional :: tol !! Wavelet-response coefficient pruning tolerance forwarded to psi_j.
      integer, intent(in), optional :: oplength !! Maximum total number of autocorrelation-wavelet values; default is 10^7.
      real(dp), allocatable :: matrix(:,:)
      type(wt_vector_t), allocatable :: psi(:)
      integer :: row
      integer :: max_length
      integer :: offset

      psi = psi_j(j, filter_number, family, tol, oplength)
      if (size(psi) == 0) then
         allocate(matrix(0, 0))
         return
      end if
      max_length = size(psi(size(psi))%values)
      allocate(matrix(size(psi), max_length), source=0.0_dp)
      do row = 1, size(psi)
         offset = (max_length - size(psi(row)%values)) / 2
         matrix(row, offset + 1:offset + size(psi(row)%values)) = psi(row)%values
      end do
   end function psi_j_mat

   function ipndacw(j, filter_number, family, tol) result(matrix)
      integer, intent(in) :: j !! Negative matrix order; result has order -j by -j.
      real(dp), intent(in), optional :: filter_number !! Wavethresh real-filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh real-filter family; default is DaubLeAsymm.
      real(dp), intent(in), optional :: tol !! Wavelet-response coefficient pruning tolerance; default is 1e-100.
      real(dp), allocatable :: matrix(:,:)
      type(wt_filter_t) :: filter
      real(dp), allocatable :: psi_matrix(:,:)
      real(dp) :: fnum
      real(dp) :: twoj
      real(dp) :: twol
      real(dp) :: two2j
      real(dp) :: two2jmo
      integer :: nscale
      integer :: row
      integer :: col
      character(len=24) :: fam

      if (j >= 0) then
         allocate(matrix(0, 0))
         return
      end if
      nscale = -j
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      filter = filter_select(fnum, trim(fam))
      if (.not. filter%ok .or. filter%is_complex) then
         allocate(matrix(0, 0))
         return
      end if

      allocate(matrix(nscale, nscale), source=0.0_dp)
      if (size(filter%low) == 2) then
         do row = 1, nscale
            do col = row, nscale
               if (col == row) then
                  twoj = 2.0_dp**real(row, dp)
                  two2j = twoj * twoj
                  matrix(row, col) = (two2j + 5.0_dp) / (3.0_dp * twoj)
               else
                  two2jmo = 2.0_dp**real(2 * row - 1, dp)
                  twol = 2.0_dp**real(col, dp)
                  matrix(row, col) = (two2jmo + 1.0_dp) / twol
               end if
               matrix(col, row) = matrix(row, col)
            end do
         end do
         return
      end if

      psi_matrix = psi_j_mat(j, fnum, trim(fam), tol)
      if (size(psi_matrix, 1) /= nscale) then
         deallocate(matrix)
         allocate(matrix(0, 0))
         return
      end if
      matrix = matmul(psi_matrix, transpose(psi_matrix))
   end function ipndacw

   pure function upsample_filter(coefficients, factor) result(upsampled)
      real(dp), intent(in) :: coefficients(:) !! Compact filter coefficients to separate by factor sample positions.
      integer, intent(in) :: factor !! Positive integer spacing between retained filter coefficients.
      real(dp), allocatable :: upsampled(:)
      integer :: i

      if (factor < 1 .or. size(coefficients) == 0) then
         allocate(upsampled(0))
         return
      end if
      allocate(upsampled((size(coefficients) - 1) * factor + 1), source=0.0_dp)
      do i = 1, size(coefficients)
         upsampled(1 + (i - 1) * factor) = coefficients(i)
      end do
   end function upsample_filter

   pure function linear_convolution(left, right) result(values)
      real(dp), intent(in) :: left(:) !! First finite real sequence in the linear convolution.
      real(dp), intent(in) :: right(:) !! Second finite real sequence in the linear convolution.
      real(dp), allocatable :: values(:)
      integer :: i
      integer :: k

      if (size(left) == 0 .or. size(right) == 0) then
         allocate(values(0))
         return
      end if
      allocate(values(size(left) + size(right) - 1), source=0.0_dp)
      do i = 1, size(left)
         do k = 1, size(right)
            values(i + k - 1) = values(i + k - 1) + left(i) * right(k)
         end do
      end do
   end function linear_convolution

   pure function autocorrelation_full(values) result(correlation)
      real(dp), intent(in) :: values(:) !! Finite real wavelet response whose full discrete autocorrelation is required.
      real(dp), allocatable :: correlation(:)
      integer :: lag
      integer :: i
      integer :: n

      n = size(values)
      if (n == 0) then
         allocate(correlation(0))
         return
      end if
      allocate(correlation(2 * n - 1), source=0.0_dp)
      do lag = -(n - 1), n - 1
         do i = max(1, 1 + lag), min(n, n + lag)
            correlation(lag + n) = correlation(lag + n) + values(i) * values(i - lag)
         end do
      end do
   end function autocorrelation_full

end module wavethresh_spectrum
