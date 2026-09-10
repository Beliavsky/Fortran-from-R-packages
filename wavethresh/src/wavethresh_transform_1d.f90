! SPDX-License-Identifier: GPL-2.0-or-later
! Dense one-dimensional wavelet transforms for wavethresh 4.7.3.
module wavethresh_transform_1d
   use wavethresh_types, only : dp, wd_t, wp_t, wt_vector_t, wt_filter_t
   use wavethresh_filters, only : filter_select
   use waveslim_transform_1d, only : dwt_step, idwt_step, modwt_step
   implicit none
   private

   public :: wd, wr_wd, wst, wr_wst, av_basis, wp, wpst, conbar, wst_conbar
   public :: getpacket_wp, putpacket_wp
   public :: is_power_of_two, nlevels_from_length

contains

   pure function is_power_of_two(n) result(answer)
      integer, intent(in) :: n !! Candidate positive integer length.
      logical :: answer
      integer :: value
      if (n < 1) then
         answer = .false.
         return
      end if
      value = n
      do while (mod(value, 2) == 0 .and. value > 1)
         value = value / 2
      end do
      answer = value == 1
   end function is_power_of_two

   pure function nlevels_from_length(n) result(levels)
      integer, intent(in) :: n !! Data length; must be a positive power of two.
      integer :: levels
      integer :: value
      if (.not. is_power_of_two(n)) then
         levels = -1
         return
      end if
      levels = 0
      value = n
      do while (value > 1)
         value = value / 2
         levels = levels + 1
      end do
   end function nlevels_from_length

   function wd(data, filter_number, family, boundary, nlevels) result(object)
      real(dp), intent(in) :: data(:) !! Input real series; length must be a power of two.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: boundary !! Boundary mode; periodic is implemented in this release.
      integer, intent(in), optional :: nlevels !! Number of decomposition levels; default is the full dyadic depth.
      type(wd_t) :: object
      type(wt_vector_t), allocatable :: temp_detail(:)
      type(wt_vector_t), allocatable :: temp_scaling(:)
      real(dp), allocatable :: work(:)
      real(dp), allocatable :: detail(:)
      real(dp), allocatable :: smooth(:)
      real(dp) :: fnum
      character(len=24) :: fam
      character(len=16) :: bc
      integer :: maximum
      integer :: levels
      integer :: step
      integer :: level

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      bc = "periodic"
      if (present(boundary)) bc = boundary
      object%n_original = size(data)
      object%boundary = bc
      object%filter = filter_select(fnum, trim(fam))
      if (.not. object%filter%ok) then
         object%message = object%filter%message
         return
      end if
      if (object%filter%is_complex) then
         object%message = "wd real API requires a real wavethresh filter"
         return
      end if
      if (trim(bc) /= "periodic") then
         object%message = "this dense wd implementation currently supports periodic boundaries"
         return
      end if
      maximum = nlevels_from_length(size(data))
      if (maximum < 0) then
         object%message = "data length must be a power of two"
         return
      end if
      levels = maximum
      if (present(nlevels)) levels = nlevels
      if (levels < 0 .or. levels > maximum) then
         object%message = "nlevels is outside the valid dyadic range"
         return
      end if
      object%nlevels = levels
      object%transform_type = "wavelet"
      allocate(object%detail(0:levels - 1))
      allocate(object%scaling(0:levels))
      object%scaling(levels)%values = data
      if (levels == 0) then
         object%ok = .true.
         object%message = "ok"
         return
      end if
      allocate(temp_detail(levels))
      allocate(temp_scaling(levels))
      work = data
      do step = 1, levels
         call dwt_step(work, object%filter%high, object%filter%low, detail, smooth)
         temp_detail(step)%values = detail
         temp_scaling(step)%values = smooth
         work = smooth
      end do
      do step = 1, levels
         level = levels - step
         object%detail(level)%values = temp_detail(step)%values
         object%scaling(level)%values = temp_scaling(step)%values
      end do
      object%ok = .true.
      object%message = "ok"
   end function wd

   function wr_wd(object, start_level) result(data)
      type(wd_t), intent(in) :: object !! Decimated wavelet decomposition to reconstruct.
      integer, intent(in), optional :: start_level !! Coarsest R-style level to reconstruct from; default is zero.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: work(:)
      real(dp), allocatable :: next(:)
      integer :: first
      integer :: level

      if (.not. object%ok .or. object%nlevels < 1) then
         if (allocated(object%scaling)) then
            data = object%scaling(object%nlevels)%values
         else
            allocate(data(0))
         end if
         return
      end if
      first = 0
      if (present(start_level)) first = start_level
      if (first < 0 .or. first >= object%nlevels) then
         allocate(data(0))
         return
      end if
      work = object%scaling(first)%values
      do level = first, object%nlevels - 1
         call idwt_step(object%detail(level)%values, work, object%filter%high, object%filter%low, next)
         work = next
      end do
      data = work
   end function wr_wd

   function wst(data, filter_number, family, nlevels) result(object)
      real(dp), intent(in) :: data(:) !! Input real series; length must be a power of two.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: nlevels !! Full dyadic stationary depth; when present it must equal log2(size(data)).
      type(wd_t) :: object
      real(dp) :: fnum
      character(len=24) :: fam
      integer :: maximum
      integer :: levels
      integer :: level

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      object%n_original = size(data)
      object%boundary = "periodic"
      object%filter = filter_select(fnum, trim(fam))
      if (.not. object%filter%ok .or. object%filter%is_complex) then
         object%message = "wst requires a supported real filter"
         return
      end if
      maximum = nlevels_from_length(size(data))
      if (maximum < 0) then
         object%message = "data length must be a power of two"
         return
      end if
      levels = maximum
      if (present(nlevels)) levels = nlevels
      if (levels /= maximum) then
         object%message = "wavethresh wst uses the full dyadic decomposition depth"
         return
      end if
      object%nlevels = levels
      object%transform_type = "station"
      allocate(object%detail(0:levels - 1))
      allocate(object%scaling(0:levels))
      object%scaling(levels)%values = data
      if (levels == 0) then
         object%ok = .true.
         object%message = "ok"
         return
      end if
      do level = 0, levels - 1
         allocate(object%detail(level)%values(size(data)), source=0.0_dp)
         allocate(object%scaling(level)%values(size(data)), source=0.0_dp)
      end do
      call wst_decompose_branch(object, data, 1, size(data) / 2 + 1, levels)
      object%ok = .true.
      object%message = "ok"
   end function wst

   function wr_wst(object, start_level) result(data)
      type(wd_t), intent(in) :: object !! Stationary wavelet decomposition to invert by average-basis synthesis.
      integer, intent(in), optional :: start_level !! Lowest R-style detail level retained; default is zero.
      real(dp), allocatable :: data(:)
      integer :: first

      if (.not. object%ok .or. trim(object%transform_type) /= "station") then
         allocate(data(0))
         return
      end if
      if (object%nlevels == 0) then
         data = object%scaling(0)%values
         return
      end if
      first = 0
      if (present(start_level)) first = start_level
      if (first < 0 .or. first >= object%nlevels) then
         allocate(data(0))
         return
      end if
      data = av_basis_core(object, object%nlevels - 1, 0, 1, first)
   end function wr_wst

   function av_basis(object, level, ix1, ix2) result(answer)
      type(wd_t), intent(in) :: object !! Stationary wavelet transform supplying scaling and detail packet rows.
      integer, intent(in) :: level !! Recursive R-style packet level from zero through object%nlevels minus one.
      integer, intent(in) :: ix1 !! Zero-based left packet index at the requested level.
      integer, intent(in) :: ix2 !! Zero-based right packet index at the requested level.
      real(dp), allocatable :: answer(:)

      answer = av_basis_core(object, level, ix1, ix2, 0)
   end function av_basis

   recursive function av_basis_core(object, level, ix1, ix2, minimum_detail_level) result(answer)
      type(wd_t), intent(in) :: object !! Stationary transform supplying packet rows for recursive average synthesis.
      integer, intent(in) :: level !! Recursive R-style packet level currently being synthesized.
      integer, intent(in) :: ix1 !! Zero-based packet index used for the unshifted synthesis branch.
      integer, intent(in) :: ix2 !! Zero-based packet index used for the one-sample-shifted synthesis branch.
      integer, intent(in) :: minimum_detail_level !! Lowest detail level retained; lower-frequency detail rows are zeroed.
      real(dp), allocatable :: answer(:)
      real(dp), allocatable :: scaling(:)
      real(dp), allocatable :: detail(:)
      real(dp), allocatable :: left(:)
      real(dp), allocatable :: right(:)
      real(dp), allocatable :: partial(:)

      allocate(scaling(0), detail(0), left(0), right(0), partial(0))
      if (.not. object%ok .or. trim(object%transform_type) /= "station") then
         allocate(answer(0))
         return
      end if
      if (level < 0 .or. level >= object%nlevels .or. ix1 < 0 .or. ix2 < 0) then
         allocate(answer(0))
         return
      end if
      if (minimum_detail_level < 0 .or. minimum_detail_level >= object%nlevels) then
         allocate(answer(0))
         return
      end if

      if (level == 0) then
         scaling = wst_packet(object, 0, ix1, .true.)
         call select_wst_detail(object, 0, ix1, minimum_detail_level, detail)
         if (size(scaling) == 0 .or. size(detail) == 0) then
            allocate(answer(0))
            return
         end if
         call wst_synthesis_step(scaling, detail, object%filter%low, left)
         scaling = wst_packet(object, 0, ix2, .true.)
         call select_wst_detail(object, 0, ix2, minimum_detail_level, detail)
         if (size(scaling) == 0 .or. size(detail) == 0) then
            allocate(answer(0))
            return
         end if
         call wst_synthesis_step(scaling, detail, object%filter%low, partial)
         right = rotate_right_one(partial)
      else
         scaling = av_basis_core(object, level - 1, 2 * ix1, 2 * ix1 + 1, minimum_detail_level)
         call select_wst_detail(object, level, ix1, minimum_detail_level, detail)
         if (size(scaling) == 0 .or. size(detail) == 0) then
            allocate(answer(0))
            return
         end if
         call wst_synthesis_step(scaling, detail, object%filter%low, left)
         scaling = av_basis_core(object, level - 1, 2 * ix2, 2 * ix2 + 1, minimum_detail_level)
         call select_wst_detail(object, level, ix2, minimum_detail_level, detail)
         if (size(scaling) == 0 .or. size(detail) == 0) then
            allocate(answer(0))
            return
         end if
         call wst_synthesis_step(scaling, detail, object%filter%low, partial)
         right = rotate_right_one(partial)
      end if
      if (size(left) /= size(right)) then
         allocate(answer(0))
         return
      end if
      answer = 0.5_dp * (left + right)
   end function av_basis_core

   subroutine select_wst_detail(object, level, index, minimum_detail_level, detail)
      type(wd_t), intent(in) :: object !! Stationary transform containing the packed detail row to select.
      integer, intent(in) :: level !! R-style detail level whose packet is requested.
      integer, intent(in) :: index !! Zero-based packet index within the selected detail row.
      integer, intent(in) :: minimum_detail_level !! Lowest retained detail level; lower rows return zeros.
      real(dp), allocatable, intent(out) :: detail(:) !! Selected packet, or a same-size zero packet when suppressed.
      real(dp), allocatable :: scaling(:)

      if (level >= minimum_detail_level) then
         detail = wst_packet(object, level, index, .false.)
      else
         scaling = wst_packet(object, level, index, .true.)
         allocate(detail(size(scaling)), source=0.0_dp)
      end if
   end subroutine select_wst_detail

   recursive subroutine wst_decompose_branch(object, book, first_one, first_two, level)
      type(wd_t), intent(inout) :: object !! Stationary transform object receiving one recursive branch of coefficients.
      real(dp), intent(in) :: book(:) !! Current smooth packet to split using unshifted and one-sample-shifted DWTs.
      integer, intent(in) :: first_one !! One-based first storage position for the unshifted child packet.
      integer, intent(in) :: first_two !! One-based first storage position for the shifted child packet.
      integer, intent(in) :: level !! Current one-based decomposition level; output is stored in row level minus one.
      real(dp), allocatable :: detail_one(:)
      real(dp), allocatable :: smooth_one(:)
      real(dp), allocatable :: detail_two(:)
      real(dp), allocatable :: smooth_two(:)
      real(dp), allocatable :: shifted(:)
      integer :: half

      if (level <= 0 .or. size(book) < 2) return
      call wst_analysis_step(book, object%filter%low, detail_one, smooth_one)
      half = size(smooth_one)
      object%detail(level - 1)%values(first_one:first_one + half - 1) = detail_one
      object%scaling(level - 1)%values(first_one:first_one + half - 1) = smooth_one
      allocate(shifted(size(book)))
      shifted(1:size(book) - 1) = book(2:size(book))
      shifted(size(book)) = book(1)
      call wst_analysis_step(shifted, object%filter%low, detail_two, smooth_two)
      object%detail(level - 1)%values(first_two:first_two + half - 1) = detail_two
      object%scaling(level - 1)%values(first_two:first_two + half - 1) = smooth_two
      if (half > 1) then
         call wst_decompose_branch(object, smooth_one, first_one, first_one + half / 2, level - 1)
         call wst_decompose_branch(object, smooth_two, first_two, first_two + half / 2, level - 1)
      end if
   end subroutine wst_decompose_branch

   pure subroutine wst_analysis_step(input, low, detail, smooth)
      !! Applies the periodic analysis convolutions used by upstream wavepackst.
      real(dp), intent(in) :: input(:) !! Even-length scaling packet to decompose.
      real(dp), intent(in) :: low(:) !! Upstream wavethresh low-pass filter coefficients.
      real(dp), allocatable, intent(out) :: detail(:) !! Decimated mother-wavelet coefficients.
      real(dp), allocatable, intent(out) :: smooth(:) !! Decimated father-wavelet coefficients.
      integer :: coefficient
      integer :: filter_index
      integer :: input_index
      real(dp) :: detail_sign

      allocate(detail(size(input) / 2), source=0.0_dp)
      allocate(smooth(size(input) / 2), source=0.0_dp)
      do coefficient = 0, size(detail) - 1
         do filter_index = 0, size(low) - 1
            input_index = modulo(2 * coefficient + filter_index, size(input)) + 1
            smooth(coefficient + 1) = smooth(coefficient + 1) + low(filter_index + 1) * input(input_index)
            input_index = modulo(2 * coefficient + 1 - filter_index, size(input)) + 1
            detail_sign = merge(-1.0_dp, 1.0_dp, modulo(filter_index, 2) == 0)
            detail(coefficient + 1) = detail(coefficient + 1) + &
               detail_sign * low(filter_index + 1) * input(input_index)
         end do
      end do
   end subroutine wst_analysis_step

   pure subroutine wst_synthesis_step(smooth, detail, low, output)
      !! Applies the periodic synthesis convolution used by upstream av.basis.
      real(dp), intent(in) :: smooth(:) !! Father-wavelet coefficients at one packet node.
      real(dp), intent(in) :: detail(:) !! Conforming mother-wavelet coefficients.
      real(dp), intent(in) :: low(:) !! Upstream wavethresh low-pass filter coefficients.
      real(dp), allocatable, intent(out) :: output(:) !! Reconstructed packet at the next finer scale.
      integer :: coefficient
      integer :: filter_index
      integer :: output_index
      real(dp) :: detail_sign

      allocate(output(2 * size(smooth)), source=0.0_dp)
      do coefficient = 0, size(smooth) - 1
         do filter_index = 0, size(low) - 1
            output_index = modulo(2 * coefficient + filter_index, size(output)) + 1
            output(output_index) = output(output_index) + low(filter_index + 1) * smooth(coefficient + 1)
            output_index = modulo(2 * coefficient + 1 - filter_index, size(output)) + 1
            detail_sign = merge(-1.0_dp, 1.0_dp, modulo(filter_index, 2) == 0)
            output(output_index) = output(output_index) + &
               detail_sign * low(filter_index + 1) * detail(coefficient + 1)
         end do
      end do
   end subroutine wst_synthesis_step

   function wst_packet(object, level, index, scaling_packet) result(values)
      type(wd_t), intent(in) :: object !! Stationary transform containing the requested packed coefficient row.
      integer, intent(in) :: level !! R-style row level; packet length is 2**level.
      integer, intent(in) :: index !! Zero-based packet index within the selected row.
      logical, intent(in) :: scaling_packet !! Select scaling coefficients when true and detail coefficients when false.
      real(dp), allocatable :: values(:)
      integer :: packet_length
      integer :: first
      integer :: last

      if (.not. object%ok .or. trim(object%transform_type) /= "station") then
         allocate(values(0))
         return
      end if
      if (level < 0 .or. level >= object%nlevels) then
         allocate(values(0))
         return
      end if
      packet_length = 2**level
      if (index < 0 .or. index >= object%n_original / packet_length) then
         allocate(values(0))
         return
      end if
      first = index * packet_length + 1
      last = first + packet_length - 1
      if (scaling_packet) then
         values = object%scaling(level)%values(first:last)
      else
         values = object%detail(level)%values(first:last)
      end if
   end function wst_packet

   pure function rotate_right_one(values) result(rotated)
      real(dp), intent(in) :: values(:) !! Vector to rotate cyclically right by one sample.
      real(dp), allocatable :: rotated(:)
      integer :: n

      n = size(values)
      allocate(rotated(n))
      if (n == 0) return
      rotated(1) = values(n)
      if (n > 1) rotated(2:n) = values(1:n - 1)
   end function rotate_right_one

   function wp(data, filter_number, family, nlevels) result(object)
      real(dp), intent(in) :: data(:) !! Input real series for the decimated wavelet-packet tree.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: nlevels !! Packet-tree depth; default is full dyadic depth.
      type(wp_t) :: object
      real(dp), allocatable :: detail(:)
      real(dp), allocatable :: smooth(:)
      real(dp) :: fnum
      character(len=24) :: fam
      integer :: maximum
      integer :: levels
      integer :: depth
      integer :: node

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      object%n_original = size(data)
      object%filter = filter_select(fnum, trim(fam))
      if (.not. object%filter%ok .or. object%filter%is_complex) then
         object%message = "wp requires a supported real filter"
         return
      end if
      maximum = nlevels_from_length(size(data))
      if (maximum < 0) then
         object%message = "data length must be a power of two"
         return
      end if
      levels = maximum
      if (present(nlevels)) levels = nlevels
      if (levels < 0 .or. levels > maximum) then
         object%message = "nlevels is outside the valid dyadic range"
         return
      end if
      object%nlevels = levels
      object%stationary = .false.
      allocate(object%packet(0:levels, 0:2**levels - 1))
      object%packet(0, 0)%values = data
      do depth = 1, levels
         do node = 0, 2**(depth - 1) - 1
            call dwt_step(object%packet(depth - 1, node)%values, object%filter%high, &
               object%filter%low, detail, smooth)
            object%packet(depth, 2 * node)%values = smooth
            object%packet(depth, 2 * node + 1)%values = detail
         end do
      end do
      object%ok = .true.
      object%message = "ok"
   end function wp

   function wpst(data, filter_number, family, nlevels) result(object)
      real(dp), intent(in) :: data(:) !! Input real series for the stationary wavelet-packet tree.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: nlevels !! Packet-tree depth; default is full dyadic depth.
      type(wp_t) :: object
      real(dp), allocatable :: detail(:)
      real(dp), allocatable :: smooth(:)
      real(dp), allocatable :: high(:)
      real(dp), allocatable :: low(:)
      real(dp) :: fnum
      character(len=24) :: fam
      integer :: maximum
      integer :: levels
      integer :: depth
      integer :: node

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      object%n_original = size(data)
      object%filter = filter_select(fnum, trim(fam))
      if (.not. object%filter%ok .or. object%filter%is_complex) then
         object%message = "wpst requires a supported real filter"
         return
      end if
      maximum = nlevels_from_length(size(data))
      if (maximum < 0) then
         object%message = "data length must be a power of two"
         return
      end if
      levels = maximum
      if (present(nlevels)) levels = nlevels
      if (levels < 0 .or. levels > maximum) then
         object%message = "nlevels is outside the valid dyadic range"
         return
      end if
      object%nlevels = levels
      object%stationary = .true.
      allocate(object%packet(0:levels, 0:2**levels - 1))
      object%packet(0, 0)%values = data
      high = object%filter%high / sqrt(2.0_dp)
      low = object%filter%low / sqrt(2.0_dp)
      do depth = 1, levels
         do node = 0, 2**(depth - 1) - 1
            call modwt_step(object%packet(depth - 1, node)%values, depth, high, low, detail, smooth)
            object%packet(depth, 2 * node)%values = smooth
            object%packet(depth, 2 * node + 1)%values = detail
         end do
      end do
      object%ok = .true.
      object%message = "ok"
   end function wpst

   function conbar(c_in, d_in, filter) result(c_out)
      real(dp), intent(in) :: c_in(:) !! Coarse scaling coefficients to synthesize at the next finer resolution.
      real(dp), intent(in) :: d_in(:) !! Detail coefficients paired elementwise with c_in at the same resolution.
      type(wt_filter_t), intent(in) :: filter !! Real wavelet filter descriptor supplying synthesis low/high coefficients.
      real(dp), allocatable :: c_out(:)

      if (.not. filter%ok .or. filter%is_complex .or. size(c_in) /= size(d_in)) then
         allocate(c_out(0))
         return
      end if
      call idwt_step(d_in, c_in, filter%high, filter%low, c_out)
   end function conbar

   pure function wst_conbar(c_in, d_in, filter) result(c_out)
      !! Reconstructs one stationary packet using upstream wavethresh convolution conventions.
      real(dp), intent(in) :: c_in(:) !! Coarse scaling coefficients to synthesize.
      real(dp), intent(in) :: d_in(:) !! Conforming detail coefficients to synthesize.
      type(wt_filter_t), intent(in) :: filter !! Real wavethresh filter supplying low-pass coefficients.
      real(dp), allocatable :: c_out(:)

      if (.not. filter%ok .or. filter%is_complex .or. size(c_in) /= size(d_in)) then
         allocate(c_out(0))
         return
      end if
      call wst_synthesis_step(c_in, d_in, filter%low, c_out)
   end function wst_conbar

   function getpacket_wp(object, level, index) result(values)
      type(wp_t), intent(in) :: object !! Wavelet-packet decomposition containing the requested node.
      integer, intent(in) :: level !! Packet depth, from zero at the root to object%nlevels.
      integer, intent(in) :: index !! Zero-based node index at the selected packet depth.
      real(dp), allocatable :: values(:)
      if (.not. object%ok .or. level < 0 .or. level > object%nlevels) then
         allocate(values(0))
      else if (index < 0 .or. index >= 2**level) then
         allocate(values(0))
      else
         values = object%packet(level, index)%values
      end if
   end function getpacket_wp

   subroutine putpacket_wp(object, level, index, values, ok)
      type(wp_t), intent(inout) :: object !! Wavelet-packet object whose node is to be replaced.
      integer, intent(in) :: level !! Packet depth, from zero at the root to object%nlevels.
      integer, intent(in) :: index !! Zero-based node index at the selected packet depth.
      real(dp), intent(in) :: values(:) !! Replacement packet coefficients; shape must match the node.
      logical, intent(out) :: ok !! True when the requested packet existed and its shape matched.
      ok = .false.
      if (.not. object%ok .or. level < 0 .or. level > object%nlevels) return
      if (index < 0 .or. index >= 2**level) return
      if (.not. allocated(object%packet(level, index)%values)) return
      if (size(values) /= size(object%packet(level, index)%values)) return
      object%packet(level, index)%values = values
      ok = .true.
   end subroutine putpacket_wp

end module wavethresh_transform_1d
