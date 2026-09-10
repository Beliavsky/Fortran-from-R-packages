! SPDX-License-Identifier: GPL-2.0-or-later
! Typed coefficient access/update utilities for wavethresh objects.
module wavethresh_access
   use wavethresh_types, only : dp, wd_t, wp_t, mwd_t, imwd_t, imwd_level_t, wd3d_t
   use wavethresh_transform_1d, only : getpacket_wp, putpacket_wp
   use wavethresh_stats, only : logabs
   implicit none
   private

   public :: access_d_wd, access_c_wd, put_d_wd, put_c_wd
   public :: access_d_mwd, access_c_mwd, put_d_mwd, put_c_mwd
   public :: access_d_imwd, put_d_imwd, access_d_wd3d, put_d_wd3d
   public :: getpacket, putpacket, access_d_wp, access_d_wpst, getpacket_wpst, put_d_wp_level, getpacket_r, putpacket_r
   public :: getpacket_wst, putpacket_wst, rm_det
   public :: getpacket_wst2d, putpacket_wst2d
   public :: nlevels_wd, nlevels_wp, nlevels_mwd, nlevels_imwd, nlevels_wd3d
   public :: nullevels_wd, nullevels_imwd
   public :: convert_wd, compress_imwd, uncompress_imwd

contains

   function access_d_wd(object, level) result(values)
      type(wd_t), intent(in) :: object !! Wavelet or stationary decomposition containing detail coefficients.
      integer, intent(in) :: level !! R-style zero-based level, zero denoting the coarsest detail scale.
      real(dp), allocatable :: values(:)
      if (.not. object%ok .or. level < 0 .or. level >= object%nlevels) then
         allocate(values(0))
      else
         values = object%detail(level)%values
      end if
   end function access_d_wd

   function access_c_wd(object, level) result(values)
      type(wd_t), intent(in) :: object !! Wavelet or stationary decomposition containing scaling coefficients.
      integer, intent(in) :: level !! R-style zero-based scaling level; object%nlevels denotes the original data.
      real(dp), allocatable :: values(:)
      if (.not. object%ok .or. level < 0 .or. level > object%nlevels) then
         allocate(values(0))
      else
         values = object%scaling(level)%values
      end if
   end function access_c_wd

   subroutine put_d_wd(object, level, values, ok)
      type(wd_t), intent(inout) :: object !! Wavelet object whose selected detail level is to be replaced.
      integer, intent(in) :: level !! R-style zero-based detail level.
      real(dp), intent(in) :: values(:) !! Replacement detail coefficients; length must match the stored level.
      logical, intent(out) :: ok !! True when the level existed and the replacement length matched.
      ok = .false.
      if (.not. object%ok .or. level < 0 .or. level >= object%nlevels) return
      if (size(values) /= size(object%detail(level)%values)) return
      object%detail(level)%values = values
      ok = .true.
   end subroutine put_d_wd

   subroutine put_c_wd(object, level, values, ok)
      type(wd_t), intent(inout) :: object !! Wavelet object whose selected scaling level is to be replaced.
      integer, intent(in) :: level !! R-style zero-based scaling level.
      real(dp), intent(in) :: values(:) !! Replacement scaling coefficients; length must match the stored level.
      logical, intent(out) :: ok !! True when the level existed and the replacement length matched.
      ok = .false.
      if (.not. object%ok .or. level < 0 .or. level > object%nlevels) return
      if (size(values) /= size(object%scaling(level)%values)) return
      object%scaling(level)%values = values
      ok = .true.
   end subroutine put_c_wd

   function access_d_mwd(object, level) result(values)
      type(mwd_t), intent(in) :: object !! Multiple-wavelet decomposition containing vector-valued detail coefficients.
      integer, intent(in) :: level !! R-style zero-based detail level, from zero through nlevels minus one.
      real(dp), allocatable :: values(:,:)

      if (.not. object%ok .or. level < 0 .or. level >= object%nlevels) then
         allocate(values(0, 0))
      else if (.not. allocated(object%detail(level)%values)) then
         allocate(values(0, 0))
      else
         values = object%detail(level)%values
      end if
   end function access_d_mwd

   function access_c_mwd(object, level) result(values)
      type(mwd_t), intent(in) :: object !! Multiple-wavelet decomposition containing vector-valued scaling coefficients.
      integer, intent(in) :: level !! R-style zero-based scaling level, from zero through nlevels.
      real(dp), allocatable :: values(:,:)

      if (.not. object%ok .or. level < 0 .or. level > object%nlevels) then
         allocate(values(0, 0))
      else if (.not. allocated(object%scaling(level)%values)) then
         allocate(values(0, 0))
      else
         values = object%scaling(level)%values
      end if
   end function access_c_mwd

   subroutine put_d_mwd(object, level, values, ok)
      type(mwd_t), intent(inout) :: object !! Multiple-wavelet decomposition whose detail matrix is to be replaced.
      integer, intent(in) :: level !! R-style zero-based detail level, from zero through nlevels minus one.
      real(dp), intent(in) :: values(:,:) !! Replacement detail matrix; both dimensions must match the stored level.
      logical, intent(out) :: ok !! True when the requested detail level exists and the replacement shape matches.

      ok = .false.
      if (.not. object%ok .or. level < 0 .or. level >= object%nlevels) return
      if (.not. allocated(object%detail(level)%values)) return
      if (any(shape(values) /= shape(object%detail(level)%values))) return
      object%detail(level)%values = values
      ok = .true.
   end subroutine put_d_mwd

   subroutine put_c_mwd(object, level, values, ok)
      type(mwd_t), intent(inout) :: object !! Multiple-wavelet decomposition whose scaling matrix is to be replaced.
      integer, intent(in) :: level !! R-style zero-based scaling level, from zero through nlevels.
      real(dp), intent(in) :: values(:,:) !! Replacement scaling matrix; both dimensions must match the stored level.
      logical, intent(out) :: ok !! True when the requested scaling level exists and the replacement shape matches.

      ok = .false.
      if (.not. object%ok .or. level < 0 .or. level > object%nlevels) return
      if (.not. allocated(object%scaling(level)%values)) return
      if (any(shape(values) /= shape(object%scaling(level)%values))) return
      object%scaling(level)%values = values
      ok = .true.
   end subroutine put_c_mwd

   function access_d_imwd(object, level) result(details)
      type(imwd_t), intent(in) :: object !! Two-dimensional wavelet decomposition containing three detail orientations.
      integer, intent(in) :: level !! R-style zero-based 2-D detail level.
      type(imwd_level_t) :: details
      if (.not. object%ok .or. level < 0 .or. level >= object%nlevels) return
      details = object%level(level)
   end function access_d_imwd

   subroutine put_d_imwd(object, level, details, ok)
      type(imwd_t), intent(inout) :: object !! Two-dimensional wavelet object whose detail triplet is replaced.
      integer, intent(in) :: level !! R-style zero-based 2-D detail level.
      type(imwd_level_t), intent(in) :: details !! Replacement LH, HL, and HH coefficient matrices.
      logical, intent(out) :: ok !! True when all three replacement matrices have the expected shape.
      ok = .false.
      if (.not. object%ok .or. level < 0 .or. level >= object%nlevels) return
      if (.not. allocated(details%lh) .or. .not. allocated(details%hl) .or. .not. allocated(details%hh)) return
      if (any(shape(details%lh) /= shape(object%level(level)%lh))) return
      if (any(shape(details%hl) /= shape(object%level(level)%hl))) return
      if (any(shape(details%hh) /= shape(object%level(level)%hh))) return
      if (allocated(details%smooth)) then
         if (.not. allocated(object%level(level)%smooth)) return
         if (any(shape(details%smooth) /= shape(object%level(level)%smooth))) return
      end if
      object%level(level)%lh = details%lh
      object%level(level)%hl = details%hl
      object%level(level)%hh = details%hh
      if (allocated(details%smooth)) then
         object%level(level)%smooth = details%smooth
      end if
      ok = .true.
   end subroutine put_d_imwd

   function access_d_wd3d(object, level) result(bands)
      type(wd3d_t), intent(in) :: object !! Three-dimensional wavelet decomposition containing seven detail bands.
      integer, intent(in) :: level !! R-style zero-based 3-D detail level.
      real(dp), allocatable :: bands(:,:,:,:)
      if (.not. object%ok .or. level < 0 .or. level >= object%nlevels) then
         allocate(bands(0, 0, 0, 0))
      else
         bands = object%level(level)%band
      end if
   end function access_d_wd3d

   subroutine put_d_wd3d(object, level, bands, ok)
      type(wd3d_t), intent(inout) :: object !! Three-dimensional wavelet object whose selected detail bands are replaced.
      integer, intent(in) :: level !! R-style zero-based 3-D detail level.
      real(dp), intent(in) :: bands(:,:,:,:) !! Replacement seven-band coefficient array.
      logical, intent(out) :: ok !! True when the level exists and the array shape matches.
      ok = .false.
      if (.not. object%ok .or. level < 0 .or. level >= object%nlevels) return
      if (any(shape(bands) /= shape(object%level(level)%band))) return
      object%level(level)%band = bands
      ok = .true.
   end subroutine put_d_wd3d


   function access_d_wp(object, level) result(values)
      type(wp_t), intent(in) :: object !! Decimated wavelet-packet object whose R-style level is requested.
      integer, intent(in) :: level !! R-style packet level; zero is the finest level and nlevels-1 is the coarsest detail row.
      real(dp), allocatable :: values(:)
      integer :: depth
      integer :: node
      integer :: packet_length
      integer :: position
      if (.not. object%ok .or. object%stationary .or. level < 0 .or. level >= object%nlevels) then
         allocate(values(0))
         return
      end if
      depth = object%nlevels - level
      packet_length = size(object%packet(depth, 0)%values)
      allocate(values(object%n_original))
      position = 1
      do node = 0, 2**depth - 1
         values(position:position + packet_length - 1) = object%packet(depth, node)%values
         position = position + packet_length
      end do
   end function access_d_wp

   subroutine put_d_wp_level(object, level, values, ok)
      type(wp_t), intent(inout) :: object !! Decimated packet tree whose complete R-style coefficient row is replaced.
      integer, intent(in) :: level !! R-style zero-based packet level; zero is the finest stored row.
      real(dp), intent(in) :: values(:) !! Replacement row containing exactly object%n_original coefficients.
      logical, intent(out) :: ok !! True when the level and row length are valid and every packet is replaced.
      integer :: depth
      integer :: packet_length
      integer :: node
      integer :: position
      logical :: packet_ok

      ok = .false.
      if (.not. object%ok .or. object%stationary) return
      if (level < 0 .or. level >= object%nlevels) return
      if (size(values) /= object%n_original) return
      depth = object%nlevels - level
      packet_length = size(object%packet(depth, 0)%values)
      position = 1
      do node = 0, 2**depth - 1
         call putpacket_wp(object, depth, node, values(position:position + packet_length - 1), packet_ok)
         if (.not. packet_ok) return
         position = position + packet_length
      end do
      ok = .true.
   end subroutine put_d_wp_level


   function access_d_wpst(object, level, index) result(values)
      type(wp_t), intent(in) :: object !! Stationary packet tree containing the requested nondecimated packet coefficients.
      integer, intent(in) :: level !! R-style zero-based level; zero is finest and nlevels minus one is coarsest.
      integer, intent(in) :: index !! Zero-based packet index; valid range is zero through 2**(nlevels-level)-1.
      real(dp), allocatable :: values(:)
      integer :: depth

      if (.not. object%ok .or. .not. object%stationary .or. level < 0 .or. level >= object%nlevels) then
         allocate(values(0))
         return
      end if
      depth = object%nlevels - level
      if (index < 0 .or. index >= 2**depth) then
         allocate(values(0))
         return
      end if
      values = getpacket_wp(object, depth, index)
   end function access_d_wpst

   function getpacket_wpst(object, level, index) result(values)
      type(wp_t), intent(in) :: object !! Stationary packet tree from which one raw upstream-style packet is requested.
      integer, intent(in) :: level !! R-style raw-packet level; zero is finest and nlevels is the original-data root.
      integer, intent(in) :: index !! Zero-based raw packet index in the range zero through 4**(nlevels-level)-1.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: woven(:)
      integer :: depth
      integer :: digit
      integer :: digit_position
      integer :: high_bit
      integer :: ordinary_index
      integer :: phase_index
      integer :: packet_length
      integer :: packet_position
      integer :: packets_per_basis
      integer :: work_index

      if (.not. object%ok .or. .not. object%stationary) then
         allocate(values(0))
         return
      end if
      if (level < 0 .or. level > object%nlevels .or. index < 0) then
         allocate(values(0))
         return
      end if
      if (level == object%nlevels) then
         if (index /= 0) then
            allocate(values(0))
         else
            values = object%packet(0, 0)%values
         end if
         return
      end if

      depth = object%nlevels - level
      work_index = index
      ordinary_index = 0
      phase_index = 0
      do digit_position = 0, depth - 1
         digit = mod(work_index, 4)
         work_index = work_index / 4
         ordinary_index = ordinary_index + mod(digit, 2) * 2**digit_position
         high_bit = digit / 2
         phase_index = phase_index + high_bit * 2**(depth - 1 - digit_position)
      end do
      if (work_index /= 0) then
         allocate(values(0))
         return
      end if

      woven = access_d_wpst(object, level, ordinary_index)
      if (size(woven) /= object%n_original) then
         allocate(values(0))
         return
      end if
      packet_length = 2**level
      packets_per_basis = 2**depth
      allocate(values(packet_length))
      do packet_position = 0, packet_length - 1
         values(packet_position + 1) = woven(packet_position * packets_per_basis + phase_index + 1)
      end do
   end function getpacket_wpst



   function getpacket_wst(object, level, index, coefficient_type, aspect) result(values)
      type(wd_t), intent(in) :: object !! Stationary wavelet transform whose selected packet-sized row segment is requested.
      integer, intent(in) :: level !! R-style zero-based stationary level; packet length is 2**level.
      integer, intent(in) :: index !! Zero-based packet segment index within the selected stationary coefficient row.
      character(len=*), intent(in), optional :: coefficient_type !! Row selector D or C; default is D.
      character(len=*), intent(in), optional :: aspect !! Elementwise aspect, Identity or logabs; default is Identity.
      real(dp), allocatable :: values(:)
      character(len=16) :: ctype
      character(len=16) :: transform
      integer :: packet_length
      integer :: first
      integer :: last

      ctype = "D"
      if (present(coefficient_type)) ctype = coefficient_type
      transform = "Identity"
      if (present(aspect)) transform = aspect
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
      select case (trim(ctype))
      case ("D")
         values = object%detail(level)%values(first:last)
      case ("C")
         values = object%scaling(level)%values(first:last)
      case default
         allocate(values(0))
         return
      end select
      select case (trim(transform))
      case ("Identity")
      case ("logabs")
         values = logabs(values)
      case default
         deallocate(values)
         allocate(values(0))
      end select
   end function getpacket_wst

   subroutine putpacket_wst(object, level, index, packet, ok)
      type(wd_t), intent(inout) :: object !! Stationary wavelet transform whose detail-row segment is to be replaced.
      integer, intent(in) :: level !! R-style zero-based stationary detail level; packet length is 2**level.
      integer, intent(in) :: index !! Zero-based packet segment index within the selected detail row.
      real(dp), intent(in) :: packet(:) !! Replacement coefficients; length must equal 2**level.
      logical, intent(out) :: ok !! True when the stationary level, packet address, and replacement length are valid.
      integer :: packet_length
      integer :: first
      integer :: last

      ok = .false.
      if (.not. object%ok .or. trim(object%transform_type) /= "station") return
      if (level < 0 .or. level >= object%nlevels) return
      packet_length = 2**level
      if (size(packet) /= packet_length) return
      if (index < 0 .or. index >= object%n_original / packet_length) return
      first = index * packet_length + 1
      last = first + packet_length - 1
      object%detail(level)%values(first:last) = packet
      ok = .true.
   end subroutine putpacket_wst

   pure function getpacket_wst2d(object, level, index, coefficient_type) result(packet)
      !! Extracts one square packet from a stationary two-dimensional transform.
      type(imwd_t), intent(in) :: object !! Stationary two-dimensional wavelet transform.
      integer, intent(in) :: level !! R-style level; the packet side length is two raised to this level.
      character(len=*), intent(in) :: index !! Base-four packet path of length object%nlevels minus level.
      character(len=*), intent(in), optional :: coefficient_type !! Band selector S, H, V, or D; default is S.
      real(dp), allocatable :: packet(:,:)
      character(len=1) :: band
      integer :: x
      integer :: y
      integer :: side
      logical :: valid

      band = "S"
      if (present(coefficient_type)) then
         if (len_trim(coefficient_type) < 1) then
            allocate(packet(0, 0))
            return
         end if
         band = coefficient_type(1:1)
      end if
      call wst2d_packet_coordinates(object, level, index, x, y, side, valid)
      if (.not. valid) then
         allocate(packet(0, 0))
         return
      end if
      select case (band)
      case ("S")
         packet = object%level(level)%smooth(x + 1:x + side, y + 1:y + side)
      case ("H")
         packet = object%level(level)%lh(x + 1:x + side, y + 1:y + side)
      case ("V")
         packet = object%level(level)%hl(x + 1:x + side, y + 1:y + side)
      case ("D")
         packet = object%level(level)%hh(x + 1:x + side, y + 1:y + side)
      case default
         allocate(packet(0, 0))
      end select
   end function getpacket_wst2d

   pure subroutine putpacket_wst2d(object, level, index, coefficient_type, packet, ok)
      !! Replaces one square packet in a stationary two-dimensional transform.
      type(imwd_t), intent(inout) :: object !! Stationary two-dimensional wavelet transform to modify.
      integer, intent(in) :: level !! R-style level; the packet side length is two raised to this level.
      character(len=*), intent(in) :: index !! Base-four packet path of length object%nlevels minus level.
      character(len=*), intent(in) :: coefficient_type !! Band selector S, H, V, or D.
      real(dp), intent(in) :: packet(:,:) !! Square replacement packet of side length two raised to level.
      logical, intent(out) :: ok !! True when the address, band, and packet shape are valid.
      integer :: x
      integer :: y
      integer :: side
      logical :: valid

      ok = .false.
      if (len_trim(coefficient_type) < 1) return
      call wst2d_packet_coordinates(object, level, index, x, y, side, valid)
      if (.not. valid .or. any(shape(packet) /= [side, side])) return
      select case (coefficient_type(1:1))
      case ("S")
         object%level(level)%smooth(x + 1:x + side, y + 1:y + side) = packet
         if (level == 0 .and. allocated(object%smooth)) then
            object%smooth(x + 1:x + side, y + 1:y + side) = packet
         end if
      case ("H")
         object%level(level)%lh(x + 1:x + side, y + 1:y + side) = packet
      case ("V")
         object%level(level)%hl(x + 1:x + side, y + 1:y + side) = packet
      case ("D")
         object%level(level)%hh(x + 1:x + side, y + 1:y + side) = packet
      case default
         return
      end select
      ok = .true.
   end subroutine putpacket_wst2d

   pure subroutine wst2d_packet_coordinates(object, level, index, x, y, side, ok)
      !! Converts an upstream base-four stationary packet path to zero-based matrix coordinates.
      type(imwd_t), intent(in) :: object !! Stationary two-dimensional wavelet transform being addressed.
      integer, intent(in) :: level !! R-style packet level.
      character(len=*), intent(in) :: index !! Base-four packet path, with the finest digit at the right.
      integer, intent(out) :: x !! Zero-based first matrix coordinate.
      integer, intent(out) :: y !! Zero-based second matrix coordinate.
      integer, intent(out) :: side !! Packet side length.
      logical, intent(out) :: ok !! True when the object, level, and path are valid.
      integer :: digit
      integer :: digit_position
      integer :: place
      integer :: ios

      x = 0
      y = 0
      side = 0
      ok = .false.
      if (.not. object%ok .or. .not. object%stationary) return
      if (level < 0 .or. level >= object%nlevels) return
      if (len_trim(index) /= object%nlevels - level) return
      side = 2**level
      place = side
      do digit_position = len_trim(index), 1, -1
         read(index(digit_position:digit_position), *, iostat=ios) digit
         if (ios /= 0 .or. digit < 0 .or. digit > 3) return
         ! Upstream interleaves S/H/V/D blocks and therefore advances by
         ! 2*place.  Separate typed band matrices remove that interleave.
         x = x + modulo(digit, 2) * place
         y = y + (digit / 2) * place
         place = 2 * place
      end do
      if (x + side > object%nrow_original .or. y + side > object%ncol_original) return
      ok = .true.
   end subroutine wst2d_packet_coordinates

   pure function rm_det(object) result(cleaned)
      type(wd_t), intent(in) :: object !! Typed one-dimensional transform whose detail coefficients are to be removed.
      type(wd_t) :: cleaned
      integer :: level

      cleaned = object
      if (.not. cleaned%ok) return
      do level = 0, cleaned%nlevels - 1
         cleaned%detail(level)%values = 0.0_dp
      end do
   end function rm_det

   function getpacket_r(object, level, index) result(values)
      type(wp_t), intent(in) :: object !! Decimated packet tree addressed with wavethresh R level conventions.
      integer, intent(in) :: level !! R-style packet level; zero contains scalar packets and nlevels contains the root.
      integer, intent(in) :: index !! Zero-based packet index at the selected R-style level.
      real(dp), allocatable :: values(:)
      integer :: depth
      if (.not. object%ok .or. object%stationary .or. level < 0 .or. level > object%nlevels) then
         allocate(values(0))
         return
      end if
      depth = object%nlevels - level
      values = getpacket_wp(object, depth, index)
   end function getpacket_r

   subroutine putpacket_r(object, level, index, values, ok)
      type(wp_t), intent(inout) :: object !! Decimated packet tree modified using wavethresh R level conventions.
      integer, intent(in) :: level !! R-style packet level; zero contains scalar packets and nlevels contains the root.
      integer, intent(in) :: index !! Zero-based packet index at the selected R-style level.
      real(dp), intent(in) :: values(:) !! Replacement packet coefficients.
      logical, intent(out) :: ok !! True when the R-style packet address and replacement shape are valid.
      integer :: depth
      ok = .false.
      if (.not. object%ok .or. object%stationary .or. level < 0 .or. level > object%nlevels) return
      depth = object%nlevels - level
      call putpacket_wp(object, depth, index, values, ok)
   end subroutine putpacket_r

   function getpacket(object, level, index) result(values)
      type(wp_t), intent(in) :: object !! Decimated or stationary wavelet-packet object.
      integer, intent(in) :: level !! Zero-based packet-tree depth.
      integer, intent(in) :: index !! Zero-based node number at the selected depth.
      real(dp), allocatable :: values(:)
      values = getpacket_wp(object, level, index)
   end function getpacket

   subroutine putpacket(object, level, index, values, ok)
      type(wp_t), intent(inout) :: object !! Decimated or stationary wavelet-packet object to modify.
      integer, intent(in) :: level !! Zero-based packet-tree depth.
      integer, intent(in) :: index !! Zero-based node number at the selected depth.
      real(dp), intent(in) :: values(:) !! Replacement packet coefficients.
      logical, intent(out) :: ok !! True when the node and replacement length are valid.
      call putpacket_wp(object, level, index, values, ok)
   end subroutine putpacket

   pure function nlevels_wd(object) result(levels)
      type(wd_t), intent(in) :: object !! One-dimensional wavelet object whose depth is requested.
      integer :: levels
      levels = object%nlevels
   end function nlevels_wd

   pure function nlevels_wp(object) result(levels)
      type(wp_t), intent(in) :: object !! Wavelet-packet object whose tree depth is requested.
      integer :: levels
      levels = object%nlevels
   end function nlevels_wp

   pure function nlevels_mwd(object) result(levels)
      type(mwd_t), intent(in) :: object !! Multiple-wavelet object whose decomposition depth is requested.
      integer :: levels
      levels = object%nlevels
   end function nlevels_mwd

   pure function nlevels_imwd(object) result(levels)
      type(imwd_t), intent(in) :: object !! Two-dimensional wavelet object whose depth is requested.
      integer :: levels
      levels = object%nlevels
   end function nlevels_imwd

   pure function nlevels_wd3d(object) result(levels)
      type(wd3d_t), intent(in) :: object !! Three-dimensional wavelet object whose depth is requested.
      integer :: levels
      levels = object%nlevels
   end function nlevels_wd3d

   subroutine nullevels_wd(object, levels)
      type(wd_t), intent(inout) :: object !! One-dimensional transform whose selected details are zeroed.
      integer, intent(in) :: levels(:) !! R-style zero-based detail levels to set to zero.
      integer :: i
      integer :: level
      do i = 1, size(levels)
         level = levels(i)
         if (level >= 0 .and. level < object%nlevels) object%detail(level)%values = 0.0_dp
      end do
   end subroutine nullevels_wd

   subroutine nullevels_imwd(object, levels)
      type(imwd_t), intent(inout) :: object !! Two-dimensional transform whose selected detail triplets are zeroed.
      integer, intent(in) :: levels(:) !! R-style zero-based 2-D detail levels to set to zero.
      integer :: i
      integer :: level
      do i = 1, size(levels)
         level = levels(i)
         if (level < 0 .or. level >= object%nlevels) cycle
         object%level(level)%lh = 0.0_dp
         object%level(level)%hl = 0.0_dp
         object%level(level)%hh = 0.0_dp
      end do
   end subroutine nullevels_imwd

   pure function convert_wd(object, transform_type) result(converted)
      type(wd_t), intent(in) :: object !! Unified dense wavethresh object to relabel between wd and wst representations.
      character(len=*), intent(in) :: transform_type !! Target label, wavelet or station; coefficients are not recomputed.
      type(wd_t) :: converted
      converted = object
      converted%transform_type = transform_type
   end function convert_wd

   pure function compress_imwd(object) result(compressed)
      type(imwd_t), intent(in) :: object !! Dense 2-D object; periodic storage has no redundant boundary coefficients.
      type(imwd_t) :: compressed
      compressed = object
   end function compress_imwd

   pure function uncompress_imwd(object) result(uncompressed)
      type(imwd_t), intent(in) :: object !! Dense 2-D object to return in uncompressed periodic representation.
      type(imwd_t) :: uncompressed
      uncompressed = object
   end function uncompress_imwd

end module wavethresh_access
