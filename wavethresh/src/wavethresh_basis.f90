! SPDX-License-Identifier: GPL-2.0-or-later
! Support and transform-matrix utilities translated from wavethresh 4.7.3.
module wavethresh_basis
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   use wavethresh_types, only : dp, support_t, wd_t, wp_t, node_vector_t
   use wavethresh_transform_1d, only : wd, conbar, is_power_of_two, wst_conbar
   use wavethresh_stats, only : shannon_entropy
   implicit none
   private

   public :: wavelet_support, gen_w, auto_basis, inv_basis_wp, manove_wst, inv_basis_wst, manove_wp

contains

   pure function wavelet_support(filter_number, family, m, n) result(bounds)
      real(dp), intent(in), optional :: filter_number !! Daubechies filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Supported family: DaubExPhase or DaubLeAsymm.
      integer, intent(in), optional :: m !! Upstream dilation index before its internal +1 adjustment; default is zero.
      integer, intent(in), optional :: n !! Integer translation index; default is zero.
      type(support_t) :: bounds
      real(dp) :: fnum
      integer :: scale_index
      integer :: shift
      real(dp) :: a
      real(dp) :: b
      character(len=24) :: fam

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      scale_index = 1
      if (present(m)) scale_index = m + 1
      shift = 0
      if (present(n)) shift = n

      if (trim(fam) /= "DaubExPhase" .and. trim(fam) /= "DaubLeAsymm") then
         bounds%message = "support is implemented for DaubExPhase and DaubLeAsymm"
         return
      end if
      a = -(fnum - 1.0_dp)
      b = fnum
      bounds%left = 2.0_dp**scale_index * (a + real(shift, dp))
      bounds%right = 2.0_dp**scale_index * (b + real(shift, dp))
      bounds%psi_left = a
      bounds%psi_right = b
      bounds%phi_left = 0.0_dp
      bounds%phi_right = 2.0_dp * fnum - 1.0_dp
      bounds%ok = .true.
      bounds%message = "ok"
   end function wavelet_support

   function gen_w(n, filter_number, family, boundary) result(matrix)
      integer, intent(in), optional :: n !! Transform-matrix order; default is 8 and must be a power of two.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: boundary !! Boundary mode; periodic is supported.
      real(dp), allocatable :: matrix(:,:)
      real(dp), allocatable :: basis(:)
      integer :: order
      integer :: i
      integer :: level
      integer :: position
      real(dp) :: fnum
      character(len=24) :: fam
      character(len=16) :: bc

      order = 8
      if (present(n)) order = n
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      bc = "periodic"
      if (present(boundary)) bc = boundary
      if (.not. is_power_of_two(order) .or. order < 2 .or. trim(bc) /= "periodic") then
         allocate(matrix(0, 0))
         return
      end if
      allocate(matrix(order, order), basis(order))
      matrix = 0.0_dp
      do i = 1, order
         basis = 0.0_dp
         basis(i) = 1.0_dp
         associate (object => wd(basis, filter_number=fnum, family=trim(fam), boundary=trim(bc)))
            if (.not. object%ok) then
               deallocate(matrix)
               allocate(matrix(0, 0))
               return
            end if
            position = 1
            matrix(i, position) = object%scaling(0)%values(1)
            position = position + 1
            do level = 0, object%nlevels - 1
               matrix(i, position:position + size(object%detail(level)%values) - 1) = &
                  object%detail(level)%values
               position = position + size(object%detail(level)%values)
            end do
         end associate
      end do
   end function gen_w

   pure function auto_basis(object, zilchtol) result(best)
      type(wp_t), intent(in) :: object !! Decimated wavelet-packet tree from which the minimum-entropy basis is selected.
      real(dp), intent(in), optional :: zilchtol !! Squared-energy tolerance passed to Shannon entropy; default is 1e-8.
      type(wp_t) :: best
      real(dp) :: e
      real(dp) :: e1
      real(dp) :: e2
      real(dp) :: nan_value
      real(dp) :: tolerance
      integer :: depth_child
      integer :: depth_parent
      integer :: j
      integer :: npackets
      integer :: r_level

      best = object
      if (.not. object%ok .or. object%stationary) then
         best%ok = .false.
         best%message = "AutoBasis requires a valid decimated wavelet-packet tree"
         return
      end if
      tolerance = 1.0e-8_dp
      if (present(zilchtol)) tolerance = zilchtol
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)

      do r_level = 1, object%nlevels - 1
         depth_parent = object%nlevels - r_level
         depth_child = depth_parent + 1
         npackets = 2**depth_parent
         do j = 0, npackets - 1
            if (ieee_is_nan(best%packet(depth_child, 2 * j)%values(1)) .or. &
               ieee_is_nan(best%packet(depth_child, 2 * j + 1)%values(1))) then
               best%packet(depth_parent, j)%values = nan_value
               cycle
            end if

            e1 = shannon_entropy(best%packet(depth_child, 2 * j)%values, tolerance)
            e2 = shannon_entropy(best%packet(depth_child, 2 * j + 1)%values, tolerance)
            e = shannon_entropy(best%packet(depth_parent, j)%values, tolerance)
            if (e < e1 + e2 .or. &
               (.not. ieee_is_finite(e) .and. .not. ieee_is_finite(e1) .and. .not. ieee_is_finite(e2))) then
               best%packet(depth_child, 2 * j)%values = nan_value
               best%packet(depth_child, 2 * j + 1)%values = nan_value
            else
               best%packet(depth_parent, j)%values = nan_value
            end if
         end do
      end do
      best%message = "ok"
   end function auto_basis

   function inv_basis_wp(object) result(data)
      type(wp_t), intent(in) :: object !! Packet tree whose finite nodes mark a valid basis, as produced by auto_basis.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: detail(:)
      real(dp), allocatable :: smooth(:)

      if (.not. object%ok .or. object%stationary) then
         allocate(data(0))
         return
      end if
      if (object%nlevels == 0) then
         data = object%packet(0, 0)%values
         return
      end if

      smooth = inv_basis_node(object, 1, 0)
      detail = inv_basis_node(object, 1, 1)
      if (size(smooth) == 0 .or. size(detail) == 0 .or. size(smooth) /= size(detail)) then
         allocate(data(0))
         return
      end if
      data = conbar(smooth, detail, object%filter)
   end function inv_basis_wp

   recursive function inv_basis_node(object, depth, node) result(values)
      type(wp_t), intent(in) :: object !! Packet tree containing the selected-basis marker state.
      integer, intent(in) :: depth !! Internal packet-tree depth, with zero at the original-series root.
      integer, intent(in) :: node !! Zero-based packet node index at the requested internal depth.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: detail(:)
      real(dp), allocatable :: smooth(:)

      if (depth < 0 .or. depth > object%nlevels .or. node < 0 .or. node >= 2**depth) then
         allocate(values(0))
         return
      end if
      if (.not. allocated(object%packet(depth, node)%values)) then
         allocate(values(0))
         return
      end if
      if (.not. ieee_is_nan(object%packet(depth, node)%values(1))) then
         values = object%packet(depth, node)%values
         return
      end if
      if (depth >= object%nlevels) then
         allocate(values(0))
         return
      end if

      smooth = inv_basis_node(object, depth + 1, 2 * node)
      detail = inv_basis_node(object, depth + 1, 2 * node + 1)
      if (size(smooth) == 0 .or. size(detail) == 0 .or. size(smooth) /= size(detail)) then
         allocate(values(0))
         return
      end if
      values = conbar(smooth, detail, object%filter)
   end function inv_basis_node


   function manove_wst(object, zilchtol) result(node_vector)
      type(wd_t), intent(in) :: object !! Stationary wavelet transform for minimum-entropy node selection.
      real(dp), intent(in), optional :: zilchtol !! Squared-energy tolerance used by Shannon entropy; default is 1e-8.
      type(node_vector_t) :: node_vector
      real(dp), allocatable :: packet(:)
      real(dp), allocatable :: scaling(:)
      real(dp), allocatable :: combined(:)
      real(dp) :: mother_entropy
      real(dp) :: left_entropy
      real(dp) :: right_entropy
      real(dp) :: tolerance
      integer :: level
      integer :: node
      integer :: left_node
      integer :: right_node
      integer :: nupper

      if (.not. object%ok .or. trim(object%transform_type) /= "station" .or. object%nlevels < 1) then
         node_vector%message = "MaNoVe.wst requires a valid stationary wavelet transform"
         return
      end if
      tolerance = 1.0e-8_dp
      if (present(zilchtol)) tolerance = zilchtol
      allocate(node_vector%node(object%nlevels))

      do level = 0, object%nlevels - 1
         nupper = 2**(object%nlevels - level - 1)
         allocate(node_vector%node(level + 1)%upperctrl(nupper), source="S")
         allocate(node_vector%node(level + 1)%upperl(nupper), source=0.0_dp)
         do node = 0, nupper - 1
            if (level + 1 == object%nlevels) then
               scaling = object%scaling(object%nlevels)%values
            else
               scaling = basis_wst_packet(object, level + 1, node, .true.)
            end if
            if (size(scaling) == 0) then
               node_vector%message = "MaNoVe.wst could not access a mother scaling packet"
               return
            end if
            mother_entropy = shannon_entropy(scaling, tolerance)
            left_node = 2 * node
            right_node = left_node + 1

            if (level == 0) then
               packet = basis_wst_packet(object, 0, left_node, .false.)
               scaling = basis_wst_packet(object, 0, left_node, .true.)
               if (size(packet) == 0 .or. size(scaling) == 0) then
                  node_vector%message = "MaNoVe.wst could not access a left daughter packet"
                  return
               end if
               allocate(combined(size(packet) + size(scaling)))
               combined(1:size(packet)) = packet
               combined(size(packet) + 1:size(combined)) = scaling
               left_entropy = shannon_entropy(combined, tolerance)
               deallocate(combined)

               packet = basis_wst_packet(object, 0, right_node, .false.)
               scaling = basis_wst_packet(object, 0, right_node, .true.)
               if (size(packet) == 0 .or. size(scaling) == 0) then
                  node_vector%message = "MaNoVe.wst could not access a right daughter packet"
                  return
               end if
               allocate(combined(size(packet) + size(scaling)))
               combined(1:size(packet)) = packet
               combined(size(packet) + 1:size(combined)) = scaling
               right_entropy = shannon_entropy(combined, tolerance)
               deallocate(combined)
            else
               packet = basis_wst_packet(object, level, left_node, .false.)
               if (size(packet) == 0) then
                  node_vector%message = "MaNoVe.wst could not access a left detail packet"
                  return
               end if
               left_entropy = shannon_entropy(packet, tolerance) + &
                  node_vector%node(level)%upperl(left_node + 1)
               packet = basis_wst_packet(object, level, right_node, .false.)
               if (size(packet) == 0) then
                  node_vector%message = "MaNoVe.wst could not access a right detail packet"
                  return
               end if
               right_entropy = shannon_entropy(packet, tolerance) + &
                  node_vector%node(level)%upperl(right_node + 1)
            end if

            if (mother_entropy < left_entropy .and. mother_entropy < right_entropy) then
               node_vector%node(level + 1)%upperctrl(node + 1) = "S"
               node_vector%node(level + 1)%upperl(node + 1) = mother_entropy
            else if (left_entropy < right_entropy) then
               node_vector%node(level + 1)%upperctrl(node + 1) = "L"
               node_vector%node(level + 1)%upperl(node + 1) = left_entropy
            else
               node_vector%node(level + 1)%upperctrl(node + 1) = "R"
               node_vector%node(level + 1)%upperl(node + 1) = right_entropy
            end if
         end do
      end do
      node_vector%nlevels = object%nlevels
      node_vector%ok = .true.
      node_vector%message = "ok"
   end function manove_wst

   function inv_basis_wst(object, node_vector) result(data)
      type(wd_t), intent(in) :: object !! Stationary wavelet transform supplying coefficients for one selected basis path.
      type(node_vector_t), intent(in) :: node_vector !! S/L/R node decisions such as those produced by MaNoVe.wst or numtonv.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: work(:)
      real(dp), allocatable :: detail(:)
      character(len=1), allocatable :: action(:)
      integer, allocatable :: packet_index(:)
      integer :: active_index
      integer :: tree_level
      integer :: reconstruction_step
      integer :: coefficient_level
      integer :: path_position
      integer :: nsteps

      if (.not. object%ok .or. trim(object%transform_type) /= "station") then
         allocate(data(0))
         return
      end if
      if (.not. node_vector%ok .or. node_vector%nlevels /= object%nlevels) then
         allocate(data(0))
         return
      end if
      if (.not. allocated(node_vector%node) .or. size(node_vector%node) /= object%nlevels) then
         allocate(data(0))
         return
      end if

      allocate(action(object%nlevels), source="S")
      allocate(packet_index(object%nlevels), source=0)
      active_index = 0
      nsteps = 0
      do tree_level = object%nlevels, 1, -1
         if (.not. allocated(node_vector%node(tree_level)%upperctrl)) then
            allocate(data(0))
            return
         end if
         if (active_index < 0 .or. active_index >= size(node_vector%node(tree_level)%upperctrl)) then
            allocate(data(0))
            return
         end if
         select case (node_vector%node(tree_level)%upperctrl(active_index + 1))
         case ("S")
            exit
         case ("L")
            active_index = 2 * active_index
         case ("R")
            active_index = 2 * active_index + 1
         case default
            allocate(data(0))
            return
         end select
         nsteps = nsteps + 1
      end do

      if (nsteps == 0) then
         data = object%scaling(object%nlevels)%values
         return
      end if

      ! Re-read the path controls because the active packet index changes after each decision.
      active_index = 0
      do path_position = 1, nsteps
         tree_level = object%nlevels - path_position + 1
         action(path_position) = node_vector%node(tree_level)%upperctrl(active_index + 1)
         if (action(path_position) == "L") then
            active_index = 2 * active_index
         else if (action(path_position) == "R") then
            active_index = 2 * active_index + 1
         else
            allocate(data(0))
            return
         end if
         packet_index(path_position) = active_index
      end do

      coefficient_level = object%nlevels - nsteps
      work = basis_wst_packet(object, coefficient_level, packet_index(nsteps), .true.)
      if (size(work) == 0) then
         allocate(data(0))
         return
      end if
      do reconstruction_step = 1, nsteps
         path_position = nsteps - reconstruction_step + 1
         coefficient_level = object%nlevels - nsteps + reconstruction_step - 1
         detail = basis_wst_packet(object, coefficient_level, packet_index(path_position), .false.)
         if (size(detail) == 0 .or. size(detail) /= size(work)) then
            allocate(data(0))
            return
         end if
         work = wst_conbar(work, detail, object%filter)
         if (action(path_position) == "R") work = rotate_right_local(work)
      end do
      data = work
   end function inv_basis_wst

   function manove_wp(object, zilchtol) result(node_vector)
      type(wp_t), intent(in) :: object !! Decimated wavelet-packet tree for minimum-entropy top/bottom basis decisions.
      real(dp), intent(in), optional :: zilchtol !! Squared-energy tolerance used by Shannon entropy; default is 1e-8.
      type(node_vector_t) :: node_vector
      real(dp), allocatable :: packet(:)
      real(dp) :: mother_entropy
      real(dp) :: daughter_entropy
      real(dp) :: tolerance
      integer :: level
      integer :: node
      integer :: depth
      integer :: nupper

      if (.not. object%ok .or. object%stationary .or. object%nlevels < 1) then
         node_vector%message = "MaNoVe.wp requires a valid decimated wavelet-packet tree"
         return
      end if
      tolerance = 1.0e-8_dp
      if (present(zilchtol)) tolerance = zilchtol
      allocate(node_vector%node(object%nlevels))

      do level = 0, object%nlevels - 1
         nupper = 2**(object%nlevels - level - 1)
         allocate(node_vector%node(level + 1)%upperctrl(nupper), source="T")
         allocate(node_vector%node(level + 1)%upperl(nupper), source=0.0_dp)
         depth = object%nlevels - level - 1
         do node = 0, nupper - 1
            packet = object%packet(depth, node)%values
            mother_entropy = shannon_entropy(packet, tolerance)
            if (level == 0) then
               daughter_entropy = shannon_entropy(object%packet(depth + 1, 2 * node)%values, tolerance) + &
                  shannon_entropy(object%packet(depth + 1, 2 * node + 1)%values, tolerance)
            else
               daughter_entropy = node_vector%node(level)%upperl(2 * node + 1) + &
                  node_vector%node(level)%upperl(2 * node + 2)
            end if
            if (mother_entropy < daughter_entropy) then
               node_vector%node(level + 1)%upperctrl(node + 1) = "T"
               node_vector%node(level + 1)%upperl(node + 1) = mother_entropy
            else
               node_vector%node(level + 1)%upperctrl(node + 1) = "B"
               node_vector%node(level + 1)%upperl(node + 1) = daughter_entropy
            end if
         end do
      end do
      node_vector%nlevels = object%nlevels
      node_vector%ok = .true.
      node_vector%message = "ok"
   end function manove_wp

   function basis_wst_packet(object, level, index, scaling_packet) result(values)
      type(wd_t), intent(in) :: object !! Stationary transform containing the requested packed coefficient row.
      integer, intent(in) :: level !! R-style coefficient level, from zero through nlevels minus one.
      integer, intent(in) :: index !! Zero-based packet index within the requested row.
      logical, intent(in) :: scaling_packet !! True selects scaling coefficients; false selects detail coefficients.
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
   end function basis_wst_packet

   pure function rotate_right_local(values) result(rotated)
      real(dp), intent(in) :: values(:) !! Vector to rotate cyclically right by one sample after an R branch.
      real(dp), allocatable :: rotated(:)
      integer :: n

      n = size(values)
      allocate(rotated(n))
      if (n == 0) return
      rotated(1) = values(n)
      if (n > 1) rotated(2:n) = values(1:n - 1)
   end function rotate_right_local

end module wavethresh_basis
