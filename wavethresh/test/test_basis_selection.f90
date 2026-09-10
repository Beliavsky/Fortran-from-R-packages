program test_basis_selection
   use wavethresh_types, only : dp, wd_t, wp_t, node_vector_t, wpst_regression_t, wpst_matrix_t
   use wavethresh_transform_1d, only : wst, wp, wpst
   use wavethresh_basis, only : manove_wst, inv_basis_wst, manove_wp
   use wavethresh_access, only : access_d_wpst, getpacket_wpst
   use wavethresh_stats, only : c2to4
   use wavethresh_utilities, only : numtonv, makewpst_ro, wpst_regr
   implicit none

   real(dp) :: x(8)
   real(dp) :: zeroes(8)
   real(dp) :: response(8)
   real(dp), allocatable :: reconstructed(:)
   real(dp), allocatable :: woven(:)
   real(dp), allocatable :: rebuilt(:)
   real(dp), allocatable :: raw(:)
   integer, allocatable :: primary(:)
   integer, allocatable :: old_primary(:)
   type(wd_t) :: stationary
   type(wd_t) :: stationary_zero
   type(wp_t) :: packets
   type(wp_t) :: packets_zero
   type(wp_t) :: stationary_packets
   type(node_vector_t) :: node_vector
   type(wpst_regression_t) :: model
   type(wpst_matrix_t) :: features
   integer :: code
   integer :: depth
   integer :: em
   integer :: i
   integer :: level
   integer :: nold
   integer :: ordinary_index
   integer :: packet_length
   integer :: packet_position
   integer :: phase
   integer :: packets_per_basis

   x = [1.0_dp, -2.0_dp, 3.0_dp, 4.0_dp, -1.0_dp, 0.5_dp, 2.5_dp, -3.0_dp]
   zeroes = 0.0_dp
   response = [0.0_dp, 1.0_dp, 0.5_dp, 2.0_dp, -1.0_dp, 0.25_dp, 1.5_dp, -0.5_dp]

   stationary = wst(x, filter_number=1.0_dp, family="DaubExPhase")
   if (.not. stationary%ok) error stop "wst failed in best-basis test"
   do code = 0, 2**stationary%nlevels - 1
      node_vector = numtonv(code, stationary%nlevels)
      if (.not. node_vector%ok) error stop "numtonv failed in InvBasis.wst test"
      reconstructed = inv_basis_wst(stationary, node_vector)
      if (size(reconstructed) /= size(x)) error stop "InvBasis.wst returned the wrong length"
      if (maxval(abs(reconstructed - x)) > 1.0e-12_dp) error stop "InvBasis.wst path reconstruction failed"
   end do

   node_vector = manove_wst(stationary)
   if (.not. node_vector%ok) error stop "MaNoVe.wst failed"
   reconstructed = inv_basis_wst(stationary, node_vector)
   if (size(reconstructed) /= size(x)) error stop "MaNoVe.wst/InvBasis.wst returned the wrong length"
   if (maxval(abs(reconstructed - x)) > 1.0e-12_dp) error stop "MaNoVe.wst basis reconstruction failed"

   stationary_zero = wst(zeroes, filter_number=1.0_dp, family="DaubExPhase")
   node_vector = manove_wst(stationary_zero)
   if (.not. node_vector%ok) error stop "MaNoVe.wst failed on zero data"
   do level = 1, node_vector%nlevels
      if (any(node_vector%node(level)%upperctrl /= "R")) &
         error stop "MaNoVe.wst does not match the upstream C tie convention"
   end do

   packets = wp(x, filter_number=1.0_dp, family="DaubExPhase")
   node_vector = manove_wp(packets)
   if (.not. node_vector%ok) error stop "MaNoVe.wp failed"
   do level = 1, node_vector%nlevels
      if (any(node_vector%node(level)%upperctrl /= "T" .and. node_vector%node(level)%upperctrl /= "B")) &
         error stop "MaNoVe.wp returned an invalid control code"
   end do

   packets_zero = wp(zeroes, filter_number=1.0_dp, family="DaubExPhase")
   node_vector = manove_wp(packets_zero)
   if (.not. node_vector%ok) error stop "MaNoVe.wp failed on zero data"
   do level = 1, node_vector%nlevels
      if (any(node_vector%node(level)%upperctrl /= "B")) &
         error stop "MaNoVe.wp does not match the upstream C tie convention"
   end do

   stationary_packets = wpst(x, filter_number=1.0_dp, family="DaubExPhase")
   if (.not. stationary_packets%ok) error stop "wpst failed in raw-packet test"
   do level = 0, stationary_packets%nlevels - 1
      depth = stationary_packets%nlevels - level
      packets_per_basis = 2**depth
      packet_length = 2**level
      do ordinary_index = 0, packets_per_basis - 1
         allocate(primary(1))
         primary(1) = c2to4(ordinary_index)
         do i = level, stationary_packets%nlevels - 1
            em = 2**(2 * stationary_packets%nlevels - 2 * i - 1)
            nold = size(primary)
            old_primary = primary
            deallocate(primary)
            allocate(primary(2 * nold))
            primary(1:nold) = old_primary
            primary(nold + 1:2 * nold) = em + old_primary
         end do
         woven = access_d_wpst(stationary_packets, level, ordinary_index)
         allocate(rebuilt(size(woven)), source=0.0_dp)
         do packet_position = 1, packet_length
            do phase = 1, packets_per_basis
               raw = getpacket_wpst(stationary_packets, level, primary(phase))
               if (size(raw) /= packet_length) error stop "getpacket.wpst returned the wrong packet length"
               rebuilt((packet_position - 1) * packets_per_basis + phase) = raw(packet_position)
            end do
         end do
         if (maxval(abs(rebuilt - woven)) > 1.0e-13_dp) error stop "getpacket.wpst interweaving mismatch"
         deallocate(primary, rebuilt)
      end do
   end do
   associate (root_packet => getpacket_wpst(stationary_packets, stationary_packets%nlevels, 0))
      if (maxval(abs(root_packet - x)) > 1.0e-13_dp) error stop "getpacket.wpst root packet mismatch"
   end associate
   associate (invalid_packet => getpacket_wpst(stationary_packets, stationary_packets%nlevels, 1))
      if (size(invalid_packet) /= 0) error stop "getpacket.wpst accepted an invalid root packet index"
   end associate

   model = makewpst_ro(x, response, filter_number=1.0_dp, family="DaubExPhase", &
      transform_name="identity", percentage=50.0_dp)
   if (.not. model%ok) error stop "makewpstRO failed"
   if (size(model%matrix, 1) /= 8 .or. size(model%matrix, 2) /= 4) error stop "makewpstRO matrix shape mismatch"
   features = wpst_regr(x, model%level, model%packet_index, filter_number=model%filter%filter_number, &
      family=trim(model%filter%family), transform_name=trim(model%transform))
   if (.not. features%ok) error stop "wpstREGR failed for the makewpstRO selection"
   if (maxval(abs(features%matrix - model%matrix)) > 1.0e-13_dp) error stop "makewpstRO/wpstREGR mismatch"

   print '(a)', "basis-selection tests passed"
end program test_basis_selection
