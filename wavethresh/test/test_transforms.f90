! SPDX-License-Identifier: GPL-2.0-or-later
program test_transforms
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use wavethresh
   implicit none
   real(dp), parameter :: tol = 5.0e-12_dp
   real(dp) :: x(8)
   real(dp) :: haar_input(4)
   real(dp), allocatable :: y(:)
   real(dp), allocatable :: packet_values(:)
   real(dp), allocatable :: row(:)
   real(dp), allocatable :: transform_matrix(:,:)
   real(dp) :: identity(8, 8)
   type(wd_t) :: object
   type(wd_t) :: stationary
   type(wd_t) :: haar_stationary
   type(wd_t) :: detail_removed
   type(wp_t) :: packets
   type(wp_t) :: best_packets
   type(support_t) :: bounds
   integer :: i
   integer :: depth
   integer :: node
   integer :: nan_packets
   logical :: update_ok

   x = [1.0_dp, -2.0_dp, 3.0_dp, 0.5_dp, -1.5_dp, 4.0_dp, 2.0_dp, -0.25_dp]
   object = wd(x, filter_number=1.0_dp, family="DaubExPhase")
   call assert_true(object%ok, "wd status")
   y = wr_wd(object)
   call assert_close(maxval(abs(y - x)), 0.0_dp, tol, "wd round trip")
   y = conbar(object%scaling(0)%values, object%detail(0)%values, object%filter)
   call assert_close(maxval(abs(y - object%scaling(1)%values)), 0.0_dp, tol, "conbar synthesis")

   stationary = wst(x, filter_number=1.0_dp, family="DaubExPhase")
   call assert_true(stationary%ok, "wst status")
   y = wr_wst(stationary)
   call assert_close(maxval(abs(y - x)), 0.0_dp, tol, "wst round trip")
   haar_input = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   haar_stationary = wst(haar_input, filter_number=1.0_dp, family="DaubExPhase")
   call assert_close(maxval(abs(haar_stationary%detail(0)%values - [-2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp])), &
      0.0_dp, tol, "wst upstream wavepackst finest Haar row")
   call assert_close(maxval(abs(haar_stationary%detail(1)%values - &
      [-1.0_dp, -1.0_dp, -1.0_dp, 3.0_dp] / sqrt(2.0_dp))), 0.0_dp, tol, &
      "wst upstream wavepackst coarsest Haar row")
   call assert_close(maxval(abs(haar_stationary%scaling(0)%values - 5.0_dp)), 0.0_dp, tol, &
      "wst upstream wavepackst scaling row")
   packet_values = getpacket_wst(stationary, 1, 1)
   call assert_close(maxval(abs(packet_values - stationary%detail(1)%values(3:4))), 0.0_dp, tol, &
      "getpacket.wst detail")
   packet_values = getpacket_wst(stationary, 1, 0, coefficient_type="C")
   call assert_close(maxval(abs(packet_values - stationary%scaling(1)%values(1:2))), 0.0_dp, tol, &
      "getpacket.wst scaling")
   packet_values = getpacket_wst(stationary, 1, 0, aspect="logabs")
   call assert_close(maxval(abs(packet_values - logabs(stationary%detail(1)%values(1:2)))), 0.0_dp, tol, &
      "getpacket.wst logabs")
   call putpacket_wst(stationary, 1, 1, [11.0_dp, 12.0_dp], update_ok)
   call assert_true(update_ok, "putpacket.wst status")
   packet_values = getpacket_wst(stationary, 1, 1)
   call assert_close(maxval(abs(packet_values - [11.0_dp, 12.0_dp])), 0.0_dp, tol, "putpacket.wst values")
   detail_removed = rm_det(stationary)
   do i = 0, detail_removed%nlevels - 1
      call assert_close(maxval(abs(detail_removed%detail(i)%values)), 0.0_dp, tol, "rm.det details")
   end do

   packets = wp(x, filter_number=1.0_dp, family="DaubExPhase")
   call assert_true(packets%ok, "wp status")
   packet_values = getpacket_r(packets, packets%nlevels, 0)
   call assert_close(maxval(abs(packet_values - x)), 0.0_dp, tol, "R-style packet root")
   row = access_d_wp(packets, 0)
   call assert_true(size(row) == size(x), "packet level size")
   call assert_close(sum(row**2), sum(x**2), 2.0e-11_dp, "packet energy")

   best_packets = auto_basis(packets)
   call assert_true(best_packets%ok, "AutoBasis status")
   nan_packets = 0
   do depth = 1, best_packets%nlevels
      do node = 0, 2**depth - 1
         if (ieee_is_nan(best_packets%packet(depth, node)%values(1))) nan_packets = nan_packets + 1
      end do
   end do
   call assert_true(nan_packets > 0, "AutoBasis marks excluded packets")
   y = inv_basis_wp(best_packets)
   call assert_true(size(y) == size(x), "InvBasis.wp result length")
   call assert_close(maxval(abs(y - x)), 0.0_dp, 2.0e-11_dp, "AutoBasis/InvBasis.wp round trip")

   call put_d_wp_level(packets, 0, row(size(row):1:-1), update_ok)
   call assert_true(update_ok, "putD.wp status")
   packet_values = access_d_wp(packets, 0)
   call assert_close(maxval(abs(packet_values - row(size(row):1:-1))), 0.0_dp, tol, "putD.wp row replacement")

   transform_matrix = gen_w(8, 1.0_dp, "DaubExPhase")
   call assert_true(all(shape(transform_matrix) == [8, 8]), "GenW shape")
   identity = 0.0_dp
   do i = 1, 8
      identity(i, i) = 1.0_dp
   end do
   call assert_close(maxval(abs(matmul(transform_matrix, transpose(transform_matrix)) - identity)), &
      0.0_dp, 2.0e-11_dp, "GenW orthogonality")

   bounds = wavelet_support(2.0_dp, "DaubExPhase", 0, 0)
   call assert_true(bounds%ok, "support status")
   call assert_close(bounds%left, -2.0_dp, tol, "support left")
   call assert_close(bounds%right, 4.0_dp, tol, "support right")
   call assert_close(bounds%phi_right, 3.0_dp, tol, "support phi right")

   print *, "test_transforms: PASS"
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must hold.
      character(len=*), intent(in) :: message !! Diagnostic label for a failed assertion.
      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Value produced by the translation.
      real(dp), intent(in) :: expected !! Deterministic reference value.
      real(dp), intent(in) :: tolerance !! Maximum accepted absolute error.
      character(len=*), intent(in) :: message !! Diagnostic label for a failed assertion.
      if (abs(actual - expected) > tolerance) then
         print *, "FAIL: ", trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close
end program test_transforms
