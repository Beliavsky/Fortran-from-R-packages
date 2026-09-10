! SPDX-License-Identifier: GPL-2.0-or-later
program test_utilities
   use wavethresh
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp) :: v(4)
   real(dp) :: t(3)
   real(dp) :: y(3)
   real(dp) :: response(8)
   real(dp) :: matrix(8, 4)
   real(dp), allocatable :: rotated(:)
   real(dp), allocatable :: dwwt(:)
   real(dp), allocatable :: direct_packet(:)
   integer, allocatable :: grot(:)
   integer, allocatable :: arrvec(:,:)
   type(wd_t) :: zero_wst
   type(first_last_t) :: fl
   type(first_last_t) :: fl_dh
   type(grid_data_t) :: grid
   type(scaling_function_t) :: sf
   type(basis_selection_t) :: selection
   type(wp_t) :: stationary_packets
   type(wpst_matrix_t) :: packet_matrix
   type(wpst_matrix_t) :: discrimination_matrix
   type(wpst_matrix_t) :: selected_packets
   type(node_vector_t) :: node_vector
   type(wt_filter_t) :: yates
   type(multiple_filter_t) :: multiple_filter
   type(first_last_t) :: multiple_fl
   real(dp), allocatable :: av_reconstruction(:)
   real(dp), allocatable :: direct_reconstruction(:)
   type(wd_t) :: av_object
   integer :: i
   integer :: group_labels(8)

   v = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   rotated = guyrot(v, 1)
   call assert_vector_close(rotated, [4.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], tol, "guyrot")
   rotated = rotateback(v)
   call assert_vector_close(rotated, [4.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], tol, "rotateback")

   grot = compgrot(4, 1.0_dp, "DaubExPhase")
   call assert_true(all(grot == [0, 1, 3, 7]), "compgrot Haar")
   grot = compgrot(4, 4.0_dp, "DaubExPhase")
   call assert_true(all(grot == [2, 6, 15, 31]), "compgrot general")

   arrvec = getarrvec(3)
   call assert_true(all(shape(arrvec) == [8, 2]), "getarrvec shape")
   call assert_true(all(arrvec(:, 1) == [1, 5, 2, 6, 3, 7, 4, 8]), "getarrvec level 1")
   call assert_true(all(arrvec(:, 2) == [1, 5, 3, 7, 2, 6, 4, 8]), "getarrvec level 2")
   node_vector = numtonv(5, 3)
   call assert_true(node_vector%ok, "numtonv status")
   call assert_true(all(node_vector%node(1)%upperctrl == "R"), "numtonv first controls")
   call assert_true(all(node_vector%node(2)%upperctrl == "L"), "numtonv second controls")
   call assert_true(all(node_vector%node(3)%upperctrl == "R"), "numtonv third controls")
   dwwt = make_dwwt(3, 1.0_dp, "DaubExPhase")
   call assert_vector_close(dwwt, [1.0_dp, 1.0_dp, 1.0_dp], 2.0e-11_dp, "make.dwwt Haar")

   zero_wst = cns(8)
   call assert_true(zero_wst%ok, "cns status")
   do i = 0, zero_wst%nlevels - 1
      call assert_close(maxval(abs(zero_wst%detail(i)%values)), 0.0_dp, tol, "cns detail")
   end do

   fl = first_last(2, 8)
   call assert_true(fl%ok, "first.last status")
   call assert_true(fl%ntotal == 15 .and. fl%ntotal_detail == 7, "first.last totals")
   call assert_true(all(fl%scaling(:, 2) == [0, 1, 3, 7]), "first.last scaling last")
   call assert_true(all(fl%scaling(:, 3) == [14, 12, 8, 0]), "first.last scaling offsets")
   call assert_true(all(fl%detail(:, 3) == [6, 4, 0]), "first.last detail offsets")
   fl_dh = first_last_dh(4, 8, boundary="zero")
   call assert_true(fl_dh%ok, "first.last.dh zero status")
   call assert_true(fl_dh%ntotal == 20 .and. fl_dh%ntotal_detail == 12, "first.last.dh zero totals")
   call assert_true(all(fl_dh%scaling(:, 1) == [-2, -2, -1, 0]), "first.last.dh zero first C")
   call assert_true(all(fl_dh%scaling(:, 3) == [17, 13, 8, 0]), "first.last.dh zero C offsets")
   call assert_true(all(fl_dh%detail(:, 1) == [-1, -1, 0]), "first.last.dh zero first D")
   call assert_true(all(fl_dh%detail(:, 3) == [9, 5, 0]), "first.last.dh zero D offsets")

   t = [0.1_dp, 0.9_dp, 0.4_dp]
   y = [1.0_dp, 9.0_dp, 4.0_dp]
   grid = makegrid(t, y, 4)
   call assert_true(grid%ok, "makegrid status")
   call assert_vector_close(grid%grid_t, [0.125_dp, 0.375_dp, 0.625_dp, 0.875_dp], tol, "makegrid t")
   call assert_vector_close(grid%grid_y, [1.25_dp, 3.75_dp, 6.25_dp, 8.75_dp], tol, "makegrid y")

   sf = scaling_function(1.0_dp, "DaubExPhase", 64, 20)
   call assert_true(sf%ok, "ScalingFunction status")
   call assert_close(minval(sf%y), 1.0_dp, 2.0e-12_dp, "ScalingFunction Haar minimum")
   call assert_close(maxval(sf%y), 1.0_dp, 2.0e-12_dp, "ScalingFunction Haar maximum")

   call assert_close(wvmoments(1.0_dp, "DaubExPhase", 0, .true.), 1.0_dp, 2.0e-8_dp, "wvmoments Haar scaling zeroth moment")
   call assert_close(wvmoments(1.0_dp, "DaubExPhase", 0, .false.), 0.0_dp, 2.0e-8_dp, "wvmoments Haar wavelet zeroth moment")

   do i = 1, 8
      response(i) = real(i, dp)
   end do
   matrix(:, 1) = 1.0_dp
   matrix(:, 2) = response
   matrix(:, 3) = -response
   matrix(:, 4) = [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
   selection = bestm(matrix, response, 25.0_dp)
   call assert_true(selection%ok, "bestm status")
   call assert_true(all(selection%index == [2, 3]), "bestm selected columns")
   selection = best_1d_cols(matrix, response, 0.9_dp)
   call assert_true(selection%ok, "Best1DCols status")
   call assert_true(all(selection%index == [2, 3]), "Best1DCols selected columns")

   stationary_packets = wpst(response, 1.0_dp, "DaubExPhase")
   direct_packet = access_d_wpst(stationary_packets, 0, 3)
   call assert_vector_close(direct_packet, getpacket(stationary_packets, 3, 3), tol, "accessD.wpst")
   packet_matrix = wpst2m(stationary_packets)
   call assert_true(packet_matrix%ok, "wpst2m status")
   call assert_true(all(shape(packet_matrix%matrix) == [8, 14]), "wpst2m dimensions")
   call assert_true(all(packet_matrix%level(1:8) == 0), "wpst2m finest R-style levels")
   call assert_true(all(packet_matrix%packet_index(1:8) == [(i - 1, i=1, 8)]), "wpst2m packet metadata")
   group_labels = [1, 1, 1, 1, 2, 2, 2, 2]
   discrimination_matrix = wpst2discr(stationary_packets, group_labels)
   call assert_true(discrimination_matrix%ok, "wpst2discr status")
   call assert_true(all(discrimination_matrix%groups == group_labels), "wpst2discr groups")
   call assert_close(maxval(abs(exp(discrimination_matrix%matrix) - packet_matrix%matrix**2)), &
      0.0_dp, 2.0e-11_dp, "wpst2discr log-square transform")
   selected_packets = wpst_regr(response, [0, 1, 2], [3, 1, 0], 1.0_dp, "DaubExPhase")
   call assert_true(selected_packets%ok, "wpstREGR status")
   call assert_vector_close(selected_packets%matrix(:, 1), packet_matrix%matrix(:, 4), tol, "wpstREGR level 0")
   call assert_vector_close(selected_packets%matrix(:, 2), packet_matrix%matrix(:, 10), tol, "wpstREGR level 1")
   call assert_vector_close(selected_packets%matrix(:, 3), packet_matrix%matrix(:, 13), tol, "wpstREGR level 2")

   multiple_filter = mfilter_select("Geronimo")
   call assert_true(multiple_filter%ok, "mfilter.select Geronimo status")
   call assert_true(multiple_filter%nphi == 2 .and. multiple_filter%npsi == 2, "mfilter.select Geronimo ranks")
   call assert_true(size(multiple_filter%h) == 16 .and. size(multiple_filter%g) == 16, &
      "mfilter.select Geronimo lengths")
   call assert_close(multiple_filter%h(2), 0.8_dp, tol, "mfilter.select Geronimo H")
   multiple_filter = mfilter_select("Donovan3")
   call assert_true(multiple_filter%ok, "mfilter.select Donovan3 status")
   call assert_true(multiple_filter%nphi == 3 .and. size(multiple_filter%h) == 36, "mfilter.select Donovan3 rank")
   call assert_close(multiple_filter%h(19), 1.0_dp / sqrt(2.0_dp), tol, "mfilter.select Donovan3 H19")
   call assert_close(multiple_filter%g(28), 13.0_dp / 22.0_dp, tol, "mfilter.select Donovan3 G28")
   multiple_fl = mfirst_last(4, 2, 3)
   call assert_true(multiple_fl%ok, "mfirst.last status")
   call assert_true(multiple_fl%ntotal == 13 .and. multiple_fl%ntotal_detail == 4, "mfirst.last totals")
   call assert_true(all(multiple_fl%scaling(:, 2) == [0, 2, 8]), "mfirst.last scaling last")
   call assert_true(all(multiple_fl%scaling(:, 3) == [12, 9, 0]), "mfirst.last scaling offsets")
   call assert_true(all(multiple_fl%detail(:, 3) == [3, 0]), "mfirst.last detail offsets")

   yates = filter_select(1.0_dp, "Yates")
   call assert_true(yates%ok, "Yates filter status")
   call assert_true(size(yates%low) == 2, "Yates filter length")
   call assert_vector_close(yates%low, [-1.0_dp / sqrt(2.0_dp), 1.0_dp / sqrt(2.0_dp)], tol, "Yates coefficients")

   av_object = wst(response, filter_number=1.0_dp, family="DaubExPhase")
   av_reconstruction = av_basis(av_object, av_object%nlevels - 1, 0, 1)
   direct_reconstruction = wr_wst(av_object)
   call assert_true(size(av_reconstruction) == size(response), "av.basis reconstruction size")
   call assert_true(maxval(abs(av_reconstruction - direct_reconstruction)) < 2.0e-11_dp, &
      "av.basis agrees with average stationary reconstruction")

   print *, "test_utilities: PASS"

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must evaluate true for the test to pass.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Scalar value produced by the translated implementation.
      real(dp), intent(in) :: expected !! Deterministic scalar reference value.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference from the reference.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      if (abs(actual - expected) > tolerance) then
         print *, "FAIL: ", trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_vector_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:) !! Vector produced by the translated implementation.
      real(dp), intent(in) :: expected(:) !! Deterministic vector of reference values.
      real(dp), intent(in) :: tolerance !! Maximum permitted elementwise absolute error.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      if (size(actual) /= size(expected)) then
         print *, "FAIL size: ", trim(message)
         error stop 1
      end if
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tolerance) then
            print *, "FAIL: ", trim(message)
            error stop 1
         end if
      end if
   end subroutine assert_vector_close

end program test_utilities
