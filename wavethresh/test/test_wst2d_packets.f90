program test_wst2d_packets
   use wavethresh, only : dp, imwd_t, getpacket_wst2d, putpacket_wst2d
   implicit none

   real(dp), parameter :: expected_s(4) = [1013.0_dp, 1023.0_dp, 1014.0_dp, 1024.0_dp]
   real(dp), parameter :: expected_h(4) = [2013.0_dp, 2023.0_dp, 2014.0_dp, 2024.0_dp]
   real(dp), parameter :: expected_v(4) = [3013.0_dp, 3023.0_dp, 3014.0_dp, 3024.0_dp]
   real(dp), parameter :: expected_d(4) = [4013.0_dp, 4023.0_dp, 4014.0_dp, 4024.0_dp]
   type(imwd_t) :: object
   real(dp), allocatable :: packet(:,:)
   real(dp) :: replacement(2, 2)
   integer :: i
   integer :: j
   integer :: level
   logical :: ok

   object%nrow_original = 8
   object%ncol_original = 8
   object%nlevels = 3
   object%stationary = .true.
   object%ok = .true.
   allocate(object%level(0:2))
   do level = 0, 2
      allocate(object%level(level)%smooth(8, 8), object%level(level)%lh(8, 8), source=0.0_dp)
      allocate(object%level(level)%hl(8, 8), object%level(level)%hh(8, 8), source=0.0_dp)
   end do
   do j = 1, 8
      do i = 1, 8
         object%level(1)%smooth(i, j) = 1000.0_dp + 10.0_dp * i + j
         object%level(1)%lh(i, j) = 2000.0_dp + 10.0_dp * i + j
         object%level(1)%hl(i, j) = 3000.0_dp + 10.0_dp * i + j
         object%level(1)%hh(i, j) = 4000.0_dp + 10.0_dp * i + j
      end do
   end do
   packet = getpacket_wst2d(object, 1, "02", "S")
   call check_close(maxval(abs(reshape(packet, [4]) - expected_s)), 0.0_dp, 2.0e-12_dp, "smooth packet")
   packet = getpacket_wst2d(object, 1, "02", "H")
   call check_close(maxval(abs(reshape(packet, [4]) - expected_h)), 0.0_dp, 2.0e-12_dp, "horizontal packet")
   packet = getpacket_wst2d(object, 1, "02", "V")
   call check_close(maxval(abs(reshape(packet, [4]) - expected_v)), 0.0_dp, 2.0e-12_dp, "vertical packet")
   packet = getpacket_wst2d(object, 1, "02", "D")
   call check_close(maxval(abs(reshape(packet, [4]) - expected_d)), 0.0_dp, 2.0e-12_dp, "diagonal packet")

   replacement = reshape([11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp], [2, 2])
   call putpacket_wst2d(object, 1, "02", "D", replacement, ok)
   call check(ok, "putpacket.wst2D status")
   packet = getpacket_wst2d(object, 1, "02", "D")
   call check_close(maxval(abs(packet - replacement)), 0.0_dp, 0.0_dp, "put/get packet round trip")
   packet = getpacket_wst2d(object, 1, "05", "D")
   call check(size(packet) == 0, "invalid base-four index")
   print *, "test_wst2d_packets: PASS"

contains

   subroutine check(condition, label)
      !! Stops the test program when a logical assertion fails.
      logical, intent(in) :: condition !! Assertion value.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      if (.not. condition) then
         print *, "FAIL: ", trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(value, expected, tolerance, label)
      !! Stops the test program when two scalar values differ beyond tolerance.
      real(dp), intent(in) :: value !! Computed value.
      real(dp), intent(in) :: expected !! Expected value.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      if (abs(value - expected) > tolerance) then
         print *, "FAIL: ", trim(label), value, expected
         error stop 1
      end if
   end subroutine check_close

end program test_wst2d_packets
