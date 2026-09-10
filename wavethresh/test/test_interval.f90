program test_interval
   use wavethresh, only : dp, interval_wavelet_t, wd_int, wr_int
   implicit none

   integer, parameter :: reference_indices(15) = [1, 2, 4, 5, 9, 16, 17, 29, 32, 33, 37, 61, 62, 63, 64]
   real(dp), parameter :: reference_plain(15) = [ &
      3.0383927749102151_dp, -1.1824566987315204_dp, 2.1065993902182143_dp, &
      -1.5126473822564019_dp, 0.071623230852758124_dp, -0.44103523454078108_dp, &
      -0.033901870124388778_dp, 0.0068489872566106502_dp, 0.0041723646381482543_dp, &
      0.010424870815724435_dp, -0.00043088930314025969_dp, 0.0001261915301887101_dp, &
      -0.0020605148155520836_dp, -0.010666171285137931_dp, -0.01975662358007834_dp ]
   real(dp), parameter :: reference_preconditioned(15) = [ &
      3.9418876254973654_dp, -1.1833172232131934_dp, 2.7121566694380976_dp, &
      -0.24534258791784325_dp, 0.20323421570597566_dp, -0.34148692557424748_dp, &
      -0.0045300745956600719_dp, 0.0069318741406596551_dp, 0.0012177046668137015_dp, &
      -0.000076522596083011818_dp, -0.00043088930314025969_dp, -0.00013191240430333873_dp, &
      -0.00026372438046361191_dp, 0.0003113260009528062_dp, 0.000074076212791771359_dp ]
   type(interval_wavelet_t) :: object
   real(dp), allocatable :: data(:)
   real(dp), allocatable :: reconstructed(:)
   integer :: order
   integer :: i

   allocate(data(64))
   do i = 1, size(data)
      data(i) = sin(real(i, dp) * 0.17_dp) + real(i, dp) / 100.0_dp
   end do

   object = wd_int(data, 4, 2, .false.)
   call check(object%ok, "wd.int unpreconditioned status")
   call check(all(object%filters_used == [4, 4, 2, 1]), "wd.int adaptive filter history")
   call check_close(maxval(abs(object%transformed(reference_indices) - reference_plain)), 0.0_dp, 2.0e-11_dp, &
      "wd.int upstream unpreconditioned coefficients")
   reconstructed = wr_int(object)
   call check_close(maxval(abs(reconstructed - data)), 0.0_dp, 5.0e-9_dp, "wr.int unpreconditioned reconstruction")

   object = wd_int(data, 4, 2, .true.)
   call check(object%ok, "wd.int preconditioned status")
   call check_close(maxval(abs(object%transformed(reference_indices) - reference_preconditioned)), 0.0_dp, 2.0e-11_dp, &
      "wd.int upstream preconditioned coefficients")
   reconstructed = wr_int(object)
   call check_close(maxval(abs(reconstructed - data)), 0.0_dp, 5.0e-9_dp, "wr.int preconditioned reconstruction")

   deallocate(data)
   allocate(data(256))
   do i = 1, size(data)
      data(i) = cos(real(i, dp) * 0.031_dp) + 0.2_dp * sin(real(i, dp) * 0.19_dp)
   end do
   do order = 1, 8
      object = wd_int(data, order, 3, .true.)
      call check(object%ok, "wd.int all filter orders")
      reconstructed = wr_int(object)
      call check_close(maxval(abs(reconstructed - data)), 0.0_dp, 2.0e-7_dp, "wr.int all filter orders")
   end do

   object = wd_int(data(:250), 4, 2)
   call check(.not. object%ok, "wd.int rejects non-power-of-two data")
   print *, "test_interval: PASS"

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

end program test_interval
