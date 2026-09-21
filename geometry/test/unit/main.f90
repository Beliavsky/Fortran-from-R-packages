program test_unit
   use geometry
   implicit none

   real(dp), allocatable :: a(:, :), b(:, :), v(:), sorted(:, :), unique_x(:, :), dist(:)
   real(dp) :: xpoly(4), ypoly(4), data(8)
   integer, allocatable :: ord(:)
   integer :: dims(3), idx(2, 3), i
   logical :: ok

   allocate(a(3, 3))
   a = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                1.0_dp, 1.0_dp, 1.0_dp, &
                2.0_dp, 2.0_dp, 2.0_dp], [3, 3], order=[2, 1])
   call cart2sph_points(a, b)
   call assert_close(b(2, 1), 0.25_dp*acos(-1.0_dp), 1.0e-12_dp, "cart2sph theta")
   call assert_close(b(2, 3), sqrt(3.0_dp), 1.0e-12_dp, "cart2sph radius")
   call sph2cart_points(b, a)
   call assert_close(a(3, 1), 2.0_dp, 1.0e-12_dp, "sph2cart")

   xpoly = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
   ypoly = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
   call assert_close(polyarea_vector(xpoly, ypoly), 1.0_dp, 1.0e-12_dp, "polyarea")

   deallocate(a)
   allocate(a(4, 3))
   a = reshape([3.0_dp, 1.0_dp, 2.0_dp, &
                2.0_dp, 2.0_dp, 1.0_dp, &
                3.0_dp, 1.0_dp, 2.0_dp, &
                1.0_dp, 4.0_dp, 0.0_dp], [4, 3], order=[2, 1])
   call matsort_rows(a, sorted)
   call assert_close(sorted(1, 1), 1.0_dp, 0.0_dp, "matsort first")
   call assert_close(sorted(1, 3), 3.0_dp, 0.0_dp, "matsort last")
   call unique_rows(a, .false., unique_x)
   call assert_true(size(unique_x, 1) == 3, "Unique row count")
   call matorder_rows(a, ord)
   call assert_true(ord(1) == 4, "matorder first row")

   data = [(real(i, dp), i=1, 8)]
   dims = [2, 2, 2]
   idx(1, :) = [1, 1, 1]
   idx(2, :) = [2, 2, 2]
   call entry_value_real(data, dims, idx, v, ok)
   call assert_true(ok, "entry.value status")
   call assert_close(v(1), 1.0_dp, 0.0_dp, "entry.value first")
   call assert_close(v(2), 8.0_dp, 0.0_dp, "entry.value last")
   v = [10.0_dp, 20.0_dp]
   call entry_set_real(data, dims, idx, v, ok)
   call assert_true(ok, "entry.value replacement status")
   call assert_close(data(1), 10.0_dp, 0.0_dp, "entry.value replacement")

   deallocate(a)
   allocate(a(3, 2))
   a = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], [3, 2], order=[2, 1])
   call mesh_dcircle(a, 1.0_dp, dist)
   call assert_close(dist(1), -1.0_dp, 1.0e-12_dp, "mesh.dcircle inside")
   call assert_close(dist(3), 1.0_dp, 1.0e-12_dp, "mesh.dcircle outside")
   call mesh_drectangle(a, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, dist)
   call assert_close(dist(1), -1.0_dp, 1.0e-12_dp, "mesh.drectangle inside")
   call assert_close(dist(3), 1.0_dp, 1.0e-12_dp, "mesh.drectangle outside")

   print '(a)', 'All geometry unit tests passed'

contains

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: expected !! Expected scalar reference value.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
      character(*), intent(in) :: message !! Test label printed when the assertion fails.
      if (abs(actual - expected) > tolerance) then
         write(*, '(a,2(1x,es24.16))') trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Boolean condition that must evaluate true.
      character(*), intent(in) :: message !! Test label printed when the assertion fails.
      if (.not. condition) then
         write(*, '(a)') trim(message)
         error stop 1
      end if
   end subroutine assert_true

end program test_unit
