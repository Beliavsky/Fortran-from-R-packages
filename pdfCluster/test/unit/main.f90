program test_pdfcluster_unit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use pdfcluster, only : dp, adjusted_rand_index, connected_components
   use pdfcluster, only : graph_unidimensional, h_norm, hprop2f, kepdf_estimate, valley_measure
   implicit none

   real(dp), allocatable :: estimate(:)
   real(dp), allocatable :: h(:)
   real(dp), allocatable :: hx(:, :)
   real(dp) :: area
   real(dp) :: expected
   real(dp) :: x1(3, 1)
   real(dp) :: x2(2, 1)
   real(dp) :: eval0(1, 1)
   real(dp) :: fval(3)
   logical, allocatable :: adjacency(:, :)
   logical :: active(4)
   logical :: ok
   integer, allocatable :: comp(:)
   integer :: ncomp

   x1(:, 1) = [-1.0_dp, 0.0_dp, 1.0_dp]
   call h_norm(x1, h)
   expected = (4.0_dp/9.0_dp)**(1.0_dp/5.0_dp)
   call assert_close(h(1), expected, 1.0e-13_dp, "h.norm normal-reference rule")

   x2(:, 1) = [-1.0_dp, 1.0_dp]
   eval0(1, 1) = 0.0_dp
   call kepdf_estimate(x2, eval0, h=[1.0_dp], kernel="gaussian", estimate=estimate, success=ok)
   call assert_true(ok, "Gaussian kepdf succeeds")
   expected = 0.3989423_dp*exp(-0.5_dp)
   call assert_close(estimate(1), expected, 1.0e-13_dp, "Gaussian kepdf value")

   call hprop2f(x2, [1.0_dp], hx, alpha=0.5_dp, kernel="gaussian", success=ok)
   call assert_true(ok, "hprop2f succeeds")
   call assert_close(hx(1, 1), 1.0_dp, 1.0e-13_dp, "symmetric adaptive bandwidth row 1")
   call assert_close(hx(2, 1), 1.0_dp, 1.0e-13_dp, "symmetric adaptive bandwidth row 2")

   fval = [2.0_dp, 1.0_dp, 2.0_dp]
   call valley_measure(fval, area)
   call assert_close(area, 1.0_dp/6.0_dp, 1.0e-13_dp, "single valley area")

   call graph_unidimensional([3.0_dp, 1.0_dp, 2.0_dp], adjacency)
   call assert_true(adjacency(2, 3) .and. adjacency(3, 1), "unidimensional nearest-neighbour chain")
   call assert_true(.not. adjacency(1, 2), "unidimensional graph omits non-neighbour edge")

   deallocate(adjacency)
   allocate(adjacency(4, 4))
   adjacency = .false.
   adjacency(1, 1) = .true.
   adjacency(2, 2) = .true.
   adjacency(3, 3) = .true.
   adjacency(4, 4) = .true.
   adjacency(1, 2) = .true.
   adjacency(2, 1) = .true.
   adjacency(3, 4) = .true.
   adjacency(4, 3) = .true.
   active = .true.
   call connected_components(adjacency, active, comp, ncomp)
   call assert_true(ncomp == 2, "connected-components count")
   call assert_true(comp(1) == comp(2) .and. comp(3) == comp(4), "connected-components membership")
   call assert_true(comp(1) /= comp(3), "connected-components separation")

   call assert_close(adjusted_rand_index([1, 1, 2, 2], [1, 1, 2, 2]), 1.0_dp, 1.0e-13_dp, &
                     "adjusted Rand identical partitions")
   call assert_close(adjusted_rand_index([1, 1, 2, 2], [1, 2, 1, 2]), -0.5_dp, 1.0e-13_dp, &
                     "adjusted Rand crossed partitions")
   call assert_true(ieee_is_nan(adjusted_rand_index([1, 1, 1], [1, 1, 1])), &
                    "adjusted Rand degenerate denominator follows R NaN behavior")

   print *, "All pdfCluster unit tests passed."

contains

   subroutine assert_close(actual, target, tol, message)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: target !! Expected scalar reference value.
      real(dp), intent(in) :: tol !! Maximum permitted absolute error.
      character(len=*), intent(in) :: message !! Assertion description printed on failure.

      if (abs(actual - target) > tol) then
         print *, "FAIL:", trim(message), actual, target
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Boolean condition that must evaluate true.
      character(len=*), intent(in) :: message !! Assertion description printed on failure.

      if (.not. condition) then
         print *, "FAIL:", trim(message)
         error stop 1
      end if
   end subroutine assert_true

end program test_pdfcluster_unit
