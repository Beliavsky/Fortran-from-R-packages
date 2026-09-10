! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program test_fossa
   use rssa_kinds, only : dp
   use rssa_oblique, only : decompose_ossa, fossa
   use rssa_types, only : rssa_not_supported, rssa_success, ssa_result
   implicit none
   type(ssa_result) :: fit, unsupported, adjusted, filtered_only
   real(dp) :: expected(3, 5), metric(2, 2), rotation(2, 2), transformed(2, 2)
   real(dp), parameter :: tol = 2.0e-11_dp

   call make_fixture(fit)
   expected = matmul(fit%u * spread(fit%sigma, 1, size(fit%u, 1)), transpose(fit%v))

   adjusted = fossa(fit, gamma=2.0_dp)
   call require(adjusted%info == rssa_success, "finite-gamma FOSSA status")
   call require(maxval(abs(adjusted%sigma - 1.0_dp)) < tol, "FOSSA unit oblique sigma")
   call require(maxval(abs(matmul(adjusted%u, transpose(adjusted%v)) - expected)) < tol, &
                "FOSSA preserves selected trajectory matrix")

   metric = reshape([5.0_dp, -4.0_dp, -4.0_dp, 9.0_dp], [2, 2])
   rotation = matmul(transpose(fit%v), adjusted%v)
   transformed = matmul(transpose(rotation), matmul(metric, rotation))
   call require(abs(transformed(1, 2)) < tol .and. abs(transformed(2, 1)) < tol, &
                "FOSSA diagonalizes filter-adjusted metric")
   call require(abs(transformed(1, 1) - 11.47213595499958_dp) < tol, "FOSSA leading metric eigenvalue")
   call require(abs(transformed(2, 2) - 2.527864045000421_dp) < tol, "FOSSA trailing metric eigenvalue")

   filtered_only = fossa(fit)
   call require(filtered_only%info == rssa_success, "default infinite-gamma FOSSA status")
   rotation = matmul(transpose(fit%v), filtered_only%v)
   metric = reshape([1.0_dp, -1.0_dp, -1.0_dp, 2.0_dp], [2, 2])
   transformed = matmul(transpose(rotation), matmul(metric, rotation))
   call require(abs(transformed(1, 1) - 2.618033988749895_dp) < tol, "filter-only leading eigenvalue")
   call require(abs(transformed(2, 2) - 0.3819660112501051_dp) < tol, "filter-only trailing eigenvalue")

   unsupported = decompose_ossa(adjusted)
   call require(unsupported%info == rssa_not_supported, "OSSA continuation remains unsupported like upstream")
   call require(maxval(abs(unsupported%series - adjusted%series)) < tol, "unsupported continuation preserves object")

   print '(a)', 'test_fossa: PASS'
contains
   subroutine make_fixture(object)
      type(ssa_result), intent(out) :: object !! Deterministic two-component SSA fixture with orthonormal U and V.

      allocate(object%series(7), object%trajectory(3, 5), object%sigma(2), object%u(3, 2), object%v(5, 2))
      object%series = [1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp]
      object%trajectory = 0.0_dp
      object%sigma = [4.0_dp, 2.0_dp]
      object%u = 0.0_dp
      object%u(1, 1) = 1.0_dp
      object%u(2, 2) = 1.0_dp
      object%v = 0.0_dp
      object%v(1, 1) = 1.0_dp
      object%v(2, 2) = 1.0_dp
      object%window = 3
      object%info = rssa_success
   end subroutine make_fixture

   subroutine require(condition, message)
      logical, intent(in) :: condition !! Assertion condition that must hold.
      character(len=*), intent(in) :: message !! Failure description printed before termination.

      if (.not. condition) then
         print '(a)', 'test_fossa: FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine require
end program test_fossa
