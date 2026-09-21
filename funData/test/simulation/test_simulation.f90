program test_simulation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use funData_api
   implicit none

   type(fun_data) :: basis
   type(fun_data) :: sparse
   type(fun_data) :: noisy
   type(sim_fun_result) :: sim
   type(sim_multi_result) :: msim
   type(real_vector) :: grids(2)
   type(basis_spec) :: specs(2)
   real(dp), allocatable :: vals(:)
   real(dp), allocatable :: norms(:)
   integer, allocatable :: seed(:)
   integer :: seed_size
   logical :: ok

   call eigenvalues(4, "linear", vals, ok)
   call assert_true(ok, "linear eigenvalues")
   call assert_close_vec(vals, [1.0_dp, 0.75_dp, 0.5_dp, 0.25_dp], 1.0e-12_dp, "linear eigenvalues values")
   call eigenvalues(3, "exponential", vals, ok)
   call assert_true(ok, "exponential eigenvalues")
   call assert_close(vals(2), exp(-0.5_dp), 1.0e-12_dp, "exponential eigenvalue")

   call eigenfunctions(linspace(0.0_dp, 1.0_dp, 101), 4, "Fourier", basis, ok = ok)
   call assert_true(ok, "Fourier basis")
   call norm_fun_data(basis, norms, ok = ok)
   call assert_true(ok, "Fourier norms")
   call assert_close_vec(norms, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 1.0e-10_dp, "Fourier orthonormal norms")

   call eigenfunctions(linspace(0.0_dp, 1.0_dp, 101), 4, "FourierLin", basis, ok = ok)
   call assert_true(ok, "FourierLin basis")
   call norm_fun_data(basis, norms, ok = ok)
   call assert_close(norms(4), 1.00153_dp, 3.0e-5_dp, "FourierLin discrete norm")

   call random_seed(size = seed_size)
   allocate(seed(seed_size))
   seed = 20260920
   call random_seed(put = seed)
   grids(1)%values = linspace(0.0_dp, 1.0_dp, 21)
   call sim_fun_data(grids(1:1), [4], [character(len=16) :: "Fourier"], "linear", 3, sim, ok = ok)
   call assert_true(ok, "simFunData 1D")
   call assert_true(all(sim%sim_data%dims == [3, 21]), "simFunData dimensions")
   call assert_true(nobs_fun_data(sim%true_funs) == 4, "simFunData true basis count")

   call sparsify_fun_data(sim%sim_data, 5, 5, sparse, ok)
   call assert_true(ok, "sparsify")
   call assert_true(count(.not. ieee_is_nan(sparse%x)) == 15, "sparsify retained values")
   call add_error_fun_data(sim%sim_data, 0.25_dp, noisy, ok)
   call assert_true(ok, "addError")
   call assert_true(maxval(abs(noisy%x - sim%sim_data%x)) > 0.0_dp, "addError changes data")

   grids(1)%values = linspace(0.0_dp, 1.0_dp, 11)
   grids(2)%values = linspace(-1.0_dp, 1.0_dp, 9)
   call sim_fun_data(grids, [2, 3], [character(len=16) :: "Poly", "Fourier"], "linear", 2, sim, ok = ok)
   call assert_true(ok, "simFunData 2D")
   call assert_true(all(sim%true_funs%dims == [6, 11, 9]), "2D true basis dimensions")
   call assert_true(all(sim%sim_data%dims == [2, 11, 9]), "2D simulated dimensions")

   allocate(specs(1)%argvals(1), specs(2)%argvals(1))
   specs(1)%argvals(1)%values = linspace(0.0_dp, 1.0_dp, 13)
   specs(2)%argvals(1)%values = linspace(-0.5_dp, 0.5_dp, 9)
   allocate(specs(1)%m(1), specs(2)%m(1))
   specs(1)%m = [3]
   specs(2)%m = [3]
   allocate(specs(1)%basis_type(1), specs(2)%basis_type(1))
   specs(1)%basis_type = [character(len=16) :: "Poly"]
   specs(2)%basis_type = [character(len=16) :: "Fourier"]
   seed = 20260920
   call random_seed(put = seed)
   call sim_multi_fun_data("weighted", specs, "linear", 4, msim, ok = ok)
   call assert_true(ok, "simMultiFunData weighted")
   call assert_true(size(msim%sim_data%components) == 2, "simMulti component count")
   call assert_true(nobs_fun_data(msim%sim_data%components(1)) == 4, "simMulti observation count")
   call assert_true(nobs_fun_data(msim%true_funs%components(1)) == 3, "simMulti basis count")

   print '(a)', "All funData simulation tests passed"

contains

   pure function linspace(a, b, n) result(x)
      real(dp), intent(in) :: a !! Lower endpoint included in the generated regular grid.
      real(dp), intent(in) :: b !! Upper endpoint included in the generated regular grid.
      integer, intent(in) :: n !! Number of grid points, required to be at least two.
      real(dp), allocatable :: x(:)
      integer :: i

      allocate(x(n))
      do i = 1, n
         x(i) = a + (b - a) * real(i - 1, dp) / real(n - 1, dp)
      end do
   end function linspace

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must evaluate true for the test to continue.
      character(len=*), intent(in) :: message !! Human-readable label printed when the assertion fails.
      if (.not. condition) then
         print '(a)', "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tol, message)
      real(dp), intent(in) :: actual !! Computed scalar result being checked.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tol !! Maximum accepted absolute difference.
      character(len=*), intent(in) :: message !! Human-readable assertion label.
      if (abs(actual - expected) > tol) then
         print '(a,2es24.15)', "FAIL: " // trim(message) // " actual/expected ", actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_close_vec(actual, expected, tol, message)
      real(dp), intent(in) :: actual(:) !! Computed vector result being checked.
      real(dp), intent(in) :: expected(:) !! Reference vector with the same shape as actual.
      real(dp), intent(in) :: tol !! Maximum accepted absolute elementwise difference.
      character(len=*), intent(in) :: message !! Human-readable assertion label.
      if (size(actual) /= size(expected)) then
         print '(a)', "FAIL size: " // trim(message)
         error stop 1
      end if
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tol) then
            print '(a)', "FAIL values: " // trim(message)
            error stop 1
         end if
      end if
   end subroutine assert_close_vec

end program test_simulation
