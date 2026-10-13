! SPDX-License-Identifier: GPL-2.0-or-later
program test_oracle
   use diptest_api
   implicit none
   integer :: unit, count_cases, case_number, n, check_hulls, low, high, ng, nl, nfg, nfl, j, ios, status
   integer, allocatable :: indices(:)
   real(dp), allocatable :: x(:), uniforms(:, :)
   real(dp) :: expected_dip, expected_zero, expected_p, actual_p
   character(100) :: label
   type(dip_result) :: fit
   type(dip_test_result) :: test
   logical :: asymptotic

   open(newunit=unit, file='test/data/oracle_cases.txt', status='old', action='read')
   read(unit, *) count_cases
   do case_number = 1, count_cases
      read(unit, '(a)') label
      read(unit, *) n, check_hulls
      read(unit, *) expected_dip, expected_zero, expected_p
      read(unit, *) low, high
      allocate(x(n))
      read(unit, *) x
      call dip_compute(x, fit)
      call require(fit%status == dip_ok, trim(label) // ' status')
      call near(fit%statistic, expected_dip, 3.0e-14_dp, trim(label) // ' statistic')
      call near(dip(x, .true.), expected_zero, 3.0e-14_dp, trim(label) // ' min.is.0')
      call dip_test(x, test)
      call require(test%status == dip_ok, trim(label) // ' test status')
      call near(test%p_value, expected_p, 3.0e-12_dp, trim(label) // ' p-value')
      call require(all(fit%x(2:) >= fit%x(:n - 1)), trim(label) // ' sorted data')
      if (check_hulls == 1) then
         call require(fit%low == low .and. fit%high == high, trim(label) // ' modal endpoints')
         allocate(indices(n))
         read(unit, *) indices
         call require(all(indices == fit%mn), trim(label) // ' minorant predecessors')
         read(unit, *) indices
         call require(all(indices == fit%mj), trim(label) // ' majorant successors')
         deallocate(indices)
         read(unit, *) ng, nl, nfg, nfl
         call check_nodes(fit%gcm, ng, unit, trim(label) // ' modal GCM')
         call check_nodes(fit%lcm, nl, unit, trim(label) // ' modal LCM')
         call check_nodes(fit%full_gcm, nfg, unit, trim(label) // ' full GCM')
         call check_nodes(fit%full_lcm, nfl, unit, trim(label) // ' full LCM')
      end if
      deallocate(x)
   end do
   close(unit)
   print '(a,i0,a)', 'Passed ', count_cases, ' actual-R data/statistic/p-value fixtures.'

   open(newunit=unit, file='test/data/pvalue_cases.txt', status='old', action='read')
   j = 0
   do
      read(unit, *, iostat=ios) n, expected_dip, expected_p
      if (ios < 0) exit
      call require(ios == 0, 'table fixture parse')
      call dip_p_value(expected_dip, n, actual_p, status, asymptotic)
      call require(status == dip_ok, 'table status')
      call near(actual_p, expected_p, 8.0e-12_dp, 'R table interpolation')
      call require(asymptotic .eqv. (n > 72000), 'asymptotic flag')
      j = j + 1
   end do
   close(unit)
   print '(a,i0,a)', 'Passed ', j, ' actual-R table interpolation fixtures.'

   open(newunit=unit, file='test/data/monte_carlo.txt', status='old', action='read')
   read(unit, *) n, j
   allocate(x(n), uniforms(n, j))
   read(unit, *) expected_dip, expected_p
   read(unit, *) x
   read(unit, *) uniforms
   close(unit)
   call dip_test(x, test, uniforms)
   call require(test%status == dip_ok .and. test%simulated .and. test%replicates == j, 'MC metadata')
   call near(test%statistic, expected_dip, 3.0e-14_dp, 'MC observed dip')
   call near(test%p_value, expected_p, 1.0e-15_dp, 'exact caller-uniform MC p-value')
   print '(a)', 'Passed deterministic 257-replicate actual-R Monte Carlo fixture.'

contains

   subroutine require(condition, message)
      logical, intent(in) :: condition !! Assertion that must hold for the current reference case.
      character(*), intent(in) :: message !! Human-readable name of the failing reference operation.
      if (.not. condition) then
         print '(a)', 'FAIL: ' // message
         error stop 1
      end if
   end subroutine require

   subroutine near(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Result from the Fortran translation.
      real(dp), intent(in) :: expected !! Independent stored result from the actual upstream R package.
      real(dp), intent(in) :: tolerance !! Absolute maximum allowed numerical discrepancy.
      character(*), intent(in) :: message !! Human-readable operation name on failure.
      if (.not. (abs(actual - expected) <= tolerance)) then
         print '(a,2es25.16)', 'FAIL: ' // message // ': ', actual, expected
         error stop 1
      end if
   end subroutine near

   subroutine check_nodes(actual, count_nodes, unit, message)
      integer, intent(in) :: actual(:) !! Computed ordered 1-based hull node indices.
      integer, intent(in) :: count_nodes !! Expected hull node count from the R fixture.
      integer, intent(in) :: unit !! Open fixture unit whose next record is the expected hull.
      character(*), intent(in) :: message !! Name of the hull being checked.
      integer :: expected(count_nodes)
      read(unit, *) expected
      call require(size(actual) == count_nodes, message // ' length')
      call require(all(actual == expected), message // ' indices')
   end subroutine check_nodes
end program test_oracle
