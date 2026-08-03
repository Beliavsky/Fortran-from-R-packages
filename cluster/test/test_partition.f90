! SPDX-License-Identifier: GPL-2.0-or-later
program test_partition
  use cluster, only: dp, partition_result, daisy, daisy_mixed, pam, clara, fanny, &
    variable_numeric, variable_binary_asymmetric, variable_nominal, cluster_success
  implicit none

  real(dp) :: x(6,2), xm(4,3)
  real(dp), allocatable :: d(:, :)
  integer :: types(3), status
  type(partition_result) :: result
  character(len=:), allocatable :: message

  x = reshape([0.0_dp, 0.2_dp, -0.1_dp, 5.0_dp, 5.2_dp, 4.9_dp, &
               0.0_dp, -0.1_dp, 0.2_dp, 5.1_dp, 4.8_dp, 5.2_dp], [6,2])

  call daisy(x, d, status=status, message=message)
  call check(status == cluster_success, 'numeric daisy status')
  call check(maxval(abs(d-transpose(d))) < 1.0e-12_dp, 'daisy symmetry')
  call check(maxval(abs([(d(status,status), status=1,6)])) < 1.0e-12_dp, 'daisy diagonal')

  call pam(x, 2, result)
  call check(result%ok(), 'PAM status')
  call check(all(result%clustering(1:3) == result%clustering(1)), 'PAM first group')
  call check(all(result%clustering(4:6) == result%clustering(4)), 'PAM second group')
  call check(result%clustering(1) /= result%clustering(4), 'PAM separation')

  call clara(x, 2, result, samples=4, sample_size=5, seed=17)
  call check(result%ok(), 'CLARA status')
  call check(result%clustering(1) /= result%clustering(4), 'CLARA separation')

  call fanny(x, 2, result, membership_exponent=2.0_dp)
  call check(result%ok(), 'FANNY status')
  call check(maxval(abs(sum(result%membership, dim=2)-1.0_dp)) < 1.0e-10_dp, &
    'FANNY membership sums')
  call check(result%clustering(1) /= result%clustering(4), 'FANNY separation')

  xm = reshape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
                0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
                1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [4,3])
  types = [variable_numeric, variable_binary_asymmetric, variable_nominal]
  call daisy_mixed(xm, types, d, status=status, message=message)
  call check(status == cluster_success, 'mixed daisy status')
  call check(d(1,2) > 0.0_dp .and. d(1,4) > d(1,2), 'mixed Gower ordering')

  print '(a)', 'test_partition: PASS'

contains

  subroutine check(condition, text)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: text
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//text
      error stop 1
    end if
  end subroutine check

end program test_partition
