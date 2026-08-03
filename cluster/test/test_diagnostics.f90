! SPDX-License-Identifier: GPL-2.0-or-later
module test_gap_callback_mod
  use cluster, only: dp, partition_result, pam, cluster_success
  implicit none
contains
  subroutine pam_callback(x, k, labels, status)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: k
    integer, allocatable, intent(out) :: labels(:)
    integer, intent(out) :: status
    type(partition_result) :: result
    call pam(x, k, result)
    if (result%ok()) then
      labels = result%clustering
      status = cluster_success
    else
      allocate(labels(0))
      status = result%status
    end if
  end subroutine pam_callback
end module test_gap_callback_mod

program test_diagnostics
  use cluster, only: dp, silhouette_result, gap_result, daisy, silhouette, medoids, &
    meanabsdev, size_diss, lower_to_upper_tri_inds, upper_to_lower_tri_inds, &
    clus_gap, max_se, cluster_success
  use test_gap_callback_mod, only: pam_callback
  implicit none

  real(dp) :: x(8,2), objective
  real(dp), allocatable :: d(:, :)
  integer :: labels(8), status
  integer, allocatable :: med(:), a(:), b(:)
  type(silhouette_result) :: sil
  type(gap_result) :: gap
  character(len=:), allocatable :: message

  x = reshape([0.0_dp,0.1_dp,-0.1_dp,0.2_dp,5.0_dp,5.1_dp,4.9_dp,5.2_dp, &
               0.0_dp,-0.1_dp,0.2_dp,0.1_dp,5.0_dp,5.2_dp,4.8_dp,5.1_dp], [8,2])
  labels = [1,1,1,1,2,2,2,2]
  call daisy(x, d, status=status, message=message)
  call check(status == cluster_success, 'daisy status')
  call silhouette(labels, d, sil)
  call check(sil%ok(), 'silhouette status')
  call check(sil%average_width > 0.8_dp, 'silhouette separation')

  call medoids(labels, d, med, objective, status)
  call check(status == cluster_success .and. size(med) == 2, 'medoids')
  call check(objective >= 0.0_dp, 'medoid objective')
  call check(abs(meanabsdev([1.0_dp,2.0_dp,3.0_dp])-2.0_dp/3.0_dp) < 1.0e-12_dp, &
    'meanabsdev')
  call check(size_diss(10) == 5 .and. size_diss(11) == -1, 'sizeDiss')
  call lower_to_upper_tri_inds(5, a)
  call upper_to_lower_tri_inds(5, b)
  call check(size(a) == 10 .and. size(b) == 10, 'triangle index sizes')
  call check(all(a(b) == [(status,status=1,10)]), 'triangle index inverse')

  call clus_gap(x, 3, 4, pam_callback, gap, seed=91)
  call check(gap%ok(), 'gap status')
  call check(size(gap%gap) == 3 .and. gap%selected_k >= 1 .and. gap%selected_k <= 3, &
    'gap dimensions')
  call check(max_se([1.0_dp,2.0_dp,1.5_dp], [0.1_dp,0.1_dp,0.1_dp], 'firstmax') == 2, &
    'maxSE firstmax')

  print '(a)', 'test_diagnostics: PASS'

contains
  subroutine check(condition, text)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: text
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//text
      error stop 1
    end if
  end subroutine check
end program test_diagnostics
