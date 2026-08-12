program test_discrete_lp
  use adagio
  implicit none
  type(assignment_result) :: ar
  type(change_result) :: cr
  type(setcover_result) :: sr
  type(mknapsack_result) :: mr
  type(knapsack_result) :: kr
  type(subset_result) :: ss
  type(binpack_result) :: br
  real(dp) :: c(3,3)
  integer :: sets(4,4), sv(8)

  c = reshape([4._dp,2._dp,3._dp, 1._dp,0._dp,2._dp, 3._dp,5._dp,2._dp], [3,3])
  ar = assignment(c)
  call check(ar%err == 0, 'assignment status')
  call check(abs(ar%value - 5._dp) < 1e-10_dp, 'assignment objective')
  call check(all(ar%perm == [2,1,3]), 'assignment permutation')

  ar = assignment(c, maximize=.true.)
  call check(abs(ar%value - 11._dp) < 1e-10_dp, 'max assignment objective')

  cr = change_making([2,5,10,50,100], 999)
  call check(cr%feasible, 'change feasible')
  call check(cr%count == 17, 'change count')
  call check(sum(cr%solution * [2,5,10,50,100]) == 999, 'change reconstruction')

  sets = reshape([1,0,0,1, 1,1,0,0, 0,1,1,0, 0,0,1,1], [4,4])
  sr = setcover(sets)
  call check(sr%feasible, 'set cover feasible')
  call check(abs(sr%objective - 2._dp) < 1e-10_dp, 'set cover objective')
  call check(all(sum(sets(sr%sets,:), dim=1) >= 1), 'set cover coverage')

  mr = mknapsack([40._dp,60._dp,30._dp,40._dp,20._dp,5._dp], &
                 [110._dp,150._dp,70._dp,80._dp,30._dp,5._dp], [85._dp,65._dp])
  call check(abs(mr%value - 345._dp) < 1e-8_dp, 'multiple knapsack objective')
  call check(sum(pack([40._dp,60._dp,30._dp,40._dp,20._dp,5._dp], mr%ksack==1)) <= 85._dp+1e-10_dp, &
             'multiple knapsack capacity 1')
  call check(sum(pack([40._dp,60._dp,30._dp,40._dp,20._dp,5._dp], mr%ksack==2)) <= 65._dp+1e-10_dp, &
             'multiple knapsack capacity 2')

  kr = knapsack([2,20,20,30,40,30,60,10], &
                [15._dp,100._dp,90._dp,60._dp,40._dp,15._dp,10._dp,1._dp], 102)
  call check(abs(kr%profit - 280._dp) < 1e-10_dp, 'knapsack objective')
  call check(kr%capacity == 102, 'knapsack capacity')

  sv = [15,13,11,9,8,7,6,5]
  ss = subsetsum(sv, 31, 'dynamic')
  call check(ss%val == 31 .and. ss%found, 'subset sum dynamic')
  call check(sum(sv(ss%inds)) == 31, 'subset sum reconstruction')
  call check(sss_test(sv, 31) == 31, 'subset sum test')

  br = bpp_approx([100._dp,99._dp,89._dp,88._dp,87._dp,75._dp,67._dp,65._dp, &
                   65._dp,57._dp,57._dp,49._dp,47._dp,31._dp,27._dp,18._dp, &
                   13._dp,9._dp,8._dp,1._dp], 100._dp, 'firstfit')
  call check(br%nbins == 12, 'bin packing first fit')

  print *, 'test_discrete_lp: PASS'
contains
  subroutine check(ok, msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if (.not. ok) then
       print *, 'FAIL: ', trim(msg)
       error stop 1
    end if
  end subroutine check
end program test_discrete_lp
