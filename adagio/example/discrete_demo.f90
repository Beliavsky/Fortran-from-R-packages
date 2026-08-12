program discrete_demo
  use adagio
  implicit none
  type(mknapsack_result) :: r
  r = mknapsack([40._dp,60._dp,30._dp,40._dp,20._dp,5._dp], &
                [110._dp,150._dp,70._dp,80._dp,30._dp,5._dp], [85._dp,65._dp])
  print '(a,f8.2)', 'multiple-knapsack value = ', r%value
  print '(a,*(i0,1x))', 'knapsack assignment     = ', r%ksack
end program discrete_demo
