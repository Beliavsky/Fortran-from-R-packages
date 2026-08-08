program test_stats
  use survival
  implicit none
  real(dp)::time(6),risk(6)
  integer::status(6),group(6)
  type(concordance_result)::cc
  type(survdiff_result)::sd
  time=[1._dp,2._dp,3._dp,4._dp,5._dp,6._dp]
  status=[1,1,1,1,0,0]
  group=[1,1,1,2,2,2]
  risk=[6._dp,5._dp,4._dp,3._dp,2._dp,1._dp]
  call concordance_right(time,status,risk,cc)
  if(cc%cindex<0.99_dp) error stop 'concordance'
  call survdiff(time,status,group,2,sd)
  if(sd%chisq<=0._dp) error stop 'survdiff'
  print *, 'test_stats PASS', cc%cindex,sd%chisq
end program
