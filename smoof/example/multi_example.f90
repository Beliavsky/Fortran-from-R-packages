program multi_example
  use smoof_kinds, only : dp
  use smoof_multi, only : dtlz2
  implicit none
  real(dp)::x(7),f(3)
  x=0.5_dp
  call dtlz2(x,3,f)
  print '(a,3f14.8)','DTLZ2 = ',f
end program multi_example
