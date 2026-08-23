program example
  use pgnorm
  use pgnorm_special, only: dp
  implicit none
  real(dp) :: x(10)
  call rpgnorm(size(x),x,p=3.0_dp,mean=1.0_dp,scale=2.0_dp)
  print '(10(f9.4,1x))',x
end program example
