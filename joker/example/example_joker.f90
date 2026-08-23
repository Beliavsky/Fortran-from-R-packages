
program example_joker
  use joker
  implicit none
  real(dp) :: x(5)
  type(estimate2) :: fit
  integer :: i
  do i=1,size(x)
    x(i)=rgamma_j(3.0_dp,2.0_dp)
  end do
  fit=egamma_mle(x)
  print '(a,2f12.6)',"gamma MLE shape, scale: ",fit%p1,fit%p2
end program
