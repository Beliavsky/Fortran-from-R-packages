program test_utils
  use ga
  implicit none
  integer :: b(8), g(8), bb(8)
  real(dp) :: x(3), lo(3), up(3), pm, prob(4), fx(4), hist(5)
  call decimal2binary(173,b)
  if(binary2decimal(b)/=173) error stop "binary conversion"
  call binary2gray(b,g)
  call gray2binary(g,bb)
  if(any(bb/=b)) error stop "gray conversion"
  x=[-2.0_dp,0.5_dp,3.0_dp];lo=0.0_dp;up=1.0_dp
  call repair_solution(x,lo,up)
  if(any(x<lo).or.any(x>up)) error stop "repair solution"
  x=[-4.2_dp,0.5_dp,5.7_dp]
  call reflect_solution(x,lo,up)
  if(any(x<lo).or.any(x>up)) error stop "reflect solution"
  pm=ga_pmutation(1,100,p0=0.5_dp,p=0.01_dp,t=50.0_dp)
  if(abs(pm-0.5_dp)>1.0e-14_dp) error stop "adaptive mutation initial"
  fx=[1.0_dp,4.0_dp,2.0_dp,3.0_dp]
  call optim_probsel(fx,prob)
  if(abs(sum(prob)-1.0_dp)>1.0e-14_dp) error stop "optim probability normalization"
  if(maxloc(prob,dim=1)/=2) error stop "optim probability ordering"
  hist=[1.0_dp,2.0_dp,1.0_dp,2.0_dp,2.0_dp]
  if(garun(hist)/=3) error stop "garun package semantics"
  print *, "test_utils: PASS"
end program test_utils
