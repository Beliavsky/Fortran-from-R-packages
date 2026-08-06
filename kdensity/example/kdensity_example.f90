program kdensity_example
  use kdensity
  implicit none
  real(dp)::x(12)=[0.05_dp,0.08_dp,0.12_dp,0.16_dp,0.20_dp,0.25_dp,0.31_dp,0.38_dp,0.47_dp,0.58_dp,0.70_dp,0.82_dp]
  real(dp)::grid(5)=[0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp]
  type(kdensity_options)::opt
  type(kdensity_fit)::fit
  integer::i
  opt%kernel='beta';opt%start='beta';opt%bandwidth='HS';opt%support=[0.0_dp,1.0_dp];opt%support_supplied=.true.
  fit=fit_kdensity(x,opt)
  if(fit%status/=0)error stop trim(fit%message)
  print '(a,f10.6)', 'bandwidth = ',fit%bw
  print '(a)', ' x        density'
  do i=1,size(grid);print '(f6.2,2x,f12.7)',grid(i),fit%pdf(grid(i));enddo
end program
