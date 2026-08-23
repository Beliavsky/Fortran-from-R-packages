program demo
  use l1pack
  implicit none
  real(dp) :: x(5,2), y(5)
  type(lad_result) :: fit

  x(:,1)=1.0_dp
  x(:,2)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  y=[1.0_dp,3.0_dp,5.0_dp,7.0_dp,100.0_dp]

  call lad_fit_br(x,y,fit)
  print '(a,2f12.6)', 'LAD coefficients: ',fit%coefficients
  print '(a,f12.6)', 'Sum absolute deviations: ',fit%sad
end program demo
