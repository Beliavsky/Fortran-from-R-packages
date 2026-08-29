program basic
  use matrixNormal
  implicit none
  real(dp)::m(2,2),u(2,2),v(2,2),a(2,2),lower(2,2),upper(2,2)
  type(probability_result)::pr
  m=0.0_dp
  u=reshape([1.0_dp,0.25_dp,0.25_dp,1.5_dp],[2,2])
  v=reshape([2.0_dp,0.3_dp,0.3_dp,1.0_dp],[2,2])
  a=rmatnorm(m,u,v,1234)
  print '(a,es14.6)', 'log density = ',dmatnorm(a,m,u,v)
  lower=-1.0_dp
  upper=1.0_dp
  pr=pmatnorm(lower,upper,m,u,v)
  print '(a,f10.6,a,es10.2)', 'rectangle probability = ',pr%value,' +/- ',pr%error
end program basic
