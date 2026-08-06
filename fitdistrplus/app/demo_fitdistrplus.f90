program demo_fitdistrplus
  use fitdistrplus
  implicit none
  type(distribution_model) :: dist
  type(fit_result) :: mle, qme
  type(gof_result) :: gof
  type(descriptive_result) :: description
  type(fit_control) :: control
  real(dp), allocatable :: x(:)
  integer :: i

  allocate(x(100))
  call make_weibull(dist)
  do i = 1, size(x)
    x(i) = dist%quantile((real(i,dp)-0.5_dp)/real(size(x),dp), [1.8_dp, 2.5_dp])
  end do

  control%tolerance = 1.0e-9_dp
  call fitdist_auto(x, dist, method_mle, mle, control)
  call qmedist(x, dist, [1.5_dp, 2.0_dp], [0.25_dp, 0.75_dp], qme, control)
  call gofstat(x, dist, mle, gof)
  call descdist(x, description)

  write(*,'(a,2(1x,f10.5))') "Weibull MLE shape/scale:", mle%estimate
  write(*,'(a,2(1x,f10.5))') "Weibull QME shape/scale:", qme%estimate
  write(*,'(a,3(1x,f10.5))') "KS, CvM, AD:", gof%ks, gof%cvm, gof%ad
  write(*,'(a,2(1x,f10.5))') "Sample skewness/kurtosis:", &
    description%skewness, description%kurtosis
end program demo_fitdistrplus
