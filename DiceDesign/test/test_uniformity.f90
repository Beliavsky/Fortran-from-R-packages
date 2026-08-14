program test_uniformity
  use dicedesign, only : dp, unif_test_statistic, unif_test_quantile, rss2d_result, rss3d_result, rss2d, rss3d, &
    sumof2uniforms_cdf, sumof3uniforms_cdf
  implicit none
  real(dp) :: x(4), d2(6,2), d3(6,3)
  type(rss2d_result) :: r2
  type(rss3d_result) :: r3

  x=[0.1_dp,0.4_dp,0.6_dp,0.9_dp]
  call assert_close(unif_test_statistic(x,'greenwood'),0.24_dp,1e-14_dp,'greenwood')
  call assert_close(unif_test_statistic(x,'qm'),0.42_dp,1e-14_dp,'qm')
  call assert_close(unif_test_statistic(x,'ks'),0.15_dp,1e-14_dp,'ks')
  call assert_close(unif_test_statistic(x,'V'),0.30_dp,1e-14_dp,'V')
  call assert_close(unif_test_statistic(x,'cvm'),0.023333333333333334_dp,1e-14_dp,'cvm')
  call assert_close(unif_test_quantile('greenwood',10,0.05_dp),0.2404_dp,1e-14_dp,'greenwood q')
  call assert_close(unif_test_quantile('ks',10,0.05_dp),0.409398349356582_dp,1e-14_dp,'ks q')
  if (unif_test_statistic(x,'ks',transform_spacings=.true.) < 0.0_dp) error stop 'spacing transform'
  call assert_close(sumof2uniforms_cdf(0.2_dp,0.6_dp,0.8_dp), &
    0.62500000000000011_dp,2e-15_dp,'2-uniform CDF')
  call assert_close(sumof3uniforms_cdf(0.3_dp,0.4_dp,0.5_dp,0.7_dp), &
    0.69568452380952372_dp,2e-15_dp,'3-uniform CDF')

  d2(1,:)=[0.1_dp,0.2_dp]
  d2(2,:)=[0.3_dp,0.8_dp]
  d2(3,:)=[0.7_dp,0.4_dp]
  d2(4,:)=[0.9_dp,0.9_dp]
  d2(5,:)=[0.2_dp,0.5_dp]
  d2(6,:)=[0.6_dp,0.1_dp]
  call rss2d(d2,[0.0_dp,0.0_dp],[1.0_dp,1.0_dp],r2,gof_test_type='ks',n_angle=24)
  if (r2%worst_case(1)/=1 .or. r2%worst_case(2)/=2) error stop 'rss2d worst pair'
  call assert_close(sum(r2%worst_dir**2),1.0_dp,1e-13_dp,'rss2d direction norm')
  if (size(r2%stat)/=24) error stop 'rss2d stat length'
  call assert_close(maxval(r2%stat),0.28038475772933680_dp,2e-14_dp,'rss2d max')
  call assert_close(r2%worst_dir(1),-0.96592582628906820_dp,2e-14_dp,'rss2d dir x')
  call assert_close(r2%worst_dir(2),0.25881904510252102_dp,2e-14_dp,'rss2d dir y')

  d3(:,1)=d2(:,1)
  d3(:,2)=d2(:,2)
  d3(:,3)=[0.4_dp,0.7_dp,0.2_dp,0.8_dp,0.95_dp,0.05_dp]
  call rss3d(d3,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],r3,gof_test_type='ks',n_angle=12)
  if (any(r3%worst_case/=[1,2,3])) error stop 'rss3d worst triplet'
  call assert_close(sum(r3%worst_dir**2),1.0_dp,1e-13_dp,'rss3d direction norm')
  if (any(shape(r3%stat)/=[12,12])) error stop 'rss3d stat shape'
  call assert_close(maxval(r3%stat),0.37471408890160707_dp,2e-14_dp,'rss3d max')
  call assert_close(r3%worst_dir(1),-0.34008636888545951_dp,2e-14_dp,'rss3d dir x')
  call assert_close(r3%worst_dir(2),0.58904686987122756_dp,2e-14_dp,'rss3d dir y')
  call assert_close(r3%worst_dir(3),-0.73305187182982634_dp,2e-14_dp,'rss3d dir z')

  print *, 'test_uniformity: PASS'
contains
  subroutine assert_close(a,b,tol,name)
    real(dp), intent(in) :: a,b,tol
    character(len=*), intent(in) :: name
    if (abs(a-b)>tol) then
      print *, trim(name),a,b,abs(a-b)
      error stop 'assert_close failed'
    end if
  end subroutine assert_close
end program test_uniformity
