! SPDX-License-Identifier: GPL-2.0-only
program test_boundary_support
  use ks, only: dp, beta_kernel2_pdf, boundary_kde_pdf, copula_density_empirical, &
                convex_hull_2d, support_mask, balloon_kde_2d, kcde_model, fit_kcde, kcde_eval
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(dp) :: x(6,2), eval(3,2), H(2,2), lo(2), hi(2), f(3), c(3), hloc(3)
  real(dp) :: square(5,2), density(4), level, xcdf(4,2), pcdf(3,2), hcdf(2,2), vals(3)
  real(dp), allocatable :: hull(:,:)
  logical :: inside(4)
  type(kcde_model) :: cm

  if(.not.ieee_is_finite(beta_kernel2_pdf(0.05_dp,0.10_dp,0.15_dp))) error stop 'beta kernel finite'
  if(beta_kernel2_pdf(0.05_dp,0.10_dp,0.15_dp)<=0.0_dp) error stop 'beta kernel positive'

  x(1,:)=[0.05_dp,0.10_dp]; x(2,:)=[0.20_dp,0.35_dp]; x(3,:)=[0.40_dp,0.80_dp]
  x(4,:)=[0.65_dp,0.55_dp]; x(5,:)=[0.85_dp,0.25_dp]; x(6,:)=[0.95_dp,0.90_dp]
  eval(1,:)=[0.15_dp,0.20_dp]; eval(2,:)=[0.50_dp,0.50_dp]; eval(3,:)=[0.85_dp,0.75_dp]
  H=0.0_dp; H(1,1)=0.02_dp; H(2,2)=0.025_dp
  lo=0.0_dp; hi=1.0_dp
  f=boundary_kde_pdf(x,H,eval,lo,hi)
  c=copula_density_empirical(x,H,eval)
  if(any(.not.ieee_is_finite(f)).or.any(f<0.0_dp)) error stop 'boundary kde'
  if(any(.not.ieee_is_finite(c)).or.any(c<0.0_dp)) error stop 'copula density'

  square(1,:)=[0.0_dp,0.0_dp]; square(2,:)=[1.0_dp,0.0_dp]; square(3,:)=[1.0_dp,1.0_dp]
  square(4,:)=[0.0_dp,1.0_dp]; square(5,:)=[0.5_dp,0.5_dp]
  call convex_hull_2d(square,hull)
  if(size(hull,1)/=4.or.size(hull,2)/=2) error stop 'convex hull'

  density=[0.5_dp,0.3_dp,0.1_dp,0.1_dp]
  call support_mask(density,0.5_dp,0.8_dp,inside,level)
  if(count(inside)/=2.or.abs(level-0.3_dp)>1.0e-14_dp) error stop 'support mask'

  call balloon_kde_2d(x,H,eval,f,hloc)
  if(any(.not.ieee_is_finite(f)).or.any(f<0.0_dp).or.any(hloc<=0.0_dp)) error stop 'balloon kde'

  xcdf(1,:)=[-1.0_dp,0.2_dp]; xcdf(2,:)=[0.3_dp,-0.5_dp]
  xcdf(3,:)=[1.1_dp,0.8_dp]; xcdf(4,:)=[-0.2_dp,1.0_dp]
  pcdf(1,:)=[0.0_dp,0.0_dp]; pcdf(2,:)=[0.7_dp,0.4_dp]; pcdf(3,:)=[-0.5_dp,1.2_dp]
  hcdf=reshape([0.7_dp,0.2_dp,0.2_dp,0.5_dp],[2,2])
  call fit_kcde(xcdf,cm,H=hcdf)
  call kcde_eval(cm,pcdf,vals)
  if(maxval(abs(vals-[0.19258760166067265_dp,0.3911206905968508_dp, &
                       0.2875329219991539_dp]))>3.0e-5_dp) error stop '2d kcde reference'

  print *, 'test_boundary_support: PASS'
end program test_boundary_support
