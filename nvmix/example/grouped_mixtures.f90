! SPDX-License-Identifier: GPL-3.0-or-later
program grouped_mixtures
  use nvmix
  implicit none
  type(nvmix_model) :: model
  type(sample_result) :: sample
  real(dp), allocatable :: correlation(:,:)
  real(dp) :: loc(3),scale(3,3)
  integer :: groupings(3),families(2)
  loc=0.0_dp
  scale=reshape([1.0_dp,0.4_dp,0.2_dp,0.4_dp,1.0_dp,0.3_dp,0.2_dp,0.3_dp,1.0_dp],[3,3])
  groupings=[1,1,2]; families=mix_inverse_gamma
  model=make_grouped_model(loc,scale,groupings,families,[6.0_dp,14.0_dp])
  correlation=corgnvmix(model)
  print '(a)'; print '(a)','implied correlation matrix:'
  print '(3f12.6)',transpose(correlation)
  sample=rgnvmix(4,model,2026_i8)
  print '(a)'; print '(a)','sample:'
  print '(3f12.6)',transpose(sample%x)
end program
