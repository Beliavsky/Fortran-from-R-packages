! SPDX-License-Identifier: GPL-3.0-or-later
program nvmix_demo
  use nvmix
  implicit none
  type(nvmix_model) :: model
  type(sample_result) :: sample
  real(dp) :: loc(2),scale(2,2)
  loc=0.0_dp
  scale=reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
  model=make_nvmix_model(loc,scale,mix_inverse_gamma,7.0_dp)
  print '(a,f12.8)','density at the origin: ',dnvmix([0.0_dp,0.0_dp],model)
  sample=rnvmix(5,model,12345_i8)
  print '(a)'; print '(a)','five Student-t observations:'
  print '(2f14.6)',transpose(sample%x)
end program
