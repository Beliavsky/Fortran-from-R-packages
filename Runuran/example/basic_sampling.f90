! SPDX-License-Identifier: GPL-2.0-or-later
program basic_sampling
  use runuran
  implicit none
  type(rng_state)::rng
  type(continuous_distribution)::d
  type(unuran_generator)::g
  real(dp)::x(10)
  call rng_seed(rng,20260816_i8)
  d=udgamma(2.5_dp,1.2_dp)
  g=tdr_new(d)
  call g%sample_n(rng,x)
  print '(a,10f10.4)','gamma sample: ',x
end program
