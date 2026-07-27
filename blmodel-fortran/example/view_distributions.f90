! SPDX-License-Identifier: GPL-3.0-only
program view_distributions
  use blmodel, only : dp, observ_normal, observ_powerexp, observ_student_t
  implicit none

  real(dp) :: points(2, 3), q(2), covariance(2, 2), params(1)
  real(dp), allocatable :: density(:)
  integer :: info

  points = reshape([0.0_dp, 0.1_dp, 0.2_dp, -0.2_dp, -0.1_dp, 0.3_dp], [2, 3])
  q = [0.05_dp, -0.05_dp]
  covariance = reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], [2, 2])

  call observ_normal(points, q, covariance, density=density, info=info)
  write(*, '(a,3(1x,f12.8))') 'normal:', density
  params = 5.0_dp
  call observ_student_t(points, q, covariance, params, density, info)
  write(*, '(a,3(1x,f12.8))') 'student t:', density
  params = 0.6_dp
  call observ_powerexp(points, q, covariance, params, density, info)
  write(*, '(a,3(1x,f12.8))') 'power exponential:', density
end program view_distributions
