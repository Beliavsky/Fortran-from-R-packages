! SPDX-License-Identifier: GPL-3.0-only
program test_distributions
  use blmodel, only : dp, observ_normal, observ_powerexp, observ_student_t
  use test_support, only : assert_true, assert_vector_close
  implicit none

  real(dp) :: points(2, 3), q(2), covariance(2, 2), parameter(1)
  real(dp), allocatable :: density(:)
  integer :: info

  points = reshape([0.0_dp, 0.1_dp, 0.2_dp, -0.2_dp, -0.1_dp, 0.3_dp], [2, 3])
  q = [0.05_dp, -0.05_dp]
  covariance = reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], [2, 2])

  call observ_normal(points, q, covariance, density=density, info=info)
  call assert_true(info == 0, 'normal density status')
  call assert_vector_close(density, [2.24224194063972_dp, 1.66109368480264_dp, &
    0.860994962392189_dp], 2.0e-14_dp, 2.0e-14_dp, 'normal density')

  parameter = 5.0_dp
  call observ_student_t(points, q, covariance, parameter, density, info)
  call assert_true(info == 0, 'student density status')
  call assert_vector_close(density, [3.00216253208603_dp, 1.69037061542154_dp, &
    0.620514366658504_dp], 3.0e-14_dp, 3.0e-14_dp, 'student density')

  parameter = 0.6_dp
  call observ_powerexp(points, q, covariance, parameter, density, info)
  call assert_true(info == 0, 'power exponential density status')
  call assert_vector_close(density, [2.71455970198338_dp, 1.54973407153831_dp, &
    0.658679565119626_dp], 5.0e-14_dp, 5.0e-14_dp, 'power exponential density')

  covariance(2, 2) = -1.0_dp
  call observ_normal(points, q, covariance, density=density, info=info)
  call assert_true(info /= 0, 'non-positive-definite covariance rejected')

  write(*, '(a)') 'test_distributions: PASS'
end program test_distributions
