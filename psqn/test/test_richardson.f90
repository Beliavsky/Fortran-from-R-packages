! SPDX-License-Identifier: Apache-2.0
program test_richardson
  use psqn, only : dp, richardson_vector_derivative
  implicit none

  real(dp) :: res(2), truth(2), tol
  real(dp), parameter :: x = 1.5_dp

  tol = epsilon(1.0_dp)**(3.0_dp/5.0_dp)
  truth = [2.0_dp * exp(2.0_dp*x), 3.0_dp * cos(3.0_dp*x)]

  call richardson_vector_derivative(fun, x, res, eps=1.0e-4_dp, scale=2.0_dp, tol=tol, order=6)
  if (any(abs(res - truth) > 10.0_dp * abs(truth) * tol)) error stop "Richardson mismatch scale=2"

  call richardson_vector_derivative(fun, x, res, eps=1.0e-4_dp, scale=4.0_dp, tol=tol, order=6)
  if (any(abs(res - truth) > 10.0_dp * abs(truth) * tol)) error stop "Richardson mismatch scale=4"

  print *, "test_richardson: PASS"

contains

  subroutine fun(z, out)
    real(dp), intent(in) :: z
    real(dp), intent(out) :: out(:)
    out(1) = exp(2.0_dp*z)
    out(2) = sin(3.0_dp*z)
  end subroutine fun

end program test_richardson
