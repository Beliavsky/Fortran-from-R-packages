program change_points_example
  use ecp_api, only: dp, e_cp3o, cp_result
  implicit none

  real(dp) :: x(12, 1)
  type(cp_result) :: fit

  x(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 4.0_dp, &
    4.0_dp, 4.0_dp, -2.0_dp, -2.0_dp, -2.0_dp, -2.0_dp]
  fit = e_cp3o(x, k=2, minsize=3, alpha=1.0_dp)

  write (*, '(a,*(i0,1x))') 'Two-change-point path: ', fit%cp_loc(2, 1:2)
end program change_points_example
