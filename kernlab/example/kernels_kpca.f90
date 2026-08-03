! SPDX-License-Identifier: GPL-2.0-only
program kernels_kpca
  use kernlab
  implicit none
  real(dp) :: x(6,2)
  type(kernel_spec) :: kernel
  type(kpca_result) :: model
  integer :: i

  x = reshape([ -2.0_dp,-1.0_dp, -1.0_dp,-2.0_dp, -1.5_dp,-1.5_dp, &
                 1.0_dp, 2.0_dp,  2.0_dp, 1.0_dp,  1.5_dp, 1.5_dp ], &
               [6,2], order=[2,1])
  kernel = rbfdot(0.5_dp)
  call kpca(x, kernel, model, features=2)
  print '(a)', 'First two kernel principal-component scores:'
  do i=1,size(x,1)
    print '(i3,2f12.6)', i, model%rotated(i,:)
  end do
end program kernels_kpca
