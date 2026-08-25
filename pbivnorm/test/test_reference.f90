! SPDX-License-Identifier: GPL-2.0-or-later
program test_reference
   use pbivnorm_mod, only : dp, pbivnorm
   implicit none
   integer :: i
   real(dp), parameter :: tol = 3.0e-12_dp
   real(dp), parameter :: x(8) = [-2.0_dp, -1.0_dp, -0.2_dp, 0.0_dp, 0.3_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   real(dp), parameter :: y(8) = [ 1.5_dp, -0.5_dp,  0.7_dp, 0.0_dp, 1.2_dp, 0.2_dp, 1.0_dp, 2.5_dp]
   real(dp), parameter :: r(8) = [-0.9_dp, -0.5_dp, 0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.95_dp, -0.8_dp]
   real(dp), parameter :: ref(8) = [ &
      0.00246554521850179259_dp, &
      0.0124469043718115208_dp, &
      0.318936433219385362_dp, &
      0.290215311627583072_dp, &
      0.586685782040582637_dp, &
      0.563653435471129782_dp, &
      0.841336147032871073_dp, &
      0.992440436642593737_dp ]

   do i = 1, size(x)
      if (abs(pbivnorm(x(i), y(i), r(i)) - ref(i)) > tol) then
         print '(a,i0,2(1x,es24.15))', 'reference mismatch ', i, pbivnorm(x(i), y(i), r(i)), ref(i)
         error stop 1
      end if
   end do
   print '(a)', 'test_reference: PASS'
end program test_reference
