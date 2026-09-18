program basic_dccmidas
   use dccmidas
   implicit none

   integer, parameter :: nt = 40
   real(dp) :: res(2, nt)
   real(dp) :: dt(2, 2, nt)
   type(dcc_matrices) :: matrices
   integer :: t
   integer :: status

   dt = 0.0_dp
   do t = 1, nt
      res(1, t) = sin(0.17_dp * real(t, dp))
      res(2, t) = 0.5_dp * res(1, t) + cos(0.11_dp * real(t, dp))
      dt(1, 1, t) = 0.01_dp
      dt(2, 2, t) = 0.012_dp
   end do

   call dcc_mat_est([0.03_dp, 0.92_dp], res, dt, matrices, status)
   if (status /= DCCMIDAS_SUCCESS) error stop 'DCC example failed'
   print '(a,f10.6)', 'Final conditional correlation: ', matrices%r(1, 2, nt)
end program basic_dccmidas
