! SPDX-License-Identifier: GPL-3.0-or-later
program information_demo
   use infotheo, only : dp, entropy, mutinformation_matrix, condinformation, &
      multiinformation, interinformation, natstobits, discretize, infotheo_ok
   implicit none
   integer :: x(4, 3), status
   integer, allocatable :: bins(:, :)
   real(dp), allocatable :: mi(:, :)
   real(dp) :: measurements(4, 2)

   x(:, 1) = [0, 0, 1, 1]
   x(:, 2) = [0, 1, 0, 1]
   x(:, 3) = [0, 1, 1, 0]
   print '(a,f8.4)', 'XOR joint entropy (bits): ', natstobits(entropy(x))
   print '(a,f8.4)', 'Total correlation (bits): ', natstobits(multiinformation(x))
   print '(a,f8.4)', 'Co-information (bits): ', natstobits(interinformation(x))
   print '(a,f8.4)', 'I(X;Y|XOR) (bits): ', natstobits(condinformation(x(:, 1:1), x(:, 2:2), x(:, 3:3)))
   mi = mutinformation_matrix(x)
   print '(a,3f8.4)', 'Pairwise MI first row (bits): ', natstobits(mi(1, :))
   measurements(:, 1) = [0.1_dp, 0.2_dp, 0.8_dp, 0.9_dp]
   measurements(:, 2) = [10.0_dp, 30.0_dp, 20.0_dp, 40.0_dp]
   call discretize(measurements, bins, 'equalfreq', 2, status)
   if (status /= infotheo_ok) error stop 'Discretization failed'
   print '(a,f8.4)', 'Discretized joint entropy (nats): ', entropy(bins)
end program information_demo
