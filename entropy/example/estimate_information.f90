! SPDX-License-Identifier: GPL-3.0-or-later
program estimate_information
   use entropy_module, only: dp, entropy, freqs_shrink, shrinkage_intensity, mi_empirical, &
                             gstatindep, test_result, discretize
   implicit none
   real(dp) :: counts(4), joint(2, 3)
   integer, allocatable :: bins(:)
   integer :: status
   type(test_result) :: independence
   counts = [8.0_dp, 3.0_dp, 1.0_dp, 0.0_dp]
   print *, 'Empirical entropy (bits):', entropy(counts, unit='log2')
   print *, 'Shrinkage entropy (bits):', entropy(counts, method='shrink', unit='log2')
   print *, 'Estimated shrinkage:', shrinkage_intensity(counts)
   print *, 'Shrunk probabilities:', freqs_shrink(counts)
   joint = reshape([5.0_dp, 1.0_dp, 2.0_dp, 6.0_dp, 3.0_dp, 4.0_dp], [2, 3])
   print *, 'Mutual information (bits):', mi_empirical(joint, 'log2')
   independence = gstatindep(joint)
   print *, 'G statistic, df, upper-tail p:', independence%stat, independence%df, independence%pval
   call discretize([0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp], 4, bins, status=status)
   if (status /= 0) error stop 'Invalid histogram range'
   print *, 'Right-closed histogram:', bins
end program estimate_information
