! SPDX-License-Identifier: GPL-2.0-or-later
! Public umbrella module for the wavethresh computational translation.
module wavethresh
   use wavethresh_types
   use wavethresh_filters
   use wavethresh_multiwavelet
   use wavethresh_transform_1d
   use wavethresh_transform_nd
   use wavethresh_access
   use wavethresh_stats
   use wavethresh_threshold
   use wavethresh_threshold_extra
   use wavethresh_bayes
   use wavethresh_fourier
   use wavethresh_signals
   use wavethresh_density
   use wavethresh_interval
   use wavethresh_irregular
   use wavethresh_complex_threshold
   use wavethresh_spectrum
   use wavethresh_basis
   use wavethresh_cv
   use wavethresh_utilities
   implicit none
   public
end module wavethresh
