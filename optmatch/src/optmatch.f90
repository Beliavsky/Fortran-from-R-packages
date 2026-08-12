! SPDX-License-Identifier: MIT
module optmatch
   use optmatch_kinds, only : dp, optmatch_inf
   use optmatch_types
   use optmatch_stats, only : mean_value, sample_variance, median_value, mad_value
   use optmatch_distance
   use optmatch_matching
   use optmatch_utilities
   use optmatch_feasibility
   implicit none
   public
end module optmatch
