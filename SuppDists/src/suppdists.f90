module suppdists
   use suppdists_kinds, only : dp, i8
   use suppdists_stats, only : dist_stats, sample_moments, moments
   use suppdists_inverse_gaussian
   use suppdists_johnson
   use suppdists_ghyper
   use suppdists_rank
   use suppdists_pearson
   use suppdists_max_fratio
   implicit none
   public
end module suppdists
