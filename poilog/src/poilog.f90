! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
module poilog
   use poilog_kinds, only : dp
   use poilog_distribution, only : dpoilog, dpoilog_vec, dbipoilog, dbipoilog_vec
   use poilog_rng, only : poilog_seed, rpoilog, rbipoilog
   use poilog_mle, only : poilog_fit, bipoilog_fit, poilog_mle_fit, bipoilog_mle_fit
   implicit none
   public
end module poilog
