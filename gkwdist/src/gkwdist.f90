! SPDX-License-Identifier: MIT
module gkwdist
   use gkwdist_kinds, only : dp
   use gkwdist_core, only : fam_gkw,fam_bkw,fam_kkw,fam_ekw,fam_mc,fam_kw,fam_beta, &
      family_from_name,family_npar,family_name
   use gkwdist_distributions
   use gkwdist_startvalues, only : gkwgetstartvalues, theoretical_moment
   use gkwdist_rng, only : seed_rng
   implicit none
   public
end module gkwdist
