! Public facade for the modern Fortran translation of the R package ordinal.
! Copyright (C) 2011-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal
   use ordinal_kinds, only : dp
   use ordinal_links
   use ordinal_thresholds
   use ordinal_clm
   use ordinal_clmm
   use ordinal_clmm_general
   use ordinal_profile
   use ordinal_rank
   implicit none
   public
end module ordinal
