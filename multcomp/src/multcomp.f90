! SPDX-License-Identifier: GPL-2.0-only
module multcomp
  use multcomp_kinds
  use multcomp_types
  use multcomp_parm
  use multcomp_contrasts
  use multcomp_glht
  use multcomp_cld
  use multcomp_helpers
  use multcomp_mmm
  use mvtnorm, only : probability_control, genz_bretz, tvpack, miwa
  implicit none
  public
end module multcomp
