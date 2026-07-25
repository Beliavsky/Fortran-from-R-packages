! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fcopulae
  use fcopulae_kinds, only : dp, i8, pi
  use fcopulae_archimedean
  use fcopulae_elliptical
  use fcopulae_extreme_value
  use fcopulae_empirical
  use fcopulae_utils, only : pseudo_observations, kendall_tau_sample, spearman_rho_sample
  implicit none
  public
end module fcopulae
