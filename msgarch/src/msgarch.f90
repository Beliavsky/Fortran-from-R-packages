! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch
   use msgarch_kinds
   use msgarch_rng
   use msgarch_special
   use msgarch_distributions
   use msgarch_types
   use msgarch_models
   use msgarch_filter
   use msgarch_simulation
   use msgarch_mapping
   use msgarch_estimation
   use msgarch_forecast
   use msgarch_risk
   use msgarch_hmm
   use msgarch_posterior
   implicit none
   public
end module msgarch
