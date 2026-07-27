! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch
   use rmgarch_kinds
   use rmgarch_math
   use rmgarch_rng
   use rmgarch_distributions
   use rmgarch_fft
   use rmgarch_types
   use rmgarch_univariate
   use rmgarch_dcc
   use rmgarch_fdcc
   use rmgarch_rolling
   use rmgarch_scenario
   use rmgarch_copula
   use rmgarch_cgarch
   use rmgarch_ica
   use rmgarch_gogarch
   use rmgarch_varx
   use rmgarch_diagnostics
   use rmgarch_model
   implicit none
   public
end module rmgarch
