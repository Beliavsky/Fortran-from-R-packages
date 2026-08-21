! locfit-fortran; GPL-2-or-later. See LICENSE_NOTICE and upstream/locfit-R.
! Public convenience facade for locfit-fortran.
module locfit
  use locfit_kinds
  use locfit_constants
  use locfit_kernels
  use locfit_math
  use locfit_families
  use locfit_basis
  use locfit_core
  use locfit_density
  use locfit_bandwidth
  use locfit_interpolation
  use locfit_diagnostics
  use locfit_robust
  implicit none
  public
end module locfit
