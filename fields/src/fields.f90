! SPDX-License-Identifier: GPL-2.0-or-later
! Public umbrella module for the modern Fortran translation of fields.
module fields
use fields_kinds
use fields_distance
use fields_covariance
use fields_polynomial
use fields_spline1d
use fields_kriging
use fields_tps
use fields_sparse_kriging
use fields_fast_tps
use fields_stats
use fields_variogram
use fields_grid
use fields_geometry
use fields_fft_grid
use fields_spatial_process
use fields_simulation
use fields_quantile
use spam_types, only: csr_matrix, spam_chol
implicit none
public
end module fields
