! SPDX-License-Identifier: GPL-2.0-or-later
module splines
   use splines_kinds, only : dp
   use splines_core, only : b_spline_t, poly_spline_t, spline_design, &
      spline_basis_nonzero, linear_interp, fit_interpolating_spline, &
      fit_periodic_spline, to_polynomial_spline, inverse_monotone_spline, &
      type7_quantile
   use splines_basis, only : bs_basis, natural_spline_basis
   implicit none
   public
end module splines
