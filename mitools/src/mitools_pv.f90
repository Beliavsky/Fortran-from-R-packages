! SPDX-License-Identifier: GPL-2.0-only
! Computational translation of CRAN mitools 2.4 by Thomas Lumley.
! Fortran translation and modifications: 2026-08-30.
module mitools_pv
   use r_kinds, only : dp
   use mitools_types, only : mitools_invalid_index, mitools_invalid_shape, mitools_success
   implicit none
   private

   public :: pv_materialize, pv_select

contains

   pure subroutine pv_select(plausible_values, replicate, selected, status)
      real(dp), intent(in) :: plausible_values(:, :, :) !! Plausible values, shape (n_row, n_variable, n_replicate).
      integer, intent(in) :: replicate !! One-based plausible-value replicate index in 1:n_replicate.
      real(dp), allocatable, intent(out) :: selected(:, :) !! Values for all mapped variables at the selected replicate.
      integer, intent(out) :: status !! Zero on success; nonzero for an invalid replicate or empty replicate dimension.

      if (size(plausible_values, 3) < 1) then
         status = mitools_invalid_shape
         return
      end if
      if (replicate < 1 .or. replicate > size(plausible_values, 3)) then
         status = mitools_invalid_index
         return
      end if
      selected = plausible_values(:, :, replicate)
      status = mitools_success
   end subroutine pv_select

   pure subroutine pv_materialize(data, plausible_values, replicate, materialized, status)
      real(dp), intent(in) :: data(:, :) !! Base numeric data matrix retained for every plausible-value analysis.
      real(dp), intent(in) :: plausible_values(:, :, :) !! Plausible values to append, shape (n_row, n_variable, n_replicate).
      integer, intent(in) :: replicate !! One-based plausible-value replicate index in 1:n_replicate.
      real(dp), allocatable, intent(out) :: materialized(:, :) !! Base data with selected plausible-value columns appended.
      integer, intent(out) :: status !! Zero on success; nonzero for row mismatch or invalid replicate index.
      integer :: nbase
      integer :: npv

      if (size(data, 1) /= size(plausible_values, 1)) then
         status = mitools_invalid_shape
         return
      end if
      if (size(plausible_values, 3) < 1) then
         status = mitools_invalid_shape
         return
      end if
      if (replicate < 1 .or. replicate > size(plausible_values, 3)) then
         status = mitools_invalid_index
         return
      end if

      nbase = size(data, 2)
      npv = size(plausible_values, 2)
      allocate(materialized(size(data, 1), nbase + npv))
      materialized(:, 1:nbase) = data
      materialized(:, nbase + 1:nbase + npv) = plausible_values(:, :, replicate)
      status = mitools_success
   end subroutine pv_materialize

end module mitools_pv
