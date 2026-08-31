! SPDX-License-Identifier: GPL-2.0-only
! Computational translation of CRAN mitools 2.4 by Thomas Lumley.
! Fortran translation and modifications: 2026-08-30.
module mitools_imputation
   use r_kinds, only : dp
   use mitools_types, only : imputation_list, mitools_invalid_index, mitools_invalid_shape, mitools_success
   implicit none
   private

   public :: imputation_cbind, imputation_dimensions, imputation_get
   public :: imputation_list_from_array, imputation_rbind

contains

   pure subroutine imputation_list_from_array(values, list, status)
      real(dp), intent(in) :: values(:, :, :) !! Numeric imputed data, shape (n_row, n_column, n_imputation).
      type(imputation_list), intent(out) :: list !! Numeric imputation-list container holding a copy of all datasets.
      integer, intent(out) :: status !! Zero on success; nonzero when no imputation datasets are supplied.

      if (size(values, 3) < 1) then
         status = mitools_invalid_shape
         return
      end if
      list%values = values
      status = mitools_success
   end subroutine imputation_list_from_array

   pure subroutine imputation_dimensions(list, nrow, ncol, nimp, status)
      type(imputation_list), intent(in) :: list !! Numeric imputation-list container to query.
      integer, intent(out) :: nrow !! Number of rows in each imputed dataset.
      integer, intent(out) :: ncol !! Number of numeric columns in each imputed dataset.
      integer, intent(out) :: nimp !! Number of imputed datasets.
      integer, intent(out) :: status !! Zero on success; nonzero when the container is uninitialized.

      if (.not. allocated(list%values)) then
         nrow = 0
         ncol = 0
         nimp = 0
         status = mitools_invalid_shape
         return
      end if
      nrow = size(list%values, 1)
      ncol = size(list%values, 2)
      nimp = size(list%values, 3)
      status = mitools_success
   end subroutine imputation_dimensions

   pure subroutine imputation_get(list, imputation, dataset, status)
      type(imputation_list), intent(in) :: list !! Numeric imputation-list container from which one dataset is extracted.
      integer, intent(in) :: imputation !! One-based imputation index in the inclusive range 1:n_imputation.
      real(dp), allocatable, intent(out) :: dataset(:, :) !! Selected imputed dataset, shape (n_row, n_column).
      integer, intent(out) :: status !! Zero on success; nonzero for an uninitialized list or invalid index.

      if (.not. allocated(list%values)) then
         status = mitools_invalid_shape
         return
      end if
      if (imputation < 1 .or. imputation > size(list%values, 3)) then
         status = mitools_invalid_index
         return
      end if
      dataset = list%values(:, :, imputation)
      status = mitools_success
   end subroutine imputation_get

   pure subroutine imputation_rbind(first, second, combined, status)
      type(imputation_list), intent(in) :: first !! First imputation list, supplying the upper rows in every dataset.
      type(imputation_list), intent(in) :: second !! Second imputation list, supplying appended rows in every dataset.
      type(imputation_list), intent(out) :: combined !! Row-bound list with the same columns and imputation count.
      integer, intent(out) :: status !! Zero on success; nonzero when columns or imputation counts differ.
      integer :: first_rows
      integer :: ncol
      integer :: nimp
      integer :: second_rows

      if (.not. allocated(first%values) .or. .not. allocated(second%values)) then
         status = mitools_invalid_shape
         return
      end if
      if (size(first%values, 2) /= size(second%values, 2) .or. &
          size(first%values, 3) /= size(second%values, 3)) then
         status = mitools_invalid_shape
         return
      end if

      first_rows = size(first%values, 1)
      second_rows = size(second%values, 1)
      ncol = size(first%values, 2)
      nimp = size(first%values, 3)
      allocate(combined%values(first_rows + second_rows, ncol, nimp))
      combined%values(1:first_rows, :, :) = first%values
      combined%values(first_rows + 1:first_rows + second_rows, :, :) = second%values
      status = mitools_success
   end subroutine imputation_rbind

   pure subroutine imputation_cbind(first, second, combined, status)
      type(imputation_list), intent(in) :: first !! First imputation list, supplying the left columns in every dataset.
      type(imputation_list), intent(in) :: second !! Second imputation list, supplying appended columns in every dataset.
      type(imputation_list), intent(out) :: combined !! Column-bound list with the same rows and imputation count.
      integer, intent(out) :: status !! Zero on success; nonzero when row or imputation counts differ.
      integer :: first_cols
      integer :: nimp
      integer :: nrow
      integer :: second_cols

      if (.not. allocated(first%values) .or. .not. allocated(second%values)) then
         status = mitools_invalid_shape
         return
      end if
      if (size(first%values, 1) /= size(second%values, 1) .or. &
          size(first%values, 3) /= size(second%values, 3)) then
         status = mitools_invalid_shape
         return
      end if

      nrow = size(first%values, 1)
      first_cols = size(first%values, 2)
      second_cols = size(second%values, 2)
      nimp = size(first%values, 3)
      allocate(combined%values(nrow, first_cols + second_cols, nimp))
      combined%values(:, 1:first_cols, :) = first%values
      combined%values(:, first_cols + 1:first_cols + second_cols, :) = second%values
      status = mitools_success
   end subroutine imputation_cbind

end module mitools_imputation
