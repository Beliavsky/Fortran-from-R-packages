! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_types
   use matrix_kinds, only : dp
   implicit none
   private

   type, public :: coo_matrix
      integer :: nrow = 0
      integer :: ncol = 0
      integer, allocatable :: row_ind(:)
      integer, allocatable :: col_ind(:)
      real(dp), allocatable :: values(:)
   contains
      procedure :: nnz => coo_nnz
      procedure :: valid => coo_valid
   end type coo_matrix

   type, public :: sparse_vector
      integer :: n = 0
      integer, allocatable :: index(:)
      real(dp), allocatable :: values(:)
   contains
      procedure :: nnz => svec_nnz
      procedure :: valid => svec_valid
   end type sparse_vector

contains

   integer function coo_nnz(self) result(n)
      class(coo_matrix), intent(in) :: self
      if (allocated(self%values)) then
         n = size(self%values)
      else
         n = 0
      end if
   end function coo_nnz

   logical function coo_valid(self) result(ok)
      class(coo_matrix), intent(in) :: self
      integer :: n
      ok = self%nrow >= 0 .and. self%ncol >= 0
      if (.not. allocated(self%row_ind) .or. .not. allocated(self%col_ind) .or. &
          .not. allocated(self%values)) then
         ok = .false.
         return
      end if
      n = size(self%values)
      if (size(self%row_ind) /= n .or. size(self%col_ind) /= n) then
         ok = .false.
         return
      end if
      if (n > 0) then
         if (any(self%row_ind < 1) .or. any(self%row_ind > self%nrow) .or. &
             any(self%col_ind < 1) .or. any(self%col_ind > self%ncol)) ok = .false.
      end if
   end function coo_valid

   integer function svec_nnz(self) result(nz)
      class(sparse_vector), intent(in) :: self
      if (allocated(self%values)) then
         nz = size(self%values)
      else
         nz = 0
      end if
   end function svec_nnz

   logical function svec_valid(self) result(ok)
      class(sparse_vector), intent(in) :: self
      integer :: nz
      ok = self%n >= 0
      if (.not. allocated(self%index) .or. .not. allocated(self%values)) then
         ok = .false.
         return
      end if
      nz = size(self%values)
      if (size(self%index) /= nz) then
         ok = .false.
         return
      end if
      if (nz > 0) then
         if (any(self%index < 1) .or. any(self%index > self%n)) ok = .false.
      end if
   end function svec_valid

end module matrixextra_types
