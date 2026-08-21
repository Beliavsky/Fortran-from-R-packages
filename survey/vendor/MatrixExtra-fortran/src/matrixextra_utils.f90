! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_utils
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix, csr_from_triplet
   use matrixextra_types, only : coo_matrix, sparse_vector
   use matrixextra_conversions, only : coo_from_csr, csr_from_coo
   implicit none
   private
   public :: csr_remove_zeros, coo_remove_zeros, sparse_vector_remove_zeros
   public :: csr_sort_indices, coo_sort_indices, sparse_vector_sort_indices
   public :: csr_check, coo_check, sparse_vector_check
   public :: csr_filter, csr_map, empty_sparse

   abstract interface
      logical function sparse_predicate(i,j,x)
         import dp
         integer, intent(in) :: i,j
         real(dp), intent(in) :: x
      end function sparse_predicate
      real(dp) function sparse_map_fun(i,j,x)
         import dp
         integer, intent(in) :: i,j
         real(dp), intent(in) :: x
      end function sparse_map_fun
   end interface

contains

   subroutine csr_remove_zeros(a,tol)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in), optional :: tol
      type(coo_matrix) :: c
      real(dp) :: eps
      integer :: n, p, info
      logical, allocatable :: keep(:)
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp,tol)
      call coo_from_csr(a,c)
      n = c%nnz()
      allocate(keep(n))
      keep = abs(c%values) > eps
      p = count(keep)
      c%row_ind = pack(c%row_ind,keep)
      c%col_ind = pack(c%col_ind,keep)
      c%values = pack(c%values,keep)
      call csr_from_coo(c,a,info)
   end subroutine csr_remove_zeros

   subroutine coo_remove_zeros(a,tol)
      type(coo_matrix), intent(inout) :: a
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      logical, allocatable :: keep(:)
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp,tol)
      allocate(keep(a%nnz()))
      keep = abs(a%values) > eps
      a%row_ind = pack(a%row_ind,keep)
      a%col_ind = pack(a%col_ind,keep)
      a%values = pack(a%values,keep)
   end subroutine coo_remove_zeros

   subroutine sparse_vector_remove_zeros(v,tol)
      type(sparse_vector), intent(inout) :: v
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      logical, allocatable :: keep(:)
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp,tol)
      allocate(keep(v%nnz()))
      keep = abs(v%values) > eps
      v%index = pack(v%index,keep)
      v%values = pack(v%values,keep)
   end subroutine sparse_vector_remove_zeros

   subroutine csr_sort_indices(a)
      type(csr_matrix), intent(inout) :: a
      type(coo_matrix) :: c
      integer :: info
      call coo_from_csr(a,c)
      call csr_from_coo(c,a,info)
   end subroutine csr_sort_indices

   subroutine coo_sort_indices(a)
      type(coo_matrix), intent(inout) :: a
      integer :: gap, k, p, tr, tc
      real(dp) :: tv
      gap = a%nnz()/2
      do while (gap > 0)
         do k = gap+1, a%nnz()
            tr=a%row_ind(k); tc=a%col_ind(k); tv=a%values(k)
            p=k
            do while (p > gap)
               if (a%row_ind(p-gap) < tr) exit
               if (a%row_ind(p-gap) == tr .and. a%col_ind(p-gap) <= tc) exit
               a%row_ind(p)=a%row_ind(p-gap)
               a%col_ind(p)=a%col_ind(p-gap)
               a%values(p)=a%values(p-gap)
               p=p-gap
            end do
            a%row_ind(p)=tr; a%col_ind(p)=tc; a%values(p)=tv
         end do
         gap=gap/2
      end do
   end subroutine coo_sort_indices

   subroutine sparse_vector_sort_indices(v)
      type(sparse_vector), intent(inout) :: v
      integer :: gap,k,p,ti
      real(dp) :: tv
      gap=v%nnz()/2
      do while (gap>0)
         do k=gap+1,v%nnz()
            ti=v%index(k); tv=v%values(k); p=k
            do while (p>gap .and. v%index(p-gap)>ti)
               v%index(p)=v%index(p-gap); v%values(p)=v%values(p-gap); p=p-gap
            end do
            v%index(p)=ti; v%values(p)=tv
         end do
         gap=gap/2
      end do
   end subroutine sparse_vector_sort_indices

   logical function csr_check(a) result(ok)
      type(csr_matrix), intent(in) :: a
      ok = a%valid()
   end function csr_check

   logical function coo_check(a) result(ok)
      type(coo_matrix), intent(in) :: a
      ok = a%valid()
   end function coo_check

   logical function sparse_vector_check(v) result(ok)
      type(sparse_vector), intent(in) :: v
      ok = v%valid()
   end function sparse_vector_check

   subroutine csr_filter(a,predicate,b)
      type(csr_matrix), intent(in) :: a
      procedure(sparse_predicate) :: predicate
      type(csr_matrix), intent(out) :: b
      type(coo_matrix) :: c
      integer :: k, info
      logical, allocatable :: keep(:)
      call coo_from_csr(a,c)
      allocate(keep(c%nnz()))
      do k=1,c%nnz()
         keep(k)=predicate(c%row_ind(k),c%col_ind(k),c%values(k))
      end do
      c%row_ind=pack(c%row_ind,keep)
      c%col_ind=pack(c%col_ind,keep)
      c%values=pack(c%values,keep)
      call csr_from_coo(c,b,info)
   end subroutine csr_filter

   subroutine csr_map(a,fun,b,drop_zeros)
      type(csr_matrix), intent(in) :: a
      procedure(sparse_map_fun) :: fun
      type(csr_matrix), intent(out) :: b
      logical, intent(in), optional :: drop_zeros
      type(coo_matrix) :: c
      integer :: k, info
      logical :: drop
      drop=.true.
      if (present(drop_zeros)) drop=drop_zeros
      call coo_from_csr(a,c)
      do k=1,c%nnz()
         c%values(k)=fun(c%row_ind(k),c%col_ind(k),c%values(k))
      end do
      if (drop) call coo_remove_zeros(c)
      call csr_from_coo(c,b,info)
   end subroutine csr_map

   subroutine empty_sparse(nrow,ncol,a)
      integer, intent(in) :: nrow,ncol
      type(csr_matrix), intent(out) :: a
      integer :: info
      call csr_from_triplet(max(0,nrow),max(0,ncol),[integer ::],[integer ::], &
         [real(dp) ::],a,info)
   end subroutine empty_sparse

end module matrixextra_utils
