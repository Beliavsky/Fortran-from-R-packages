! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on MarkowitzR, copyright 2014-2020 Steven E. Pav.
module markowitzr_linalg
   use markowitzr_kinds, only: dp
   implicit none
   private

   public :: column_mean, sample_covariance, covariance_of_mean
   public :: invert_matrix, symmetric_vech, symmetric_ivech
   public :: kronecker_product, duplication_matrix, identity_matrix
   public :: lower_vector_indices, vech_index, matrix_rank
   public :: symmetrize_matrix

contains

   pure function column_mean(x) result(mu)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: mu(size(x,2))

      if (size(x,1) > 0) then
         mu = sum(x,dim=1)/real(size(x,1),dp)
      else
         mu = 0.0_dp
      end if
   end function column_mean

   pure function sample_covariance(x) result(sigma)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: sigma(size(x,2),size(x,2))
      real(dp) :: mu(size(x,2)), d(size(x,2))
      integer :: i, n

      n = size(x,1)
      sigma = 0.0_dp
      if (n <= 1) return
      mu = column_mean(x)
      do i = 1, n
         d = x(i,:)-mu
         sigma = sigma+spread(d,2,size(d))*spread(d,1,size(d))
      end do
      sigma = sigma/real(n-1,dp)
   end function sample_covariance

   pure function covariance_of_mean(x) result(sigma)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: sigma(size(x,2),size(x,2))

      if (size(x,1) > 0) then
         sigma = sample_covariance(x)/real(size(x,1),dp)
      else
         sigma = 0.0_dp
      end if
   end function covariance_of_mean

   subroutine invert_matrix(a, ainv, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: ainv(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: aug(:, :), row_tmp(:)
      real(dp) :: pivot, factor, scale
      integer :: n, i, k, pivot_row

      n = size(a,1)
      status = 0
      ainv = 0.0_dp
      if (size(a,2) /= n .or. any(shape(ainv) /= [n,n])) then
         status = 1
         return
      end if
      if (n == 0) then
         status = 2
         return
      end if

      allocate(aug(n,2*n),row_tmp(2*n))
      aug = 0.0_dp
      aug(:,1:n) = a
      do i = 1, n
         aug(i,n+i) = 1.0_dp
      end do
      scale = max(1.0_dp,maxval(abs(a)))

      do k = 1, n
         pivot_row = k-1+maxloc(abs(aug(k:n,k)),dim=1)
         pivot = aug(pivot_row,k)
         if (abs(pivot) <= 1000.0_dp*epsilon(1.0_dp)*scale) then
            status = 2
            return
         end if
         if (pivot_row /= k) then
            row_tmp = aug(k,:)
            aug(k,:) = aug(pivot_row,:)
            aug(pivot_row,:) = row_tmp
         end if
         aug(k,:) = aug(k,:)/aug(k,k)
         do i = 1, n
            if (i == k) cycle
            factor = aug(i,k)
            aug(i,:) = aug(i,:)-factor*aug(k,:)
         end do
      end do
      ainv = aug(:,n+1:2*n)
   end subroutine invert_matrix

   pure function symmetric_vech(a) result(v)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: v(size(a,1)*(size(a,1)+1)/2)
      integer :: i, j, k, n

      n = size(a,1)
      k = 0
      do j = 1, n
         do i = j, n
            k = k+1
            v(k) = a(i,j)
         end do
      end do
   end function symmetric_vech

   function symmetric_ivech(v, status) result(a)
      real(dp), intent(in) :: v(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: a(:, :)
      integer :: n, i, j, k

      n = int((sqrt(1.0_dp+8.0_dp*real(size(v),dp))-1.0_dp)/2.0_dp)
      if (n*(n+1)/2 /= size(v)) then
         allocate(a(0,0))
         if (present(status)) status = 1
         return
      end if
      allocate(a(n,n))
      a = 0.0_dp
      k = 0
      do j = 1, n
         do i = j, n
            k = k+1
            a(i,j) = v(k)
            a(j,i) = v(k)
         end do
      end do
      if (present(status)) status = 0
   end function symmetric_ivech

   pure function kronecker_product(a, b) result(kron)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp) :: kron(size(a,1)*size(b,1),size(a,2)*size(b,2))
      integer :: i, j, nr, nc

      nr = size(b,1)
      nc = size(b,2)
      do j = 1, size(a,2)
         do i = 1, size(a,1)
            kron((i-1)*nr+1:i*nr,(j-1)*nc+1:j*nc) = a(i,j)*b
         end do
      end do
   end function kronecker_product

   pure function duplication_matrix(n) result(d)
      integer, intent(in) :: n
      real(dp) :: d(n*n,n*(n+1)/2)
      integer :: i, j, k

      d = 0.0_dp
      k = 0
      do j = 1, n
         do i = j, n
            k = k+1
            d(i+(j-1)*n,k) = 1.0_dp
            if (i /= j) d(j+(i-1)*n,k) = 1.0_dp
         end do
      end do
   end function duplication_matrix

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   pure function lower_vector_indices(n) result(indices)
      integer, intent(in) :: n
      integer :: indices(n*(n+1)/2)
      integer :: i, j, k

      k = 0
      do j = 1, n
         do i = j, n
            k = k+1
            indices(k) = i+(j-1)*n
         end do
      end do
   end function lower_vector_indices

   pure integer function vech_index(n, row_index, col_index) result(index_value)
      integer, intent(in) :: n, row_index, col_index
      integer :: i, j, k

      index_value = 0
      if (row_index < col_index .or. row_index > n .or. col_index < 1) return
      k = 0
      do j = 1, n
         do i = j, n
            k = k+1
            if (i == row_index .and. j == col_index) then
               index_value = k
               return
            end if
         end do
      end do
   end function vech_index

   function matrix_rank(a, tolerance) result(rank_value)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in), optional :: tolerance
      integer :: rank_value
      real(dp), allocatable :: work(:, :), row_tmp(:)
      real(dp) :: tol, pivot, factor
      integer :: m, n, i, j, pivot_row, active_row

      m = size(a,1)
      n = size(a,2)
      allocate(work(m,n),row_tmp(n))
      work = a
      tol = 1000.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
      if (present(tolerance)) tol = tolerance
      rank_value = 0
      active_row = 1
      do j = 1, n
         if (active_row > m) exit
         pivot_row = active_row-1+maxloc(abs(work(active_row:m,j)),dim=1)
         pivot = work(pivot_row,j)
         if (abs(pivot) <= tol) cycle
         if (pivot_row /= active_row) then
            row_tmp = work(active_row,:)
            work(active_row,:) = work(pivot_row,:)
            work(pivot_row,:) = row_tmp
         end if
         work(active_row,:) = work(active_row,:)/work(active_row,j)
         do i = active_row+1, m
            factor = work(i,j)
            work(i,:) = work(i,:)-factor*work(active_row,:)
         end do
         rank_value = rank_value+1
         active_row = active_row+1
      end do
   end function matrix_rank

   subroutine symmetrize_matrix(a)
      real(dp), intent(inout) :: a(:, :)
      real(dp) :: value
      integer :: i, j

      do j = 1, size(a,2)
         do i = j+1, size(a,1)
            value = 0.5_dp*(a(i,j)+a(j,i))
            a(i,j) = value
            a(j,i) = value
         end do
      end do
   end subroutine symmetrize_matrix

end module markowitzr_linalg
