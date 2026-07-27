! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on SharpeR, copyright 2012-2025 Steven E. Pav.
module sharper_linalg
   use sharper_kinds, only: dp
   implicit none
   private

   public :: column_mean, sample_covariance, covariance_of_mean
   public :: solve_linear, invert_matrix, quadratic_form
   public :: outer_quadratic_form, symmetric_vech, symmetric_ivech
   public :: kronecker_product, identity_matrix, matrix_rank

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
      real(dp) :: mu(size(x,2))
      real(dp) :: d(size(x,2))
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

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in) :: b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: aug(:, :), rowtmp(:)
      real(dp) :: pivot, factor
      integer :: n, i, k, p

      n = size(a,1)
      status = 0
      x = 0.0_dp
      if (size(a,2) /= n .or. size(b) /= n .or. size(x) /= n) then
         status = 1
         return
      end if
      allocate(aug(n,n+1),rowtmp(n+1))
      aug(:,1:n) = a
      aug(:,n+1) = b
      do k = 1, n
         p = k-1+maxloc(abs(aug(k:n,k)),dim=1)
         pivot = aug(p,k)
         if (abs(pivot) <= 100.0_dp*tiny(1.0_dp)) then
            status = 2
            return
         end if
         if (p /= k) then
            rowtmp = aug(k,:)
            aug(k,:) = aug(p,:)
            aug(p,:) = rowtmp
         end if
         aug(k,:) = aug(k,:)/aug(k,k)
         do i = 1, n
            if (i == k) cycle
            factor = aug(i,k)
            if (abs(factor) > 0.0_dp) aug(i,:) = aug(i,:)-factor*aug(k,:)
         end do
      end do
      x = aug(:,n+1)
   end subroutine solve_linear

   subroutine invert_matrix(a, ainv, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: ainv(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: e(:), col(:)
      integer :: n, j, local_status
      n = size(a,1)
      status = 0
      ainv = 0.0_dp
      if (size(a,2) /= n .or. any(shape(ainv) /= [n,n])) then
         status = 1
         return
      end if
      allocate(e(n),col(n))
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_linear(a,e,col,local_status)
         if (local_status /= 0) then
            status = local_status
            return
         end if
         ainv(:,j) = col
      end do
   end subroutine invert_matrix

   pure function quadratic_form(x, a) result(v)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: v
      v = dot_product(x,matmul(a,x))
   end function quadratic_form

   pure function outer_quadratic_form(h, a) result(v)
      real(dp), intent(in) :: h(:, :), a(:, :)
      real(dp) :: v(size(h,1),size(h,1))
      v = matmul(h,matmul(a,transpose(h)))
   end function outer_quadratic_form

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

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   function matrix_rank(a, tolerance) result(rank_value)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in), optional :: tolerance
      integer :: rank_value
      real(dp), allocatable :: work(:, :), rowtmp(:)
      real(dp) :: tol, pivot, factor
      integer :: m, n, i, j, p, row
      m = size(a,1)
      n = size(a,2)
      allocate(work(m,n),rowtmp(n))
      work = a
      tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
      if (present(tolerance)) tol = tolerance
      rank_value = 0
      row = 1
      do j = 1, n
         if (row > m) exit
         p = row-1+maxloc(abs(work(row:m,j)),dim=1)
         pivot = work(p,j)
         if (abs(pivot) <= tol) cycle
         if (p /= row) then
            rowtmp = work(row,:)
            work(row,:) = work(p,:)
            work(p,:) = rowtmp
         end if
         work(row,:) = work(row,:)/work(row,j)
         do i = row+1, m
            factor = work(i,j)
            work(i,:) = work(i,:)-factor*work(row,:)
         end do
         rank_value = rank_value+1
         row = row+1
      end do
   end function matrix_rank

end module sharper_linalg
