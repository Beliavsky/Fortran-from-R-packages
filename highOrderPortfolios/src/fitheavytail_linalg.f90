! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_linalg
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use fitheavytail_kinds, only: dp
   use fitheavytail_status, only: ht_success, ht_invalid_argument, ht_singular_matrix
   implicit none
   private

   public :: clean_complete_rows, column_mean, sample_covariance, weighted_covariance
   public :: inverse_matrix, solve_linear, logdet_spd, quadratic_forms
   public :: spatial_median, symmetric_eigen, identity_matrix, trace_matrix
   public :: outer_product, standardize_columns, factor_decomposition
   public :: frobenius_norm, is_symmetric_positive_definite, symmetrize

   interface solve_linear
      module procedure solve_linear_vector
      module procedure solve_linear_matrix
   end interface solve_linear

contains

   subroutine clean_complete_rows(x, clean, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: clean(:,:)
      integer, intent(out), optional :: status
      integer :: i, j, k, n_keep
      logical :: good

      n_keep = 0
      do i = 1, size(x,1)
         good = .true.
         do j = 1, size(x,2)
            if (ieee_is_nan(x(i,j))) then
               good = .false.
               exit
            end if
         end do
         if (good) n_keep = n_keep + 1
      end do
      allocate(clean(n_keep, size(x,2)))
      k = 0
      do i = 1, size(x,1)
         good = .true.
         do j = 1, size(x,2)
            if (ieee_is_nan(x(i,j))) then
               good = .false.
               exit
            end if
         end do
         if (good) then
            k = k + 1
            clean(k,:) = x(i,:)
         end if
      end do
      if (present(status)) status = ht_success
   end subroutine clean_complete_rows

   function column_mean(x) result(mu)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: mu(size(x,2))
      if (size(x,1) > 0) then
         mu = sum(x, dim=1)/real(size(x,1), dp)
      else
         mu = 0.0_dp
      end if
   end function column_mean

   function sample_covariance(x, center) result(cov)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: center(:)
      real(dp) :: cov(size(x,2), size(x,2))
      real(dp) :: mu(size(x,2)), row(size(x,2))
      integer :: i, n

      n = size(x,1)
      if (present(center)) then
         mu = center
      else
         mu = column_mean(x)
      end if
      cov = 0.0_dp
      do i = 1, n
         row = x(i,:) - mu
         cov = cov + outer_product(row, row)
      end do
      if (n > 1) cov = cov/real(n-1, dp)
   end function sample_covariance

   function weighted_covariance(xc, weights, divisor) result(cov)
      real(dp), intent(in) :: xc(:,:), weights(:)
      real(dp), intent(in), optional :: divisor
      real(dp) :: cov(size(xc,2), size(xc,2))
      real(dp) :: d
      integer :: i

      cov = 0.0_dp
      do i = 1, size(xc,1)
         cov = cov + weights(i)*outer_product(xc(i,:), xc(i,:))
      end do
      if (present(divisor)) then
         d = divisor
      else
         d = real(size(xc,1), dp)
      end if
      if (abs(d) > tiny(1.0_dp)) cov = cov/d
      cov = symmetrize(cov)
   end function weighted_covariance

   pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a), size(b))
      integer :: j
      do j = 1, size(b)
         c(:,j) = a*b(j)
      end do
   end function outer_product

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   pure function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, min(size(a,1), size(a,2))
         value = value + a(i,i)
      end do
   end function trace_matrix

   pure function symmetrize(a) result(s)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: s(size(a,1),size(a,2))
      s = 0.5_dp*(a + transpose(a))
   end function symmetrize

   function frobenius_norm(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      value = sqrt(sum(a*a))
   end function frobenius_norm

   subroutine solve_linear_vector(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: aug(:,:)
      real(dp) :: pivot, factor, temp
      integer :: n, i, j, k, p, istat

      n = size(a,1)
      istat = ht_success
      if (size(a,2) /= n .or. size(b) /= n .or. size(x) /= n) then
         x = 0.0_dp
         istat = ht_invalid_argument
         if (present(status)) status = istat
         return
      end if
      allocate(aug(n,n+1))
      aug(:,1:n) = a
      aug(:,n+1) = b
      do k = 1, n
         p = k
         pivot = abs(aug(k,k))
         do i = k+1, n
            if (abs(aug(i,k)) > pivot) then
               pivot = abs(aug(i,k))
               p = i
            end if
         end do
         if (pivot <= epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))) then
            x = 0.0_dp
            istat = ht_singular_matrix
            if (present(status)) status = istat
            return
         end if
         if (p /= k) then
            do j = k, n+1
               temp = aug(k,j)
               aug(k,j) = aug(p,j)
               aug(p,j) = temp
            end do
         end if
         do i = k+1, n
            factor = aug(i,k)/aug(k,k)
            aug(i,k:n+1) = aug(i,k:n+1) - factor*aug(k,k:n+1)
         end do
      end do
      x = 0.0_dp
      do i = n, 1, -1
         x(i) = (aug(i,n+1) - dot_product(aug(i,i+1:n), x(i+1:n)))/aug(i,i)
      end do
      if (present(status)) status = istat
   end subroutine solve_linear_vector

   subroutine solve_linear_matrix(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), intent(out) :: x(:,:)
      integer, intent(out), optional :: status
      integer :: j, istat, current

      current = ht_success
      if (size(x,1) /= size(a,1) .or. size(x,2) /= size(b,2)) then
         x = 0.0_dp
         current = ht_invalid_argument
      else
         do j = 1, size(b,2)
            call solve_linear_vector(a, b(:,j), x(:,j), istat)
            if (istat /= ht_success) then
               current = istat
               exit
            end if
         end do
      end if
      if (present(status)) status = current
   end subroutine solve_linear_matrix

   subroutine inverse_matrix(a, inva, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: inva(:,:)
      integer, intent(out), optional :: status
      integer :: istat
      call solve_linear_matrix(a, identity_matrix(size(a,1)), inva, istat)
      if (present(status)) status = istat
   end subroutine inverse_matrix

   function logdet_spd(a, status) result(value)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp), allocatable :: l(:,:)
      real(dp) :: s
      integer :: n, i, j, k, istat

      n = size(a,1)
      allocate(l(n,n))
      l = 0.0_dp
      istat = ht_success
      do i = 1, n
         do j = 1, i
            s = a(i,j)
            do k = 1, j-1
               s = s - l(i,k)*l(j,k)
            end do
            if (i == j) then
               if (s <= 0.0_dp) then
                  value = -huge(1.0_dp)
                  istat = ht_singular_matrix
                  if (present(status)) status = istat
                  return
               end if
               l(i,j) = sqrt(s)
            else
               l(i,j) = s/l(j,j)
            end if
         end do
      end do
      value = 0.0_dp
      do i = 1, n
         value = value + 2.0_dp*log(l(i,i))
      end do
      if (present(status)) status = istat
   end function logdet_spd

   function is_symmetric_positive_definite(a) result(ok)
      real(dp), intent(in) :: a(:,:)
      logical :: ok
      integer :: status
      real(dp) :: dummy
      dummy = logdet_spd(symmetrize(a), status)
      ok = status == ht_success .and. dummy > -huge(1.0_dp)
   end function is_symmetric_positive_definite

   subroutine quadratic_forms(xc, sigma, q, status)
      real(dp), intent(in) :: xc(:,:), sigma(:,:)
      real(dp), intent(out) :: q(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: invs(:,:)
      integer :: i, istat

      allocate(invs(size(sigma,1),size(sigma,2)))
      call inverse_matrix(sigma, invs, istat)
      if (istat /= ht_success) then
         q = huge(1.0_dp)
      else
         do i = 1, size(xc,1)
            q(i) = dot_product(xc(i,:), matmul(invs, xc(i,:)))
         end do
      end if
      if (present(status)) status = istat
   end subroutine quadratic_forms

   subroutine spatial_median(x, median, max_iter, tol, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: median(:)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: status
      real(dp) :: current(size(x,2)), next(size(x,2)), d, denom, tolerance
      integer :: i, iter, niter

      niter = 200
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      current = column_mean(x)
      do iter = 1, niter
         next = 0.0_dp
         denom = 0.0_dp
         do i = 1, size(x,1)
            d = sqrt(sum((x(i,:)-current)**2))
            if (d <= 1.0e-14_dp) then
               median = x(i,:)
               if (present(status)) status = ht_success
               return
            end if
            next = next + x(i,:)/d
            denom = denom + 1.0_dp/d
         end do
         next = next/denom
         if (sqrt(sum((next-current)**2)) <= tolerance*(1.0_dp + sqrt(sum(current**2)))) exit
         current = next
      end do
      median = next
      if (present(status)) status = ht_success
   end subroutine spatial_median

   subroutine symmetric_eigen(a, values, vectors, status)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: values(:), vectors(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: work(:,:)
      real(dp) :: app, aqq, apq, tau, t, c, s, temp, maxoff
      integer :: n, p, q, i, j, iter, max_iter, k, istat

      n = size(a,1)
      istat = ht_success
      if (size(a,2) /= n .or. size(values) /= n .or. any(shape(vectors) /= [n,n])) then
         values = 0.0_dp
         vectors = 0.0_dp
         istat = ht_invalid_argument
         if (present(status)) status = istat
         return
      end if
      allocate(work(n,n))
      work = symmetrize(a)
      vectors = identity_matrix(n)
      max_iter = max(50, 100*n*n)
      do iter = 1, max_iter
         maxoff = 0.0_dp
         p = 1
         q = min(2,n)
         do i = 1, n-1
            do j = i+1, n
               if (abs(work(i,j)) > maxoff) then
                  maxoff = abs(work(i,j))
                  p = i
                  q = j
               end if
            end do
         end do
         if (maxoff <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(work)))) exit
         app = work(p,p)
         aqq = work(q,q)
         apq = work(p,q)
         tau = (aqq-app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau + sqrt(1.0_dp + tau*tau))
         else
            t = -1.0_dp/(-tau + sqrt(1.0_dp + tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp+t*t)
         s = t*c
         do k = 1, n
            if (k /= p .and. k /= q) then
               temp = work(k,p)
               work(k,p) = c*temp - s*work(k,q)
               work(p,k) = work(k,p)
               work(k,q) = s*temp + c*work(k,q)
               work(q,k) = work(k,q)
            end if
         end do
         work(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
         work(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
         work(p,q) = 0.0_dp
         work(q,p) = 0.0_dp
         do k = 1, n
            temp = vectors(k,p)
            vectors(k,p) = c*temp - s*vectors(k,q)
            vectors(k,q) = s*temp + c*vectors(k,q)
         end do
      end do
      do i = 1, n
         values(i) = work(i,i)
      end do
      do i = 1, n-1
         k = i
         do j = i+1, n
            if (values(j) > values(k)) k = j
         end do
         if (k /= i) then
            temp = values(i)
            values(i) = values(k)
            values(k) = temp
            do j = 1, n
               temp = vectors(j,i)
               vectors(j,i) = vectors(j,k)
               vectors(j,k) = temp
            end do
         end if
      end do
      if (present(status)) status = istat
   end subroutine symmetric_eigen

   subroutine standardize_columns(x, z, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: z(:,:)
      integer, intent(out), optional :: status
      real(dp) :: mu(size(x,2)), sd
      integer :: j
      mu = column_mean(x)
      do j = 1, size(x,2)
         if (size(x,1) > 1) then
            sd = sqrt(sum((x(:,j)-mu(j))**2)/real(size(x,1)-1,dp))
         else
            sd = 0.0_dp
         end if
         if (sd > 0.0_dp) then
            z(:,j) = (x(:,j)-mu(j))/sd
         else
            z(:,j) = 0.0_dp
         end if
      end do
      if (present(status)) status = ht_success
   end subroutine standardize_columns

   subroutine factor_decomposition(s, factors, psi, b, sigma, status)
      real(dp), intent(in) :: s(:,:)
      integer, intent(in) :: factors
      real(dp), intent(inout) :: psi(:)
      real(dp), intent(out) :: b(:,:), sigma(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: scaled(:,:), eigval(:), eigvec(:,:)
      real(dp) :: invsqrt(size(psi)), root
      integer :: i, j, n, istat

      n = size(s,1)
      invsqrt = 1.0_dp/sqrt(max(psi, 1.0e-12_dp))
      allocate(scaled(n,n), eigval(n), eigvec(n,n))
      do j = 1, n
         do i = 1, n
            scaled(i,j) = s(i,j)*invsqrt(i)*invsqrt(j)
         end do
      end do
      call symmetric_eigen(scaled, eigval, eigvec, istat)
      b = 0.0_dp
      do j = 1, factors
         root = sqrt(max(1.0_dp, eigval(j))-1.0_dp)
         b(:,j) = sqrt(max(psi,1.0e-12_dp))*eigvec(:,j)*root
      end do
      sigma = matmul(b, transpose(b))
      do i = 1, n
         psi(i) = max(0.0_dp, s(i,i)-sigma(i,i))
         sigma(i,i) = sigma(i,i)+psi(i)
      end do
      if (present(status)) status = istat
   end subroutine factor_decomposition

end module fitheavytail_linalg
