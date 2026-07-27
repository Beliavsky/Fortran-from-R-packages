! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_linalg
   use parma_kinds, only: dp
   implicit none
   private
   public :: solve_linear, inverse_matrix, jacobi_eigen, symmetric_sqrt
   public :: nearest_positive_definite, covariance_matrix, mean_columns
   public :: project_box_budget, project_linear_equalities, vector_norm

contains

   function vector_norm(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sqrt(sum(x * x))
   end function vector_norm

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aug(:,:)
      real(dp) :: factor, pivot, tmp
      integer :: i, j, k, n, p

      n = size(b)
      info = 0
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      allocate(aug(n,n+1))
      aug(:,1:n) = a
      aug(:,n+1) = b
      do k = 1, n
         p = k
         pivot = abs(aug(k,k))
         do i = k + 1, n
            if (abs(aug(i,k)) > pivot) then
               pivot = abs(aug(i,k))
               p = i
            end if
         end do
         if (pivot <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
            info = k
            x = 0.0_dp
            return
         end if
         if (p /= k) then
            do j = k, n + 1
               tmp = aug(k,j)
               aug(k,j) = aug(p,j)
               aug(p,j) = tmp
            end do
         end if
         do i = k + 1, n
            factor = aug(i,k) / aug(k,k)
            aug(i,k:n+1) = aug(i,k:n+1) - factor * aug(k,k:n+1)
         end do
      end do
      do i = n, 1, -1
         x(i) = (aug(i,n+1) - dot_product(aug(i,i+1:n), x(i+1:n))) / aug(i,i)
      end do
   end subroutine solve_linear

   subroutine inverse_matrix(a, ainv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: e(:), x(:)
      integer :: j, n, ierr

      n = size(a,1)
      info = 0
      if (size(a,2) /= n .or. any(shape(ainv) /= [n,n])) then
         info = -1
         return
      end if
      allocate(e(n), x(n))
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_linear(a, e, x, ierr)
         if (ierr /= 0) then
            info = ierr
            return
         end if
         ainv(:,j) = x
      end do
   end subroutine inverse_matrix

   subroutine jacobi_eigen(a, values, vectors, info, tol, max_sweeps)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: values(:), vectors(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_sweeps
      real(dp), allocatable :: d(:,:)
      real(dp) :: app, aqq, apq, phi, c, s, tau, threshold
      real(dp) :: dip, diq, vip, viq
      integer :: i, p, q, n, sweep, nsweep

      n = size(a,1)
      info = 0
      if (size(a,2) /= n .or. size(values) /= n .or. any(shape(vectors) /= [n,n])) then
         info = -1
         return
      end if
      threshold = sqrt(epsilon(1.0_dp))
      if (present(tol)) threshold = tol
      nsweep = 100 * max(1,n)
      if (present(max_sweeps)) nsweep = max_sweeps
      allocate(d(n,n))
      d = 0.5_dp * (a + transpose(a))
      vectors = 0.0_dp
      do i = 1, n
         vectors(i,i) = 1.0_dp
      end do
      do sweep = 1, nsweep
         apq = 0.0_dp
         p = 1
         q = min(2,n)
         do i = 1, n - 1
            if (maxval(abs(d(i,i+1:n))) > apq) then
               q = i + maxloc(abs(d(i,i+1:n)), dim=1)
               p = i
               apq = abs(d(p,q))
            end if
         end do
         if (n <= 1 .or. apq <= threshold * max(1.0_dp, maxval(abs(d)))) exit
         app = d(p,p)
         aqq = d(q,q)
         apq = d(p,q)
         phi = 0.5_dp * atan2(2.0_dp * apq, aqq - app)
         c = cos(phi)
         s = sin(phi)
         do i = 1, n
            if (i /= p .and. i /= q) then
               dip = d(i,p)
               diq = d(i,q)
               d(i,p) = c * dip - s * diq
               d(p,i) = d(i,p)
               d(i,q) = s * dip + c * diq
               d(q,i) = d(i,q)
            end if
         end do
         tau = 2.0_dp * c * s * apq
         d(p,p) = c*c*app - tau + s*s*aqq
         d(q,q) = s*s*app + tau + c*c*aqq
         d(p,q) = 0.0_dp
         d(q,p) = 0.0_dp
         do i = 1, n
            vip = vectors(i,p)
            viq = vectors(i,q)
            vectors(i,p) = c * vip - s * viq
            vectors(i,q) = s * vip + c * viq
         end do
      end do
      if (sweep > nsweep) info = 1
      do i = 1, n
         values(i) = d(i,i)
      end do
      call sort_eigenpairs(values, vectors)
   end subroutine jacobi_eigen

   subroutine sort_eigenpairs(values, vectors)
      real(dp), intent(inout) :: values(:), vectors(:,:)
      real(dp) :: tv
      real(dp), allocatable :: col(:)
      integer :: i, j, k, n

      n = size(values)
      allocate(col(size(vectors,1)))
      do i = 1, n - 1
         k = i
         do j = i + 1, n
            if (values(j) > values(k)) k = j
         end do
         if (k /= i) then
            tv = values(i)
            values(i) = values(k)
            values(k) = tv
            col = vectors(:,i)
            vectors(:,i) = vectors(:,k)
            vectors(:,k) = col
         end if
      end do
   end subroutine sort_eigenpairs

   subroutine symmetric_sqrt(a, root, invroot, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: root(:,:)
      real(dp), intent(out), optional :: invroot(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: values(:), vectors(:,:), d(:,:), dinv(:,:)
      integer :: i, n

      n = size(a,1)
      allocate(values(n), vectors(n,n), d(n,n), dinv(n,n))
      call jacobi_eigen(a, values, vectors, info)
      if (info < 0) return
      values = max(values, 1.0e-14_dp)
      d = 0.0_dp
      dinv = 0.0_dp
      do i = 1, n
         d(i,i) = sqrt(values(i))
         dinv(i,i) = 1.0_dp / d(i,i)
      end do
      root = matmul(vectors, matmul(d, transpose(vectors)))
      if (present(invroot)) invroot = matmul(vectors, matmul(dinv, transpose(vectors)))
   end subroutine symmetric_sqrt

   subroutine nearest_positive_definite(a, out, floor_value)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: out(:,:)
      real(dp), intent(in), optional :: floor_value
      real(dp), allocatable :: values(:), vectors(:,:), d(:,:)
      real(dp) :: floorx
      integer :: i, info, n

      n = size(a,1)
      floorx = 1.0e-10_dp
      if (present(floor_value)) floorx = floor_value
      allocate(values(n), vectors(n,n), d(n,n))
      call jacobi_eigen(0.5_dp*(a+transpose(a)), values, vectors, info)
      values = max(values, floorx)
      d = 0.0_dp
      do i = 1, n
         d(i,i) = values(i)
      end do
      out = matmul(vectors, matmul(d, transpose(vectors)))
   end subroutine nearest_positive_definite

   function mean_columns(x) result(mu)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: mu(size(x,2))
      mu = sum(x, dim=1) / real(size(x,1), dp)
   end function mean_columns

   function covariance_matrix(x, center) result(cov)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: center
      real(dp) :: cov(size(x,2),size(x,2))
      real(dp), allocatable :: xc(:,:)
      real(dp) :: mu(size(x,2))
      logical :: do_center
      integer :: n

      n = size(x,1)
      do_center = .true.
      if (present(center)) do_center = center
      allocate(xc(n,size(x,2)))
      if (do_center) then
         mu = mean_columns(x)
         xc = x - spread(mu,1,n)
      else
         xc = x
      end if
      if (n > 1) then
         cov = matmul(transpose(xc), xc) / real(n - 1, dp)
      else
         cov = 0.0_dp
      end if
   end function covariance_matrix

   subroutine project_box_budget(x, lb, ub, budget, projected, info)
      real(dp), intent(in) :: x(:), lb(:), ub(:), budget
      real(dp), intent(out) :: projected(:)
      integer, intent(out) :: info
      real(dp) :: lo, hi, mid, s
      integer :: iter

      info = 0
      if (size(x) /= size(lb) .or. size(x) /= size(ub) .or. size(projected) /= size(x)) then
         info = -1
         return
      end if
      if (sum(lb) > budget + 1.0e-12_dp .or. sum(ub) < budget - 1.0e-12_dp) then
         info = 1
         projected = min(max(x,lb),ub)
         return
      end if
      lo = minval(x - ub) - abs(budget) - 1.0_dp
      hi = maxval(x - lb) + abs(budget) + 1.0_dp
      do iter = 1, 200
         mid = 0.5_dp * (lo + hi)
         projected = min(max(x - mid, lb), ub)
         s = sum(projected)
         if (abs(s - budget) <= 1.0e-13_dp * max(1.0_dp,abs(budget))) then
            lo = mid
            hi = mid
            exit
         end if
         if (s > budget) then
            lo = mid
         else
            hi = mid
         end if
      end do
      projected = min(max(x - 0.5_dp*(lo+hi), lb), ub)
   end subroutine project_box_budget

   subroutine project_linear_equalities(x, a, b, projected, info)
      real(dp), intent(in) :: x(:), a(:,:), b(:)
      real(dp), intent(out) :: projected(:)
      integer, intent(out) :: info
      real(dp), allocatable :: gram(:,:), rhs(:), lambda(:)
      integer :: m

      m = size(a,1)
      if (m == 0) then
         projected = x
         info = 0
         return
      end if
      allocate(gram(m,m), rhs(m), lambda(m))
      gram = matmul(a, transpose(a))
      rhs = matmul(a, x) - b
      call solve_linear(gram, rhs, lambda, info)
      if (info == 0) projected = x - matmul(transpose(a), lambda)
   end subroutine project_linear_equalities

end module parma_linalg
