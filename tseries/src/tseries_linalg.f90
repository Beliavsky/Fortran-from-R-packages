! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries_linalg
   use tseries_kinds, only : dp
   implicit none
   private

   public :: solve_linear
   public :: invert_matrix
   public :: least_squares
   public :: covariance_matrix
   public :: jacobi_eigen
   public :: standardize_columns

contains

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in) :: b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: aug(:, :), rowtmp(:)
      real(dp) :: pivot, factor
      integer :: n, i, k, pivot_row

      n = size(b)
      status = 0
      x = 0.0_dp
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         status = 1
         return
      end if
      allocate(aug(n,n+1), rowtmp(n+1))
      aug(:,1:n) = a
      aug(:,n+1) = b
      do k = 1, n
         pivot_row = k
         do i = k+1, n
            if (abs(aug(i,k)) > abs(aug(pivot_row,k))) pivot_row = i
         end do
         if (abs(aug(pivot_row,k)) <= 100.0_dp*epsilon(1.0_dp)) then
            status = 2
            return
         end if
         if (pivot_row /= k) then
            rowtmp = aug(k,:)
            aug(k,:) = aug(pivot_row,:)
            aug(pivot_row,:) = rowtmp
         end if
         pivot = aug(k,k)
         aug(k,:) = aug(k,:)/pivot
         do i = 1, n
            if (i == k) cycle
            factor = aug(i,k)
            if (abs(factor) > tiny(1.0_dp)) aug(i,:) = aug(i,:) - factor*aug(k,:)
         end do
      end do
      x = aug(:,n+1)
   end subroutine solve_linear

   subroutine invert_matrix(a, inverse, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: inverse(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: e(:), col(:)
      integer :: n, j, info

      n = size(a,1)
      status = 0
      inverse = 0.0_dp
      if (size(a,2) /= n .or. size(inverse,1) /= n .or. size(inverse,2) /= n) then
         status = 1
         return
      end if
      allocate(e(n), col(n))
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_linear(a, e, col, info)
         if (info /= 0) then
            status = info
            return
         end if
         inverse(:,j) = col
      end do
   end subroutine invert_matrix

   subroutine least_squares(x, y, beta, residuals, covariance, status)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: beta(:)
      real(dp), intent(out), optional :: residuals(:)
      real(dp), intent(out), optional :: covariance(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: xtx(:, :), xty(:), fitted(:), inv(:, :)
      real(dp) :: sigma2
      integer :: n, p, info

      n = size(x,1)
      p = size(x,2)
      status = 0
      beta = 0.0_dp
      if (size(y) /= n .or. size(beta) /= p) then
         status = 1
         return
      end if
      allocate(xtx(p,p), xty(p), fitted(n))
      xtx = matmul(transpose(x),x)
      xty = matmul(transpose(x),y)
      call solve_linear(xtx, xty, beta, info)
      if (info /= 0) then
         status = info
         if (present(residuals)) residuals = 0.0_dp
         if (present(covariance)) covariance = 0.0_dp
         return
      end if
      fitted = matmul(x,beta)
      if (present(residuals)) residuals = y-fitted
      if (present(covariance)) then
         if (size(covariance,1) /= p .or. size(covariance,2) /= p) then
            status = 1
            return
         end if
         allocate(inv(p,p))
         call invert_matrix(xtx,inv,info)
         if (info /= 0) then
            status = info
            covariance = 0.0_dp
            return
         end if
         if (n > p) then
            sigma2 = sum((y-fitted)**2)/real(n-p,dp)
         else
            sigma2 = 0.0_dp
         end if
         covariance = sigma2*inv
      end if
   end subroutine least_squares

   subroutine covariance_matrix(x, cov)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(out) :: cov(:, :)
      real(dp), allocatable :: centered(:, :), means(:)
      integer :: n, p, j

      n = size(x,1)
      p = size(x,2)
      cov = 0.0_dp
      if (size(cov,1) /= p .or. size(cov,2) /= p .or. n < 2) return
      allocate(centered(n,p),means(p))
      do j = 1, p
         means(j) = sum(x(:,j))/real(n,dp)
         centered(:,j) = x(:,j)-means(j)
      end do
      cov = matmul(transpose(centered),centered)/real(n-1,dp)
   end subroutine covariance_matrix

   subroutine jacobi_eigen(a, eigenvalues, eigenvectors, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: eigenvalues(:)
      real(dp), intent(out) :: eigenvectors(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: d(:, :)
      real(dp) :: app, aqq, apq, tau, t, c, s, maxoff, dip, diq, vip, viq
      integer :: n, p, q, qbest, i, iter

      n = size(a,1)
      status = 0
      if (size(a,2) /= n .or. size(eigenvalues) /= n .or. &
          size(eigenvectors,1) /= n .or. size(eigenvectors,2) /= n) then
         status = 1
         return
      end if
      allocate(d(n,n))
      d = 0.5_dp*(a+transpose(a))
      eigenvectors = 0.0_dp
      do i = 1, n
         eigenvectors(i,i)=1.0_dp
      end do
      do iter = 1, max(50,100*n*n)
         maxoff = 0.0_dp
         p = 1
         qbest = min(2,n)
         do i = 1, n-1
            do q = i+1, n
               if (abs(d(i,q)) > maxoff) then
                  maxoff = abs(d(i,q))
                  p = i
                  qbest = q
               end if
            end do
         end do
         q = qbest
         if (maxoff <= 1.0e-12_dp*max(1.0_dp,maxval(abs(d)))) exit
         app = d(p,p)
         aqq = d(q,q)
         apq = d(p,q)
         tau = (aqq-app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
         else
            t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp+t*t)
         s = t*c
         do i = 1, n
            if (i /= p .and. i /= q) then
               dip = d(i,p)
               diq = d(i,q)
               d(i,p)=c*dip-s*diq
               d(p,i)=d(i,p)
               d(i,q)=s*dip+c*diq
               d(q,i)=d(i,q)
            end if
         end do
         d(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         d(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
         d(p,q)=0.0_dp
         d(q,p)=0.0_dp
         do i = 1, n
            vip=eigenvectors(i,p)
            viq=eigenvectors(i,q)
            eigenvectors(i,p)=c*vip-s*viq
            eigenvectors(i,q)=s*vip+c*viq
         end do
      end do
      if (iter > max(50,100*n*n)) status=2
      do i=1,n
         eigenvalues(i)=d(i,i)
      end do
      call sort_eigenpairs_descending(eigenvalues,eigenvectors)
   end subroutine jacobi_eigen

   subroutine sort_eigenpairs_descending(values, vectors)
      real(dp), intent(inout) :: values(:)
      real(dp), intent(inout) :: vectors(:, :)
      real(dp) :: tv
      real(dp), allocatable :: col(:)
      integer :: i,j,k,n
      n=size(values)
      allocate(col(size(vectors,1)))
      do i=1,n-1
         k=i
         do j=i+1,n
            if(values(j)>values(k)) k=j
         end do
         if(k/=i) then
            tv=values(i); values(i)=values(k); values(k)=tv
            col=vectors(:,i); vectors(:,i)=vectors(:,k); vectors(:,k)=col
         end if
      end do
   end subroutine sort_eigenpairs_descending

   subroutine standardize_columns(x, z, means, scales)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(out) :: z(:, :)
      real(dp), intent(out), optional :: means(:), scales(:)
      real(dp) :: mu, sd
      integer :: n,p,j
      n=size(x,1); p=size(x,2)
      do j=1,p
         mu=sum(x(:,j))/real(n,dp)
         if(n>1) then
            sd=sqrt(sum((x(:,j)-mu)**2)/real(n-1,dp))
         else
            sd=1.0_dp
         end if
         if(sd<=epsilon(1.0_dp)) sd=1.0_dp
         z(:,j)=(x(:,j)-mu)/sd
         if(present(means)) means(j)=mu
         if(present(scales)) scales(j)=sd
      end do
   end subroutine standardize_columns

end module tseries_linalg
