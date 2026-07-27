! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_linalg
   use optionpricing_kinds, only : dp
   implicit none
   private
   public :: solve_linear, least_squares, symmetric_eigen, orthonormal_complete
   public :: identity_matrix, lower_ones
contains

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i=1,n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   pure function lower_ones(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i, j
      a = 0.0_dp
      do j=1,n
         do i=j,n
            a(i,j) = 1.0_dp
         end do
      end do
   end function lower_ones

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: aug(:,:)
      real(dp) :: pivot, factor, scale
      integer :: n, i, j, k, p

      n = size(b)
      status = 0
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         status = 1
         x = 0.0_dp
         return
      end if
      allocate(aug(n,n+1))
      aug(:,1:n) = a
      aug(:,n+1) = b
      scale = max(1.0_dp,maxval(abs(a)))
      do k=1,n
         p = k-1+maxloc(abs(aug(k:n,k)),dim=1)
         if (abs(aug(p,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            status = 2
            x = 0.0_dp
            return
         end if
         if (p /= k) then
            do j=k,n+1
               pivot = aug(k,j)
               aug(k,j) = aug(p,j)
               aug(p,j) = pivot
            end do
         end if
         pivot = aug(k,k)
         aug(k,k:n+1) = aug(k,k:n+1)/pivot
         do i=1,n
            if (i == k) cycle
            factor = aug(i,k)
            aug(i,k:n+1) = aug(i,k:n+1)-factor*aug(k,k:n+1)
         end do
      end do
      x = aug(:,n+1)
   end subroutine solve_linear

   subroutine least_squares(x, y, beta, status, intercept)
      ! Solve multivariate least squares using pivoted modified Gram-Schmidt.
      real(dp), intent(in) :: x(:,:), y(:,:)
      real(dp), intent(out) :: beta(:,:)
      integer, intent(out) :: status
      logical, intent(in), optional :: intercept
      real(dp), allocatable :: a(:,:), q(:,:), r(:,:), work(:), qty(:,:), bperm(:,:)
      real(dp) :: normv, tol, temp
      integer, allocatable :: piv(:)
      integer :: n, p0, p, m, i, j, k, best, rank
      logical :: with_intercept

      n = size(x,1)
      p0 = size(x,2)
      m = size(y,2)
      with_intercept = .true.
      if (present(intercept)) with_intercept = intercept
      p = p0 + merge(1,0,with_intercept)
      status = 0
      if (size(y,1) /= n .or. size(beta,1) /= p .or. size(beta,2) /= m .or. n < 1) then
         status = 1
         beta = 0.0_dp
         return
      end if
      allocate(a(n,p),q(n,p),r(p,p),work(n),qty(p,m),bperm(p,m),piv(p))
      if (with_intercept) then
         a(:,1) = 1.0_dp
         a(:,2:p) = x
      else
         a = x
      end if
      q = 0.0_dp
      r = 0.0_dp
      do j=1,p
         piv(j)=j
      end do
      tol = sqrt(epsilon(1.0_dp))*max(1.0_dp,maxval(abs(a)))
      rank = 0
      do k=1,p
         best = k
         normv = -1.0_dp
         do j=k,p
            temp = dot_product(a(:,j),a(:,j))
            if (temp > normv) then
               normv = temp
               best = j
            end if
         end do
         if (best /= k) then
            work = a(:,k); a(:,k)=a(:,best); a(:,best)=work
            i=piv(k); piv(k)=piv(best); piv(best)=i
            if (k>1) then
               work(1:k-1)=r(1:k-1,k)
               r(1:k-1,k)=r(1:k-1,best)
               r(1:k-1,best)=work(1:k-1)
            end if
         end if
         normv = sqrt(max(0.0_dp,dot_product(a(:,k),a(:,k))))
         if (normv <= tol) exit
         rank = k
         r(k,k)=normv
         q(:,k)=a(:,k)/normv
         do j=k+1,p
            r(k,j)=dot_product(q(:,k),a(:,j))
            a(:,j)=a(:,j)-r(k,j)*q(:,k)
            ! One reorthogonalization step.
            temp=dot_product(q(:,k),a(:,j))
            r(k,j)=r(k,j)+temp
            a(:,j)=a(:,j)-temp*q(:,k)
         end do
      end do
      qty=0.0_dp
      bperm=0.0_dp
      if (rank > 0) qty(1:rank,:) = matmul(transpose(q(:,1:rank)),y)
      do j=1,m
         do i=rank,1,-1
            bperm(i,j)=(qty(i,j)-dot_product(r(i,i+1:rank),bperm(i+1:rank,j)))/r(i,i)
         end do
      end do
      beta=0.0_dp
      do i=1,p
         beta(piv(i),:)=bperm(i,:)
      end do
      if (rank < p) status=2
   end subroutine least_squares

   subroutine symmetric_eigen(a, values, vectors, status)
      ! Jacobi eigensolver. Eigenvalues/vectors are sorted descending.
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: values(:), vectors(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: w(:,:)
      real(dp) :: app, aqq, apq, phi, c, s, aip, aiq, vip, viq, maxoff, tmp
      integer :: n, i, j, p, q, iter, maxiter, idx

      n=size(a,1)
      status=0
      if (size(a,2)/=n .or. size(values)/=n .or. any(shape(vectors)/=[n,n])) then
         status=1; values=0.0_dp; vectors=0.0_dp; return
      end if
      allocate(w(n,n))
      w=0.5_dp*(a+transpose(a))
      vectors=identity_matrix(n)
      maxiter=max(50,50*n*n)
      do iter=1,maxiter
         maxoff=0.0_dp; p=1; q=min(2,n)
         do j=2,n
            do i=1,j-1
               if (abs(w(i,j))>maxoff) then
                  maxoff=abs(w(i,j)); p=i; q=j
               end if
            end do
         end do
         if (maxoff <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(w)))) exit
         app=w(p,p); aqq=w(q,q); apq=w(p,q)
         phi=0.5_dp*atan2(2.0_dp*apq,aqq-app)
         c=cos(phi); s=sin(phi)
         do i=1,n
            if (i==p .or. i==q) cycle
            aip=w(i,p); aiq=w(i,q)
            w(i,p)=c*aip-s*aiq; w(p,i)=w(i,p)
            w(i,q)=s*aip+c*aiq; w(q,i)=w(i,q)
         end do
         w(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         w(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
         w(p,q)=0.0_dp; w(q,p)=0.0_dp
         do i=1,n
            vip=vectors(i,p); viq=vectors(i,q)
            vectors(i,p)=c*vip-s*viq
            vectors(i,q)=s*vip+c*viq
         end do
      end do
      if (iter>maxiter) status=2
      do i=1,n
         values(i)=w(i,i)
      end do
      do i=1,n-1
         idx=i-1+maxloc(values(i:n),dim=1)
         if (idx/=i) then
            tmp=values(i); values(i)=values(idx); values(idx)=tmp
            do j=1,n
               tmp=vectors(j,i); vectors(j,i)=vectors(j,idx); vectors(j,idx)=tmp
            end do
         end if
      end do
   end subroutine symmetric_eigen

   subroutine orthonormal_complete(first, q, status)
      real(dp), intent(in) :: first(:,:)
      real(dp), intent(out) :: q(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: v(:)
      real(dp) :: nrm
      integer :: n, k, i, j, col
      n=size(q,1); k=size(first,2); status=0
      if (size(q,2)/=n .or. size(first,1)/=n .or. k>n) then
         status=1; q=0.0_dp; return
      end if
      allocate(v(n)); q=0.0_dp; col=0
      do j=1,k
         v=first(:,j)
         do i=1,col
            v=v-dot_product(q(:,i),v)*q(:,i)
         end do
         nrm=sqrt(dot_product(v,v))
         if (nrm>sqrt(epsilon(1.0_dp))) then
            col=col+1; q(:,col)=v/nrm
         end if
      end do
      do j=1,n
         v=0.0_dp; v(j)=1.0_dp
         do i=1,col
            v=v-dot_product(q(:,i),v)*q(:,i)
         end do
         nrm=sqrt(dot_product(v,v))
         if (nrm>sqrt(epsilon(1.0_dp))) then
            col=col+1; q(:,col)=v/nrm
            if (col==n) exit
         end if
      end do
      if (col<n) status=2
   end subroutine orthonormal_complete
end module optionpricing_linalg
