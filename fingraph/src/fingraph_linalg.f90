! SPDX-License-Identifier: GPL-3.0-only
module fingraph_linalg
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fingraph_kinds, only : dp
   use fingraph_status, only : fg_ok, fg_invalid_input, fg_size_mismatch, &
      fg_singular_matrix, fg_no_convergence, fg_not_positive_definite
   implicit none
   private

   public :: solve_linear_system, inverse_matrix, symmetric_pseudoinverse
   public :: symmetric_eigen_jacobi, inverse_symmetric_positive_definite
   public :: frobenius_norm, vector_norm2, trace_matrix, identity_matrix
   public :: project_simplex, project_simplex_sum, nnqp_projected_gradient
   public :: isotonic_increasing, helmert_basis, rows_correlation
   public :: is_positive_definite, diagonal_of_triple, outer_product
   public :: relative_frobenius_error, symmetrize, sort_indices_ascending

   interface solve_linear_system
      module procedure solve_matrix_rhs
      module procedure solve_vector_rhs
   end interface solve_linear_system
contains
   pure function frobenius_norm(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      value = sqrt(sum(a*a))
   end function frobenius_norm

   pure function vector_norm2(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sqrt(dot_product(x,x))
   end function vector_norm2

   pure function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      integer :: i, n
      n = min(size(a,1), size(a,2))
      value = 0.0_dp
      do i = 1, n
         value = value + a(i,i)
      end do
   end function trace_matrix

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   pure function outer_product(x,y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x),size(y))
      integer :: j
      do j = 1, size(y)
         a(:,j) = x*y(j)
      end do
   end function outer_product

   pure subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:,:)
      a = 0.5_dp*(a + transpose(a))
   end subroutine symmetrize

   pure function relative_frobenius_error(a,b) result(value)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp) :: value, den
      den = frobenius_norm(b)
      if (den <= tiny(1.0_dp)) then
         value = frobenius_norm(a-b)
      else
         value = frobenius_norm(a-b)/den
      end if
   end function relative_frobenius_error

   subroutine solve_matrix_rhs(a,b,x,status)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: aa(:,:), bb(:,:)
      real(dp) :: pivot_value, scale, factor, temp
      integer :: n, nrhs, i, j, k, pivot, local_status

      n = size(a,1)
      nrhs = size(b,2)
      if (n < 1 .or. size(a,2) /= n .or. size(b,1) /= n .or. nrhs < 1) then
         allocate(x(0,0))
         if (present(status)) status = fg_size_mismatch
         return
      end if
      allocate(aa(n,n),bb(n,nrhs),x(n,nrhs))
      aa = a
      bb = b
      scale = max(1.0_dp,maxval(abs(aa)))
      local_status = fg_ok
      do k = 1, n
         pivot = k
         pivot_value = abs(aa(k,k))
         do i = k+1,n
            if (abs(aa(i,k)) > pivot_value) then
               pivot = i
               pivot_value = abs(aa(i,k))
            end if
         end do
         if (pivot_value <= 1000.0_dp*epsilon(1.0_dp)*scale) then
            x = 0.0_dp
            local_status = fg_singular_matrix
            if (present(status)) status = local_status
            return
         end if
         if (pivot /= k) then
            do j = 1,n
               temp = aa(k,j); aa(k,j)=aa(pivot,j); aa(pivot,j)=temp
            end do
            do j = 1,nrhs
               temp = bb(k,j); bb(k,j)=bb(pivot,j); bb(pivot,j)=temp
            end do
         end if
         do i = k+1,n
            factor = aa(i,k)/aa(k,k)
            aa(i,k)=0.0_dp
            if (k < n) aa(i,k+1:n)=aa(i,k+1:n)-factor*aa(k,k+1:n)
            bb(i,:)=bb(i,:)-factor*bb(k,:)
         end do
      end do
      x = 0.0_dp
      do i = n,1,-1
         if (i < n) then
            x(i,:)=(bb(i,:)-matmul(aa(i,i+1:n),x(i+1:n,:)))/aa(i,i)
         else
            x(i,:)=bb(i,:)/aa(i,i)
         end if
      end do
      if (present(status)) status = local_status
   end subroutine solve_matrix_rhs

   subroutine solve_vector_rhs(a,b,x,status)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: bm(:,:), xm(:,:)
      integer :: local_status
      if (size(a,1) /= size(b)) then
         allocate(x(0))
         if (present(status)) status = fg_size_mismatch
         return
      end if
      allocate(bm(size(b),1)); bm(:,1)=b
      call solve_matrix_rhs(a,bm,xm,local_status)
      if (local_status == fg_ok) then
         allocate(x(size(b))); x=xm(:,1)
      else
         allocate(x(0))
      end if
      if (present(status)) status=local_status
   end subroutine solve_vector_rhs

   subroutine inverse_matrix(a,ainv,status)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: eye(:,:)
      integer :: n, local_status
      n=size(a,1)
      if (n < 1 .or. size(a,2)/=n) then
         allocate(ainv(0,0))
         if (present(status)) status=fg_size_mismatch
         return
      end if
      allocate(eye(n,n)); eye=identity_matrix(n)
      call solve_matrix_rhs(a,eye,ainv,local_status)
      if (present(status)) status=local_status
   end subroutine inverse_matrix

   subroutine symmetric_eigen_jacobi(a,values,vectors,status,tolerance,max_iterations)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), allocatable :: work(:,:), tmp_col(:)
      real(dp) :: tol, offmax, app, aqq, apq, tau, t, c, s
      real(dp) :: aik, aqk, vip, viq, temp
      integer :: n, i, j, k, p, q, iter, maxiter, best, local_status

      n=size(a,1)
      if (n < 1 .or. size(a,2)/=n) then
         allocate(values(0),vectors(0,0))
         if (present(status)) status=fg_size_mismatch
         return
      end if
      allocate(work(n,n),values(n),vectors(n,n),tmp_col(n))
      work=0.5_dp*(a+transpose(a))
      vectors=identity_matrix(n)
      tol=1000.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(work)))
      if (present(tolerance)) tol=max(tolerance,epsilon(1.0_dp))
      maxiter=max(200,100*n*n)
      if (present(max_iterations)) maxiter=max(1,max_iterations)
      local_status=fg_no_convergence
      if (n==1) then
         values(1)=work(1,1)
         local_status=fg_ok
      else
         do iter=1,maxiter
            offmax=0.0_dp; p=1; q=2
            do i=1,n-1
               do j=i+1,n
                  if (abs(work(i,j))>offmax) then
                     offmax=abs(work(i,j)); p=i; q=j
                  end if
               end do
            end do
            if (offmax<=tol) then
               local_status=fg_ok
               exit
            end if
            app=work(p,p); aqq=work(q,q); apq=work(p,q)
            tau=(aqq-app)/(2.0_dp*apq)
            if (tau>=0.0_dp) then
               t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
            else
               t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
            end if
            c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
            do k=1,n
               if (k/=p .and. k/=q) then
                  aik=work(k,p); aqk=work(k,q)
                  work(k,p)=c*aik-s*aqk; work(p,k)=work(k,p)
                  work(k,q)=s*aik+c*aqk; work(q,k)=work(k,q)
               end if
            end do
            work(p,p)=app-t*apq; work(q,q)=aqq+t*apq
            work(p,q)=0.0_dp; work(q,p)=0.0_dp
            do k=1,n
               vip=vectors(k,p); viq=vectors(k,q)
               vectors(k,p)=c*vip-s*viq
               vectors(k,q)=s*vip+c*viq
            end do
         end do
         values=[(work(i,i),i=1,n)]
      end if
      ! Armadillo/R eigenvalues are ascending.
      do i=1,n-1
         best=i
         do j=i+1,n
            if (values(j)<values(best)) best=j
         end do
         if (best/=i) then
            temp=values(i); values(i)=values(best); values(best)=temp
            tmp_col=vectors(:,i); vectors(:,i)=vectors(:,best); vectors(:,best)=tmp_col
         end if
      end do
      if (present(status)) status=local_status
   end subroutine symmetric_eigen_jacobi

   subroutine symmetric_pseudoinverse(a,ainv,status,relative_tolerance)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: relative_tolerance
      real(dp), allocatable :: values(:),vectors(:,:),scaled(:,:)
      real(dp) :: tol, scale
      integer :: n, i, local_status
      call symmetric_eigen_jacobi(a,values,vectors,local_status)
      n=size(values)
      if (n==0) then
         allocate(ainv(0,0))
         if (present(status)) status=local_status
         return
      end if
      scale=max(1.0_dp,maxval(abs(values)))
      tol=10000.0_dp*epsilon(1.0_dp)*scale
      if (present(relative_tolerance)) tol=max(relative_tolerance*scale,epsilon(1.0_dp))
      allocate(scaled(n,n)); scaled=vectors
      do i=1,n
         if (abs(values(i))>tol) then
            scaled(:,i)=scaled(:,i)/values(i)
         else
            scaled(:,i)=0.0_dp
         end if
      end do
      allocate(ainv(n,n)); ainv=matmul(scaled,transpose(vectors)); call symmetrize(ainv)
      if (present(status)) status=local_status
   end subroutine symmetric_pseudoinverse

   subroutine inverse_symmetric_positive_definite(a,ainv,status)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: values(:),vectors(:,:),scaled(:,:)
      real(dp) :: tol
      integer :: i,n,local_status
      call symmetric_eigen_jacobi(a,values,vectors,local_status)
      n=size(values)
      if (n==0) then
         allocate(ainv(0,0)); if (present(status)) status=local_status; return
      end if
      tol=10000.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(values)))
      if (minval(values)<=tol) then
         allocate(ainv(n,n)); ainv=0.0_dp
         if (present(status)) status=fg_not_positive_definite
         return
      end if
      allocate(scaled(n,n)); scaled=vectors
      do i=1,n
         scaled(:,i)=scaled(:,i)/values(i)
      end do
      allocate(ainv(n,n)); ainv=matmul(scaled,transpose(vectors)); call symmetrize(ainv)
      if (present(status)) status=local_status
   end subroutine inverse_symmetric_positive_definite

   logical function is_positive_definite(a,tolerance)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: values(:),vectors(:,:)
      real(dp) :: tol
      integer :: status
      call symmetric_eigen_jacobi(a,values,vectors,status)
      if (size(values)==0) then
         is_positive_definite=.false.; return
      end if
      tol=1000.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(values)))
      if (present(tolerance)) tol=tolerance
      is_positive_definite=minval(values)>tol
   end function is_positive_definite

   subroutine project_simplex(v,w,status)
      real(dp), intent(in) :: v(:)
      real(dp), intent(out) :: w(:)
      integer, intent(out), optional :: status
      call project_simplex_sum(v,1.0_dp,w,status)
   end subroutine project_simplex

   subroutine project_simplex_sum(v,total,w,status)
      real(dp), intent(in) :: v(:), total
      real(dp), intent(out) :: w(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: u(:)
      real(dp) :: key,csum,theta
      integer :: n,i,j,rho
      n=size(v)
      if (n<1 .or. size(w)/=n .or. total<0.0_dp) then
         if (size(w)>0) w=0.0_dp
         if (present(status)) status=fg_invalid_input
         return
      end if
      if (total<=tiny(1.0_dp)) then
         w=0.0_dp; if (present(status)) status=fg_ok; return
      end if
      allocate(u(n)); u=v
      do i=2,n
         key=u(i); j=i-1
         do while (j>=1)
            if (u(j)>=key) exit
            u(j+1)=u(j); j=j-1
         end do
         u(j+1)=key
      end do
      csum=0.0_dp; rho=0
      do i=1,n
         csum=csum+u(i)
         if (u(i)+(total-csum)/real(i,dp)>0.0_dp) rho=i
      end do
      if (rho==0) then
         w=total/real(n,dp)
      else
         theta=(sum(u(1:rho))-total)/real(rho,dp)
         w=max(v-theta,0.0_dp)
         if (sum(w)>0.0_dp) w=w*(total/sum(w))
      end if
      if (present(status)) status=fg_ok
   end subroutine project_simplex_sum

   subroutine nnqp_projected_gradient(h,b,x,status,max_iterations,tolerance)
      real(dp), intent(in) :: h(:,:), b(:)
      real(dp), intent(inout) :: x(:)
      integer, intent(out), optional :: status
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: xnew(:),grad(:)
      real(dp) :: lips,step,tol,den
      integer :: n,k,maxiter,local_status
      n=size(b)
      if (size(h,1)/=n .or. size(h,2)/=n .or. size(x)/=n .or. n<1) then
         if (present(status)) status=fg_size_mismatch; return
      end if
      maxiter=5000; if (present(max_iterations)) maxiter=max(1,max_iterations)
      tol=1e-10_dp; if (present(tolerance)) tol=max(tolerance,epsilon(1.0_dp))
      lips=maxval(sum(abs(h),dim=2))
      if (lips<=tiny(1.0_dp)) then
         x=max(x,0.0_dp); if (present(status)) status=fg_ok; return
      end if
      step=1.0_dp/lips
      allocate(xnew(n),grad(n)); x=max(x,0.0_dp)
      local_status=fg_no_convergence
      do k=1,maxiter
         grad=matmul(h,x)-b
         xnew=max(x-step*grad,0.0_dp)
         den=max(1.0_dp,vector_norm2(x))
         if (vector_norm2(xnew-x)/den<tol) then
            x=xnew; local_status=fg_ok; exit
         end if
         x=xnew
      end do
      if (present(status)) status=local_status
   end subroutine nnqp_projected_gradient

   subroutine isotonic_increasing(y,fit)
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: fit(:)
      real(dp), allocatable :: level(:),weight(:)
      integer, allocatable :: first(:),last(:)
      integer :: n,m,i,j
      n=size(y)
      if (size(fit)/=n) return
      allocate(level(n),weight(n),first(n),last(n))
      m=0
      do i=1,n
         m=m+1; level(m)=y(i); weight(m)=1.0_dp; first(m)=i; last(m)=i
         do while (m>1)
            if (level(m-1)<=level(m)) exit
            level(m-1)=(weight(m-1)*level(m-1)+weight(m)*level(m))/(weight(m-1)+weight(m))
            weight(m-1)=weight(m-1)+weight(m)
            last(m-1)=last(m)
            m=m-1
         end do
      end do
      do i=1,m
         do j=first(i),last(i)
            fit(j)=level(i)
         end do
      end do
   end subroutine isotonic_increasing

   pure function helmert_basis(n) result(p)
      integer, intent(in) :: n
      real(dp) :: p(n,max(0,n-1))
      integer :: j
      real(dp) :: denom
      p=0.0_dp
      do j=1,n-1
         denom=sqrt(real(j*(j+1),dp))
         p(1:j,j)=1.0_dp/denom
         p(j+1,j)=-real(j,dp)/denom
      end do
   end function helmert_basis

   subroutine rows_correlation(x,c,status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: c(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: centered(:,:),sd(:)
      real(dp) :: meanv,den
      integer :: p,n,i,j
      p=size(x,1); n=size(x,2)
      if (p<1 .or. n<2) then
         allocate(c(0,0)); if (present(status)) status=fg_invalid_input; return
      end if
      allocate(centered(p,n),sd(p),c(p,p))
      do i=1,p
         meanv=sum(x(i,:))/real(n,dp)
         centered(i,:)=x(i,:)-meanv
         sd(i)=sqrt(dot_product(centered(i,:),centered(i,:))/real(n-1,dp))
      end do
      c=0.0_dp
      do i=1,p
         c(i,i)=1.0_dp
         do j=i+1,p
            den=real(n-1,dp)*sd(i)*sd(j)
            if (den>tiny(1.0_dp)) then
               c(i,j)=dot_product(centered(i,:),centered(j,:))/den
            else
               c(i,j)=0.0_dp
            end if
            c(j,i)=c(i,j)
         end do
      end do
      if (present(status)) status=fg_ok
   end subroutine rows_correlation

   pure function diagonal_of_triple(u,a) result(d)
      real(dp), intent(in) :: u(:,:), a(:,:)
      real(dp) :: d(size(u,2))
      real(dp) :: temp(size(a,1),size(u,2))
      integer :: j
      temp=matmul(a,u)
      do j=1,size(u,2)
         d(j)=dot_product(u(:,j),temp(:,j))
      end do
   end function diagonal_of_triple

   subroutine sort_indices_ascending(x,index)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: index(:)
      integer :: i,j,key
      if (size(index)/=size(x)) return
      index=[(i,i=1,size(x))]
      do i=2,size(index)
         key=index(i); j=i-1
         do while (j>=1)
            if (x(index(j))<=x(key)) exit
            index(j+1)=index(j); j=j-1
         end do
         index(j+1)=key
      end do
   end subroutine sort_indices_ascending
end module fingraph_linalg
