! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_linalg
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use pracma_kinds, only : dp, eps_dp
   use pracma_status
   use pracma_types, only : symmetric_eigen_result, linear_solve_result
   implicit none
   private

   public :: outer_product, identity_matrix, symmetrize, solve_linear
   public :: inverse_matrix, inv, mldivide, mrdivide, determinant, logdet_spd
   public :: cholesky, isposdef, nearest_spd, symmetric_eigen, eigjacobi, eig
   public :: pinv, Rank, nullspace, cond, rref, gramSchmidt, qrSolve
   public :: givens, householder, lu, lu_crout, lufact, lusys
   public :: expm, logm, sqrtm, signm, rootm, hessenberg, arnoldi, gmres
   public :: orth, subspace, trace_matrix, frobenius_norm, charpoly
   public :: linearproj, affineproj, procrustes, kabsch, nearest_symmetric

   interface solve_linear
      module procedure solve_linear_vector
      module procedure solve_linear_matrix
   end interface solve_linear

   interface inv
      module procedure inv_function
   end interface inv

   interface mldivide
      module procedure mldivide_vector
      module procedure mldivide_matrix
   end interface mldivide

contains

   pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a), size(b))
      integer :: j
      do j=1,size(b)
         c(:,j)=a*b(j)
      end do
   end function outer_product

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a=0.0_dp
      do i=1,n
         a(i,i)=1.0_dp
      end do
   end function identity_matrix

   pure function symmetrize(a) result(s)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: s(size(a,1),size(a,2))
      s=0.5_dp*(a+transpose(a))
   end function symmetrize

   pure function nearest_symmetric(a) result(s)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: s(size(a,1),size(a,2))
      s=symmetrize(a)
   end function nearest_symmetric

   pure real(dp) function trace_matrix(a)
      real(dp),intent(in)::a(:,:)
      integer::i
      trace_matrix=0.0_dp
      do i=1,min(size(a,1),size(a,2))
         trace_matrix=trace_matrix+a(i,i)
      end do
   end function trace_matrix

   pure real(dp) function frobenius_norm(a)
      real(dp),intent(in)::a(:,:)
      frobenius_norm=sqrt(sum(a*a))
   end function frobenius_norm

   subroutine solve_linear_vector(a,b,x,status)
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(out)::x(:)
      integer,intent(out),optional::status
      real(dp),allocatable::aug(:,:)
      real(dp)::pivot,factor,temp,scale
      integer::n,i,j,k,p,istat
      n=size(a,1); istat=pracma_ok
      if(size(a,2)/=n .or. size(b)/=n .or. size(x)/=n)then
         x=0.0_dp; istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      allocate(aug(n,n+1)); aug(:,1:n)=a; aug(:,n+1)=b
      scale=max(1.0_dp,maxval(abs(a)))
      do k=1,n
         p=k; pivot=abs(aug(k,k))
         do i=k+1,n
            if(abs(aug(i,k))>pivot)then
               pivot=abs(aug(i,k)); p=i
            end if
         end do
         if(pivot<=eps_dp*real(max(1,n),dp)*scale)then
            x=0.0_dp; istat=pracma_singular
            if(present(status))status=istat
            return
         end if
         if(p/=k)then
            do j=k,n+1
               temp=aug(k,j); aug(k,j)=aug(p,j); aug(p,j)=temp
            end do
         end if
         do i=k+1,n
            factor=aug(i,k)/aug(k,k)
            aug(i,k:n+1)=aug(i,k:n+1)-factor*aug(k,k:n+1)
         end do
      end do
      x=0.0_dp
      do i=n,1,-1
         if(i<n)then
            x(i)=(aug(i,n+1)-dot_product(aug(i,i+1:n),x(i+1:n)))/aug(i,i)
         else
            x(i)=aug(i,n+1)/aug(i,i)
         end if
      end do
      if(present(status))status=istat
   end subroutine solve_linear_vector

   subroutine solve_linear_matrix(a,b,x,status)
      real(dp),intent(in)::a(:,:),b(:,:)
      real(dp),intent(out)::x(:,:)
      integer,intent(out),optional::status
      integer::j,istat,current
      current=pracma_ok
      if(size(a,1)/=size(a,2) .or. size(b,1)/=size(a,1) .or. &
         size(x,1)/=size(a,1) .or. size(x,2)/=size(b,2))then
         x=0.0_dp; current=pracma_dimension_mismatch
      else
         do j=1,size(b,2)
            call solve_linear_vector(a,b(:,j),x(:,j),istat)
            if(istat/=pracma_ok)then
               current=istat; exit
            end if
         end do
      end if
      if(present(status))status=current
   end subroutine solve_linear_matrix

   subroutine inverse_matrix(a,inva,status)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::inva(:,:)
      integer,intent(out),optional::status
      integer::istat
      if(size(a,1)/=size(a,2) .or. any(shape(inva)/=shape(a)))then
         inva=0.0_dp; istat=pracma_dimension_mismatch
      else
         call solve_linear_matrix(a,identity_matrix(size(a,1)),inva,istat)
      end if
      if(present(status))status=istat
   end subroutine inverse_matrix

   function inv_function(a,status) result(inva)
      real(dp),intent(in)::a(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::inva(:,:)
      integer::istat
      allocate(inva(size(a,1),size(a,2)))
      call inverse_matrix(a,inva,istat)
      if(istat/=pracma_ok)inva=ieee_value(0.0_dp,ieee_quiet_nan)
      if(present(status))status=istat
   end function inv_function

   function mldivide_vector(a,b,status) result(x)
      real(dp),intent(in)::a(:,:),b(:)
      integer,intent(out),optional::status
      real(dp),allocatable::x(:)
      integer::istat
      if(size(a,1)==size(a,2))then
         allocate(x(size(a,2))); call solve_linear(a,b,x,istat)
      else
         x=least_squares_vector(a,b,istat)
      end if
      if(present(status))status=istat
   end function mldivide_vector

   function mldivide_matrix(a,b,status) result(x)
      real(dp),intent(in)::a(:,:),b(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::x(:,:)
      integer::istat,j,current
      allocate(x(size(a,2),size(b,2))); current=pracma_ok
      do j=1,size(b,2)
         x(:,j)=mldivide_vector(a,b(:,j),istat)
         if(istat/=pracma_ok)current=istat
      end do
      if(present(status))status=current
   end function mldivide_matrix

   function mrdivide(a,b,status) result(x)
      real(dp),intent(in)::a(:,:),b(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::x(:,:)
      integer::istat
      x=transpose(mldivide_matrix(transpose(b),transpose(a),istat))
      if(present(status))status=istat
   end function mrdivide

   function least_squares_vector(a,b,status) result(x)
      real(dp),intent(in)::a(:,:),b(:)
      integer,intent(out)::status
      real(dp),allocatable::x(:)
      real(dp),allocatable::ata(:,:),atb(:),p(:,:)
      integer::n
      n=size(a,2); allocate(x(n))
      if(size(a,1)/=size(b))then
         x=0.0_dp; status=pracma_dimension_mismatch; return
      end if
      ata=matmul(transpose(a),a)
      atb=matmul(transpose(a),b)
      call solve_linear_vector(ata,atb,x,status)
      if(status/=pracma_ok)then
         p=pinv(a)
         x=matmul(p,b)
         status=pracma_ok
      end if
   end function least_squares_vector

   function determinant(a,status) result(d)
      real(dp),intent(in)::a(:,:)
      integer,intent(out),optional::status
      real(dp)::d
      real(dp),allocatable::u(:,:)
      real(dp)::pivot,temp,factor,scale
      integer::n,i,j,k,p,sgn,istat
      n=size(a,1); istat=pracma_ok
      if(size(a,2)/=n)then
         d=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      allocate(u(n,n)); u=a; sgn=1; scale=max(1.0_dp,maxval(abs(a)))
      do k=1,n
         p=k; pivot=abs(u(k,k))
         do i=k+1,n
            if(abs(u(i,k))>pivot)then
               pivot=abs(u(i,k)); p=i
            end if
         end do
         if(pivot<=eps_dp*real(max(1,n),dp)*scale)then
            d=0.0_dp; istat=pracma_singular
            if(present(status))status=istat
            return
         end if
         if(p/=k)then
            do j=1,n
               temp=u(k,j); u(k,j)=u(p,j); u(p,j)=temp
            end do
            sgn=-sgn
         end if
         do i=k+1,n
            factor=u(i,k)/u(k,k)
            u(i,k+1:n)=u(i,k+1:n)-factor*u(k,k+1:n)
         end do
      end do
      d=real(sgn,dp)
      do i=1,n
         d=d*u(i,i)
      end do
      if(present(status))status=istat
   end function determinant

   subroutine cholesky(a,l,status)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::l(:,:)
      integer,intent(out),optional::status
      real(dp)::s
      integer::n,i,j,k,istat
      n=size(a,1); l=0.0_dp; istat=pracma_ok
      if(size(a,2)/=n .or. any(shape(l)/=[n,n]))then
         istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      do i=1,n
         do j=1,i
            s=a(i,j)
            do k=1,j-1
               s=s-l(i,k)*l(j,k)
            end do
            if(i==j)then
               if(s<=0.0_dp)then
                  istat=pracma_not_positive_definite
                  if(present(status))status=istat
                  return
               end if
               l(i,j)=sqrt(s)
            else
               l(i,j)=s/l(j,j)
            end if
         end do
      end do
      if(present(status))status=istat
   end subroutine cholesky

   function logdet_spd(a,status) result(v)
      real(dp),intent(in)::a(:,:)
      integer,intent(out),optional::status
      real(dp)::v
      real(dp),allocatable::l(:,:)
      integer::i,istat
      allocate(l(size(a,1),size(a,2)))
      call cholesky(a,l,istat)
      if(istat/=pracma_ok)then
         v=-huge(1.0_dp)
      else
         v=0.0_dp
         do i=1,size(a,1)
            v=v+2.0_dp*log(l(i,i))
         end do
      end if
      if(present(status))status=istat
   end function logdet_spd

   logical function isposdef(a,tolerance)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::values(:),vectors(:,:)
      real(dp)::tol
      integer::status,n
      n=size(a,1); tol=100.0_dp*eps_dp*max(1.0_dp,maxval(abs(a)))
      if(present(tolerance))tol=tolerance
      if(size(a,2)/=n .or. maxval(abs(a-transpose(a)))>tol)then
         isposdef=.false.; return
      end if
      allocate(values(n),vectors(n,n))
      call symmetric_eigen(a,values,vectors,status)
      isposdef=status==pracma_ok .and. minval(values)>tol
   end function isposdef

   function nearest_spd(a) result(b)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable::b(:,:)
      real(dp),allocatable::v(:,:),w(:),h(:,:),i_n(:,:)
      real(dp)::spacing_value,mineig
      integer::n,status,k
      n=size(a,1); allocate(b(n,n),v(n,n),w(n),h(n,n),i_n(n,n))
      b=symmetrize(a)
      call symmetric_eigen(b,w,v,status)
      w=max(w,0.0_dp)
      h=matmul(v,matmul(diagonal_matrix(w),transpose(v)))
      b=symmetrize(0.5_dp*(b+h))
      i_n=identity_matrix(n)
      k=1
      do while(.not.isposdef(b) .and. k<=100)
         call symmetric_eigen(b,w,v,status)
         mineig=minval(w)
         spacing_value=eps_dp*max(1.0_dp,frobenius_norm(b))
         b=b+i_n*(-mineig*real(k*k,dp)+spacing_value)
         k=k+1
      end do
   end function nearest_spd

   pure function diagonal_matrix(x) result(a)
      real(dp),intent(in)::x(:)
      real(dp)::a(size(x),size(x))
      integer::i
      a=0.0_dp
      do i=1,size(x)
         a(i,i)=x(i)
      end do
   end function diagonal_matrix

   subroutine symmetric_eigen(a,values,vectors,status,max_iter,tolerance)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::values(:),vectors(:,:)
      integer,intent(out),optional::status
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::work(:,:)
      real(dp)::app,aqq,apq,tau,t,c,s,temp,maxoff,tol
      integer::n,p,q,i,j,iter,niter,k,istat
      n=size(a,1); istat=pracma_ok
      if(size(a,2)/=n .or. size(values)/=n .or. any(shape(vectors)/=[n,n]))then
         values=0.0_dp; vectors=0.0_dp; istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      allocate(work(n,n)); work=symmetrize(a); vectors=identity_matrix(n)
      niter=max(50,100*n*n)
      if(present(max_iter))niter=max_iter
      tol=100.0_dp*eps_dp*max(1.0_dp,maxval(abs(work)))
      if(present(tolerance))tol=tolerance
      do iter=1,niter
         maxoff=0.0_dp; p=1; q=min(2,n)
         do i=1,n-1
            do j=i+1,n
               if(abs(work(i,j))>maxoff)then
                  maxoff=abs(work(i,j)); p=i; q=j
               end if
            end do
         end do
         if(maxoff<=tol)exit
         app=work(p,p); aqq=work(q,q); apq=work(p,q)
         tau=(aqq-app)/(2.0_dp*apq)
         if(tau>=0.0_dp)then
            t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
         else
            t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
         end if
         c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
         do k=1,n
            if(k/=p .and. k/=q)then
               temp=work(k,p)
               work(k,p)=c*temp-s*work(k,q); work(p,k)=work(k,p)
               work(k,q)=s*temp+c*work(k,q); work(q,k)=work(k,q)
            end if
         end do
         work(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         work(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
         work(p,q)=0.0_dp; work(q,p)=0.0_dp
         do k=1,n
            temp=vectors(k,p)
            vectors(k,p)=c*temp-s*vectors(k,q)
            vectors(k,q)=s*temp+c*vectors(k,q)
         end do
      end do
      if(iter>niter)istat=pracma_not_converged
      do i=1,n
         values(i)=work(i,i)
      end do
      call sort_eigenpairs(values,vectors)
      if(present(status))status=istat
   end subroutine symmetric_eigen

   subroutine sort_eigenpairs(values,vectors)
      real(dp),intent(inout)::values(:),vectors(:,:)
      integer::i,j,k,n
      real(dp)::temp
      n=size(values)
      do i=1,n-1
         k=i
         do j=i+1,n
            if(values(j)>values(k))k=j
         end do
         if(k/=i)then
            temp=values(i); values(i)=values(k); values(k)=temp
            do j=1,n
               temp=vectors(j,i); vectors(j,i)=vectors(j,k); vectors(j,k)=temp
            end do
         end if
      end do
   end subroutine sort_eigenpairs

   function eigjacobi(a) result(res)
      real(dp),intent(in)::a(:,:)
      type(symmetric_eigen_result)::res
      integer::n
      n=size(a,1)
      allocate(res%values(n),res%vectors(n,n))
      call symmetric_eigen(a,res%values,res%vectors,res%status)
   end function eigjacobi

   function eig(a) result(res)
      real(dp),intent(in)::a(:,:)
      type(symmetric_eigen_result)::res
      integer::n
      n=size(a,1); allocate(res%values(n),res%vectors(n,n))
      if(size(a,2)/=n .or. maxval(abs(a-transpose(a)))>1.0e-10_dp*max(1.0_dp,maxval(abs(a))))then
         res%values=ieee_value(0.0_dp,ieee_quiet_nan)
         res%vectors=0.0_dp; res%status=pracma_unsupported
      else
         call symmetric_eigen(a,res%values,res%vectors,res%status)
      end if
   end function eig

   function pinv(a,tolerance) result(p)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::p(:,:)
      real(dp),allocatable::ata(:,:),values(:),vectors(:,:),inv_values(:)
      real(dp)::tol
      integer::n,status,i
      n=size(a,2); allocate(p(n,size(a,1)),ata(n,n),values(n),vectors(n,n),inv_values(n))
      ata=matmul(transpose(a),a)
      call symmetric_eigen(ata,values,vectors,status)
      tol=max(size(a,1),size(a,2))*eps_dp*max(1.0_dp,maxval(abs(values)))
      if(present(tolerance))tol=tolerance
      inv_values=0.0_dp
      do i=1,n
         if(values(i)>tol)inv_values(i)=1.0_dp/values(i)
      end do
      p=matmul(vectors,matmul(diagonal_matrix(inv_values),matmul(transpose(vectors),transpose(a))))
   end function pinv

   integer function Rank(a,tolerance)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::values(:),vectors(:,:),ata(:,:)
      real(dp)::tol
      integer::n,status
      n=size(a,2); allocate(ata(n,n),values(n),vectors(n,n))
      ata=matmul(transpose(a),a)
      call symmetric_eigen(ata,values,vectors,status)
      tol=max(size(a,1),size(a,2))*eps_dp*max(1.0_dp,maxval(abs(values)))
      if(present(tolerance))tol=tolerance
      Rank=count(values>tol)
   end function Rank

   function nullspace(a,tolerance) result(z)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::z(:,:)
      real(dp),allocatable::values(:),vectors(:,:),ata(:,:)
      real(dp)::tol
      integer::n,status,r,i,k
      n=size(a,2); allocate(ata(n,n),values(n),vectors(n,n))
      ata=matmul(transpose(a),a)
      call symmetric_eigen(ata,values,vectors,status)
      tol=max(size(a,1),size(a,2))*eps_dp*max(1.0_dp,maxval(abs(values)))
      if(present(tolerance))tol=tolerance
      r=count(values>tol); allocate(z(n,n-r)); k=0
      do i=1,n
         if(values(i)<=tol)then
            k=k+1; z(:,k)=vectors(:,i)
         end if
      end do
   end function nullspace

   real(dp) function cond(a,p)
      real(dp),intent(in)::a(:,:)
      character(len=*),intent(in),optional::p
      real(dp),allocatable::values(:),vectors(:,:),ata(:,:)
      integer::n,status
      if(present(p))then
         if(trim(p)=='1')then
            cond=matrix_norm1(a)*matrix_norm1(inv_function(a)); return
         else if(trim(p)=='I' .or. trim(p)=='inf')then
            cond=matrix_norm_inf(a)*matrix_norm_inf(inv_function(a)); return
         end if
      end if
      n=size(a,2); allocate(ata(n,n),values(n),vectors(n,n))
      ata=matmul(transpose(a),a)
      call symmetric_eigen(ata,values,vectors,status)
      if(minval(values)<=eps_dp*max(1.0_dp,maxval(values)))then
         cond=huge(1.0_dp)
      else
         cond=sqrt(maxval(values)/minval(values))
      end if
   end function cond

   pure real(dp) function matrix_norm1(a)
      real(dp),intent(in)::a(:,:)
      matrix_norm1=maxval(sum(abs(a),dim=1))
   end function matrix_norm1

   pure real(dp) function matrix_norm_inf(a)
      real(dp),intent(in)::a(:,:)
      matrix_norm_inf=maxval(sum(abs(a),dim=2))
   end function matrix_norm_inf

   function rref(a,tolerance) result(r)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::r(:,:)
      real(dp)::tol,pivot,temp
      integer::row,col,m,n,i,p
      m=size(a,1); n=size(a,2); allocate(r(m,n)); r=a
      tol=100.0_dp*eps_dp*max(1.0_dp,maxval(abs(a)))
      if(present(tolerance))tol=tolerance
      row=1
      do col=1,n
         if(row>m)exit
         p=row; pivot=abs(r(row,col))
         do i=row+1,m
            if(abs(r(i,col))>pivot)then
               pivot=abs(r(i,col)); p=i
            end if
         end do
         if(pivot<=tol)cycle
         if(p/=row)then
            do i=1,n
               temp=r(row,i); r(row,i)=r(p,i); r(p,i)=temp
            end do
         end if
         r(row,:)=r(row,:)/r(row,col)
         do i=1,m
            if(i/=row)r(i,:)=r(i,:)-r(i,col)*r(row,:)
         end do
         row=row+1
      end do
      where(abs(r)<=tol)r=0.0_dp
   end function rref

   subroutine gramSchmidt(a,q,r,status)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::q(:,:),r(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::v(:)
      real(dp)::nv
      integer::m,n,i,j,istat
      m=size(a,1); n=size(a,2); q=0.0_dp; r=0.0_dp; istat=pracma_ok
      if(any(shape(q)/=[m,n]) .or. any(shape(r)/=[n,n]))then
         istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      allocate(v(m))
      do j=1,n
         v=a(:,j)
         do i=1,j-1
            r(i,j)=dot_product(q(:,i),v)
            v=v-r(i,j)*q(:,i)
         end do
         nv=sqrt(sum(v*v)); r(j,j)=nv
         if(nv<=eps_dp*max(1.0_dp,maxval(abs(a))))then
            istat=pracma_singular; exit
         end if
         q(:,j)=v/nv
      end do
      if(present(status))status=istat
   end subroutine gramSchmidt

   function qrSolve(a,b,status) result(x)
      real(dp),intent(in)::a(:,:),b(:)
      integer,intent(out),optional::status
      real(dp),allocatable::x(:),q(:,:),r(:,:),rhs(:)
      integer::m,n,i,istat
      m=size(a,1); n=size(a,2)
      allocate(x(n),q(m,n),r(n,n),rhs(n))
      call gramSchmidt(a,q,r,istat)
      if(istat==pracma_ok)then
         rhs=matmul(transpose(q),b)
         x=0.0_dp
         do i=n,1,-1
            if(i<n)then
               x(i)=(rhs(i)-dot_product(r(i,i+1:n),x(i+1:n)))/r(i,i)
            else
               x(i)=rhs(i)/r(i,i)
            end if
         end do
      else
         x=0.0_dp
      end if
      if(present(status))status=istat
   end function qrSolve

   subroutine givens(a,b,c,s)
      real(dp),intent(in)::a,b
      real(dp),intent(out)::c,s
      real(dp)::r
      if(abs(b)<=tiny(1.0_dp))then
         c=sign(1.0_dp,a); s=0.0_dp
      else if(abs(a)<=tiny(1.0_dp))then
         c=0.0_dp; s=sign(1.0_dp,b)
      else
         r=sign(sqrt(a*a+b*b),a); c=a/r; s=b/r
      end if
   end subroutine givens

   subroutine householder(x,v,beta)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::v(:),beta
      real(dp)::sigma,mu,v1
      v=x
      if(size(x)<=1)then
         beta=0.0_dp; return
      end if
      sigma=sum(x(2:)**2)
      if(sigma<=tiny(1.0_dp))then
         beta=0.0_dp; v=0.0_dp; v(1)=1.0_dp
      else
         mu=sqrt(x(1)*x(1)+sigma)
         if(x(1)<=0.0_dp)then
            v1=x(1)-mu
         else
            v1=-sigma/(x(1)+mu)
         end if
         beta=2.0_dp*v1*v1/(sigma+v1*v1)
         v=v/v1; v(1)=1.0_dp
      end if
   end subroutine householder

   subroutine lu(a,l,u,piv,status)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::l(:,:),u(:,:)
      integer,intent(out)::piv(:)
      integer,intent(out),optional::status
      real(dp)::temp,factor,pivot
      integer::n,i,j,k,p,istat,itmp
      n=size(a,1); l=identity_matrix(n); u=a; piv=[(i,i=1,n)]; istat=pracma_ok
      if(size(a,2)/=n .or. any(shape(l)/=[n,n]) .or. any(shape(u)/=[n,n]) .or. size(piv)/=n)then
         istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      do k=1,n-1
         p=k; pivot=abs(u(k,k))
         do i=k+1,n
            if(abs(u(i,k))>pivot)then
               pivot=abs(u(i,k)); p=i
            end if
         end do
         if(pivot<=eps_dp*max(1.0_dp,maxval(abs(a))))then
            istat=pracma_singular; exit
         end if
         if(p/=k)then
            do j=1,n
               temp=u(k,j); u(k,j)=u(p,j); u(p,j)=temp
            end do
            if(k>1)then
               do j=1,k-1
                  temp=l(k,j); l(k,j)=l(p,j); l(p,j)=temp
               end do
            end if
            itmp=piv(k); piv(k)=piv(p); piv(p)=itmp
         end if
         do i=k+1,n
            factor=u(i,k)/u(k,k); l(i,k)=factor
            u(i,k:n)=u(i,k:n)-factor*u(k,k:n)
         end do
      end do
      if(abs(u(n,n))<=eps_dp*max(1.0_dp,maxval(abs(a))))istat=pracma_singular
      if(present(status))status=istat
   end subroutine lu

   subroutine lu_crout(a,l,u,status)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::l(:,:),u(:,:)
      integer,intent(out),optional::status
      integer::n,i,j,k,istat
      real(dp)::s
      n=size(a,1); l=0.0_dp; u=0.0_dp; istat=pracma_ok
      do j=1,n
         u(j,j)=1.0_dp
         do i=j,n
            s=a(i,j)
            do k=1,j-1
               s=s-l(i,k)*u(k,j)
            end do
            l(i,j)=s
         end do
         if(abs(l(j,j))<=eps_dp*max(1.0_dp,maxval(abs(a))))then
            istat=pracma_singular; exit
         end if
         do i=j+1,n
            s=a(j,i)
            do k=1,j-1
               s=s-l(j,k)*u(k,i)
            end do
            u(j,i)=s/l(j,j)
         end do
      end do
      if(present(status))status=istat
   end subroutine lu_crout

   subroutine lufact(a,l,u,piv,status)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::l(:,:),u(:,:)
      integer,intent(out)::piv(:)
      integer,intent(out),optional::status
      call lu(a,l,u,piv,status)
   end subroutine lufact

   function lusys(a,b,status) result(x)
      real(dp),intent(in)::a(:,:),b(:)
      integer,intent(out),optional::status
      real(dp),allocatable::x(:)
      integer::istat
      allocate(x(size(b))); call solve_linear(a,b,x,istat)
      if(present(status))status=istat
   end function lusys

   function expm(a) result(e)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable::e(:,:),term(:,:),b(:,:)
      real(dp)::norma
      integer::n,s,k
      n=size(a,1); allocate(e(n,n),term(n,n),b(n,n))
      if(size(a,2)/=n)then
         e=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      norma=matrix_norm1(a)
      if(norma<=0.5_dp)then
         s=0
      else
         s=max(0,ceiling(log(norma/0.5_dp)/log(2.0_dp)))
      end if
      b=a/(2.0_dp**s); e=identity_matrix(n); term=e
      do k=1,80
         term=matmul(term,b)/real(k,dp)
         e=e+term
         if(frobenius_norm(term)<=eps_dp*max(1.0_dp,frobenius_norm(e)))exit
      end do
      do k=1,s
         e=matmul(e,e)
      end do
   end function expm

   function sqrtm(a,status) result(s)
      real(dp),intent(in)::a(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::s(:,:),v(:,:),w(:)
      integer::n,istat
      n=size(a,1); allocate(s(n,n),v(n,n),w(n))
      call symmetric_eigen(a,w,v,istat)
      if(istat/=pracma_ok .or. minval(w)<-100.0_dp*eps_dp*max(1.0_dp,maxval(abs(w))))then
         s=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_not_positive_definite
      else
         w=sqrt(max(w,0.0_dp)); s=matmul(v,matmul(diagonal_matrix(w),transpose(v)))
      end if
      if(present(status))status=istat
   end function sqrtm

   function logm(a,status) result(lg)
      real(dp),intent(in)::a(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::lg(:,:),v(:,:),w(:)
      integer::n,istat
      n=size(a,1); allocate(lg(n,n),v(n,n),w(n))
      call symmetric_eigen(a,w,v,istat)
      if(istat/=pracma_ok .or. minval(w)<=0.0_dp)then
         lg=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_invalid_argument
      else
         w=log(w); lg=matmul(v,matmul(diagonal_matrix(w),transpose(v)))
      end if
      if(present(status))status=istat
   end function logm

   function signm(a,status) result(sg)
      real(dp),intent(in)::a(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::sg(:,:),v(:,:),w(:)
      integer::n,istat,i
      n=size(a,1); allocate(sg(n,n),v(n,n),w(n))
      call symmetric_eigen(a,w,v,istat)
      do i=1,n
         if(w(i)>0.0_dp)then
            w(i)=1.0_dp
         else if(w(i)<0.0_dp)then
            w(i)=-1.0_dp
         else
            w(i)=0.0_dp
         end if
      end do
      sg=matmul(v,matmul(diagonal_matrix(w),transpose(v)))
      if(present(status))status=istat
   end function signm

   function rootm(a,power,status) result(r)
      real(dp),intent(in)::a(:,:)
      integer,intent(in)::power
      integer,intent(out),optional::status
      real(dp),allocatable::r(:,:),v(:,:),w(:)
      integer::n,istat
      n=size(a,1); allocate(r(n,n),v(n,n),w(n))
      call symmetric_eigen(a,w,v,istat)
      if(power<=0 .or. (modulo(power,2)==0 .and. minval(w)<0.0_dp))then
         r=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_invalid_argument
      else
         w=sign(1.0_dp,w)*abs(w)**(1.0_dp/real(power,dp))
         r=matmul(v,matmul(diagonal_matrix(w),transpose(v)))
      end if
      if(present(status))status=istat
   end function rootm

   function hessenberg(a) result(h)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable::h(:,:),v(:)
      real(dp)::beta
      integer::n,k
      n=size(a,1); allocate(h(n,n)); h=a
      do k=1,n-2
         allocate(v(n-k))
         call householder(h(k+1:n,k),v,beta)
         h(k+1:n,k:n)=h(k+1:n,k:n)-beta*outer_product(v,matmul(transpose(h(k+1:n,k:n)),v))
         h(1:n,k+1:n)=h(1:n,k+1:n)-beta*outer_product(matmul(h(1:n,k+1:n),v),v)
         deallocate(v)
      end do
   end function hessenberg

   subroutine arnoldi(a,b,k,q,h,status)
      real(dp),intent(in)::a(:,:),b(:)
      integer,intent(in)::k
      real(dp),intent(out)::q(:,:),h(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::v(:)
      real(dp)::nv
      integer::j,i,istat,n
      n=size(a,1); q=0.0_dp; h=0.0_dp; istat=pracma_ok
      if(size(a,2)/=n .or. size(b)/=n .or. size(q,1)/=n .or. size(q,2)<k+1 .or. &
         size(h,1)<k+1 .or. size(h,2)<k)then
         istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      nv=sqrt(sum(b*b))
      if(nv<=tiny(1.0_dp))then
         istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      q(:,1)=b/nv; allocate(v(n))
      do j=1,k
         v=matmul(a,q(:,j))
         do i=1,j
            h(i,j)=dot_product(q(:,i),v); v=v-h(i,j)*q(:,i)
         end do
         h(j+1,j)=sqrt(sum(v*v))
         if(h(j+1,j)<=eps_dp*max(1.0_dp,frobenius_norm(a)))exit
         q(:,j+1)=v/h(j+1,j)
      end do
      if(present(status))status=istat
   end subroutine arnoldi

   function gmres(a,b,x0,tolerance,max_iter,restart,status) result(x)
      real(dp),intent(in)::a(:,:),b(:),x0(:)
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iter,restart
      integer,intent(out),optional::status
      real(dp),allocatable::x(:),r(:),q(:,:),h(:,:),y(:),e1(:),candidate(:)
      real(dp)::tol,beta
      integer::n,m,niter,outer,j,istat
      n=size(b); m=min(30,n); if(present(restart))m=min(max(1,restart),n)
      niter=100; if(present(max_iter))niter=max_iter
      tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      allocate(x(n),r(n),q(n,m+1),h(m+1,m),e1(m+1),candidate(n)); x=x0
      istat=pracma_not_converged
      do outer=1,niter
         r=b-matmul(a,x); beta=sqrt(sum(r*r))
         if(beta<=tol*max(1.0_dp,sqrt(sum(b*b))))then
            istat=pracma_ok; exit
         end if
         call arnoldi(a,r,m,q,h)
         do j=1,m
            e1=0.0_dp; e1(1)=beta
            y=least_squares_vector(h(1:j+1,1:j),e1(1:j+1),istat)
            candidate=x+matmul(q(:,1:j),y)
            if(sqrt(sum((b-matmul(a,candidate))**2))<=tol*max(1.0_dp,sqrt(sum(b*b))))then
               x=candidate; istat=pracma_ok; exit
            end if
         end do
         if(istat==pracma_ok)exit
         x=candidate
      end do
      if(present(status))status=istat
   end function gmres

   function orth(a,tolerance) result(q)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(in),optional::tolerance
      real(dp),allocatable::q(:,:),values(:),vectors(:,:),aat(:,:)
      real(dp)::tol
      integer::m,status,r,i,k
      m=size(a,1); allocate(aat(m,m),values(m),vectors(m,m))
      aat=matmul(a,transpose(a)); call symmetric_eigen(aat,values,vectors,status)
      tol=max(size(a,1),size(a,2))*eps_dp*max(1.0_dp,maxval(values))
      if(present(tolerance))tol=tolerance
      r=count(values>tol); allocate(q(m,r)); k=0
      do i=1,m
         if(values(i)>tol)then
            k=k+1; q(:,k)=vectors(:,i)
         end if
      end do
   end function orth

   real(dp) function subspace(a,b)
      real(dp),intent(in)::a(:,:),b(:,:)
      real(dp),allocatable::qa(:,:),qb(:,:),m(:,:),values(:),vectors(:,:)
      integer::n,status
      qa=orth(a); qb=orth(b); m=matmul(transpose(qa),qb)
      n=size(m,2); allocate(values(n),vectors(n,n))
      call symmetric_eigen(matmul(transpose(m),m),values,vectors,status)
      subspace=acos(min(1.0_dp,sqrt(max(0.0_dp,minval(values)))))
   end function subspace

   function charpoly(a) result(c)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable::c(:),b(:,:),i_n(:,:)
      integer::n,k
      real(dp)::ck
      n=size(a,1); allocate(c(n+1),b(n,n),i_n(n,n))
      c=0.0_dp; c(1)=1.0_dp; b=identity_matrix(n); i_n=b
      do k=1,n
         b=matmul(a,b)
         ck=-trace_matrix(b)/real(k,dp)
         c(k+1)=ck
         b=b+ck*i_n
      end do
   end function charpoly

   function linearproj(a) result(p)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable::p(:,:)
      p=matmul(a,pinv(a))
   end function linearproj

   function affineproj(points) result(p)
      real(dp),intent(in)::points(:,:)
      real(dp),allocatable::p(:,:)
      real(dp),allocatable::center(:),centered(:,:),q(:,:)
      integer::i
      allocate(center(size(points,2))); center=sum(points,dim=1)/real(size(points,1),dp)
      allocate(centered(size(points,1),size(points,2)))
      do i=1,size(points,1)
         centered(i,:)=points(i,:)-center
      end do
      q=orth(transpose(centered)); p=matmul(q,transpose(q))
   end function affineproj

   subroutine kabsch(x,y,r,status)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp),intent(out)::r(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::h(:,:),u(:,:),v(:,:),w(:)
      integer::n,istat
      n=size(x,2); allocate(h(n,n),u(n,n),v(n,n),w(n))
      h=matmul(transpose(x),y)
      call symmetric_eigen(matmul(transpose(h),h),w,v,istat)
      w=sqrt(max(w,0.0_dp))
      u=matmul(h,v)
      call normalize_columns(u)
      r=matmul(v,transpose(u))
      if(determinant(r)<0.0_dp)then
         v(:,n)=-v(:,n); r=matmul(v,transpose(u))
      end if
      if(present(status))status=istat
   end subroutine kabsch

   subroutine normalize_columns(a)
      real(dp),intent(inout)::a(:,:)
      integer::j
      real(dp)::nrm
      do j=1,size(a,2)
         nrm=sqrt(sum(a(:,j)**2))
         if(nrm>tiny(1.0_dp))a(:,j)=a(:,j)/nrm
      end do
   end subroutine normalize_columns

   subroutine procrustes(x,y,r,scale,translation,residual,status)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp),intent(out)::r(:,:),scale,translation(:),residual
      integer,intent(out),optional::status
      real(dp),allocatable::xc(:,:),yc(:,:),mx(:),my(:)
      integer::i,istat
      allocate(mx(size(x,2)),my(size(y,2)),xc(size(x,1),size(x,2)),yc(size(y,1),size(y,2)))
      if(any(shape(x)/=shape(y)) .or. size(translation)/=size(x,2))then
         r=0.0_dp; scale=0.0_dp; translation=0.0_dp; residual=huge(1.0_dp)
         istat=pracma_dimension_mismatch
         if(present(status))status=istat
         return
      end if
      mx=sum(x,dim=1)/real(size(x,1),dp); my=sum(y,dim=1)/real(size(y,1),dp)
      do i=1,size(x,1)
         xc(i,:)=x(i,:)-mx; yc(i,:)=y(i,:)-my
      end do
      call kabsch(xc,yc,r,istat)
      scale=sum(matmul(xc,r)*yc)/max(tiny(1.0_dp),sum(xc*xc))
      translation=my-scale*matmul(mx,r)
      residual=sqrt(sum((scale*matmul(x,r)+spread(translation,1,size(x,1))-y)**2))
      if(present(status))status=istat
   end subroutine procrustes

end module pracma_linalg
