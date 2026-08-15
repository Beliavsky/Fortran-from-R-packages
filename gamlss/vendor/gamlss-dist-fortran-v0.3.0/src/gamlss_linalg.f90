! SPDX-License-Identifier: GPL-3.0-only
module gamlss_linalg
   use gamlss_kinds, only : dp
   implicit none
   private
   public :: solve_linear, invert_matrix, weighted_least_squares
   public :: penalized_weighted_least_squares, cholesky_factor, matrix_rank
contains

   subroutine solve_linear(a,b,x,status)
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),allocatable,intent(out)::x(:)
      integer,intent(out)::status
      real(dp),allocatable::aa(:,:),bb(:),row(:)
      real(dp)::fac,piv
      integer::n,i,j,k,ip
      n=size(b)
      status=0
      if(size(a,1)/=n .or. size(a,2)/=n)then
         allocate(x(0))
      status=1
      return
      end if
      aa=a
      bb=b
      do k=1,n-1
         ip=k-1+maxloc(abs(aa(k:n,k)),dim=1)
         piv=abs(aa(ip,k))
         if(piv<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa))))then
            allocate(x(0))
      status=2
      return
         end if
         if(ip/=k)then
            allocate(row(n))
      row=aa(k,:)
      aa(k,:)=aa(ip,:)
      aa(ip,:)=row
            deallocate(row)
            piv=bb(k)
      bb(k)=bb(ip)
      bb(ip)=piv
         end if
         do i=k+1,n
            fac=aa(i,k)/aa(k,k)
            aa(i,k)=0.0_dp
            aa(i,k+1:n)=aa(i,k+1:n)-fac*aa(k,k+1:n)
            bb(i)=bb(i)-fac*bb(k)
         end do
      end do
      if(abs(aa(n,n))<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa))))then
         allocate(x(0))
      status=2
      return
      end if
      allocate(x(n))
      x=0.0_dp
      do i=n,1,-1
         if(i<n)then
            x(i)=(bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
         else
            x(i)=bb(i)/aa(i,i)
         end if
      end do
   end subroutine solve_linear

   subroutine invert_matrix(a,ainv,status)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable,intent(out)::ainv(:,:)
      integer,intent(out)::status
      real(dp),allocatable::e(:),col(:)
      integer::n,j,istat
      n=size(a,1)
      status=0
      if(size(a,2)/=n)then
      allocate(ainv(0,0))
      status=1
      return
      end if
      allocate(ainv(n,n),e(n))
      ainv=0.0_dp
      do j=1,n
         e=0.0_dp
      e(j)=1.0_dp
         call solve_linear(a,e,col,istat)
         if(istat/=0)then
      deallocate(ainv)
      allocate(ainv(0,0))
      status=istat
      return
      end if
         ainv(:,j)=col
      end do
   end subroutine invert_matrix

   subroutine weighted_least_squares(x,y,w,beta,cov,status)
      real(dp),intent(in)::x(:,:),y(:),w(:)
      real(dp),allocatable,intent(out)::beta(:),cov(:,:)
      integer,intent(out)::status
      real(dp),allocatable::xtwx(:,:),xtwy(:)
      integer::n,p,i
      n=size(x,1)
      p=size(x,2)
      status=0
      if(size(y)/=n .or. size(w)/=n)then
         allocate(beta(0),cov(0,0))
      status=1
      return
      end if
      allocate(xtwx(p,p),xtwy(p))
      xtwx=0.0_dp
      xtwy=0.0_dp
      do i=1,n
         xtwx=xtwx+w(i)*outer(x(i,:),x(i,:))
         xtwy=xtwy+w(i)*x(i,:)*y(i)
      end do
      call solve_linear(xtwx,xtwy,beta,status)
      if(status/=0)then
      allocate(cov(0,0))
      return
      end if
      call invert_matrix(xtwx,cov,status)
   end subroutine weighted_least_squares

   subroutine penalized_weighted_least_squares(x,y,w,penalty,beta,cov,status)
      real(dp),intent(in)::x(:,:),y(:),w(:),penalty(:,:)
      real(dp),allocatable,intent(out)::beta(:),cov(:,:)
      integer,intent(out)::status
      real(dp),allocatable::xtwx(:,:),xtwy(:)
      integer::n,p,i
      n=size(x,1)
      p=size(x,2)
      status=0
      if(size(y)/=n.or.size(w)/=n.or.size(penalty,1)/=p.or.size(penalty,2)/=p)then
         allocate(beta(0),cov(0,0))
      status=1
      return
      end if
      allocate(xtwx(p,p),xtwy(p))
      xtwx=penalty
      xtwy=0.0_dp
      do i=1,n
         xtwx=xtwx+w(i)*outer(x(i,:),x(i,:))
         xtwy=xtwy+w(i)*x(i,:)*y(i)
      end do
      call solve_linear(xtwx,xtwy,beta,status)
      if(status/=0)then
      allocate(cov(0,0))
      return
      end if
      call invert_matrix(xtwx,cov,status)
   end subroutine penalized_weighted_least_squares

   pure function outer(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::i
      do i=1,size(a)
      c(i,:)=a(i)*b
      end do
   end function outer

   subroutine cholesky_factor(a,l,status)
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable,intent(out)::l(:,:)
      integer,intent(out)::status
      integer::n,i,j,k
      real(dp)::s
      n=size(a,1)
      status=0
      if(size(a,2)/=n)then
      allocate(l(0,0))
      status=1
      return
      end if
      allocate(l(n,n))
      l=0.0_dp
      do i=1,n
         do j=1,i
            s=a(i,j)
            do k=1,j-1
      s=s-l(i,k)*l(j,k)
      end do
            if(i==j)then
               if(s<=0.0_dp)then
      status=2
      return
      end if
               l(i,j)=sqrt(s)
            else
               l(i,j)=s/l(j,j)
            end if
         end do
      end do
   end subroutine cholesky_factor

   integer function matrix_rank(a,tol) result(r)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(in),optional::tol
      real(dp),allocatable::m(:,:),tmp(:)
      real(dp)::t,piv,fac
      integer::nr,nc,i,j,k,ip
      nr=size(a,1)
      nc=size(a,2)
      m=a
      r=0
      t=max(nr,nc)*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
      if(present(tol))t=tol
      k=1
      do j=1,nc
         if(k>nr)exit
         ip=k-1+maxloc(abs(m(k:nr,j)),dim=1)
         piv=abs(m(ip,j))
         if(piv<=t)cycle
         if(ip/=k)then
            allocate(tmp(nc))
      tmp=m(k,:)
      m(k,:)=m(ip,:)
      m(ip,:)=tmp
      deallocate(tmp)
         end if
         do i=k+1,nr
            fac=m(i,j)/m(k,j)
      m(i,j:nc)=m(i,j:nc)-fac*m(k,j:nc)
         end do
         r=r+1
      k=k+1
      end do
   end function matrix_rank

end module gamlss_linalg
