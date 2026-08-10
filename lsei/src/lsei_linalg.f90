! SPDX-License-Identifier: GPL-2.0-or-later
module lsei_linalg
   use lsei_kinds, only : dp
   implicit none
   private
   public :: symmetric_eigen, pseudoinverse, null_space, dense_solve
   public :: least_squares_pivoted, identity_matrix, matrix_rank
contains
   function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp), allocatable :: a(:,:)
      integer :: i
      allocate(a(n,n)); a = 0.0_dp
      do i=1,n; a(i,i)=1.0_dp; end do
   end function identity_matrix

   subroutine dense_solve(a,b,x,info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: m(:,:), rhs(:), row(:)
      real(dp) :: piv, fac, tmp, scale
      integer :: n,i,j,k,p
      n=size(a,1); info=0; x=0.0_dp
      if(size(a,2)/=n .or. size(b)/=n .or. size(x)/=n) then; info=-1; return; end if
      if(n==0) return
      allocate(m(n,n),rhs(n),row(n)); m=a; rhs=b
      scale=max(1.0_dp,maxval(abs(m)))
      do k=1,n-1
         p=k; piv=abs(m(k,k))
         do i=k+1,n
            if(abs(m(i,k))>piv) then; piv=abs(m(i,k)); p=i; end if
         end do
         if(piv <= 100.0_dp*epsilon(1.0_dp)*scale) then; info=k; return; end if
         if(p/=k) then
            row=m(k,:); m(k,:)=m(p,:); m(p,:)=row
            tmp=rhs(k); rhs(k)=rhs(p); rhs(p)=tmp
         end if
         do i=k+1,n
            fac=m(i,k)/m(k,k); m(i,k)=0.0_dp
            do j=k+1,n; m(i,j)=m(i,j)-fac*m(k,j); end do
            rhs(i)=rhs(i)-fac*rhs(k)
         end do
      end do
      if(abs(m(n,n)) <= 100.0_dp*epsilon(1.0_dp)*scale) then; info=n; return; end if
      x(n)=rhs(n)/m(n,n)
      do i=n-1,1,-1
         x(i)=(rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
      end do
   end subroutine dense_solve

   subroutine symmetric_eigen(a,eval,evec,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: eval(:), evec(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: d(:,:)
      real(dp) :: app,aqq,apq,tau,t,c,s,dip,diq,vip,viq,off,eps,tmp
      integer :: n,p,q,i,sweep,imax
      n=size(a,1); info=0
      if(size(a,2)/=n .or. size(eval)/=n .or. size(evec,1)/=n .or. size(evec,2)/=n) then
         info=-1; return
      end if
      allocate(d(n,n)); d=0.5_dp*(a+transpose(a)); evec=0.0_dp
      do i=1,n; evec(i,i)=1.0_dp; end do
      eps=100.0_dp*epsilon(1.0_dp)
      do sweep=1,max(50,30*n*n)
         off=0.0_dp
         do p=1,n-1; do q=p+1,n; off=max(off,abs(d(p,q))); end do; end do
         if(off <= eps*max(1.0_dp,maxval(abs(d)))) exit
         do p=1,n-1
            do q=p+1,n
               apq=d(p,q)
               if(abs(apq) <= eps*max(1.0_dp,abs(d(p,p))+abs(d(q,q)))) cycle
               app=d(p,p); aqq=d(q,q); tau=(aqq-app)/(2.0_dp*apq)
               if(tau>=0.0_dp) then
                  t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
               else
                  t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
               end if
               c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
               do i=1,n
                  if(i==p .or. i==q) cycle
                  dip=d(i,p); diq=d(i,q)
                  d(i,p)=c*dip-s*diq; d(p,i)=d(i,p)
                  d(i,q)=s*dip+c*diq; d(q,i)=d(i,q)
               end do
               d(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
               d(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
               d(p,q)=0.0_dp; d(q,p)=0.0_dp
               do i=1,n
                  vip=evec(i,p); viq=evec(i,q)
                  evec(i,p)=c*vip-s*viq; evec(i,q)=s*vip+c*viq
               end do
            end do
         end do
      end do
      do i=1,n; eval(i)=d(i,i); end do
      do i=1,n-1
         imax=i
         do p=i+1,n; if(eval(p)>eval(imax)) imax=p; end do
         if(imax/=i) then
            tmp=eval(i); eval(i)=eval(imax); eval(imax)=tmp
            do q=1,n
               tmp=evec(q,i); evec(q,i)=evec(q,imax); evec(q,imax)=tmp
            end do
         end if
      end do
   end subroutine symmetric_eigen

   subroutine pseudoinverse(a,ap,rank,tol,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ap(:,:)
      integer, intent(out) :: rank
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: info
      real(dp), allocatable :: ata(:,:),ev(:),v(:,:),av(:)
      real(dp) :: th,smax,rtol
      integer :: m,n,i,st
      m=size(a,1); n=size(a,2); ap=0.0_dp; rank=0; st=0
      if(size(ap,1)/=n .or. size(ap,2)/=m) then; st=-1; if(present(info)) info=st; return; end if
      if(n==0 .or. m==0) then; if(present(info)) info=0; return; end if
      allocate(ata(n,n),ev(n),v(n,n),av(m)); ata=matmul(transpose(a),a)
      call symmetric_eigen(ata,ev,v,st)
      smax=sqrt(max(0.0_dp,ev(1))); rtol=1.0e-14_dp; if(present(tol)) rtol=tol
      th=(rtol*max(1.0_dp,smax))**2
      do i=1,n
         if(ev(i)>th) then
            av=matmul(a,v(:,i)); ap=ap+spread(v(:,i)/ev(i),2,m)*spread(av,1,n); rank=rank+1
         end if
      end do
      if(present(info)) info=st
   end subroutine pseudoinverse

   integer function matrix_rank(a,tol) result(r)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: ap(:,:)
      integer :: st
      allocate(ap(size(a,2),size(a,1)))
      if(present(tol)) then; call pseudoinverse(a,ap,r,tol,st); else; call pseudoinverse(a,ap,r,info=st); end if
   end function matrix_rank

   subroutine null_space(a,z,rank,tol,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: z(:,:)
      integer, intent(out) :: rank
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: info
      real(dp), allocatable :: ata(:,:),ev(:),v(:,:)
      real(dp) :: th,smax,rtol
      integer :: n,st
      n=size(a,2); allocate(ata(n,n),ev(n),v(n,n)); ata=matmul(transpose(a),a)
      call symmetric_eigen(ata,ev,v,st); rtol=1.0e-14_dp; if(present(tol)) rtol=tol
      if(n==0) then; rank=0; allocate(z(0,0)); if(present(info)) info=st; return; end if
      smax=sqrt(max(0.0_dp,ev(1))); th=(rtol*max(1.0_dp,smax))**2; rank=count(ev>th)
      allocate(z(n,n-rank)); if(n-rank>0) z=v(:,rank+1:n)
      if(present(info)) info=st
   end subroutine null_space

   subroutine least_squares_pivoted(a,b,x,rank,rnorm,pivot,tol)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), intent(out) :: x(:,:)
      integer, intent(out) :: rank
      real(dp), intent(out) :: rnorm(:)
      integer, intent(out) :: pivot(:)
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: ap(:,:)
      integer :: n,j,st
      n=size(a,2); allocate(ap(n,size(a,1)))
      if(present(tol)) then; call pseudoinverse(a,ap,rank,tol,st); else; call pseudoinverse(a,ap,rank,info=st); end if
      x=matmul(ap,b); pivot=[(j,j=1,n)]
      do j=1,size(b,2); rnorm(j)=norm2(matmul(a,x(:,j))-b(:,j)); end do
   end subroutine least_squares_pivoted
end module lsei_linalg
