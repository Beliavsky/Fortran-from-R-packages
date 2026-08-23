module directional_linalg
   use directional_kinds, only : dp
   implicit none
   private
   public :: solve_linear, symmetric_eigen, sort_eigen_desc, det3, svd3
contains
   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      integer, intent(out), optional :: info
      real(dp) :: aa(size(b),size(b)), bb(size(b)), row(size(b)), t, f
      integer :: n, i, j, k, piv, ierr
      n=size(b); aa=a; bb=b; ierr=0
      do k=1,n-1
         piv=k
         do i=k+1,n
            if(abs(aa(i,k))>abs(aa(piv,k))) piv=i
         end do
         if(abs(aa(piv,k))<=epsilon(1.0_dp)) then; ierr=k; exit; end if
         if(piv/=k) then
            row=aa(k,:); aa(k,:)=aa(piv,:); aa(piv,:)=row
            t=bb(k); bb(k)=bb(piv); bb(piv)=t
         end if
         do i=k+1,n
            f=aa(i,k)/aa(k,k)
            aa(i,k:n)=aa(i,k:n)-f*aa(k,k:n)
            bb(i)=bb(i)-f*bb(k)
         end do
      end do
      if(ierr==0 .and. abs(aa(n,n))<=epsilon(1.0_dp)) ierr=n
      if(ierr==0) then
         do i=n,1,-1
            t=bb(i)
            do j=i+1,n; t=t-aa(i,j)*x(j); end do
            x(i)=t/aa(i,i)
         end do
      else
         x=0.0_dp
      end if
      if(present(info)) info=ierr
   end subroutine

   subroutine symmetric_eigen(a, eval, evec)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: eval(size(a,1)), evec(size(a,1),size(a,1))
      real(dp) :: d(size(a,1),size(a,1)), app,aqq,apq,tau,t,c,s,dip,diq,vip,viq,maxoff
      integer :: n,i,j,p,q,it
      n=size(a,1); d=0.5_dp*(a+transpose(a)); evec=0.0_dp
      do i=1,n; evec(i,i)=1.0_dp; end do
      do it=1,100*n*n
         maxoff=0.0_dp; p=1; q=min(2,n)
         do i=1,n-1; do j=i+1,n
            if(abs(d(i,j))>maxoff) then; maxoff=abs(d(i,j)); p=i; q=j; end if
         end do; end do
         if(maxoff<1.0e-13_dp) exit
         app=d(p,p); aqq=d(q,q); apq=d(p,q)
         tau=(aqq-app)/(2.0_dp*apq)
         if(tau>=0) then; t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
         else; t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau)); end if
         c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
         do i=1,n
            if(i/=p .and. i/=q) then
               dip=d(i,p); diq=d(i,q)
               d(i,p)=c*dip-s*diq; d(p,i)=d(i,p)
               d(i,q)=s*dip+c*diq; d(q,i)=d(i,q)
            end if
         end do
         d(p,p)=c*c*app-2*c*s*apq+s*s*aqq
         d(q,q)=s*s*app+2*c*s*apq+c*c*aqq
         d(p,q)=0.0_dp; d(q,p)=0.0_dp
         do i=1,n
            vip=evec(i,p); viq=evec(i,q)
            evec(i,p)=c*vip-s*viq; evec(i,q)=s*vip+c*viq
         end do
      end do
      do i=1,n; eval(i)=d(i,i); end do
   end subroutine

   subroutine sort_eigen_desc(eval,evec)
      real(dp),intent(inout)::eval(:),evec(:,:)
      integer::i,j,k,n; real(dp)::t,col(size(evec,1))
      n=size(eval)
      do i=1,n-1
         k=i
         do j=i+1,n; if(eval(j)>eval(k)) k=j; end do
         if(k/=i) then
            t=eval(i);eval(i)=eval(k);eval(k)=t
            col=evec(:,i);evec(:,i)=evec(:,k);evec(:,k)=col
         end if
      end do
   end subroutine

   pure real(dp) function det3(a) result(d)
      real(dp),intent(in)::a(3,3)
      d=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
   end function

   subroutine svd3(a,u,s,v)
      real(dp),intent(in)::a(3,3)
      real(dp),intent(out)::u(3,3),s(3),v(3,3)
      real(dp)::ata(3,3),eval(3),q(3),nrm
      integer::j
      ata=matmul(transpose(a),a)
      call symmetric_eigen(ata,eval,v); call sort_eigen_desc(eval,v)
      s=sqrt(max(eval,0.0_dp)); u=0.0_dp
      do j=1,3
         if(s(j)>1e-12_dp) then; u(:,j)=matmul(a,v(:,j))/s(j); end if
      end do
      ! Modified Gram-Schmidt to stabilize U.
      nrm=sqrt(sum(u(:,1)**2)); if(nrm>0)u(:,1)=u(:,1)/nrm
      u(:,2)=u(:,2)-dot_product(u(:,1),u(:,2))*u(:,1); nrm=sqrt(sum(u(:,2)**2)); if(nrm>0)u(:,2)=u(:,2)/nrm
      q=[u(1,1)*u(2,2)-u(2,1)*u(1,2),u(3,1)*u(1,2)-u(1,1)*u(3,2),u(2,1)*u(3,2)-u(3,1)*u(2,2)]
      u(:,3)=q/max(sqrt(sum(q*q)),tiny(1.0_dp))
   end subroutine
end module directional_linalg
