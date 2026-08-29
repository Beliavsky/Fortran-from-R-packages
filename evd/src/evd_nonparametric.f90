! SPDX-License-Identifier: GPL-3.0-only
module evd_nonparametric
   use r_compat, only : dp
   implicit none
   private
   public :: empirical_to_exponential, abvnonpar, amvnonpar, qcbvnonpar_exp
contains

subroutine empirical_to_exponential(data, out)
   real(dp), intent(in) :: data(:,:)
   real(dp), intent(out) :: out(size(data,1),size(data,2))
   integer :: n, d, i, j, k, ntie
   real(dp) :: rsum
   n = size(data,1)
   d = size(data,2)
   do j=1,d
      do i=1,n
         rsum = 1.0_dp
         ntie = 1
         do k=1,n
            if (k == i) cycle
            if (data(k,j) < data(i,j)) then
               rsum = rsum + 1.0_dp
            else if (data(k,j) == data(i,j)) then
               rsum = rsum + 0.5_dp
               ntie = ntie + 1
            end if
         end do
         ! rsum is the average rank under ties.
         out(i,j) = -log(rsum / real(n+1,dp))
      end do
   end do
end subroutine empirical_to_exponential

pure function abvnonpar(x, data, method, madj, k) result(a)
   real(dp), intent(in) :: x(:), data(:,:)
   character(len=*), intent(in), optional :: method
   integer, intent(in), optional :: madj, k
   real(dp) :: a(size(x))
   character(len=16) :: meth
   integer :: i, j, n, adj, kk
   real(dp) :: d1(size(data,1)), d2(size(data,1)), s1, s2, v, rr(size(data,1)), rrk
   real(dp) :: a0, a1, tmp(size(data,1))
   n = size(data,1)
   if (size(data,2) /= 2) then
      a = huge(1.0_dp)
      return
   end if
   meth = 'cfg'
   if (present(method)) meth = adjustl(method)
   adj = 0
   if (present(madj)) adj = madj
   kk = max(1,n/4)
   if (present(k)) kk = k
   d1=data(:,1)
   d2=data(:,2)
   s1=sum(d1)
   s2=sum(d2)
   select case(trim(meth))
   case('cfg')
      do i=1,size(x)
         v=0.0_dp
         do j=1,n
            v=v+log(max((1.0_dp-x(i))*d1(j),x(i)*d2(j)))
         end do
         a(i)=exp((v-(1.0_dp-x(i))*sum(log(d1))-x(i)*sum(log(d2)))/real(n,dp))
         a(i)=min(1.0_dp,max(a(i),x(i),1.0_dp-x(i)))
      end do
   case('pickands')
      if (adj==2) then
         d1=d1/(s1/real(n,dp))
         d2=d2/(s2/real(n,dp))
      end if
      do i=1,size(x)
         v=0.0_dp
         do j=1,n
            if (x(i)<=0.0_dp) then
               v=v+d2(j)
            else if (x(i)>=1.0_dp) then
               v=v+d1(j)
            else
               v=v+min(d1(j)/x(i),d2(j)/(1.0_dp-x(i)))
            end if
         end do
         if (adj==1) v=v-x(i)*s1-(1.0_dp-x(i))*s2+real(n,dp)
         a(i)=real(n,dp)/v
         a(i)=min(1.0_dp,max(a(i),x(i),1.0_dp-x(i)))
      end do
   case('tdo')
      do i=1,size(x)
         v=0.0_dp
         do j=1,n
            v=v+min(x(i)/(1.0_dp+real(n,dp)*d1(j)), &
                    (1.0_dp-x(i))/(1.0_dp+real(n,dp)*d2(j)))
         end do
         a(i)=1.0_dp-v/(1.0_dp+log(real(n,dp)))
         a(i)=min(1.0_dp,max(a(i),x(i),1.0_dp-x(i)))
      end do
   case('pot')
      rr=1.0_dp/d1+1.0_dp/d2
      call kth_descending(rr,kk+1,rrk)
      do i=1,size(x)
         v=0.0_dp
         do j=1,n
            if (rr(j)>rrk) v=v+max(x(i)/(d1(j)*rr(j)),(1.0_dp-x(i))/(d2(j)*rr(j)))
         end do
         a(i)=2.0_dp*v/real(kk,dp)
      end do
      a0=0.0_dp
      a1=0.0_dp
      do j=1,n
         if (rr(j)>rrk) then
            a0=a0+1.0_dp/(d2(j)*rr(j))
            a1=a1+1.0_dp/(d1(j)*rr(j))
         end if
      end do
      a0=2.0_dp*a0/real(kk,dp)
      a1=2.0_dp*a1/real(kk,dp)
      do i=1,size(x)
         a(i)=a(i)+1.0_dp-(1.0_dp-x(i))*a0-x(i)*a1
         a(i)=min(1.0_dp,max(a(i),x(i),1.0_dp-x(i)))
      end do
   case default
      a=huge(1.0_dp)
   end select
contains
   pure subroutine kth_descending(z, kth, value)
      real(dp),intent(in)::z(:)
      integer,intent(in)::kth
      real(dp),intent(out)::value
      real(dp)::w(size(z)),t
      integer::ii,jj
      w=z
      do ii=2,size(w)
         t=w(ii)
         jj=ii-1
         do while(jj>=1)
            if(w(jj)>=t) exit
            w(jj+1)=w(jj)
            jj=jj-1
         end do
         w(jj+1)=t
      end do
      value=w(min(max(1,kth),size(w)))
   end subroutine kth_descending
end function abvnonpar

pure function amvnonpar(x, data, madj) result(a)
   real(dp), intent(in) :: x(:,:), data(:,:)
   integer, intent(in), optional :: madj
   real(dp) :: a(size(x,1))
   real(dp) :: z(size(data,1),size(data,2)), csum(size(data,2)), rs, v
   integer :: i,j,k,n,d,adj
   n=size(data,1)
   d=size(data,2)
   z=data
   adj=0
   if(present(madj)) adj=madj
   csum=sum(data,dim=1)
   if(adj==2) then
      do k=1,d
         z(:,k)=real(n,dp)*z(:,k)/csum(k)
      end do
   end if
   do i=1,size(x,1)
      rs=sum(x(i,:))
      if(rs<=0.0_dp) then
      a(i)=huge(1.0_dp)
      cycle
      end if
      v=0.0_dp
      do j=1,n
         v=v+minval(z(j,:)/(x(i,:)/rs),mask=x(i,:)>0.0_dp)
      end do
      if(adj==1) v=v-sum((x(i,:)/rs)*csum)+real(n,dp)
      a(i)=real(n,dp)/v
      a(i)=min(1.0_dp,max(a(i),maxval(x(i,:)/rs)))
   end do
end function amvnonpar

subroutine qcbvnonpar_exp(p, x, ax, q1, q2, mint)
   real(dp), intent(in) :: p(:), x(:), ax(:)
   real(dp), intent(out) :: q1(size(x),size(p)), q2(size(x),size(p))
   real(dp), intent(in), optional :: mint
   real(dp) :: mm, lp
   integer :: i,j
   mm=1.0_dp
   if(present(mint)) mm=mint
   do j=1,size(p)
      lp=log(p(j)**mm)
      do i=1,size(x)
         q1(i,j)=-x(i)*lp/ax(i)
         q2(i,j)=-(1.0_dp-x(i))*lp/ax(i)
      end do
   end do
end subroutine qcbvnonpar_exp
end module evd_nonparametric
