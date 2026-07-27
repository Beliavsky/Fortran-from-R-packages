! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_optimize
   use fhmm_kinds, only: dp
   implicit none
   private
   public :: nelder_mead

   abstract interface
      function objective_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface

contains

   subroutine nelder_mead(f,x0,max_iterations,x_tolerance,f_tolerance,xbest,fbest,iterations,ok)
      procedure(objective_function) :: f
      real(dp),intent(in)::x0(:)
      integer,intent(in)::max_iterations
      real(dp),intent(in)::x_tolerance,f_tolerance
      real(dp),intent(out)::xbest(:),fbest
      integer,intent(out)::iterations
      logical,intent(out)::ok
      real(dp),allocatable::simplex(:,:),values(:),centroid(:),xr(:),xe(:),xc(:),tmpx(:)
      real(dp)::fr,fe,fc,spreadx,spreadf,step
      integer::n,i
      real(dp),parameter::alpha=1.0_dp,gamma=2.0_dp,rho=0.5_dp,sigma=0.5_dp
      n=size(x0);allocate(simplex(n,n+1),values(n+1),centroid(n),xr(n),xe(n),xc(n),tmpx(n))
      simplex(:,1)=x0
      do i=1,n
         simplex(:,i+1)=x0
         step=0.05_dp*max(abs(x0(i)),1.0_dp)
         simplex(i,i+1)=simplex(i,i+1)+step
      end do
      do i=1,n+1;values(i)=safe_eval(f,simplex(:,i));end do
      ok=.false.
      do iterations=1,max_iterations
         call sort_simplex(simplex,values)
         spreadf=maxval(abs(values-values(1)))
         spreadx=0.0_dp
         do i=2,n+1;spreadx=max(spreadx,maxval(abs(simplex(:,i)-simplex(:,1))));end do
         if(spreadf<=f_tolerance .and. spreadx<=x_tolerance)then;ok=.true.;exit;end if
         centroid=sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr=centroid+alpha*(centroid-simplex(:,n+1));fr=safe_eval(f,xr)
         if(fr<values(1))then
            xe=centroid+gamma*(xr-centroid);fe=safe_eval(f,xe)
            if(fe<fr)then;simplex(:,n+1)=xe;values(n+1)=fe
            else;simplex(:,n+1)=xr;values(n+1)=fr;end if
         else if(fr<values(n))then
            simplex(:,n+1)=xr;values(n+1)=fr
         else
            if(fr<values(n+1))then
               xc=centroid+rho*(xr-centroid)
            else
               xc=centroid-rho*(centroid-simplex(:,n+1))
            end if
            fc=safe_eval(f,xc)
            if(fc<min(fr,values(n+1)))then
               simplex(:,n+1)=xc;values(n+1)=fc
            else
               do i=2,n+1
                  simplex(:,i)=simplex(:,1)+sigma*(simplex(:,i)-simplex(:,1))
                  values(i)=safe_eval(f,simplex(:,i))
               end do
            end if
         end if
      end do
      call sort_simplex(simplex,values);xbest=simplex(:,1);fbest=values(1)
      if(iterations>max_iterations)iterations=max_iterations
   end subroutine nelder_mead

   real(dp) function safe_eval(f,x) result(value)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      procedure(objective_function)::f
      real(dp),intent(in)::x(:)
      value=f(x)
      if(.not.ieee_is_finite(value))value=huge(1.0_dp)/100.0_dp
   end function safe_eval

   subroutine sort_simplex(simplex,values)
      real(dp),intent(inout)::simplex(:,:),values(:)
      real(dp),allocatable::col(:)
      real(dp)::v
      integer::i,j,n
      n=size(values);allocate(col(size(simplex,1)))
      do i=2,n
         v=values(i);col=simplex(:,i);j=i-1
         do while(j>=1)
            if(values(j)<=v)exit
            values(j+1)=values(j);simplex(:,j+1)=simplex(:,j);j=j-1
         end do
         values(j+1)=v;simplex(:,j+1)=col
      end do
   end subroutine sort_simplex

end module fhmm_optimize
