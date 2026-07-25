! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_optimizer
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use msgarch_kinds, only : dp
   implicit none
   private
   type,public::optimizer_result
      real(dp),allocatable::x(:)
      real(dp)::value=huge(1.0_dp)
      integer::iterations=0,evaluations=0
      logical::converged=.false.
   end type optimizer_result
   abstract interface
      function objective_function(x) result(f)
         import dp
         real(dp),intent(in)::x(:)
         real(dp)::f
      end function objective_function
   end interface
   public::nelder_mead
contains
   function nelder_mead(fun,x0,step,tolerance,max_iterations,lower,upper) result(result)
      procedure(objective_function)::fun
      real(dp),intent(in)::x0(:)
      real(dp),intent(in),optional::step(:),tolerance,lower(:),upper(:)
      integer,intent(in),optional::max_iterations
      type(optimizer_result)::result
      real(dp),allocatable::simplex(:,:),value(:),centroid(:),xr(:),xe(:),xc(:),st(:)
      real(dp)::tol,fr,fe,fc,spread_f,spread_x
      integer::n,maxit,i,j
      n=size(x0);maxit=1000;if(present(max_iterations))maxit=max_iterations
      tol=1.0e-7_dp;if(present(tolerance))tol=tolerance
      allocate(simplex(n,n+1),value(n+1),centroid(n),xr(n),xe(n),xc(n),st(n))
      st=0.08_dp*max(1.0_dp,abs(x0));if(present(step))st=step
      simplex(:,1)=x0;call apply_bounds(simplex(:,1),lower,upper)
      do j=2,n+1
         simplex(:,j)=simplex(:,1);simplex(j-1,j)=simplex(j-1,j)+st(j-1);call apply_bounds(simplex(:,j),lower,upper)
      end do
      do j=1,n+1;value(j)=safe_eval(fun,simplex(:,j));end do
      result%evaluations=n+1
      do i=1,maxit
         call sort_simplex(simplex,value)
         spread_f=maxval(abs(value-value(1)))
         spread_x=0.0_dp
         do j=2,n+1;spread_x=max(spread_x,maxval(abs(simplex(:,j)-simplex(:,1))));end do
         if(spread_f<=tol*(1.0_dp+abs(value(1))).and.spread_x<=sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,1)))))then
            result%converged=.true.;exit
         end if
         centroid=sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr=centroid+(centroid-simplex(:,n+1));call apply_bounds(xr,lower,upper)
         fr=safe_eval(fun,xr);result%evaluations=result%evaluations+1
         if(fr<value(1))then
            xe=centroid+2.0_dp*(xr-centroid);call apply_bounds(xe,lower,upper)
            fe=safe_eval(fun,xe);result%evaluations=result%evaluations+1
            if(fe<fr)then;simplex(:,n+1)=xe;value(n+1)=fe;else;simplex(:,n+1)=xr;value(n+1)=fr;end if
         else if(fr<value(n))then
            simplex(:,n+1)=xr;value(n+1)=fr
         else
            if(fr<value(n+1))then;xc=centroid+0.5_dp*(xr-centroid);else;xc=centroid+0.5_dp*(simplex(:,n+1)-centroid);end if
            call apply_bounds(xc,lower,upper);fc=safe_eval(fun,xc);result%evaluations=result%evaluations+1
            if(fc<min(fr,value(n+1)))then
               simplex(:,n+1)=xc;value(n+1)=fc
            else
               do j=2,n+1
                  simplex(:,j)=simplex(:,1)+0.5_dp*(simplex(:,j)-simplex(:,1));call apply_bounds(simplex(:,j),lower,upper)
                  value(j)=safe_eval(fun,simplex(:,j))
               end do
               result%evaluations=result%evaluations+n
            end if
         end if
         result%iterations=i
      end do
      call sort_simplex(simplex,value);allocate(result%x(n));result%x=simplex(:,1);result%value=value(1)
      if(.not.result%converged)result%iterations=maxit
   end function nelder_mead

   function safe_eval(fun,x) result(f)
      procedure(objective_function)::fun
      real(dp),intent(in)::x(:)
      real(dp)::f
      f=fun(x);if(.not.ieee_is_finite(f))f=huge(1.0_dp)/100.0_dp
   end function safe_eval

   subroutine apply_bounds(x,lower,upper)
      real(dp),intent(inout)::x(:)
      real(dp),intent(in),optional::lower(:),upper(:)
      if(present(lower))x=max(x,lower)
      if(present(upper))x=min(x,upper)
   end subroutine apply_bounds

   subroutine sort_simplex(x,f)
      real(dp),intent(inout)::x(:,:),f(:)
      real(dp)::tf,temp(size(x,1))
      integer::i,j,k
      do i=1,size(f)-1
         k=i;do j=i+1,size(f);if(f(j)<f(k))k=j;end do
         if(k/=i)then;tf=f(i);f(i)=f(k);f(k)=tf;temp=x(:,i);x(:,i)=x(:,k);x(:,k)=temp;end if
      end do
   end subroutine sort_simplex
end module msgarch_optimizer
