! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_nonlinear
   use robustbase_kinds, only: dp
   use robustbase_scale, only: mad_scale
   use robustbase_psi, only: huber_weight, tukey_weight
   use robustbase_linalg, only: least_squares
   implicit none
   private
   public :: nonlinear_model, robust_nls_result, robust_nls_fit
   abstract interface
      subroutine nonlinear_model(theta,x,yhat)
         import dp
         real(dp),intent(in)::theta(:),x(:,:)
         real(dp),intent(out)::yhat(:)
      end subroutine nonlinear_model
   end interface
   type :: robust_nls_result
      real(dp),allocatable::parameters(:),fitted(:),residuals(:),weights(:)
      real(dp)::scale=0.0_dp,objective=0.0_dp
      integer::iterations=0
      logical::converged=.false.
   end type
contains
   subroutine robust_nls_fit(model,x,y,start,result,psi,tuning,max_iter,tol)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),start(:)
      type(robust_nls_result),intent(out)::result
      character(len=*),intent(in),optional::psi
      real(dp),intent(in),optional::tuning,tol
      integer,intent(in),optional::max_iter
      character(len=16)::ps
      integer::n,p,j,it,mi,info
      real(dp)::c,tt,h,delta
      real(dp),allocatable::theta(:),newtheta(:),fit(:),fitp(:),r(:),w(:),jac(:,:),step(:)
      n=size(y);p=size(start);ps='huber';if(present(psi))ps=adjustl(psi);c=1.345_dp;if(present(tuning))c=tuning;tt=1.0e-7_dp;if(present(tol))tt=tol;mi=100;if(present(max_iter))mi=max_iter
      allocate(theta(p),newtheta(p),fit(n),fitp(n),r(n),w(n),jac(n,p),step(p));theta=start
      do it=1,mi
         call model(theta,x,fit);r=y-fit
         result%scale=mad_scale(r);if(result%scale<=1.0e-14_dp)result%scale=sqrt(sum(r*r)/real(max(1,n-p),dp))
         select case(trim(ps));case('tukey');w=tukey_weight(r/max(result%scale,1.0e-14_dp),c);case default;w=huber_weight(r/max(result%scale,1.0e-14_dp),c);end select
         do j=1,p
            newtheta=theta;h=sqrt(epsilon(1.0_dp))*(1.0_dp+abs(theta(j)));newtheta(j)=newtheta(j)+h;call model(newtheta,x,fitp);jac(:,j)=(fitp-fit)/h
         end do
         call least_squares(jac*spread(sqrt(w),2,p),r*sqrt(w),step,info);if(info/=0)exit
         newtheta=theta+step;delta=maxval(abs(step));theta=newtheta
         if(delta<=tt*(1.0_dp+maxval(abs(theta))))exit
      end do
      call model(theta,x,fit);r=y-fit
      allocate(result%parameters(p),result%fitted(n),result%residuals(n),result%weights(n));result%parameters=theta;result%fitted=fit;result%residuals=r;result%weights=w;result%objective=sum(w*r*r);result%iterations=it;result%converged=(info==0 .and. it<=mi)
   end subroutine robust_nls_fit
end module robustbase_nonlinear
