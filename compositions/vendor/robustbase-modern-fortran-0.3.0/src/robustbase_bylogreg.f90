! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_bylogreg
   use robustbase_kinds, only: dp
   use robustbase_probability, only: normal_cdf
   use robustbase_linalg, only: invert_symmetric
   use robustbase_regression, only: robust_regression_result, robust_glm_fit
   use robustbase_detmcd, only: detmcd_result, cov_detmcd
   implicit none
   private
   public :: by_logistic_result, by_logistic_fit, by_phi, by_phi_derivative, by_phi_second_derivative

   type :: by_logistic_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: standard_errors(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: objective = 0.0_dp
      integer :: iterations = 0
      logical :: converged = .false.
   end type by_logistic_result
contains
   subroutine by_logistic_fit(x,y,result,const,start,max_iter,tol)
      real(dp),intent(in)::x(:,:),y(:)
      type(by_logistic_result),intent(out)::result
      real(dp),intent(in),optional::const,start(:),tol
      integer,intent(in),optional::max_iter
      type(robust_regression_result)::initial
      type(detmcd_result)::mcd_start
      real(dp),allocatable::beta(:),xi(:),newxi(:),gradient(:),direction(:),scores(:),stscores(:), &
                              matm(:,:),ifsqr(:,:),minv(:,:),outerm(:,:),eta(:),x_start(:,:),y_start(:)
      logical,allocatable::start_mask(:)
      real(dp)::c,tt,obj,newobj,sigma,norm_beta,norm_direction,step
      integer::n,p,it,mi,info,i,half
      logical::accepted

      n=size(x,1);p=size(x,2)
      if(size(y)/=n .or. p<1 .or. any(y<0.0_dp) .or. any(y>1.0_dp)) error stop "by_logistic_fit: invalid input"
      c=0.5_dp
      if(present(const)) c=const
      if(c<=0.0_dp) error stop "by_logistic_fit: const must be positive"
      mi=1000
      if(present(max_iter)) mi=max(1,max_iter)
      tt=1.0e-7_dp
      if(present(tol)) tt=max(tol,epsilon(1.0_dp))
      allocate(beta(p),xi(p),newxi(p),gradient(p),direction(p),scores(n),stscores(n),matm(p,p),ifsqr(p,p),minv(p,p),outerm(p,p),eta(n))
      if(present(start)) then
         if(size(start)/=p) error stop "by_logistic_fit: start size mismatch"
         beta=start
      else
         if(p>1 .and. all(abs(x(:,1)-1.0_dp)<=10.0_dp*epsilon(1.0_dp))) then
            call cov_detmcd(x(:,2:p),mcd_start,alpha=0.75_dp,save_orderings=.false.)
            allocate(start_mask(n))
            start_mask=mcd_start%estimate%weights
            if(count(start_mask)>=p+1) then
               call pack_rows(x,y,start_mask,x_start,y_start)
               call robust_glm_fit(x_start,y_start,'binomial',initial,max_iter=100)
            else
               call robust_glm_fit(x,y,'binomial',initial,max_iter=100)
            end if
         else
            call robust_glm_fit(x,y,'binomial',initial,max_iter=100)
         end if
         beta=initial%coefficients
      end if
      norm_beta=sqrt(dot_product(beta,beta))
      if(norm_beta<=1.0e-12_dp) then
         beta=0.0_dp
         beta(1)=1.0e-3_dp
         norm_beta=1.0e-3_dp
      end if
      sigma=1.0_dp/norm_beta
      xi=beta*sigma
      stscores=matmul(x,xi)
      obj=objective_scores(stscores/sigma,y,c)
      result%converged=.false.
      do it=1,mi
         sigma=optimize_sigma(stscores,y,c,sigma)
         scores=stscores/sigma
         do i=1,p
            gradient(i)=sum(by_phi_derivative(scores,y,c)*x(:,i))/real(n,dp)
         end do
         direction=-gradient+dot_product(gradient,xi)*xi
         norm_direction=sqrt(dot_product(direction,direction))
         if(norm_direction<=tt) then
            result%converged=.true.
            exit
         end if
         direction=direction/norm_direction
         step=1.0_dp
         accepted=.false.
         do half=0,20
            newxi=xi+step*direction
            newxi=newxi/sqrt(dot_product(newxi,newxi))
            newobj=objective_scores(matmul(x,newxi)/sigma,y,c)
            if(newobj<=obj-1.0e-12_dp) then
               accepted=.true.
               exit
            end if
            step=0.5_dp*step
         end do
         if(.not.accepted) then
            result%converged=.true.
            exit
         end if
         xi=newxi
         stscores=matmul(x,xi)
         if(abs(obj-newobj)<=tt*(1.0_dp+abs(obj))) then
            obj=newobj
            result%converged=.true.
            exit
         end if
         obj=newobj
      end do
      beta=xi/sigma
      allocate(result%coefficients(p),result%covariance(p,p),result%standard_errors(p),result%fitted(n),result%residuals(n))
      result%coefficients=beta
      result%objective=objective_only(x,y,beta,c)
      result%iterations=min(it,mi)
      eta=matmul(x,beta)
      result%fitted=logistic(eta)
      result%residuals=y-result%fitted
      matm=0.0_dp
      ifsqr=0.0_dp
      do i=1,n
         outerm=outer_product(x(i,:),x(i,:))
         matm=matm+by_phi_second_derivative(eta(i),y(i),c)*outerm
         ifsqr=ifsqr+by_phi_derivative(eta(i),y(i),c)**2*outerm
      end do
      matm=matm/real(n,dp)
      ifsqr=ifsqr/real(n,dp)
      call invert_symmetric(matm,minv,info,ridge=1.0e-10_dp)
      if(info==0) then
         result%covariance=matmul(matmul(minv,ifsqr),minv)/real(n,dp)
      else
         result%covariance=0.0_dp
         result%converged=.false.
      end if
      do i=1,p
         result%standard_errors(i)=sqrt(max(result%covariance(i,i),0.0_dp))
      end do
   end subroutine by_logistic_fit

   function objective_scores(scores,y,c) result(value)
      real(dp),intent(in)::scores(:),y(:),c
      real(dp)::value
      integer::i
      value=0.0_dp
      do i=1,size(y)
         value=value+by_phi(scores(i),y(i),c)
      end do
      value=value/real(size(y),dp)
   end function objective_scores

   function optimize_sigma(stscores,y,c,current) result(sigma)
      real(dp),intent(in)::stscores(:),y(:),c,current
      real(dp)::sigma,left,right,x1,x2,f1,f2,ratio
      integer::i
      ratio=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
      left=log(max(current/4.0_dp,1.0e-8_dp))
      right=log(min(current*4.0_dp,1.0e8_dp))
      x1=right-ratio*(right-left)
      x2=left+ratio*(right-left)
      f1=objective_scores(stscores/exp(x1),y,c)
      f2=objective_scores(stscores/exp(x2),y,c)
      do i=1,80
         if(f1<=f2) then
            right=x2
            x2=x1
            f2=f1
            x1=right-ratio*(right-left)
            f1=objective_scores(stscores/exp(x1),y,c)
         else
            left=x1
            x1=x2
            f1=f2
            x2=left+ratio*(right-left)
            f2=objective_scores(stscores/exp(x2),y,c)
         end if
      end do
      sigma=exp(0.5_dp*(left+right))
   end function optimize_sigma

   subroutine pack_rows(x,y,mask,x_out,y_out)
      real(dp),intent(in)::x(:,:),y(:)
      logical,intent(in)::mask(:)
      real(dp),allocatable,intent(out)::x_out(:,:),y_out(:)
      integer::j
      allocate(x_out(count(mask),size(x,2)),y_out(count(mask)))
      do j=1,size(x,2)
         x_out(:,j)=pack(x(:,j),mask)
      end do
      y_out=pack(y,mask)
   end subroutine pack_rows

   function objective_only(x,y,beta,c) result(objective)
      real(dp),intent(in)::x(:,:),y(:),beta(:),c
      real(dp)::objective
      real(dp),allocatable::eta(:)
      integer::i
      eta=matmul(x,beta)
      objective=0.0_dp
      do i=1,size(y)
         objective=objective+by_phi(eta(i),y(i),c)
      end do
      objective=objective/real(size(y),dp)
   end function objective_only

   elemental function by_phi(s,y,c) result(value)
      real(dp),intent(in)::s,y,c
      real(dp)::value,dev
      dev=by_deviance(s,y)
      value=by_rho(dev,c)+gby_positive(s,c)+gby_negative(s,c)
   end function by_phi

   elemental function by_phi_derivative(s,y,c) result(value)
      real(dp),intent(in)::s,y,c
      real(dp)::value,f,ds,dev,g1,g2
      f=logistic_scalar(s)
      ds=f*(1.0_dp-f)
      dev=by_deviance(s,y)
      g1=by_deviance(s,1.0_dp)
      g2=by_deviance(-s,1.0_dp)
      value=-by_psi(dev,c)*(y-f)+ds*(by_psi(g1,c)-by_psi(g2,c))
   end function by_phi_derivative

   elemental function by_phi_second_derivative(s,y,c) result(value)
      real(dp),intent(in)::s,y,c
      real(dp)::value,f,ds,dev,g1,g2
      f=logistic_scalar(s)
      ds=f*(1.0_dp-f)
      dev=by_deviance(s,y)
      g1=by_deviance(s,1.0_dp)
      g2=by_deviance(-s,1.0_dp)
      value=by_psi_derivative(dev,c)*(f-y)**2+ds*by_psi(dev,c)
      value=value+ds*(1.0_dp-2.0_dp*f)*(by_psi(g1,c)-by_psi(g2,c))
      value=value-ds*(by_psi_derivative(g1,c)*(1.0_dp-f)+by_psi_derivative(g2,c)*f)
   end function by_phi_second_derivative

   elemental function by_deviance(s,y) result(value)
      real(dp),intent(in)::s,y
      real(dp)::value,a
      a=abs(s)
      value=log(1.0_dp+exp(-a))
      if((y-0.5_dp)*s<0.0_dp) value=value+a
   end function by_deviance

   elemental function by_rho(t,c) result(value)
      real(dp),intent(in)::t,c
      real(dp)::value,ec,st
      ec=exp(-sqrt(c))
      if(t<=c) then
         value=t*ec
      else
         st=sqrt(max(t,0.0_dp))
         value=ec*(2.0_dp+2.0_dp*sqrt(c)+c)-2.0_dp*exp(-st)*(1.0_dp+st)
      end if
   end function by_rho

   elemental function by_psi(t,c) result(value)
      real(dp),intent(in)::t,c
      real(dp)::value
      if(t<=c) then
         value=exp(-sqrt(c))
      else
         value=exp(-sqrt(max(t,0.0_dp)))
      end if
   end function by_psi

   elemental function by_psi_derivative(t,c) result(value)
      real(dp),intent(in)::t,c
      real(dp)::value,st
      if(t<=c .or. t<=0.0_dp) then
         value=0.0_dp
      else
         st=sqrt(t)
         value=-exp(-st)/(2.0_dp*st)
      end if
   end function by_psi_derivative

   elemental function gby_positive(s,c) result(value)
      real(dp),intent(in)::s,c
      real(dp)::value,f,u,ef,threshold
      ef=exp(0.25_dp)*sqrt(acos(-1.0_dp))
      f=logistic_scalar(s)
      u=-log(max(f,tiny(1.0_dp)))
      threshold=-log(exp(c)-1.0_dp)
      if(s<=threshold) then
         value=ef*(normal_cdf(sqrt(2.0_dp)*(0.5_dp+sqrt(u)))-1.0_dp)+f*exp(-sqrt(u))
      else
         value=f*exp(-sqrt(c))+ef*(normal_cdf(sqrt(2.0_dp)*(0.5_dp+sqrt(c)))-1.0_dp)
      end if
   end function gby_positive

   elemental function gby_negative(s,c) result(value)
      real(dp),intent(in)::s,c
      real(dp)::value,f,u,ef,threshold
      ef=exp(0.25_dp)*sqrt(acos(-1.0_dp))
      f=logistic_scalar(-s)
      u=-log(max(f,tiny(1.0_dp)))
      threshold=log(exp(c)-1.0_dp)
      if(s>=threshold) then
         value=ef*(normal_cdf(sqrt(2.0_dp)*(0.5_dp+sqrt(u)))-1.0_dp)+f*exp(-sqrt(u))
      else
         value=f*exp(-sqrt(c))+ef*(normal_cdf(sqrt(2.0_dp)*(0.5_dp+sqrt(c)))-1.0_dp)
      end if
   end function gby_negative

   elemental function logistic_scalar(x) result(value)
      real(dp),intent(in)::x
      real(dp)::value
      if(x>=0.0_dp) then
         value=1.0_dp/(1.0_dp+exp(-min(x,700.0_dp)))
      else
         value=exp(max(x,-700.0_dp))/(1.0_dp+exp(max(x,-700.0_dp)))
      end if
   end function logistic_scalar

   elemental function logistic(x) result(value)
      real(dp),intent(in)::x
      real(dp)::value
      value=logistic_scalar(x)
   end function logistic

   pure function outer_product(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::j
      do j=1,size(b)
         c(:,j)=a*b(j)
      end do
   end function outer_product
end module robustbase_bylogreg
