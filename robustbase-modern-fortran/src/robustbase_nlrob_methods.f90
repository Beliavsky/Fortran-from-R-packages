! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_nlrob_methods
   use robustbase_kinds, only: dp, huge_penalty
   use robustbase_nonlinear, only: nonlinear_model
   use robustbase_lmrob, only: robust_m_scale
   use robustbase_psi, only: tukey_weight, tukey_psi
   use robustbase_linalg, only: least_squares, invert_symmetric
   use robustbase_sort, only: sort_real
   implicit none
   private
   public :: nlrob_method_result, nlrob_mm_fit, nlrob_tau_fit, nlrob_cm_fit, nlrob_mtl_fit

   type :: nlrob_method_result
      real(dp),allocatable :: parameters(:)
      real(dp),allocatable :: fitted(:)
      real(dp),allocatable :: residuals(:)
      real(dp),allocatable :: weights(:)
      real(dp),allocatable :: covariance(:,:)
      real(dp),allocatable :: standard_errors(:)
      real(dp) :: scale=0.0_dp
      real(dp) :: tau_scale=0.0_dp
      real(dp) :: objective=huge_penalty
      integer :: iterations=0
      integer :: evaluations=0
      character(len=8) :: method='MM'
      logical :: converged=.false.
   end type nlrob_method_result
contains
   subroutine nlrob_mm_fit(model,x,y,start,result,lower,upper,n_starts,tuning_chi,tuning_psi,max_iter,tol)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),start(:)
      type(nlrob_method_result),intent(out)::result
      real(dp),intent(in),optional::lower(:),upper(:),tuning_chi,tuning_psi,tol
      integer,intent(in),optional::n_starts,max_iter
      integer::n,p,ns,mi,s,it,info,i,best_it
      real(dp)::cchi,cpsi,tt,scale,best_scale,u
      real(dp),allocatable::theta(:),candidate(:),best(:),fit(:),res(:),w(:),jac(:,:),step(:),lo(:),hi(:)
      n=size(y);p=size(start)
      if(size(x,1)/=n .or. p<1)error stop 'nlrob_mm_fit: invalid dimensions'
      ns=40;if(present(n_starts))ns=max(1,n_starts)
      mi=100;if(present(max_iter))mi=max(1,max_iter)
      cchi=1.54764_dp;if(present(tuning_chi))cchi=tuning_chi
      cpsi=4.685061_dp;if(present(tuning_psi))cpsi=tuning_psi
      tt=1.0e-7_dp;if(present(tol))tt=max(tol,epsilon(1.0_dp))
      allocate(theta(p),candidate(p),best(p),fit(n),res(n),w(n),jac(n,p),step(p),lo(p),hi(p))
      call set_bounds(start,lower,upper,lo,hi)
      best_scale=huge_penalty;best=start;best_it=0
      do s=1,ns
         if(s==1)then
            theta=min(hi,max(lo,start))
         else
            call random_number(candidate)
            do i=1,p
               if(is_finite_bound(lo(i),hi(i)))then
                  theta(i)=lo(i)+candidate(i)*(hi(i)-lo(i))
               else
                  call random_number(u)
                  theta(i)=start(i)+(2.0_dp*u-1.0_dp)*(1.0_dp+abs(start(i)))
                  theta(i)=min(hi(i),max(lo(i),theta(i)))
               end if
            end do
         end if
         call s_refine(model,x,y,theta,lo,hi,cchi,mi,tt,scale,it,info)
         result%evaluations=result%evaluations+it*(p+1)
         if(info==0 .and. scale<best_scale)then
            best_scale=scale;best=theta;best_it=it
         end if
      end do
      if(best_scale>=huge_penalty*0.5_dp)error stop 'nlrob_mm_fit: no valid start'
      theta=best
      call m_refine(model,x,y,theta,lo,hi,max(best_scale,1.0e-14_dp),cpsi,mi,tt,it,info)
      result%evaluations=result%evaluations+it*(p+1)
      call fill_result(model,x,y,theta,max(best_scale,1.0e-14_dp),cpsi,'MM',best_it+it,info==0,result)
   end subroutine nlrob_mm_fit

   subroutine nlrob_tau_fit(model,x,y,start,lower,upper,result,tuning_scale,tuning_tau,b_scale,b_tau,max_iter,tol)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),start(:),lower(:),upper(:)
      type(nlrob_method_result),intent(out)::result
      real(dp),intent(in),optional::tuning_scale,tuning_tau,b_scale,b_tau,tol
      integer,intent(in),optional::max_iter
      real(dp)::c1,c2,b1,b2,tt,scale,tau
      real(dp),allocatable::theta(:)
      integer::mi,it,evals
      logical::conv
      c1=1.55_dp;if(present(tuning_scale))c1=tuning_scale
      c2=6.04_dp;if(present(tuning_tau))c2=tuning_tau
      b1=0.20_dp;if(present(b_scale))b1=b_scale
      b2=0.46_dp;if(present(b_tau))b2=b_tau
      mi=500;if(present(max_iter))mi=max(20,max_iter)
      tt=1.0e-7_dp;if(present(tol))tt=max(tol,epsilon(1.0_dp))
      theta=start
      call nelder_mead_model(model,x,y,theta,lower,upper,'tau',c1,c2,b1,b2,0.0_dp,0,mi,tt,it,evals,conv)
      call criterion_stats(model,x,y,theta,'tau',c1,c2,b1,b2,0.0_dp,0,result%objective,scale,tau)
      call fill_result(model,x,y,theta,max(scale,1.0e-14_dp),c2,'tau',it,conv,result)
      result%tau_scale=tau;result%evaluations=evals;result%objective=tau*tau
   end subroutine nlrob_tau_fit

   subroutine nlrob_cm_fit(model,x,y,start,lower,upper,result,tuning_scale,tuning_m,b_scale,lambda,max_iter,tol)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),start(:),lower(:),upper(:)
      type(nlrob_method_result),intent(out)::result
      real(dp),intent(in),optional::tuning_scale,tuning_m,b_scale,lambda,tol
      integer,intent(in),optional::max_iter
      real(dp)::c1,c2,b1,lam,tt,scale,tau,obj
      real(dp),allocatable::theta(:)
      integer::mi,it,evals
      logical::conv
      c1=1.54764_dp;if(present(tuning_scale))c1=tuning_scale
      c2=4.685061_dp;if(present(tuning_m))c2=tuning_m
      b1=0.5_dp;if(present(b_scale))b1=b_scale
      lam=1.0_dp;if(present(lambda))lam=lambda
      mi=500;if(present(max_iter))mi=max(20,max_iter)
      tt=1.0e-7_dp;if(present(tol))tt=max(tol,epsilon(1.0_dp))
      theta=start
      call nelder_mead_model(model,x,y,theta,lower,upper,'CM',c1,c2,b1,1.0_dp,lam,0,mi,tt,it,evals,conv)
      call criterion_stats(model,x,y,theta,'CM',c1,c2,b1,1.0_dp,lam,0,obj,scale,tau)
      call fill_result(model,x,y,theta,max(scale,1.0e-14_dp),c2,'CM',it,conv,result)
      result%evaluations=evals;result%objective=obj
   end subroutine nlrob_cm_fit

   subroutine nlrob_mtl_fit(model,x,y,start,lower,upper,result,alpha,max_iter,tol)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),start(:),lower(:),upper(:)
      type(nlrob_method_result),intent(out)::result
      real(dp),intent(in),optional::alpha,tol
      integer,intent(in),optional::max_iter
      real(dp)::a,tt,scale,tau,obj
      real(dp),allocatable::theta(:)
      integer::mi,it,evals,h
      logical::conv
      a=0.75_dp;if(present(alpha))a=alpha
      if(a<=0.5_dp .or. a>1.0_dp)error stop 'nlrob_mtl_fit: alpha must be in (0.5,1]'
      h=max(size(start)+1,min(size(y),int(floor(a*real(size(y),dp)))))
      mi=500;if(present(max_iter))mi=max(20,max_iter)
      tt=1.0e-7_dp;if(present(tol))tt=max(tol,epsilon(1.0_dp))
      theta=start
      call nelder_mead_model(model,x,y,theta,lower,upper,'MTL',1.54764_dp,4.685061_dp,0.5_dp,1.0_dp,0.0_dp,h,mi,tt,it,evals,conv)
      call criterion_stats(model,x,y,theta,'MTL',1.54764_dp,4.685061_dp,0.5_dp,1.0_dp,0.0_dp,h,obj,scale,tau)
      call fill_result(model,x,y,theta,max(scale,1.0e-14_dp),4.685061_dp,'MTL',it,conv,result)
      result%evaluations=evals;result%objective=obj
   end subroutine nlrob_mtl_fit

   subroutine s_refine(model,x,y,theta,lower,upper,c,max_iter,tol,scale,iterations,info)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),lower(:),upper(:),c,tol
      real(dp),intent(inout)::theta(:)
      integer,intent(in)::max_iter
      real(dp),intent(out)::scale
      integer,intent(out)::iterations,info
      real(dp),allocatable::fit(:),res(:),w(:),jac(:,:),step(:)
      real(dp)::delta
      allocate(fit(size(y)),res(size(y)),w(size(y)),jac(size(y),size(theta)),step(size(theta)))
      info=0
      do iterations=1,max_iter
         call model(theta,x,fit);res=y-fit;scale=robust_m_scale(res,tuning=c,b=0.5_dp)
         if(scale<=1.0e-14_dp)return
         w=tukey_weight(res/scale,c)
         call finite_jacobian(model,theta,x,fit,jac)
         call weighted_step(jac,res,w,step,info)
         if(info/=0)return
         delta=maxval(abs(step));theta=min(upper,max(lower,theta+step))
         if(delta<=tol*(1.0_dp+maxval(abs(theta))))return
      end do
   end subroutine s_refine

   subroutine m_refine(model,x,y,theta,lower,upper,scale,c,max_iter,tol,iterations,info)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),lower(:),upper(:),scale,c,tol
      real(dp),intent(inout)::theta(:)
      integer,intent(in)::max_iter
      integer,intent(out)::iterations,info
      real(dp),allocatable::fit(:),res(:),w(:),jac(:,:),step(:)
      real(dp)::delta
      allocate(fit(size(y)),res(size(y)),w(size(y)),jac(size(y),size(theta)),step(size(theta)))
      info=0
      do iterations=1,max_iter
         call model(theta,x,fit);res=y-fit;w=tukey_weight(res/scale,c)
         call finite_jacobian(model,theta,x,fit,jac)
         call weighted_step(jac,res,w,step,info)
         if(info/=0)return
         delta=maxval(abs(step));theta=min(upper,max(lower,theta+step))
         if(delta<=tol*(1.0_dp+maxval(abs(theta))))return
      end do
   end subroutine m_refine

   subroutine fill_result(model,x,y,theta,scale,c,method,iterations,converged,result)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),theta(:),scale,c
      character(len=*),intent(in)::method
      integer,intent(in)::iterations
      logical,intent(in)::converged
      type(nlrob_method_result),intent(inout)::result
      real(dp),allocatable::fit(:),res(:),w(:),jac(:,:),a(:,:),b(:,:),ainv(:,:),xx(:,:)
      integer::n,p,i,info
      n=size(y);p=size(theta)
      allocate(fit(n),res(n),w(n),jac(n,p),a(p,p),b(p,p),ainv(p,p),xx(p,p))
      call model(theta,x,fit);res=y-fit;w=tukey_weight(res/max(scale,1.0e-14_dp),c)
      call finite_jacobian(model,theta,x,fit,jac)
      a=0.0_dp;b=0.0_dp
      do i=1,n
         xx=outer_product(jac(i,:),jac(i,:))
         a=a+tukey_derivative(res(i)/max(scale,1.0e-14_dp),c)*xx
         b=b+tukey_psi(res(i)/max(scale,1.0e-14_dp),c)**2*xx
      end do
      a=a/real(n,dp);b=b/real(n,dp)
      call invert_symmetric(a,ainv,info,ridge=1.0e-10_dp)
      allocate(result%parameters(p),result%fitted(n),result%residuals(n),result%weights(n),result%covariance(p,p),result%standard_errors(p))
      result%parameters=theta;result%fitted=fit;result%residuals=res;result%weights=w;result%scale=scale
      result%covariance=scale*scale*matmul(matmul(ainv,b),ainv)/real(n,dp)
      result%covariance=0.5_dp*(result%covariance+transpose(result%covariance))
      do i=1,p
         result%standard_errors(i)=sqrt(max(result%covariance(i,i),0.0_dp))
      end do
      result%iterations=iterations;result%method=method;result%converged=converged .and. info==0
      if(result%objective>=huge_penalty*0.5_dp)result%objective=sum(tukey_rho_norm(res/max(scale,1.0e-14_dp),c))
   end subroutine fill_result

   subroutine nelder_mead_model(model,x,y,theta,lower,upper,criterion,c1,c2,b1,b2,lambda,h,max_iter,tol,iterations,evaluations,converged)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),lower(:),upper(:),c1,c2,b1,b2,lambda,tol
      real(dp),intent(inout)::theta(:)
      character(len=*),intent(in)::criterion
      integer,intent(in)::h,max_iter
      integer,intent(out)::iterations,evaluations
      logical,intent(out)::converged
      integer::p,j,ilo,ihi,inhi
      real(dp),allocatable::simplex(:,:),value(:),centroid(:),xr(:),xe(:),xc(:)
      real(dp)::fr,fe,fc,spreadv,diameter,step
      p=size(theta)
      if(size(lower)/=p .or. size(upper)/=p .or. any(lower>upper))error stop 'nelder_mead_model: invalid bounds'
      allocate(simplex(p,p+1),value(p+1),centroid(p),xr(p),xe(p),xc(p))
      simplex(:,1)=min(upper,max(lower,theta))
      do j=1,p
         simplex(:,j+1)=simplex(:,1)
         step=0.05_dp*(1.0_dp+abs(theta(j)))
         if(upper(j)<huge(1.0_dp)*0.1_dp .and. lower(j)>-huge(1.0_dp)*0.1_dp)step=max(step,0.05_dp*(upper(j)-lower(j)))
         simplex(j,j+1)=min(upper(j),max(lower(j),simplex(j,j+1)+step))
         if(abs(simplex(j,j+1)-simplex(j,1))<=epsilon(1.0_dp)*(1.0_dp+abs(simplex(j,1)))) &
            simplex(j,j+1)=min(upper(j),max(lower(j),simplex(j,j+1)-2.0_dp*step))
      end do
      evaluations=0
      do j=1,p+1
         value(j)=evaluate(simplex(:,j));evaluations=evaluations+1
      end do
      converged=.false.
      do iterations=1,max_iter
         call order_vertices(value,ilo,ihi,inhi)
         spreadv=maxval(abs(value-value(ilo)))
         diameter=0.0_dp
         do j=1,p+1
            diameter=max(diameter,maxval(abs(simplex(:,j)-simplex(:,ilo))))
         end do
         if(spreadv<=tol*(1.0_dp+abs(value(ilo))) .and. diameter<=sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,ilo)))))then
            converged=.true.;exit
         end if
         centroid=0.0_dp
         do j=1,p+1
            if(j/=ihi)centroid=centroid+simplex(:,j)
         end do
         centroid=centroid/real(p,dp)
         xr=min(upper,max(lower,centroid+(centroid-simplex(:,ihi))))
         fr=evaluate(xr);evaluations=evaluations+1
         if(fr<value(ilo))then
            xe=min(upper,max(lower,centroid+2.0_dp*(xr-centroid)))
            fe=evaluate(xe);evaluations=evaluations+1
            if(fe<fr)then;simplex(:,ihi)=xe;value(ihi)=fe;else;simplex(:,ihi)=xr;value(ihi)=fr;end if
         else if(fr<value(inhi))then
            simplex(:,ihi)=xr;value(ihi)=fr
         else
            if(fr<value(ihi))then
               xc=min(upper,max(lower,centroid+0.5_dp*(xr-centroid)))
            else
               xc=min(upper,max(lower,centroid+0.5_dp*(simplex(:,ihi)-centroid)))
            end if
            fc=evaluate(xc);evaluations=evaluations+1
            if(fc<min(fr,value(ihi)))then
               simplex(:,ihi)=xc;value(ihi)=fc
            else
               do j=1,p+1
                  if(j/=ilo)then
                     simplex(:,j)=min(upper,max(lower,simplex(:,ilo)+0.5_dp*(simplex(:,j)-simplex(:,ilo))))
                     value(j)=evaluate(simplex(:,j));evaluations=evaluations+1
                  end if
               end do
            end if
         end if
      end do
      call order_vertices(value,ilo,ihi,inhi)
      theta=simplex(:,ilo)
   contains
      function evaluate(par) result(value_out)
         real(dp),intent(in)::par(:)
         real(dp)::value_out,scl,tau
         call criterion_stats(model,x,y,par,criterion,c1,c2,b1,b2,lambda,h,value_out,scl,tau)
      end function evaluate
   end subroutine nelder_mead_model

   subroutine criterion_stats(model,x,y,theta,criterion,c1,c2,b1,b2,lambda,h,objective,scale,tau)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::x(:,:),y(:),theta(:),c1,c2,b1,b2,lambda
      character(len=*),intent(in)::criterion
      integer,intent(in)::h
      real(dp),intent(out)::objective,scale,tau
      real(dp),allocatable::fit(:),res(:),sorted(:)
      allocate(fit(size(y)),res(size(y)))
      call model(theta,x,fit)
      if(any(.not.(fit<huge(1.0_dp) .and. fit>-huge(1.0_dp))))then
         objective=huge_penalty;scale=huge_penalty;tau=huge_penalty;return
      end if
      res=y-fit
      scale=robust_m_scale(res,tuning=c1,b=b1)
      if(scale<=1.0e-14_dp)then
         objective=0.0_dp;tau=0.0_dp;return
      end if
      select case(trim(criterion))
      case('tau')
         tau=scale*sqrt(sum(tukey_rho_norm(res/scale,c2))/real(size(y),dp)/max(b2,1.0e-12_dp))
         objective=tau*tau
      case('CM')
         tau=scale
         objective=sum(tukey_rho_norm(res/scale,c2))/real(size(y),dp)+lambda*log(max(scale,1.0e-14_dp))
      case('MTL')
         sorted=res*res;call sort_real(sorted)
         objective=sum(sorted(1:h));scale=sqrt(objective/real(max(1,h-size(theta)),dp));tau=scale
      case default
         objective=huge_penalty;tau=scale
      end select
   end subroutine criterion_stats

   subroutine finite_jacobian(model,theta,x,fit,jac)
      procedure(nonlinear_model)::model
      real(dp),intent(in)::theta(:),x(:,:),fit(:)
      real(dp),intent(out)::jac(:,:)
      real(dp),allocatable::tp(:),fp(:)
      real(dp)::step
      integer::j
      allocate(tp(size(theta)),fp(size(fit)))
      do j=1,size(theta)
         tp=theta;step=sqrt(epsilon(1.0_dp))*(1.0_dp+abs(theta(j)));tp(j)=tp(j)+step
         call model(tp,x,fp);jac(:,j)=(fp-fit)/step
      end do
   end subroutine finite_jacobian

   subroutine weighted_step(jac,res,w,step,info)
      real(dp),intent(in)::jac(:,:),res(:),w(:)
      real(dp),intent(out)::step(:)
      integer,intent(out)::info
      real(dp),allocatable::a(:,:),b(:)
      integer::j
      allocate(a(size(jac,1),size(jac,2)),b(size(res)))
      do j=1,size(jac,2)
         a(:,j)=jac(:,j)*sqrt(max(w,0.0_dp))
      end do
      b=res*sqrt(max(w,0.0_dp))
      call least_squares(a,b,step,info)
   end subroutine weighted_step

   subroutine set_bounds(start,lower,upper,lo,hi)
      real(dp),intent(in)::start(:)
      real(dp),intent(in),optional::lower(:),upper(:)
      real(dp),intent(out)::lo(:),hi(:)
      lo=-huge(1.0_dp);hi=huge(1.0_dp)
      if(present(lower))then
         if(size(lower)/=size(start))error stop 'set_bounds: lower size mismatch'
         lo=lower
      end if
      if(present(upper))then
         if(size(upper)/=size(start))error stop 'set_bounds: upper size mismatch'
         hi=upper
      end if
      if(any(lo>hi))error stop 'set_bounds: invalid bounds'
   end subroutine set_bounds

   logical function is_finite_bound(lo,hi) result(ok)
      real(dp),intent(in)::lo,hi
      ok=lo>-huge(1.0_dp)*0.1_dp .and. hi<huge(1.0_dp)*0.1_dp
   end function is_finite_bound

   subroutine order_vertices(value,ilo,ihi,inhi)
      real(dp),intent(in)::value(:)
      integer,intent(out)::ilo,ihi,inhi
      integer::j
      ilo=minloc(value,dim=1);ihi=maxloc(value,dim=1);inhi=ilo
      do j=1,size(value)
         if(j/=ihi)then
            if(inhi==ihi .or. value(j)>value(inhi))inhi=j
         end if
      end do
   end subroutine order_vertices

   elemental function tukey_rho_norm(x,c) result(value)
      real(dp),intent(in)::x,c
      real(dp)::value,u
      u=x/c
      if(abs(u)<1.0_dp)then
         value=1.0_dp-(1.0_dp-u*u)**3
      else
         value=1.0_dp
      end if
   end function tukey_rho_norm

   elemental function tukey_derivative(x,c) result(value)
      real(dp),intent(in)::x,c
      real(dp)::value,u
      u=x/c
      if(abs(u)<1.0_dp)then
         value=(1.0_dp-u*u)*(1.0_dp-5.0_dp*u*u)
      else
         value=0.0_dp
      end if
   end function tukey_derivative

   pure function outer_product(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::j
      do j=1,size(b)
         c(:,j)=a*b(j)
      end do
   end function outer_product

end module robustbase_nlrob_methods
