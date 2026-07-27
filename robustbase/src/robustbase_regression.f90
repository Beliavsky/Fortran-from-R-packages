! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_regression
   use robustbase_kinds, only: dp, huge_penalty
   use robustbase_sort, only: median, sort_real_with_index
   use robustbase_scale, only: mad_scale
   use robustbase_psi, only: tukey_weight, huber_weight
   use robustbase_linalg, only: least_squares, solve_linear, invert_symmetric
   implicit none
   private
   public :: robust_regression_result, lts_regression, mm_regression, robust_glm_fit
   type :: robust_regression_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp) :: scale = 0.0_dp
      real(dp) :: objective = 0.0_dp
      integer :: iterations = 0
      integer :: h = 0
      logical :: converged = .false.
   end type
contains
   subroutine weighted_ls(x,y,w,beta,info)
      real(dp),intent(in)::x(:,:),y(:),w(:)
      real(dp),intent(out)::beta(:)
      integer,intent(out)::info
      real(dp),allocatable::a(:,:),b(:)
      integer::j,n,p
      n=size(x,1);p=size(x,2)
      allocate(a(n,p),b(n))
      do j=1,p;a(:,j)=x(:,j)*sqrt(max(w,0.0_dp));end do
      b=y*sqrt(max(w,0.0_dp))
      call least_squares(a,b,beta,info)
   end subroutine weighted_ls

   subroutine lts_regression(x,y,result,alpha,n_starts,max_csteps,reweight_cutoff)
      real(dp),intent(in)::x(:,:),y(:)
      type(robust_regression_result),intent(out)::result
      real(dp),intent(in),optional::alpha,reweight_cutoff
      integer,intent(in),optional::n_starts,max_csteps
      integer::n,p,h,ns,mc,s,it,info,i,nw
      integer,allocatable::sub(:),idx(:),bestsub(:)
      real(dp),allocatable::beta(:),bestbeta(:),r(:),rs(:),w(:),xtwx(:,:),inv(:,:)
      real(dp)::a,bestobj,obj,cut
      n=size(x,1);p=size(x,2)
      if(size(y)/=n .or. n<p) error stop "lts_regression: invalid dimensions"
      a=0.75_dp;if(present(alpha))a=alpha;h=max(p,min(n,int(floor(a*real(n,dp)))));ns=200;if(present(n_starts))ns=n_starts;mc=30;if(present(max_csteps))mc=max_csteps
      allocate(sub(h),idx(n),bestsub(h),beta(p),bestbeta(p),r(n),rs(n),w(n),xtwx(p,p),inv(p,p))
      bestobj=huge_penalty;bestbeta=0.0_dp;bestsub=0
      do s=1,ns
         call random_subset(n,p,sub(1:p))
         call least_squares(x(sub(1:p),:),y(sub(1:p)),beta,info);if(info/=0)cycle
         do it=1,mc
            r=y-matmul(x,beta);rs=r*r;call sort_real_with_index(rs,idx);sub=idx(1:h)
            call least_squares(x(sub,:),y(sub),beta,info);if(info/=0)exit
         end do
         r=y-matmul(x,beta);rs=r*r;call sort_real_with_index(rs,idx);obj=sum(rs(1:h))
         if(obj<bestobj) then;bestobj=obj;bestbeta=beta;bestsub=idx(1:h);end if
      end do
      if(any(bestsub==0)) then;bestsub=[(i,i=1,h)];call least_squares(x(bestsub,:),y(bestsub),bestbeta,info);end if
      allocate(result%coefficients(p),result%residuals(n),result%fitted(n),result%weights(n),result%covariance(p,p))
      result%coefficients=bestbeta;result%fitted=matmul(x,bestbeta);result%residuals=y-result%fitted
      result%scale=1.482602218505602_dp*median(abs(result%residuals(bestsub)))
      if(result%scale<=1.0e-14_dp) result%scale=sqrt(bestobj/real(max(1,h-p),dp))
      cut=2.5_dp;if(present(reweight_cutoff))cut=reweight_cutoff
      if(result%scale>0.0_dp) then;result%weights=merge(1.0_dp,0.0_dp,abs(result%residuals)<=cut*result%scale);else;result%weights=1.0_dp;end if
      nw=count(result%weights>0.0_dp)
      if(nw>=p) then
         call weighted_ls(x,y,result%weights,result%coefficients,info)
         result%fitted=matmul(x,result%coefficients);result%residuals=y-result%fitted
      end if
      xtwx=matmul(transpose(x*spread(result%weights,2,p)),x)
      call invert_symmetric(xtwx,inv,info,ridge=1.0e-12_dp)
      result%covariance=result%scale**2*inv
      result%objective=bestobj;result%iterations=mc;result%h=h;result%converged=info==0
   contains
      subroutine random_subset(nn,kk,out)
         integer,intent(in)::nn,kk;integer,intent(out)::out(:)
         integer::ii,cand;real(dp)::u
         ii=0
         do while(ii<kk)
            call random_number(u);cand=min(nn,1+int(u*real(nn,dp)))
            if(ii==0 .or. .not.any(out(1:ii)==cand)) then;ii=ii+1;out(ii)=cand;end if
         end do
      end subroutine
   end subroutine lts_regression

   subroutine mm_regression(x,y,result,psi,tuning,alpha,n_starts,max_iter,tol)
      real(dp),intent(in)::x(:,:),y(:)
      type(robust_regression_result),intent(out)::result
      character(len=*),intent(in),optional::psi
      real(dp),intent(in),optional::tuning,alpha,tol
      integer,intent(in),optional::n_starts,max_iter
      type(robust_regression_result)::initial
      character(len=16)::ps
      real(dp)::c,tt,delta
      real(dp),allocatable::beta(:),newbeta(:),r(:),w(:),xtwx(:,:),inv(:,:)
      integer::n,p,it,mi,info,ns
      n=size(x,1);p=size(x,2);ps='tukey';if(present(psi))ps=adjustl(psi);c=4.685_dp;if(present(tuning))c=tuning;tt=1.0e-7_dp;if(present(tol))tt=tol;mi=100;if(present(max_iter))mi=max_iter;ns=100;if(present(n_starts))ns=n_starts
      if(present(alpha)) then;call lts_regression(x,y,initial,alpha=alpha,n_starts=ns);else;call lts_regression(x,y,initial,n_starts=ns);end if
      allocate(beta(p),newbeta(p),r(n),w(n),xtwx(p,p),inv(p,p));beta=initial%coefficients
      result%scale=max(initial%scale,mad_scale(y-matmul(x,beta)))
      do it=1,mi
         r=y-matmul(x,beta)
         select case(trim(ps))
         case('huber');w=huber_weight(r/max(result%scale,1.0e-14_dp),c)
         case default;w=tukey_weight(r/max(result%scale,1.0e-14_dp),c)
         end select
         call weighted_ls(x,y,w,newbeta,info);if(info/=0)exit
         delta=maxval(abs(newbeta-beta));beta=newbeta
         if(delta<=tt*(1.0_dp+maxval(abs(beta))))exit
      end do
      allocate(result%coefficients(p),result%residuals(n),result%fitted(n),result%weights(n),result%covariance(p,p))
      result%coefficients=beta;result%fitted=matmul(x,beta);result%residuals=y-result%fitted;result%weights=w
      xtwx=matmul(transpose(x*spread(w,2,p)),x);call invert_symmetric(xtwx,inv,info,ridge=1.0e-12_dp);result%covariance=result%scale**2*inv
      result%objective=sum(w*result%residuals**2);result%iterations=it;result%h=initial%h;result%converged=(info==0 .and. it<=mi)
   end subroutine mm_regression

   subroutine robust_glm_fit(x,y,family,result,tuning,max_iter,tol)
      real(dp),intent(in)::x(:,:),y(:)
      character(len=*),intent(in)::family
      type(robust_regression_result),intent(out)::result
      real(dp),intent(in),optional::tuning,tol
      integer,intent(in),optional::max_iter
      integer::n,p,it,mi,info
      real(dp)::c,tt,delta
      real(dp),allocatable::beta(:),newbeta(:),eta(:),mu(:),var(:),dmu(:),z(:),pearson(:),rw(:),w(:),xtwx(:,:),inv(:,:)
      n=size(x,1);p=size(x,2);c=1.345_dp;if(present(tuning))c=tuning;tt=1.0e-7_dp;if(present(tol))tt=tol;mi=100;if(present(max_iter))mi=max_iter
      allocate(beta(p),newbeta(p),eta(n),mu(n),var(n),dmu(n),z(n),pearson(n),rw(n),w(n),xtwx(p,p),inv(p,p));beta=0.0_dp
      if(trim(family)=='binomial') beta(1)=log(max(1.0e-6_dp,min(1.0_dp-1.0e-6_dp,sum(y)/real(n,dp)))/(1.0_dp-max(1.0e-6_dp,min(1.0_dp-1.0e-6_dp,sum(y)/real(n,dp)))))
      if(trim(family)=='poisson') beta(1)=log(max(sum(y)/real(n,dp),1.0e-6_dp))
      do it=1,mi
         eta=matmul(x,beta)
         select case(trim(family))
         case('binomial')
            mu=1.0_dp/(1.0_dp+exp(-max(-35.0_dp,min(35.0_dp,eta))));var=max(mu*(1.0_dp-mu),1.0e-10_dp);dmu=var
         case('poisson')
            mu=exp(max(-30.0_dp,min(30.0_dp,eta)));var=max(mu,1.0e-10_dp);dmu=mu
         case default
            error stop "robust_glm_fit: family must be binomial or poisson"
         end select
         pearson=(y-mu)/sqrt(var);rw=huber_weight(pearson,c);z=eta+(y-mu)/max(dmu,1.0e-10_dp);w=rw*dmu*dmu/var
         call weighted_ls(x,z,w,newbeta,info);if(info/=0)exit
         delta=maxval(abs(newbeta-beta));beta=newbeta
         if(delta<=tt*(1.0_dp+maxval(abs(beta))))exit
      end do
      allocate(result%coefficients(p),result%residuals(n),result%fitted(n),result%weights(n),result%covariance(p,p));result%coefficients=beta
      eta=matmul(x,beta)
      if(trim(family)=='binomial') then;mu=1.0_dp/(1.0_dp+exp(-max(-35.0_dp,min(35.0_dp,eta))));else;mu=exp(max(-30.0_dp,min(30.0_dp,eta)));end if
      result%fitted=mu;result%residuals=y-mu;result%weights=rw;result%scale=1.0_dp;xtwx=matmul(transpose(x*spread(w,2,p)),x);call invert_symmetric(xtwx,inv,info,ridge=1.0e-12_dp);result%covariance=inv
      result%objective=sum(rw*result%residuals**2);result%iterations=it;result%h=n;result%converged=(info==0 .and. it<=mi)
   end subroutine robust_glm_fit
end module robustbase_regression
