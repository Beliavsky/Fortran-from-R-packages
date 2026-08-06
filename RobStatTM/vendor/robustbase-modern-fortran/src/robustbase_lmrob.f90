! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_lmrob
   use robustbase_kinds, only: dp, huge_penalty
   use robustbase_linalg, only: least_squares, invert_symmetric, matrix_rank
   use robustbase_scale, only: mad_scale
   use robustbase_psi, only: tukey_weight, tukey_psi, huber_weight, huber_psi
   implicit none
   private
   public :: lmrob_result, lmrob_s_fit, lmrob_fit, lmrob_lar_fit, robust_m_scale

   type :: lmrob_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: standard_errors(:)
      real(dp), allocatable :: chain_scales(:)
      real(dp) :: scale = 0.0_dp
      real(dp) :: objective = huge_penalty
      integer :: iterations = 0
      integer :: subsamples = 0
      character(len=8) :: method = 'S'
      character(len=16) :: covariance_method = 'sandwich'
      logical :: exact_fit = .false.
      logical :: converged = .false.
   end type lmrob_result
contains
   function robust_m_scale(residuals,tuning,b,target,max_iter,tol) result(scale)
      real(dp),intent(in)::residuals(:)
      real(dp),intent(in),optional::tuning,b,target,tol
      integer,intent(in),optional::max_iter
      real(dp)::scale,c,bb,tar,tt,lo,hi,mid,flo,fhi,fmid,base
      integer::it,mi
      c=1.54764_dp
      if(present(tuning))c=tuning
      bb=0.5_dp
      if(present(b))bb=b
      tar=bb
      if(present(target))tar=target
      mi=200
      if(present(max_iter))mi=max(10,max_iter)
      tt=1.0e-10_dp
      if(present(tol))tt=max(tol,epsilon(1.0_dp))
      if(size(residuals)==0)then
         scale=0.0_dp
         return
      end if
      if(maxval(abs(residuals))<=100.0_dp*tiny(1.0_dp))then
         scale=0.0_dp
         return
      end if
      base=mad_scale(residuals)
      if(base<=1.0e-14_dp)base=sqrt(sum(residuals*residuals)/real(size(residuals),dp))
      if(base<=1.0e-14_dp)base=maxval(abs(residuals))
      lo=max(base*1.0e-8_dp,tiny(1.0_dp))
      hi=max(base,lo*10.0_dp)
      flo=mean_bisquare_rho(residuals/lo,c)-tar
      fhi=mean_bisquare_rho(residuals/hi,c)-tar
      do while(flo<0.0_dp .and. lo>tiny(1.0_dp)*10.0_dp)
         lo=lo*0.1_dp
         flo=mean_bisquare_rho(residuals/lo,c)-tar
      end do
      do while(fhi>0.0_dp .and. hi<huge(1.0_dp)*1.0e-10_dp)
         hi=2.0_dp*hi
         fhi=mean_bisquare_rho(residuals/hi,c)-tar
      end do
      if(flo<0.0_dp)then
         scale=lo
         return
      end if
      if(fhi>0.0_dp)then
         scale=hi
         return
      end if
      do it=1,mi
         mid=0.5_dp*(lo+hi)
         fmid=mean_bisquare_rho(residuals/mid,c)-tar
         if(fmid>0.0_dp)then
            lo=mid
         else
            hi=mid
         end if
         if(abs(hi-lo)<=tt*(1.0_dp+mid))exit
      end do
      scale=0.5_dp*(lo+hi)
   end function robust_m_scale

   subroutine lmrob_s_fit(x,y,result,n_resample,sampling,tuning_chi,bb,max_refine,tol)
      real(dp),intent(in)::x(:,:),y(:)
      type(lmrob_result),intent(out)::result
      integer,intent(in),optional::n_resample,max_refine
      character(len=*),intent(in),optional::sampling
      real(dp),intent(in),optional::tuning_chi,bb,tol
      integer::n,p,nr,mr,trial,info,refit,i,total_comb
      integer,allocatable::subset(:),comb(:),best_subset(:)
      real(dp),allocatable::beta(:),newbeta(:),bestbeta(:),fit(:),res(:),w(:)
      real(dp)::c,b,target_tol,s,best_scale,delta
      character(len=16)::mode
      logical::done,use_exact
      n=size(x,1);p=size(x,2)
      if(size(y)/=n .or. n<p .or. p<1)error stop 'lmrob_s_fit: invalid dimensions'
      if(matrix_rank(x)<p)error stop 'lmrob_s_fit: singular design'
      nr=500
      if(present(n_resample))nr=max(1,n_resample)
      mr=50
      if(present(max_refine))mr=max(1,max_refine)
      c=1.54764_dp
      if(present(tuning_chi))c=tuning_chi
      b=0.5_dp
      if(present(bb))b=bb
      target_tol=1.0e-7_dp
      if(present(tol))target_tol=max(tol,epsilon(1.0_dp))
      mode='nonsingular'
      if(present(sampling))mode=adjustl(sampling)
      total_comb=combination_count(n,p,1000000)
      use_exact=trim(mode)=='exact' .or. (trim(mode)=='best' .and. total_comb<=nr)
      if(trim(mode)/='exact' .and. trim(mode)/='best' .and. trim(mode)/='simple' .and. trim(mode)/='nonsingular') &
         error stop 'lmrob_s_fit: sampling must be exact, best, simple, or nonsingular'
      allocate(subset(p),comb(p),best_subset(p),beta(p),newbeta(p),bestbeta(p),fit(n),res(n),w(n))
      best_scale=huge_penalty
      bestbeta=0.0_dp
      best_subset=0
      result%subsamples=0
      if(use_exact)then
         comb=[(i,i=1,p)]
         done=.false.
         do while(.not.done)
            call evaluate_subset(comb)
            call next_combination(comb,n,done)
         end do
      else
         do trial=1,nr
            call random_subset(n,p,subset)
            if(trim(mode)=='nonsingular')then
               if(matrix_rank(x(subset,:))<p)cycle
            end if
            call evaluate_subset(subset)
         end do
      end if
      if(any(best_subset==0))then
         call least_squares(x,y,bestbeta,info)
         if(info/=0)error stop 'lmrob_s_fit: no nonsingular subsample'
         res=y-matmul(x,bestbeta)
         best_scale=robust_m_scale(res,c,b)
      end if
      allocate(result%coefficients(p),result%fitted(n),result%residuals(n),result%weights(n), &
               result%covariance(p,p),result%standard_errors(p),result%chain_scales(1))
      result%coefficients=bestbeta
      result%fitted=matmul(x,bestbeta)
      result%residuals=y-result%fitted
      result%scale=best_scale
      result%chain_scales=[best_scale]
      if(best_scale<=1.0e-14_dp)then
         result%weights=merge(1.0_dp,0.0_dp,abs(result%residuals)<=1.0e-9_dp)
         result%exact_fit=.true.
      else
         result%weights=tukey_weight(result%residuals/best_scale,c)
      end if
      call regression_covariance(x,result%residuals,result%scale,result%weights,'sandwich','tukey',c,result%covariance,info)
      do i=1,p
         result%standard_errors(i)=sqrt(max(result%covariance(i,i),0.0_dp))
      end do
      result%objective=best_scale
      result%iterations=mr
      result%method='S'
      result%covariance_method='sandwich'
      result%converged=.true.
   contains
      subroutine evaluate_subset(ind)
         integer,intent(in)::ind(:)
         call least_squares(x(ind,:),y(ind),beta,info)
         result%subsamples=result%subsamples+1
         if(info/=0)return
         do refit=1,mr
            res=y-matmul(x,beta)
            s=robust_m_scale(res,c,b)
            if(s<=1.0e-14_dp)exit
            w=tukey_weight(res/s,c)
            if(count(w>1.0e-12_dp)<p)return
            call weighted_least_squares(x,y,w,newbeta,info)
            if(info/=0)return
            delta=maxval(abs(newbeta-beta))
            beta=newbeta
            if(delta<=target_tol*(1.0_dp+maxval(abs(beta))))exit
         end do
         res=y-matmul(x,beta)
         s=robust_m_scale(res,c,b)
         if(s<best_scale)then
            best_scale=s
            bestbeta=beta
            best_subset=ind
         end if
      end subroutine evaluate_subset
   end subroutine lmrob_s_fit

   subroutine lmrob_lar_fit(x,y,result,max_iter,tol,smoothing)
      real(dp),intent(in)::x(:,:),y(:)
      type(lmrob_result),intent(out)::result
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tol,smoothing
      integer::n,p,mi,info,it,j
      real(dp)::tt,sm,scale,epsr,delta,h,f0
      real(dp),allocatable::beta(:),newbeta(:),res(:),w(:),xtx(:,:),xtxi(:,:)
      n=size(x,1);p=size(x,2)
      if(size(y)/=n .or. n<p .or. p<1)error stop 'lmrob_lar_fit: invalid dimensions'
      if(matrix_rank(x)<p)error stop 'lmrob_lar_fit: singular design'
      mi=1000
      if(present(max_iter))mi=max(1,max_iter)
      tt=1.0e-8_dp
      if(present(tol))tt=max(tol,epsilon(1.0_dp))
      sm=1.0e-4_dp
      if(present(smoothing))sm=max(smoothing,sqrt(epsilon(1.0_dp)))
      allocate(beta(p),newbeta(p),res(n),w(n),xtx(p,p),xtxi(p,p))
      call least_squares(x,y,beta,info)
      if(info/=0)error stop 'lmrob_lar_fit: initial least-squares fit failed'
      do it=1,mi
         res=y-matmul(x,beta)
         scale=max(mad_scale(res),sqrt(sum(res*res)/real(max(1,n),dp)))
         epsr=max(sm*max(scale,1.0_dp),sqrt(epsilon(1.0_dp)))
         w=1.0_dp/sqrt(res*res+epsr*epsr)
         call weighted_least_squares(x,y,w,newbeta,info)
         if(info/=0)exit
         delta=maxval(abs(newbeta-beta))
         beta=newbeta
         if(delta<=tt*(1.0_dp+maxval(abs(beta))))exit
      end do
      res=y-matmul(x,beta)
      scale=mad_scale(res)
      if(scale<=1.0e-14_dp)scale=sqrt(sum(res*res)/real(max(1,n-p),dp))
      epsr=max(sm*max(scale,1.0_dp),sqrt(epsilon(1.0_dp)))
      w=1.0_dp/sqrt(res*res+epsr*epsr)
      xtx=matmul(transpose(x),x)
      call invert_symmetric(xtx,xtxi,info,ridge=1.0e-12_dp)
      h=max(1.06_dp*max(scale,epsr)*real(n,dp)**(-0.2_dp),epsr)
      f0=real(count(abs(res)<=h),dp)/(2.0_dp*real(n,dp)*h)
      f0=max(f0,1.0_dp/(10.0_dp*max(scale,epsr)))
      allocate(result%coefficients(p),result%fitted(n),result%residuals(n),result%weights(n), &
               result%covariance(p,p),result%standard_errors(p),result%chain_scales(1))
      result%coefficients=beta
      result%fitted=matmul(x,beta)
      result%residuals=y-result%fitted
      result%weights=w
      result%scale=scale
      result%chain_scales=[scale]
      result%covariance=0.25_dp*xtxi/(f0*f0)
      result%covariance=0.5_dp*(result%covariance+transpose(result%covariance))
      do j=1,p
         result%standard_errors(j)=sqrt(max(result%covariance(j,j),0.0_dp))
      end do
      result%objective=sum(abs(result%residuals))
      result%iterations=min(it,mi)
      result%subsamples=0
      result%method='LAR'
      result%covariance_method='quantile'
      result%exact_fit=maxval(abs(result%residuals))<=1.0e-12_dp
      result%converged=(info==0 .and. it<=mi)
   end subroutine lmrob_lar_fit

   subroutine lmrob_fit(x,y,result,method,n_resample,sampling,tuning_chi,tuning_psi,bb,max_iter,tol,covariance_method)
      real(dp),intent(in)::x(:,:),y(:)
      type(lmrob_result),intent(out)::result
      character(len=*),intent(in),optional::method,sampling,covariance_method
      integer,intent(in),optional::n_resample,max_iter
      real(dp),intent(in),optional::tuning_chi,tuning_psi,bb,tol
      type(lmrob_result)::sfit
      character(len=8)::meth
      character(len=16)::covm
      real(dp)::cchi,cpsi,b,tt,dscale
      integer::mi,info,p,i,total_iterations
      real(dp),allocatable::beta(:),weights(:),residuals(:),scales(:)
      meth='MM'
      if(present(method))meth=adjustl(method)
      covm='sandwich'
      if(present(covariance_method))covm=adjustl(covariance_method)
      cchi=1.54764_dp
      if(present(tuning_chi))cchi=tuning_chi
      cpsi=4.685061_dp
      if(present(tuning_psi))cpsi=tuning_psi
      b=0.5_dp
      if(present(bb))b=bb
      mi=100
      if(present(max_iter))mi=max(1,max_iter)
      tt=1.0e-7_dp
      if(present(tol))tt=max(tol,epsilon(1.0_dp))
      if(present(n_resample) .and. present(sampling))then
         call lmrob_s_fit(x,y,sfit,n_resample=n_resample,sampling=sampling,tuning_chi=cchi,bb=b,max_refine=mi,tol=tt)
      else if(present(n_resample))then
         call lmrob_s_fit(x,y,sfit,n_resample=n_resample,tuning_chi=cchi,bb=b,max_refine=mi,tol=tt)
      else if(present(sampling))then
         call lmrob_s_fit(x,y,sfit,sampling=sampling,tuning_chi=cchi,bb=b,max_refine=mi,tol=tt)
      else
         call lmrob_s_fit(x,y,sfit,tuning_chi=cchi,bb=b,max_refine=mi,tol=tt)
      end if
      p=size(x,2)
      allocate(beta(p),weights(size(y)),residuals(size(y)),scales(4))
      beta=sfit%coefficients
      residuals=y-matmul(x,beta)
      scales=0.0_dp
      scales(1)=sfit%scale
      total_iterations=sfit%iterations
      select case(trim(meth))
      case('S')
         weights=sfit%weights
      case('SM','MM')
         call fixed_scale_m_step(x,y,beta,max(sfit%scale,1.0e-14_dp),cpsi,'tukey',mi,tt,weights,info,i)
         total_iterations=total_iterations+i
         residuals=y-matmul(x,beta)
         scales(2)=sfit%scale
      case('SMDM')
         call fixed_scale_m_step(x,y,beta,max(sfit%scale,1.0e-14_dp),cpsi,'tukey',mi,tt,weights,info,i)
         total_iterations=total_iterations+i
         residuals=y-matmul(x,beta)
         scales(2)=sfit%scale
         dscale=robust_m_scale(residuals,tuning=max(2.5_dp,cchi),b=0.5_dp)
         scales(3)=dscale
         call fixed_scale_m_step(x,y,beta,max(dscale,1.0e-14_dp),cpsi,'tukey',mi,tt,weights,info,i)
         total_iterations=total_iterations+i
         residuals=y-matmul(x,beta)
         scales(4)=dscale
      case default
         error stop 'lmrob_fit: method must be S, SM, MM, or SMDM'
      end select
      allocate(result%coefficients(p),result%fitted(size(y)),result%residuals(size(y)),result%weights(size(y)), &
               result%covariance(p,p),result%standard_errors(p),result%chain_scales(count(scales>0.0_dp)))
      result%coefficients=beta
      result%fitted=matmul(x,beta)
      result%residuals=y-result%fitted
      if(trim(meth)=='S')then
         result%scale=sfit%scale
         result%weights=sfit%weights
      else if(trim(meth)=='SMDM')then
         result%scale=max(scales(4),scales(3))
         result%weights=tukey_weight(result%residuals/max(result%scale,1.0e-14_dp),cpsi)
      else
         result%scale=sfit%scale
         result%weights=tukey_weight(result%residuals/max(result%scale,1.0e-14_dp),cpsi)
      end if
      result%chain_scales=pack(scales,scales>0.0_dp)
      call regression_covariance(x,result%residuals,result%scale,result%weights,covm,'tukey',cpsi,result%covariance,info)
      do i=1,p
         result%standard_errors(i)=sqrt(max(result%covariance(i,i),0.0_dp))
      end do
      result%objective=sum(tukey_rho_normalized(result%residuals/max(result%scale,1.0e-14_dp),cpsi))
      result%iterations=total_iterations
      result%subsamples=sfit%subsamples
      result%method=meth
      result%covariance_method=covm
      result%exact_fit=sfit%exact_fit
      result%converged=(info==0 .and. sfit%converged)
   end subroutine lmrob_fit

   subroutine fixed_scale_m_step(x,y,beta,scale,tuning,psi,max_iter,tol,weights,info,iterations)
      real(dp),intent(in)::x(:,:),y(:),scale,tuning,tol
      real(dp),intent(inout)::beta(:)
      character(len=*),intent(in)::psi
      integer,intent(in)::max_iter
      real(dp),intent(out)::weights(:)
      integer,intent(out)::info,iterations
      real(dp),allocatable::res(:),newbeta(:)
      real(dp)::delta
      allocate(res(size(y)),newbeta(size(beta)))
      info=0
      do iterations=1,max_iter
         res=y-matmul(x,beta)
         select case(trim(psi))
         case('huber')
            weights=huber_weight(res/max(scale,1.0e-14_dp),tuning)
         case default
            weights=tukey_weight(res/max(scale,1.0e-14_dp),tuning)
         end select
         if(count(weights>1.0e-12_dp)<size(beta))then
            info=1
            return
         end if
         call weighted_least_squares(x,y,weights,newbeta,info)
         if(info/=0)return
         delta=maxval(abs(newbeta-beta))
         beta=newbeta
         if(delta<=tol*(1.0_dp+maxval(abs(beta))))return
      end do
   end subroutine fixed_scale_m_step

   subroutine regression_covariance(x,residuals,scale,weights,method,psi,tuning,covariance,info)
      real(dp),intent(in)::x(:,:),residuals(:),scale,weights(:),tuning
      character(len=*),intent(in)::method,psi
      real(dp),intent(out)::covariance(:,:)
      integer,intent(out)::info
      integer::n,p,i
      real(dp),allocatable::a(:,:),b(:,:),ainv(:,:),xx(:,:),u(:),psiv(:),dpsi(:),xtwx(:,:)
      real(dp)::mean_d,mean_p2
      n=size(x,1);p=size(x,2)
      allocate(a(p,p),b(p,p),ainv(p,p),xx(p,p),u(n),psiv(n),dpsi(n),xtwx(p,p))
      u=residuals/max(scale,1.0e-14_dp)
      select case(trim(psi))
      case('huber')
         psiv=huber_psi(u,tuning)
         dpsi=merge(1.0_dp,0.0_dp,abs(u)<tuning)
      case default
         psiv=tukey_psi(u,tuning)
         dpsi=tukey_psi_derivative(u,tuning)
      end select
      select case(trim(method))
      case('weighted','w')
         xtwx=matmul(transpose(x*spread(weights,2,p)),x)
         call invert_symmetric(xtwx,covariance,info,ridge=1.0e-12_dp)
         covariance=scale*scale*covariance*real(n,dp)/real(max(1,n-p),dp)
      case('avar1')
         mean_d=sum(dpsi)/real(n,dp)
         mean_p2=sum(psiv*psiv)/real(n,dp)
         xtwx=matmul(transpose(x),x)
         call invert_symmetric(xtwx,covariance,info,ridge=1.0e-12_dp)
         covariance=scale*scale*mean_p2/max(mean_d*mean_d,1.0e-14_dp)*covariance
      case default
         a=0.0_dp;b=0.0_dp
         do i=1,n
            xx=outer_product(x(i,:),x(i,:))
            a=a+dpsi(i)*xx
            b=b+psiv(i)*psiv(i)*xx
         end do
         a=a/real(n,dp);b=b/real(n,dp)
         call invert_symmetric(a,ainv,info,ridge=1.0e-12_dp)
         covariance=scale*scale*matmul(matmul(ainv,b),ainv)/real(n,dp)
      end select
      covariance=0.5_dp*(covariance+transpose(covariance))
   end subroutine regression_covariance

   subroutine weighted_least_squares(x,y,w,beta,info)
      real(dp),intent(in)::x(:,:),y(:),w(:)
      real(dp),intent(out)::beta(:)
      integer,intent(out)::info
      real(dp),allocatable::a(:,:),b(:)
      integer::j,p
      p=size(x,2)
      allocate(a(size(x,1),p),b(size(y)))
      do j=1,p
         a(:,j)=x(:,j)*sqrt(max(w,0.0_dp))
      end do
      b=y*sqrt(max(w,0.0_dp))
      call least_squares(a,b,beta,info)
   end subroutine weighted_least_squares

   elemental function tukey_rho_normalized(x,c) result(value)
      real(dp),intent(in)::x,c
      real(dp)::value,u
      u=x/c
      if(abs(u)<1.0_dp)then
         value=1.0_dp-(1.0_dp-u*u)**3
      else
         value=1.0_dp
      end if
   end function tukey_rho_normalized

   elemental function tukey_psi_derivative(x,c) result(value)
      real(dp),intent(in)::x,c
      real(dp)::value,u
      u=x/c
      if(abs(u)<1.0_dp)then
         value=(1.0_dp-u*u)*(1.0_dp-5.0_dp*u*u)
      else
         value=0.0_dp
      end if
   end function tukey_psi_derivative

   function mean_bisquare_rho(x,c) result(value)
      real(dp),intent(in)::x(:),c
      real(dp)::value
      value=sum(tukey_rho_normalized(x,c))/real(size(x),dp)
   end function mean_bisquare_rho

   pure function outer_product(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::j
      do j=1,size(b)
         c(:,j)=a*b(j)
      end do
   end function outer_product

   subroutine random_subset(n,k,subset)
      integer,intent(in)::n,k
      integer,intent(out)::subset(:)
      integer::m,candidate
      real(dp)::u
      m=0
      do while(m<k)
         call random_number(u)
         candidate=min(n,1+int(u*real(n,dp)))
         if(m==0 .or. .not.any(subset(1:m)==candidate))then
            m=m+1
            subset(m)=candidate
         end if
      end do
   end subroutine random_subset

   integer function combination_count(n,k,limit) result(value)
      integer,intent(in)::n,k,limit
      integer::i,kk
      real(dp)::v
      kk=min(k,n-k);v=1.0_dp
      do i=1,kk
         v=v*real(n-kk+i,dp)/real(i,dp)
         if(v>real(limit,dp))then
            value=limit+1
            return
         end if
      end do
      value=nint(v)
   end function combination_count

   subroutine next_combination(comb,n,done)
      integer,intent(inout)::comb(:)
      integer,intent(in)::n
      logical,intent(out)::done
      integer::i,j,k
      k=size(comb);i=k
      do while(i>=1)
         if(comb(i)/=n-k+i)exit
         i=i-1
      end do
      if(i==0)then
         done=.true.
         return
      end if
      comb(i)=comb(i)+1
      do j=i+1,k
         comb(j)=comb(j-1)+1
      end do
      done=.false.
   end subroutine next_combination
end module robustbase_lmrob
