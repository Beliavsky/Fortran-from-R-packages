! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_fast_algorithms
   use robustbase_kinds, only: dp, huge_penalty
   use robustbase_covariance, only: robust_cov_result
   use robustbase_detmcd, only: detmcd_result, cov_detmcd
   use robustbase_fastlts, only: fast_lts_result, fast_lts_regression, h_alpha_n
   use robustbase_linalg, only: invert_symmetric, symmetric_eigen, covariance_matrix, least_squares, solve_linear
   use robustbase_sort, only: sort_real_with_index
   use robustbase_probability, only: chi_square_quantile, chi_square_cdf, normal_quantile
   implicit none
   private
   public :: partitioned_mcd_result, partitioned_lts_result, fast_mcd_partitioned, fast_lts_partitioned, &
             mcd_consistency_factor, mcd_finite_sample_factor, mcd_reweighted_finite_sample_factor

   type :: partitioned_mcd_result
      type(robust_cov_result) :: estimate
      integer,allocatable :: best_subset(:)
      integer :: partitions=0
      integer :: candidates=0
      integer :: csteps=0
      real(dp) :: raw_consistency_factor=1.0_dp
      real(dp) :: reweight_consistency_factor=1.0_dp
      logical :: exact_fit=.false.
      logical :: converged=.false.
   end type partitioned_mcd_result

   type :: partitioned_lts_result
      type(fast_lts_result) :: estimate
      integer :: partitions=0
      integer :: candidates=0
      logical :: converged=.false.
   end type partitioned_lts_result
contains

   function mcd_consistency_factor(p,alpha) result(factor)
      integer,intent(in)::p
      real(dp),intent(in)::alpha
      real(dp)::factor,q
      if(p<1 .or. alpha<=0.0_dp .or. alpha>1.0_dp)error stop 'mcd_consistency_factor: invalid arguments'
      q=chi_square_quantile(alpha,real(p,dp))
      factor=alpha/max(chi_square_cdf(q,real(p+2,dp)),1.0e-14_dp)
   end function mcd_consistency_factor

   function mcd_finite_sample_factor(p,n,alpha) result(factor)
      integer,intent(in)::p,n
      real(dp),intent(in)::alpha
      real(dp)::factor,fp500,fp875,fpalpha
      real(dp)::coeff500(2,2),coeff875(2,2),y500(2),y875(2),a500(2,2),a875(2,2),sol500(2),sol875(2)
      integer::info
      if(p<1 .or. n<2 .or. alpha<0.5_dp .or. alpha>1.0_dp)error stop 'mcd_finite_sample_factor: invalid arguments'
      if(p>2)then
         coeff875=reshape([-0.455179464070565_dp,1.11192541278794_dp,-0.294241208320834_dp,1.09649329149811_dp],[2,2])
         coeff500=reshape([-1.42764571687802_dp,1.26263336932151_dp,-1.06141115981725_dp,1.28907991440387_dp],[2,2])
         y500=log(-coeff500(1,:)/real(p,dp)**coeff500(2,:))
         y875=log(-coeff875(1,:)/real(p,dp)**coeff875(2,:))
         a500(:,1)=1.0_dp;a500(:,2)=-log([2.0_dp,3.0_dp]*real(p*p,dp))
         a875(:,1)=1.0_dp;a875(:,2)=-log([2.0_dp,3.0_dp]*real(p*p,dp))
         call solve_linear(a500,y500,sol500,info);if(info/=0)error stop 'mcd_finite_sample_factor: solve failure'
         call solve_linear(a875,y875,sol875,info);if(info/=0)error stop 'mcd_finite_sample_factor: solve failure'
         fp500=1.0_dp-exp(sol500(1))/real(n,dp)**sol500(2)
         fp875=1.0_dp-exp(sol875(1))/real(n,dp)**sol875(2)
      else if(p==2)then
         fp500=1.0_dp-exp(0.673292623522027_dp)/real(n,dp)**0.691365864961895_dp
         fp875=1.0_dp-exp(0.446537815635445_dp)/real(n,dp)**1.06690782995919_dp
      else
         fp500=1.0_dp-exp(0.262024211897096_dp)/real(n,dp)**0.604756680630497_dp
         fp875=1.0_dp-exp(-0.351584646688712_dp)/real(n,dp)**1.01646567502486_dp
      end if
      if(alpha<=0.875_dp)then
         fpalpha=fp500+(fp875-fp500)/0.375_dp*(alpha-0.5_dp)
      else
         fpalpha=fp875+(1.0_dp-fp875)/0.125_dp*(alpha-0.875_dp)
      end if
      factor=1.0_dp/max(fpalpha,1.0e-12_dp)
   end function mcd_finite_sample_factor

   function mcd_reweighted_finite_sample_factor(p,n,alpha) result(factor)
      integer,intent(in)::p,n
      real(dp),intent(in)::alpha
      real(dp)::factor,fp500,fp875,fpalpha
      real(dp)::coeff500(2,2),coeff875(2,2),y500(2),y875(2),amat(2,2),sol500(2),sol875(2)
      integer::info
      if(p<1 .or. n<2 .or. alpha<0.5_dp .or. alpha>1.0_dp)error stop 'mcd_reweighted_finite_sample_factor: invalid arguments'
      if(p>2)then
         coeff875=reshape([-0.544482443573914_dp,1.25994483222292_dp,-0.343791072183285_dp,1.25159004257133_dp],[2,2])
         coeff500=reshape([-1.02842572724793_dp,1.67659883081926_dp,-0.26800273450853_dp,1.35968562893582_dp],[2,2])
         y500=log(-coeff500(1,:)/real(p,dp)**coeff500(2,:))
         y875=log(-coeff875(1,:)/real(p,dp)**coeff875(2,:))
         amat(:,1)=1.0_dp;amat(:,2)=-log([2.0_dp,3.0_dp]*real(p*p,dp))
         call solve_linear(amat,y500,sol500,info);if(info/=0)error stop 'mcd_reweighted_finite_sample_factor: solve failure'
         call solve_linear(amat,y875,sol875,info);if(info/=0)error stop 'mcd_reweighted_finite_sample_factor: solve failure'
         fp500=1.0_dp-exp(sol500(1))/real(n,dp)**sol500(2)
         fp875=1.0_dp-exp(sol875(1))/real(n,dp)**sol875(2)
      else if(p==2)then
         fp500=1.0_dp-exp(3.11101712909049_dp)/real(n,dp)**1.91401056721863_dp
         fp875=1.0_dp-exp(0.79473550581058_dp)/real(n,dp)**1.10081930350091_dp
      else
         fp500=1.0_dp-exp(1.11098143415027_dp)/real(n,dp)**1.5182890270453_dp
         fp875=1.0_dp-exp(-0.66046776772861_dp)/real(n,dp)**0.88939595831888_dp
      end if
      if(alpha<=0.875_dp)then
         fpalpha=fp500+(fp875-fp500)/0.375_dp*(alpha-0.5_dp)
      else
         fpalpha=fp875+(1.0_dp-fp875)/0.125_dp*(alpha-0.875_dp)
      end if
      factor=1.0_dp/max(fpalpha,1.0e-12_dp)
   end function mcd_reweighted_finite_sample_factor

   subroutine fast_mcd_partitioned(x,result,alpha,n_partitions,max_csteps,reweight_probability)
      real(dp),intent(in)::x(:,:)
      type(partitioned_mcd_result),intent(out)::result
      real(dp),intent(in),optional::alpha,reweight_probability
      integer,intent(in),optional::n_partitions,max_csteps
      type(detmcd_result)::group_fit,full_fit
      integer::n,p,h,ng,mc,g,start_i,end_i,m,info,steps,best_steps,i,j,nw
      integer,allocatable::perm(:),subset(:),best_subset(:),idx(:)
      real(dp),allocatable::xg(:,:),center(:),cov(:,:),best_center(:),best_cov(:,:),invcov(:,:),dist(:),sorted(:),vals(:),vecs(:,:),xw(:,:)
      logical,allocatable::mask(:)
      real(dp)::a,rwp,best_obj,obj,q,cut,raw_factor,rw_factor
      n=size(x,1);p=size(x,2)
      if(n<p+2 .or. p<1)error stop 'fast_mcd_partitioned: invalid dimensions'
      a=0.75_dp;if(present(alpha))a=alpha
      if(a<0.5_dp .or. a>1.0_dp)error stop 'fast_mcd_partitioned: alpha must be in [0.5,1]'
      ng=5;if(present(n_partitions))ng=max(1,n_partitions)
      ng=min(ng,max(1,n/(p+2)))
      mc=50;if(present(max_csteps))mc=max(1,max_csteps)
      rwp=0.975_dp;if(present(reweight_probability))rwp=reweight_probability
      h=max(p+1,min(n,h_alpha_n(a,n,p)))
      allocate(perm(n),subset(h),best_subset(h),idx(n),center(p),cov(p,p),best_center(p),best_cov(p,p),invcov(p,p),dist(n),sorted(n),vals(p),vecs(p,p),mask(n))
      call random_permutation(n,perm)
      best_obj=huge_penalty;best_subset=0;best_steps=0;result%candidates=0
      do g=1,ng
         start_i=1+(g-1)*n/ng;end_i=g*n/ng;m=end_i-start_i+1
         if(m<p+2)cycle
         allocate(xg(m,p));xg=x(perm(start_i:end_i),:)
         call cov_detmcd(xg,group_fit,alpha=a,max_csteps=mc,save_orderings=.false.)
         call full_csteps(x,group_fit%estimate%raw_center,group_fit%estimate%raw_covariance,h,mc,center,cov,subset,obj,steps,info)
         result%candidates=result%candidates+1
         if(info==0 .and. obj<best_obj)then
            best_obj=obj;best_center=center;best_cov=cov;best_subset=subset;best_steps=steps
         end if
         deallocate(xg)
      end do
      call cov_detmcd(x,full_fit,alpha=a,max_csteps=mc,save_orderings=.false.)
      call full_csteps(x,full_fit%estimate%raw_center,full_fit%estimate%raw_covariance,h,mc,center,cov,subset,obj,steps,info)
      result%candidates=result%candidates+1
      if(info==0 .and. obj<best_obj)then
         best_obj=obj;best_center=center;best_cov=cov;best_subset=subset;best_steps=steps
      end if
      if(any(best_subset==0))error stop 'fast_mcd_partitioned: no nonsingular candidate'
      q=chi_square_quantile(a,real(p,dp))
      raw_factor=mcd_consistency_factor(p,a)*mcd_finite_sample_factor(p,n,a)
      best_cov=raw_factor*best_cov
      call invert_symmetric(best_cov,invcov,info,ridge=1.0e-12_dp)
      do i=1,n
         dist(i)=dot_product(x(i,:)-best_center,matmul(invcov,x(i,:)-best_center))
      end do
      cut=chi_square_quantile(rwp,real(p,dp));mask=dist<=cut;nw=count(mask)
      allocate(result%estimate%center(p),result%estimate%covariance(p,p),result%estimate%raw_center(p),result%estimate%raw_covariance(p,p), &
               result%estimate%distances(n),result%estimate%weights(n),result%best_subset(h))
      result%estimate%raw_center=best_center;result%estimate%raw_covariance=best_cov
      result%estimate%center=best_center;result%estimate%covariance=best_cov
      result%estimate%distances=sqrt(max(dist,0.0_dp));result%estimate%weights=mask
      rw_factor=mcd_consistency_factor(p,rwp)*mcd_reweighted_finite_sample_factor(p,n,a)
      if(nw>p)then
         do j=1,p
            result%estimate%center(j)=sum(pack(x(:,j),mask))/real(nw,dp)
         end do
         allocate(xw(nw,p))
         do j=1,p;xw(:,j)=pack(x(:,j),mask);end do
         call covariance_matrix(xw,result%estimate%center,result%estimate%covariance)
         result%estimate%covariance=rw_factor*result%estimate%covariance
      end if
      call symmetric_eigen(result%estimate%raw_covariance,vals,vecs,info)
      result%exact_fit=info/=0 .or. any(vals<=1.0e-14_dp)
      result%estimate%h=h;result%estimate%iterations=best_steps;result%estimate%converged=.not.result%exact_fit
      result%best_subset=best_subset;result%partitions=ng;result%csteps=best_steps
      result%raw_consistency_factor=raw_factor;result%reweight_consistency_factor=rw_factor
      result%converged=result%estimate%converged
   end subroutine fast_mcd_partitioned

   subroutine fast_lts_partitioned(x,y,result,alpha,n_partitions,n_starts,max_csteps,reweight_cutoff)
      real(dp),intent(in)::x(:,:),y(:)
      type(partitioned_lts_result),intent(out)::result
      real(dp),intent(in),optional::alpha,reweight_cutoff
      integer,intent(in),optional::n_partitions,n_starts,max_csteps
      type(fast_lts_result)::group_fit,full_fit
      integer::n,p,h,ng,ns,mc,g,start_i,end_i,m,info,steps,nw
      integer,allocatable::perm(:),subset(:),best_subset(:)
      real(dp),allocatable::xg(:,:),yg(:),beta(:),best_beta(:),res(:),weights(:),xtwx(:,:),inv(:,:)
      real(dp)::a,cut,best_obj,obj,q,phi,consistency
      n=size(x,1);p=size(x,2)
      if(size(y)/=n .or. n<p+1)error stop 'fast_lts_partitioned: invalid dimensions'
      a=0.5_dp;if(present(alpha))a=alpha
      ng=5;if(present(n_partitions))ng=max(1,n_partitions)
      ng=min(ng,max(1,n/(p+2)))
      ns=100;if(present(n_starts))ns=max(1,n_starts)
      mc=50;if(present(max_csteps))mc=max(1,max_csteps)
      cut=2.5_dp;if(present(reweight_cutoff))cut=reweight_cutoff
      h=h_alpha_n(a,n,p)
      allocate(perm(n),subset(h),best_subset(h),beta(p),best_beta(p),res(n),weights(n),xtwx(p,p),inv(p,p))
      call random_permutation(n,perm)
      best_obj=huge_penalty;best_subset=0;result%candidates=0
      do g=1,ng
         start_i=1+(g-1)*n/ng;end_i=g*n/ng;m=end_i-start_i+1
         if(m<p+1)cycle
         allocate(xg(m,p),yg(m));xg=x(perm(start_i:end_i),:);yg=y(perm(start_i:end_i))
         call fast_lts_regression(xg,yg,group_fit,alpha=a,sampling='deterministic',n_starts=ns,max_csteps=mc)
         beta=group_fit%raw_coefficients
         call lts_csteps(x,y,beta,h,mc,subset,obj,steps,info)
         result%candidates=result%candidates+1
         if(info==0 .and. obj<best_obj)then;best_obj=obj;best_beta=beta;best_subset=subset;end if
         deallocate(xg,yg)
      end do
      call fast_lts_regression(x,y,full_fit,alpha=a,sampling='deterministic',n_starts=ns,max_csteps=mc)
      beta=full_fit%raw_coefficients
      call lts_csteps(x,y,beta,h,mc,subset,obj,steps,info)
      result%candidates=result%candidates+1
      if(info==0 .and. obj<best_obj)then;best_obj=obj;best_beta=beta;best_subset=subset;end if
      if(any(best_subset==0))error stop 'fast_lts_partitioned: no valid candidate'
      allocate(result%estimate%coefficients(p),result%estimate%raw_coefficients(p),result%estimate%residuals(n),result%estimate%raw_residuals(n), &
               result%estimate%weights(n),result%estimate%raw_weights(n),result%estimate%covariance(p,p),result%estimate%best_subset(h))
      result%estimate%raw_coefficients=best_beta;result%estimate%raw_residuals=y-matmul(x,best_beta);result%estimate%objective=best_obj
      result%estimate%best_subset=best_subset;result%estimate%h=h
      result%estimate%raw_scale=sqrt(sum(result%estimate%raw_residuals(best_subset)**2)/real(h,dp))
      q=normal_quantile(real(h+n,dp)/(2.0_dp*real(n,dp)));phi=exp(-0.5_dp*q*q)/sqrt(2.0_dp*acos(-1.0_dp))
      consistency=merge(1.0_dp/sqrt(max(1.0e-12_dp,1.0_dp-2.0_dp*real(n,dp)*q*phi/real(h,dp))),1.0_dp,q>1.0e-12_dp)
      result%estimate%raw_scale=result%estimate%raw_scale*consistency
      result%estimate%raw_weights=merge(1.0_dp,0.0_dp,abs(result%estimate%raw_residuals)<=cut*max(result%estimate%raw_scale,1.0e-14_dp))
      call weighted_ls(x,y,result%estimate%raw_weights,result%estimate%coefficients,info)
      if(info/=0)result%estimate%coefficients=best_beta
      result%estimate%residuals=y-matmul(x,result%estimate%coefficients);nw=count(result%estimate%raw_weights>0.0_dp)
      result%estimate%scale=sqrt(sum(result%estimate%raw_weights*result%estimate%residuals**2)/real(max(1,nw-p),dp))
      result%estimate%weights=merge(1.0_dp,0.0_dp,abs(result%estimate%residuals)<=cut*max(result%estimate%scale,1.0e-14_dp))
      xtwx=matmul(transpose(x*spread(result%estimate%weights,2,p)),x);call invert_symmetric(xtwx,inv,info,ridge=1.0e-12_dp)
      result%estimate%covariance=result%estimate%scale**2*inv
      result%estimate%trials=result%candidates;result%estimate%csteps=mc;result%estimate%exhaustive=.false.;result%estimate%converged=info==0
      result%partitions=ng;result%converged=result%estimate%converged
   end subroutine fast_lts_partitioned

   subroutine full_csteps(x,start_center,start_cov,h,max_steps,center,cov,subset,objective,steps,info)
      real(dp),intent(in)::x(:,:),start_center(:),start_cov(:,:)
      integer,intent(in)::h,max_steps
      real(dp),intent(out)::center(:),cov(:,:),objective
      integer,intent(out)::subset(:),steps,info
      real(dp),allocatable::invcov(:,:),dist(:),sorted(:),vals(:),vecs(:,:)
      integer,allocatable::idx(:),previous(:)
      integer::n,p,i
      n=size(x,1);p=size(x,2)
      allocate(invcov(p,p),dist(n),sorted(n),vals(p),vecs(p,p),idx(n),previous(h))
      center=start_center;cov=start_cov;subset=0;info=0
      do steps=1,max_steps
         call invert_symmetric(cov,invcov,info,ridge=1.0e-12_dp);if(info/=0)return
         do i=1,n;dist(i)=dot_product(x(i,:)-center,matmul(invcov,x(i,:)-center));end do
         sorted=dist;call sort_real_with_index(sorted,idx);previous=subset;subset=idx(1:h)
         call subset_stats(x,subset,center,cov)
         if(all(subset==previous))exit
      end do
      call symmetric_eigen(cov,vals,vecs,info)
      if(info/=0 .or. any(vals<=1.0e-14_dp))then;info=1;objective=huge_penalty;else;objective=sum(log(vals));end if
   end subroutine full_csteps

   subroutine lts_csteps(x,y,beta,h,max_steps,subset,objective,steps,info)
      real(dp),intent(in)::x(:,:),y(:)
      real(dp),intent(inout)::beta(:)
      integer,intent(in)::h,max_steps
      integer,intent(out)::subset(:),steps,info
      real(dp),intent(out)::objective
      real(dp),allocatable::res(:),sq(:),sorted(:),newbeta(:)
      integer,allocatable::idx(:),previous(:)
      allocate(res(size(y)),sq(size(y)),sorted(size(y)),newbeta(size(beta)),idx(size(y)),previous(h))
      subset=0;info=0
      do steps=1,max_steps
         res=y-matmul(x,beta);sq=res*res;sorted=sq;call sort_real_with_index(sorted,idx);previous=subset;subset=idx(1:h)
         call least_squares(x(subset,:),y(subset),newbeta,info);if(info/=0)return
         beta=newbeta
         if(all(subset==previous))exit
      end do
      res=y-matmul(x,beta);sq=res*res;sorted=sq;call sort_real_with_index(sorted,idx);subset=idx(1:h);objective=sum(sorted(1:h))
   end subroutine lts_csteps

   subroutine subset_stats(x,subset,center,cov)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::subset(:)
      real(dp),intent(out)::center(:),cov(:,:)
      real(dp),allocatable::xs(:,:)
      integer::j
      allocate(xs(size(subset),size(x,2)));xs=x(subset,:)
      do j=1,size(x,2);center(j)=sum(xs(:,j))/real(size(subset),dp);end do
      call covariance_matrix(xs,center,cov)
   end subroutine subset_stats

   subroutine weighted_ls(x,y,w,beta,info)
      real(dp),intent(in)::x(:,:),y(:),w(:)
      real(dp),intent(out)::beta(:)
      integer,intent(out)::info
      real(dp),allocatable::a(:,:),b(:)
      integer::j
      allocate(a(size(x,1),size(x,2)),b(size(y)))
      do j=1,size(x,2);a(:,j)=x(:,j)*sqrt(max(w,0.0_dp));end do
      b=y*sqrt(max(w,0.0_dp));call least_squares(a,b,beta,info)
   end subroutine weighted_ls

   subroutine random_permutation(n,perm)
      integer,intent(in)::n
      integer,intent(out)::perm(:)
      integer::i,j,tmp
      real(dp)::u
      perm=[(i,i=1,n)]
      do i=n,2,-1
         call random_number(u);j=1+int(u*real(i,dp));j=min(i,max(1,j));tmp=perm(i);perm(i)=perm(j);perm(j)=tmp
      end do
   end subroutine random_permutation
end module robustbase_fast_algorithms
