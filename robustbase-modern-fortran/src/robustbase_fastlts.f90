! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_fastlts
   use robustbase_kinds, only: dp, huge_penalty
   use robustbase_sort, only: median, sort_real_with_index
   use robustbase_probability, only: normal_quantile
   use robustbase_linalg, only: least_squares, invert_symmetric, matrix_rank
   implicit none
   private
   public :: fast_lts_result, fast_lts_regression, h_alpha_n

   type :: fast_lts_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: raw_coefficients(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: raw_residuals(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: raw_weights(:)
      real(dp), allocatable :: covariance(:,:)
      integer, allocatable :: best_subset(:)
      real(dp) :: scale = 0.0_dp
      real(dp) :: raw_scale = 0.0_dp
      real(dp) :: objective = huge_penalty
      integer :: h = 0
      integer :: trials = 0
      integer :: csteps = 0
      logical :: exhaustive = .false.
      logical :: converged = .false.
   end type fast_lts_result
contains
   integer function h_alpha_n(alpha,n,p) result(h)
      real(dp),intent(in)::alpha
      integer,intent(in)::n,p
      integer::n2
      n2=(n+p+1)/2
      h=int(floor(real(2*n2-n,dp)+2.0_dp*real(n-n2,dp)*alpha))
      h=max(p,min(n,h))
   end function h_alpha_n

   subroutine fast_lts_regression(x,y,result,alpha,sampling,n_starts,max_subsets,max_csteps,adjust_intercept,reweight_cutoff)
      real(dp),intent(in)::x(:,:),y(:)
      type(fast_lts_result),intent(out)::result
      real(dp),intent(in),optional::alpha,reweight_cutoff
      character(len=*),intent(in),optional::sampling
      integer,intent(in),optional::n_starts,max_subsets,max_csteps
      logical,intent(in),optional::adjust_intercept
      real(dp),allocatable::beta(:),best_beta(:),residual(:),squared(:),sorted(:),weights(:),xtwx(:,:),inv(:,:)
      integer,allocatable::subset(:),best_subset(:),idx(:),comb(:),det_candidates(:,:)
      real(dp)::a,best_objective,objective,cutoff,phi,q,consistency
      character(len=16)::mode
      integer::n,p,h,ns,limit,mc,trial,info,best_steps,steps,nw,count_comb,n_det
      logical::adjust,done,use_exhaustive

      n=size(x,1);p=size(x,2)
      if(size(y)/=n .or. n<p .or. p<1) error stop "fast_lts_regression: invalid dimensions"
      if(matrix_rank(x)<p) error stop "fast_lts_regression: singular design"
      a=0.5_dp
      if(present(alpha)) a=alpha
      if(a<0.5_dp .or. a>1.0_dp) error stop "fast_lts_regression: alpha must be in [0.5,1]"
      mode='random'
      if(present(sampling)) mode=adjustl(sampling)
      ns=500
      if(present(n_starts)) ns=max(1,n_starts)
      limit=5000
      if(present(max_subsets)) limit=max(1,max_subsets)
      mc=50
      if(present(max_csteps)) mc=max(1,max_csteps)
      adjust=.false.
      if(present(adjust_intercept)) adjust=adjust_intercept
      cutoff=2.5_dp
      if(present(reweight_cutoff)) cutoff=reweight_cutoff
      h=h_alpha_n(a,n,p)

      allocate(beta(p),best_beta(p),residual(n),squared(n),sorted(n),weights(n),xtwx(p,p),inv(p,p),subset(h),best_subset(h),idx(n),comb(p))
      best_objective=huge_penalty
      best_beta=0.0_dp
      best_subset=0
      best_steps=0
      count_comb=combination_count(n,p,huge(1))
      select case(trim(mode))
      case('exact')
         use_exhaustive=.true.
      case('best')
         use_exhaustive=count_comb<=limit
      case('deterministic','random')
         use_exhaustive=.false.
      case default
         error stop "fast_lts_regression: sampling must be random, deterministic, best, or exact"
      end select

      result%trials=0
      if(use_exhaustive) then
         comb=[(trial,trial=1,p)]
         done=.false.
         do while(.not.done)
            call evaluate_start(comb)
            result%trials=result%trials+1
            call next_combination(comb,n,done)
         end do
      else
         call build_deterministic_starts(x,y,p,ns,det_candidates,n_det)
         do trial=1,n_det
            call evaluate_start(det_candidates(:,trial))
            result%trials=result%trials+1
         end do
         if(trim(mode)/='deterministic') then
            do trial=1,ns
               call random_subset(n,p,comb)
               call evaluate_start(comb)
               result%trials=result%trials+1
            end do
         end if
      end if
      if(any(best_subset==0)) error stop "fast_lts_regression: no valid subset"

      allocate(result%coefficients(p),result%raw_coefficients(p),result%residuals(n),result%raw_residuals(n), &
               result%weights(n),result%raw_weights(n),result%covariance(p,p),result%best_subset(h))
      result%raw_coefficients=best_beta
      result%raw_residuals=y-matmul(x,best_beta)
      result%objective=best_objective
      result%best_subset=best_subset
      result%h=h
      result%csteps=best_steps
      result%exhaustive=use_exhaustive
      result%raw_scale=sqrt(sum(result%raw_residuals(best_subset)**2)/real(h,dp))
      q=normal_quantile((real(h+n,dp))/(2.0_dp*real(n,dp)))
      phi=exp(-0.5_dp*q*q)/sqrt(2.0_dp*acos(-1.0_dp))
      if(q>1.0e-12_dp) then
         consistency=1.0_dp/sqrt(max(1.0e-12_dp,1.0_dp-2.0_dp*real(n,dp)*q*phi/real(h,dp)))
      else
         consistency=1.0_dp
      end if
      result%raw_scale=result%raw_scale*consistency
      if(result%raw_scale<=1.0e-14_dp) then
         result%raw_weights=merge(1.0_dp,0.0_dp,abs(result%raw_residuals)<=1.0e-7_dp)
      else
         result%raw_weights=merge(1.0_dp,0.0_dp,abs(result%raw_residuals)<=cutoff*result%raw_scale)
      end if
      call weighted_ls(x,y,result%raw_weights,result%coefficients,info)
      if(info/=0) result%coefficients=result%raw_coefficients
      result%residuals=y-matmul(x,result%coefficients)
      nw=count(result%raw_weights>0.0_dp)
      if(nw>p) then
         result%scale=sqrt(sum(result%raw_weights*result%residuals**2)/real(nw-p,dp))
      else
         result%scale=result%raw_scale
      end if
      if(result%scale<=1.0e-14_dp) then
         result%weights=result%raw_weights
      else
         result%weights=merge(1.0_dp,0.0_dp,abs(result%residuals)<=cutoff*result%scale)
      end if
      xtwx=matmul(transpose(x*spread(result%weights,2,p)),x)
      call invert_symmetric(xtwx,inv,info,ridge=1.0e-12_dp)
      result%covariance=result%scale**2*inv
      result%converged=(info==0)
   contains
      subroutine evaluate_start(candidate)
         integer,intent(in)::candidate(:)
         integer::local_info
         call least_squares(x(candidate,:),y(candidate),beta,local_info)
         if(local_info/=0) return
         call concentration_steps(x,y,beta,h,mc,adjust,subset,objective,steps,local_info)
         if(local_info/=0) return
         if(objective<best_objective) then
            best_objective=objective
            best_beta=beta
            best_subset=subset
            best_steps=steps
         end if
      end subroutine evaluate_start
   end subroutine fast_lts_regression

   subroutine concentration_steps(x,y,beta,h,max_steps,adjust_intercept,subset,objective,steps,info)
      real(dp),intent(in)::x(:,:),y(:)
      real(dp),intent(inout)::beta(:)
      integer,intent(in)::h,max_steps
      logical,intent(in)::adjust_intercept
      integer,intent(out)::subset(:),steps,info
      real(dp),intent(out)::objective
      real(dp),allocatable::residual(:),squared(:),sorted(:),newbeta(:)
      integer,allocatable::idx(:),previous(:)
      integer::n,p
      n=size(x,1);p=size(x,2)
      allocate(residual(n),squared(n),sorted(n),newbeta(p),idx(n),previous(h))
      subset=0
      info=0
      do steps=1,max_steps
         residual=y-matmul(x,beta)
         squared=residual**2
         sorted=squared
         call sort_real_with_index(sorted,idx)
         previous=subset
         subset=idx(1:h)
         call least_squares(x(subset,:),y(subset),newbeta,info)
         if(info/=0) return
         if(adjust_intercept .and. all(abs(x(:,1)-1.0_dp)<=10.0_dp*epsilon(1.0_dp))) then
            newbeta(1)=newbeta(1)+median(y(subset)-matmul(x(subset,:),newbeta))
         end if
         beta=newbeta
         if(all(subset==previous)) exit
      end do
      residual=y-matmul(x,beta)
      squared=residual**2
      sorted=squared
      call sort_real_with_index(sorted,idx)
      subset=idx(1:h)
      objective=sum(sorted(1:h))
   end subroutine concentration_steps

   subroutine build_deterministic_starts(x,y,p,maximum,candidates,n_candidates)
      real(dp),intent(in)::x(:,:),y(:)
      integer,intent(in)::p,maximum
      integer,allocatable,intent(out)::candidates(:,:)
      integer,intent(out)::n_candidates
      real(dp),allocatable::score(:),sorted_values(:)
      integer,allocatable::idx(:),candidate(:),work(:,:)
      integer::n,j,k,offset,pos
      n=size(x,1)
      allocate(score(n),sorted_values(n),idx(n),candidate(p),work(p,maximum))
      n_candidates=0
      score=y
      call emit_order(score)
      score=-y
      call emit_order(score)
      do j=1,size(x,2)
         score=x(:,j)
         call emit_order(score)
         score=-x(:,j)
         call emit_order(score)
      end do
      allocate(candidates(p,n_candidates))
      if(n_candidates>0) candidates=work(:,1:n_candidates)
   contains
      subroutine emit_order(values)
         real(dp),intent(in)::values(:)
         if(n_candidates>=maximum) return
         sorted_values=values
         call sort_real_with_index(sorted_values,idx)
         do offset=0,min(n-p,p)
            if(n_candidates>=maximum) exit
            do k=1,p
               pos=1+mod(offset+(k-1)*max(1,n/p),n)
               candidate(k)=idx(pos)
            end do
            if(all_unique(candidate)) then
               n_candidates=n_candidates+1
               work(:,n_candidates)=candidate
            end if
         end do
      end subroutine emit_order
   end subroutine build_deterministic_starts

   logical function all_unique(a) result(ok)
      integer,intent(in)::a(:)
      integer::i,j
      ok=.true.
      do i=1,size(a)-1
         do j=i+1,size(a)
            if(a(i)==a(j)) then
               ok=.false.
               return
            end if
         end do
      end do
   end function all_unique

   subroutine weighted_ls(x,y,w,beta,info)
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
   end subroutine weighted_ls

   subroutine random_subset(n,k,out)
      integer,intent(in)::n,k
      integer,intent(out)::out(:)
      real(dp)::u
      integer::i,candidate
      i=0
      do while(i<k)
         call random_number(u)
         candidate=min(n,1+int(u*real(n,dp)))
         if(i==0 .or. .not.any(out(1:i)==candidate)) then
            i=i+1
            out(i)=candidate
         end if
      end do
   end subroutine random_subset

   integer function combination_count(n,k,limit) result(value)
      integer,intent(in)::n,k,limit
      integer::i,kk
      real(dp)::v
      kk=min(k,n-k)
      v=1.0_dp
      do i=1,kk
         v=v*real(n-kk+i,dp)/real(i,dp)
         if(v>=real(limit,dp)) then
            value=limit
            return
         end if
      end do
      value=nint(v)
   end function combination_count

   subroutine next_combination(comb,n,done)
      integer,intent(inout)::comb(:)
      integer,intent(in)::n
      logical,intent(out)::done
      integer::k,i,j
      k=size(comb)
      i=k
      do while(i>=1)
         if(comb(i)/=n-k+i) exit
         i=i-1
      end do
      if(i==0) then
         done=.true.
         return
      end if
      comb(i)=comb(i)+1
      do j=i+1,k
         comb(j)=comb(j-1)+1
      end do
      done=.false.
   end subroutine next_combination
end module robustbase_fastlts
