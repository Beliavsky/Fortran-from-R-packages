! SPDX-License-Identifier: MIT
module cgnm_postprocess
   use cgnm_kinds, only : dp
   use cgnm_types, only : cgnm_result
   use cgnm_utils, only : sort_index, quantile_value
   implicit none
   private
   public :: top_indices, best_approximate_minimizers, accepted_max_ssr
   public :: accepted_indices, accepted_indices_binary, column_quantiles
contains
   subroutine top_indices(res,nout,idx)
      type(cgnm_result), intent(in) :: res
      integer, intent(in) :: nout
      integer, allocatable, intent(out) :: idx(:)
      integer, allocatable :: order(:)
      real(dp), allocatable :: r(:)
      integer :: n
      n=size(res%residual_history,1)
      allocate(r(n),order(n)); r=res%residual_history(:,size(res%residual_history,2))
      call sort_index(r,order)
      allocate(idx(min(max(0,nout),n)))
      if (size(idx)>0) idx=order(1:size(idx))
   end subroutine top_indices

   subroutine best_approximate_minimizers(res,nout,theta,idx)
      type(cgnm_result), intent(in) :: res
      integer, intent(in) :: nout
      real(dp), allocatable, intent(out) :: theta(:,:)
      integer, allocatable, intent(out), optional :: idx(:)
      integer, allocatable :: ii(:)
      call top_indices(res,nout,ii)
      allocate(theta(size(ii),size(res%theta,2)))
      if (size(ii)>0) theta=res%theta(ii,:)
      if (present(idx)) then; allocate(idx(size(ii))); idx=ii; end if
   end subroutine best_approximate_minimizers

   real(dp) function invnorm(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a1=-3.969683028665376e1_dp,a2=2.209460984245205e2_dp
      real(dp), parameter :: a3=-2.759285104469687e2_dp,a4=1.383577518672690e2_dp
      real(dp), parameter :: a5=-3.066479806614716e1_dp,a6=2.506628277459239_dp
      real(dp), parameter :: b1=-5.447609879822406e1_dp,b2=1.615858368580409e2_dp
      real(dp), parameter :: b3=-1.556989798598866e2_dp,b4=6.680131188771972e1_dp
      real(dp), parameter :: b5=-1.328068155288572e1_dp
      real(dp), parameter :: c1=-7.784894002430293e-3_dp,c2=-3.223964580411365e-1_dp
      real(dp), parameter :: c3=-2.400758277161838_dp,c4=-2.549732539343734_dp
      real(dp), parameter :: c5=4.374664141464968_dp,c6=2.938163982698783_dp
      real(dp), parameter :: d1=7.784695709041462e-3_dp,d2=3.224671290700398e-1_dp
      real(dp), parameter :: d3=2.445134137142996_dp,d4=3.754408661907416_dp
      real(dp) :: q,r
      if (p<=0.0_dp) then; x=-huge(1.0_dp); return; end if
      if (p>=1.0_dp) then; x=huge(1.0_dp); return; end if
      if (p<0.02425_dp) then
         q=sqrt(-2.0_dp*log(p))
         x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p>0.97575_dp) then
         q=sqrt(-2.0_dp*log(1.0_dp-p))
         x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else
         q=p-0.5_dp; r=q*q
         x=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/ &
           (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      end if
   end function invnorm

   real(dp) function t_quantile_approx(p,df) result(t)
      real(dp), intent(in) :: p
      integer, intent(in) :: df
      real(dp) :: z,v
      z=invnorm(p); v=real(max(1,df),dp)
      t=z+(z**3+z)/(4.0_dp*v)+(5.0_dp*z**5+16.0_dp*z**3+3.0_dp*z)/(96.0_dp*v*v)
   end function t_quantile_approx

   real(dp) function accepted_max_ssr(res,cutoff_pvalue,num_parameters,use_accepted,algorithm) result(thr)
      type(cgnm_result), intent(in) :: res
      real(dp), intent(in), optional :: cutoff_pvalue
      integer, intent(in), optional :: num_parameters,algorithm
      logical, intent(in), optional :: use_accepted
      real(dp), allocatable :: ssr(:),diffs(:)
      integer, allocatable :: ord(:)
      real(dp) :: alpha,tstat,gcrit,meanv,sdv,z
      integer :: n,i,np,alg
      logical :: usea
      n=size(res%residual_history,1); allocate(ssr(n),ord(n))
      ssr=res%residual_history(:,size(res%residual_history,2)); call sort_index(ssr,ord)
      alpha=0.05_dp; if (present(cutoff_pvalue)) alpha=cutoff_pvalue
      np=n; if (present(num_parameters)) np=min(max(1,num_parameters),n)
      usea=.true.; if (present(use_accepted)) usea=use_accepted
      alg=2; if (present(algorithm)) alg=algorithm
      if (.not.usea) then; thr=ssr(ord(np)); return; end if
      if (alg/=2 .or. n<5) then
         ! Elbow on the sorted SSR itself: stable array-level equivalent of the
         ! package's fallback screening path.
         thr=ssr(ord(max(1,min(n,(n+1)/2))))
         do i=2,n-1
            if ((ssr(ord(i+1))-ssr(ord(i))) > 4.0_dp*max(epsilon(1.0_dp), &
                 ssr(ord(i))-ssr(ord(i-1)))) then
               thr=ssr(ord(i)); exit
            end if
         end do
      else
         allocate(diffs(n-1)); diffs=0.0_dp; thr=ssr(ord(n))
         do i=2,n
            diffs(i-1)=sum(res%y(ord(i-1),:)-res%y(ord(i),:))
            if (i-1>3) then
               meanv=sum(diffs(1:i-1))/real(i-1,dp)
               sdv=sqrt(sum((diffs(1:i-1)-meanv)**2)/real(i-2,dp))
               tstat=t_quantile_approx(1.0_dp-alpha/(2.0_dp*real(i-1,dp)),i-3)
               gcrit=real(i-2,dp)/sqrt(real(i-1,dp))* &
                     sqrt(tstat*tstat/(real(i-3,dp)+tstat*tstat))
               if (sdv>0.0_dp) then
                  z=abs(meanv-diffs(i-1))/sdv
                  if (z>gcrit) then; thr=ssr(ord(i)); exit; end if
               end if
            end if
         end do
      end if
      thr=max(thr,minval(ssr)+sqrt(epsilon(1.0_dp)))
   end function accepted_max_ssr

   subroutine accepted_indices(res,idx,cutoff_pvalue,num_parameters,use_accepted,algorithm)
      type(cgnm_result), intent(in) :: res
      integer, allocatable, intent(out) :: idx(:)
      real(dp), intent(in), optional :: cutoff_pvalue
      integer, intent(in), optional :: num_parameters,algorithm
      logical, intent(in), optional :: use_accepted
      real(dp), allocatable :: ssr(:)
      real(dp) :: thr
      integer :: i,k,n
      n=size(res%residual_history,1); allocate(ssr(n))
      ssr=res%residual_history(:,size(res%residual_history,2))
      thr=accepted_max_ssr(res,cutoff_pvalue,num_parameters,use_accepted,algorithm)
      k=count(ssr<=thr); allocate(idx(k)); k=0
      do i=1,n
         if (ssr(i)<=thr) then; k=k+1; idx(k)=i; end if
      end do
   end subroutine accepted_indices

   subroutine accepted_indices_binary(res,mask,cutoff_pvalue,num_parameters,use_accepted,algorithm)
      type(cgnm_result), intent(in) :: res
      logical, allocatable, intent(out) :: mask(:)
      real(dp), intent(in), optional :: cutoff_pvalue
      integer, intent(in), optional :: num_parameters,algorithm
      logical, intent(in), optional :: use_accepted
      integer, allocatable :: idx(:)
      integer :: i
      call accepted_indices(res,idx,cutoff_pvalue,num_parameters,use_accepted,algorithm)
      allocate(mask(size(res%x,1))); mask=.false.
      do i=1,size(idx); mask(idx(i))=.true.; end do
   end subroutine accepted_indices_binary

   subroutine column_quantiles(data,prob,q)
      real(dp), intent(in) :: data(:,:),prob
      real(dp), intent(out) :: q(size(data,2))
      integer :: j
      do j=1,size(data,2); q(j)=quantile_value(data(:,j),prob); end do
   end subroutine column_quantiles
end module cgnm_postprocess
