module lavaan_bootstrap
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map
   use lavaan_fit, only : sem_fit_result, fit_ram_data
   implicit none
   private

   type, public :: bootstrap_sem_result
      real(dp), allocatable :: par(:, :), se(:), bias(:), ci_lower(:), ci_upper(:)
      integer :: n_success = 0, status = 0
   end type bootstrap_sem_result

   public :: bootstrap_ram_data

contains

   subroutine bootstrap_ram_data(template,map,data,nboot,result,seed,level)
      type(ram_model),intent(in)::template
      type(ram_free_map),intent(in)::map
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::nboot
      type(bootstrap_sem_result),intent(out)::result
      integer,intent(in),optional::seed
      real(dp),intent(in),optional::level
      type(sem_fit_result)::base,fit
      real(dp),allocatable::sample(:,:),vals(:),r(:)
      integer,allocatable::seedv(:)
      integer::n,p,k,b,i,idx,ns,m
      real(dp)::lev,plo,phi

      n=size(data,1)
      p=size(data,2)
      lev=0.95_dp
      if(present(level)) lev=level
      if(nboot<2 .or. lev<=0.0_dp .or. lev>=1.0_dp) then
      result%status=-1
      return
      end if
      if(present(seed)) then
         call random_seed(size=ns)
         allocate(seedv(ns))
         do i=1,ns
         seedv(i)=mod(abs(seed)+104729*i,2147483646)+1
         end do
         call random_seed(put=seedv)
      end if
      call fit_ram_data(template,map,data,base)
      if(.not.base%converged) then
      result%status=-2
      return
      end if
      k=size(base%par)
      allocate(result%par(nboot,k),sample(n,p),r(n))
      result%par=0.0_dp
      m=0
      do b=1,nboot
         call random_number(r)
         do i=1,n
            idx=1+int(r(i)*real(n,dp))
            if(idx>n) idx=n
            sample(i,:)=data(idx,:)
         end do
         call fit_ram_data(template,map,sample,fit)
         if(fit%converged .and. size(fit%par)==k) then
            m=m+1
            result%par(m,:)=fit%par
         end if
      end do
      result%n_success=m
      if(m<2) then
      result%status=-3
      return
      end if
      allocate(result%se(k),result%bias(k),result%ci_lower(k),result%ci_upper(k),vals(m))
      plo=(1.0_dp-lev)/2.0_dp
      phi=1.0_dp-plo
      do i=1,k
         vals=result%par(1:m,i)
         result%bias(i)=sum(vals)/real(m,dp)-base%par(i)
         result%se(i)=sqrt(sum((vals-sum(vals)/real(m,dp))**2)/real(m-1,dp))
         call sort_real(vals)
         result%ci_lower(i)=quantile_sorted(vals,plo)
         result%ci_upper(i)=quantile_sorted(vals,phi)
      end do
      if(m<nboot) result%par=result%par(1:m,:)
      result%status=0
   end subroutine bootstrap_ram_data

   subroutine sort_real(x)
      real(dp),intent(inout)::x(:)
      integer::i,j
      real(dp)::key
      do i=2,size(x)
         key=x(i)
         j=i-1
         do while(j>=1)
            if(x(j)<=key) exit
            x(j+1)=x(j)
            j=j-1
         end do
         x(j+1)=key
      end do
   end subroutine sort_real

   function quantile_sorted(x,p) result(q)
      real(dp),intent(in)::x(:),p
      real(dp)::q,h,w
      integer::lo,hi,n
      n=size(x)
      if(n==1) then
      q=x(1)
      return
      end if
      h=1.0_dp+(real(n-1,dp))*p
      lo=floor(h)
      hi=ceiling(h)
      w=h-real(lo,dp)
      q=(1.0_dp-w)*x(max(1,lo))+w*x(min(n,hi))
   end function quantile_sorted
end module lavaan_bootstrap
