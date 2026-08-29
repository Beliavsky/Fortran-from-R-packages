module lavaan_twolevel
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map
   use lavaan_fit, only : sem_fit_result, fit_ram_cov
   use lavaan_linalg, only : sample_mean_cov
   implicit none
   private

   type, public :: twolevel_sem_result
      type(sem_fit_result) :: within, between
      real(dp), allocatable :: within_cov(:, :), between_cov(:, :), cluster_mean(:, :), icc(:)
      integer :: ncluster = 0, status = 0
   end type twolevel_sem_result

   public :: fit_ram_twolevel

contains

   subroutine fit_ram_twolevel(within_template,within_map,between_template,between_map,data,cluster,result)
      type(ram_model),intent(in)::within_template,between_template
      type(ram_free_map),intent(in)::within_map,between_map
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::cluster(:)
      type(twolevel_sem_result),intent(out)::result
      integer,allocatable::labels(:),counts(:)
      real(dp),allocatable::centered(:,:),wmean(:),bmean(:),rawbcov(:,:),adjbcov(:,:)
      integer::n,p,g,i,ng
      real(dp)::inv_m

      n=size(data,1)
      p=size(data,2)
      if(size(cluster)/=n) then
      result%status=-1
      return
      end if
      call unique_labels(cluster,labels)
      ng=size(labels)
      result%ncluster=ng
      if(ng<2) then
      result%status=-2
      return
      end if
      allocate(result%cluster_mean(ng,p),counts(ng),centered(n,p))
      result%cluster_mean=0.0_dp
      counts=0
      do i=1,n
         g=label_index(labels,cluster(i))
         result%cluster_mean(g,:)=result%cluster_mean(g,:)+data(i,:)
         counts(g)=counts(g)+1
      end do
      do g=1,ng
         if(counts(g)<=0) then
         result%status=-3
         return
         end if
         result%cluster_mean(g,:)=result%cluster_mean(g,:)/real(counts(g),dp)
      end do
      do i=1,n
         g=label_index(labels,cluster(i))
         centered(i,:)=data(i,:)-result%cluster_mean(g,:)
      end do
      call sample_mean_cov(centered,wmean,result%within_cov)
      call sample_mean_cov(result%cluster_mean,bmean,rawbcov)
      inv_m=sum(1.0_dp/real(counts,dp))/real(ng,dp)
      adjbcov=rawbcov-inv_m*result%within_cov
      do i=1,p
         if(adjbcov(i,i)<1.0e-8_dp*max(rawbcov(i,i),1.0_dp)) adjbcov(i,i)=1.0e-8_dp*max(rawbcov(i,i),1.0_dp)
      end do
      result%between_cov=0.5_dp*(adjbcov+transpose(adjbcov))
      call fit_ram_cov(within_template,within_map,result%within_cov,wmean,n,result%within)
      call fit_ram_cov(between_template,between_map,result%between_cov,bmean,ng,result%between)
      allocate(result%icc(p))
      do i=1,p
         result%icc(i)=result%between_cov(i,i)/(result%between_cov(i,i)+result%within_cov(i,i))
         result%icc(i)=max(0.0_dp,min(1.0_dp,result%icc(i)))
      end do
      result%status=merge(0,1,result%within%converged .and. result%between%converged)
   end subroutine fit_ram_twolevel

   subroutine unique_labels(x,u)
      integer,intent(in)::x(:)
      integer,allocatable,intent(out)::u(:)
      integer,allocatable::tmp(:)
      integer::i,m
      allocate(tmp(size(x)))
      m=0
      do i=1,size(x)
         if(m==0 .or. .not.any(tmp(1:m)==x(i))) then
         m=m+1
         tmp(m)=x(i)
         end if
      end do
      allocate(u(m))
      u=tmp(1:m)
   end subroutine unique_labels

   integer function label_index(u,x) result(idx)
      integer,intent(in)::u(:),x
      integer::i
      idx=0
      do i=1,size(u)
      if(u(i)==x) then
      idx=i
      return
      end if
      end do
   end function label_index
end module lavaan_twolevel
