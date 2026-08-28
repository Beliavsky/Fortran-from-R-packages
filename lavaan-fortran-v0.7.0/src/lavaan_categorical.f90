module lavaan_categorical
   use lavaan_kinds, only : dp
   use lavaan_ordinal, only : ordinal_thresholds, polychoric_matrix
   use lavaan_linalg, only : inverse_general
   implicit none
   private

   type, public :: categorical_stats_result
      real(dp), allocatable :: thresholds(:), correlation(:, :), stats(:)
      real(dp), allocatable :: gamma(:, :), weight(:, :), dwls_weight(:)
      integer, allocatable :: threshold_offset(:), ncat(:)
      integer :: status=0, n_jackknife=0
   end type categorical_stats_result

   public :: categorical_wls_statistics

contains

   subroutine categorical_wls_statistics(data,result)
      integer,intent(in)::data(:,:)
      type(categorical_stats_result),intent(out)::result
      real(dp),allocatable::jk(:,:),meanjk(:),d(:),ridge(:,:),inv(:,:)
      integer,allocatable::sub(:,:)
      integer::n,p,q,m,r,i,info
      n=size(data,1)
      p=size(data,2)
      if(n<5 .or. p<2 .or. minval(data)<1) then
      result%status=-1
      return
      end if
      call compute_stats(data,result%ncat,result%threshold_offset,result%thresholds,result%correlation,result%stats,info)
      if(info/=0) then
      result%status=info
      return
      end if
      q=size(result%stats)
      allocate(jk(n,q),sub(n-1,p))
      jk=0.0_dp
      m=0
      do r=1,n
         if(r>1) sub(1:r-1,:)=data(1:r-1,:)
         if(r<n) sub(r:n-1,:)=data(r+1:n,:)
         block
            integer,allocatable::nc(:),off(:)
            real(dp),allocatable::th(:),cc(:,:),ss(:)
            call compute_stats(sub,nc,off,th,cc,ss,info)
            if(info==0 .and. size(ss)==q .and. all(nc==result%ncat)) then
               m=m+1
               jk(m,:)=ss
            end if
         end block
      end do
      result%n_jackknife=m
      if(m<max(4,q/2)) then
      result%status=-2
      return
      end if
      meanjk=sum(jk(1:m,:),dim=1)/real(m,dp)
      allocate(result%gamma(q,q))
      result%gamma=0.0_dp
      do r=1,m
         d=jk(r,:)-meanjk
         result%gamma=result%gamma+spread(d,2,q)*spread(d,1,q)
      end do
      result%gamma=real(n,dp)*real(m-1,dp)/real(m,dp)*result%gamma
      allocate(ridge(q,q))
      ridge=result%gamma
      do i=1,q
      ridge(i,i)=ridge(i,i)+1.0e-8_dp*max(1.0_dp,abs(ridge(i,i)))
      end do
      call inverse_general(ridge,inv,info)
      allocate(result%weight(q,q),result%dwls_weight(q))
      result%weight=0.0_dp
      result%dwls_weight=0.0_dp
      if(info==0) then
         result%weight=inv
      else
         do i=1,q
         if(result%gamma(i,i)>1.0e-12_dp) result%weight(i,i)=1.0_dp/result%gamma(i,i)
         end do
      end if
      do i=1,q
      result%dwls_weight(i)=result%weight(i,i)
      end do
      result%status=0
   end subroutine categorical_wls_statistics

   subroutine compute_stats(data,ncat,offset,thresholds,cor,stats,info)
      integer,intent(in)::data(:,:)
      integer,allocatable,intent(out)::ncat(:),offset(:)
      real(dp),allocatable,intent(out)::thresholds(:),cor(:,:),stats(:)
      integer,intent(out)::info
      integer::p,nth,ncor,j,i,k,pos
      integer,allocatable::counts(:)
      real(dp),allocatable::th(:)
      p=size(data,2)
      allocate(ncat(p),offset(p+1))
      offset(1)=1
      nth=0
      info=0
      do j=1,p
         ncat(j)=maxval(data(:,j))
         if(minval(data(:,j))<1 .or. ncat(j)<2) then
         info=j
         return
         end if
         nth=nth+ncat(j)-1
         offset(j+1)=nth+1
      end do
      allocate(thresholds(nth))
      pos=1
      do j=1,p
         allocate(counts(ncat(j)))
         counts=0
         do i=1,size(data,1)
         counts(data(i,j))=counts(data(i,j))+1
         end do
         if(any(counts==0)) then
         info=100+j
         return
         end if
         th=ordinal_thresholds(counts)
         thresholds(pos:pos+size(th)-1)=th
         pos=pos+size(th)
         deallocate(counts)
      end do
      call polychoric_matrix(data,cor,info)
      if(info/=0) return
      ncor=p*(p-1)/2
      allocate(stats(nth+ncor))
      stats(1:nth)=thresholds
      k=nth
      do j=1,p-1
      do i=j+1,p
      k=k+1
      stats(k)=cor(i,j)
      end do
      end do
   end subroutine compute_stats
end module lavaan_categorical
