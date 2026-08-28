module lavaan_predict
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model
   use lavaan_linalg, only : inverse_general, inverse_spd
   implicit none
   private
   public :: factor_scores_regression, ram_all_moments
contains
   subroutine ram_all_moments(model,cov_all,mean_all,info)
      type(ram_model),intent(in)::model
      real(dp),allocatable,intent(out)::cov_all(:,:),mean_all(:)
      integer,intent(out)::info
      real(dp),allocatable::ia(:,:),inv(:,:)
      integer::n,i
      n=size(model%a,1)
      allocate(ia(n,n))
      ia=-model%a
      do i=1,n
      ia(i,i)=ia(i,i)+1
      end do
      call inverse_general(ia,inv,info)
      if(info/=0) then
      allocate(cov_all(n,n),mean_all(n))
      cov_all=0
      mean_all=0
      return
      end if
      cov_all=matmul(inv,matmul(model%s,transpose(inv)))
      if(allocated(model%m)) then
      mean_all=matmul(inv,model%m)
      else
      allocate(mean_all(n))
      mean_all=0
      end if
   end subroutine ram_all_moments

   subroutine factor_scores_regression(model,x,scores,info)
      type(ram_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable,intent(out)::scores(:,:)
      integer,intent(out)::info
      real(dp),allocatable::covall(:,:),muall(:),syy(:,:),syi(:,:),gain(:,:),d(:)
      integer,allocatable::latent(:)
      logical,allocatable::isobs(:)
      integer::nall,p,q,i,j,k,n
      call ram_all_moments(model,covall,muall,info)
      if(info/=0) then
      allocate(scores(0,0))
      return
      end if
      nall=size(model%a,1)
      p=size(model%observed)
      allocate(isobs(nall))
      isobs=.false.
      isobs(model%observed)=.true.
      q=count(.not.isobs)
      allocate(latent(q))
      k=0
      do i=1,nall
      if(.not.isobs(i)) then
      k=k+1
      latent(k)=i
      end if
      end do
      allocate(syy(p,p),syi(p,p))
      do j=1,p
      do i=1,p
      syy(i,j)=covall(model%observed(i),model%observed(j))
      end do
      end do
      call inverse_spd(syy,syi,info)
      if(info/=0) then
      allocate(scores(0,0))
      return
      end if
      allocate(gain(q,p))
      do j=1,p
      do i=1,q
      gain(i,j)=0
      do k=1,p
         gain(i,j)=gain(i,j)+covall(latent(i),model%observed(k))*syi(k,j)
         end do
         end do
         end do
      n=size(x,1)
      allocate(scores(n,q),d(p))
      do i=1,n
         d=x(i,:)-muall(model%observed)
         scores(i,:)=muall(latent)+matmul(gain,d)
      end do
   end subroutine factor_scores_regression
end module lavaan_predict
