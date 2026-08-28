module lavaan_objectives
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_spd, logdet_spd, trace_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private
   public :: objective_ml, objective_gls, objective_uls, objective_wls, objective_dwls
   public :: mvn_loglik_complete, mvn_loglik_missing
contains
   function objective_ml(sigma, mu, data_cov, data_mean, meanstructure, info) result(fx)
      real(dp),intent(in)::sigma(:,:),mu(:),data_cov(:,:),data_mean(:)
      logical,intent(in)::meanstructure
      integer,intent(out)::info
      real(dp)::fx,ld,ld0
      real(dp),allocatable::si(:,:),w(:,:),d(:)
      integer::p
      call inverse_spd(sigma,si,info)
      if(info/=0) then
      fx=huge(1.0_dp)
      return
      end if
      ld=logdet_spd(sigma,info)
      ld0=logdet_spd(data_cov,info)
      p=size(sigma,1)
      if(info/=0) then
      fx=huge(1.0_dp)
      return
      end if
      allocate(w(p,p))
      w=data_cov
      if(meanstructure) then
         d=data_mean-mu
         w=w+spread(d,2,p)*spread(d,1,p)
      end if
      fx=ld+sum(w*si)-ld0-real(p,dp)
      if(fx<0 .and. abs(fx)<1e-10_dp) fx=0
   end function objective_ml

   function objective_gls(sigma,mu,data_cov,data_mean,meanstructure,info) result(fx)
      real(dp),intent(in)::sigma(:,:),mu(:),data_cov(:,:),data_mean(:)
      logical,intent(in)::meanstructure
      integer,intent(out)::info
      real(dp)::fx
      real(dp),allocatable::si(:,:),tmp(:,:),d(:)
      integer::p
      call inverse_spd(data_cov,si,info)
      if(info/=0) then
      fx=huge(1.0_dp)
      return
      end if
      tmp=matmul(si,data_cov-sigma)
      fx=0.5_dp*sum(tmp*transpose(tmp))
      if(meanstructure) then
      d=data_mean-mu
      p=size(d)
      fx=fx+sum(si*spread(d,2,p)*spread(d,1,p))
      end if
   end function objective_gls

   function objective_uls(est,obs) result(fx)
      real(dp),intent(in)::est(:),obs(:)
      real(dp)::fx
      fx=sum((obs-est)**2)
   end function objective_uls
   function objective_wls(est,obs,w) result(fx)
      real(dp),intent(in)::est(:),obs(:),w(:,:)
      real(dp)::fx
      real(dp),allocatable::d(:)
      d=obs-est
      fx=dot_product(d,matmul(w,d))
   end function objective_wls
   function objective_dwls(est,obs,wd) result(fx)
      real(dp),intent(in)::est(:),obs(:),wd(:)
      real(dp)::fx
      fx=sum((obs-est)**2*wd)
   end function objective_dwls

   function mvn_loglik_complete(x,mu,sigma,info) result(ll)
      real(dp),intent(in)::x(:,:),mu(:),sigma(:,:)
      integer,intent(out)::info
      real(dp)::ll,ld,q
      real(dp),allocatable::si(:,:),d(:)
      integer::i,p,n
      call inverse_spd(sigma,si,info)
      if(info/=0) then
      ll=-huge(1.0_dp)
      return
      end if
      ld=logdet_spd(sigma,info)
      p=size(x,2)
      n=size(x,1)
      ll=0
      do i=1,n
      d=x(i,:)-mu
      q=dot_product(d,matmul(si,d))
      ll=ll-0.5_dp*(real(p,dp)*log(2*acos(-1.0_dp))+ld+q)
      end do
   end function mvn_loglik_complete

   function mvn_loglik_missing(x,mu,sigma,info) result(ll)
      real(dp),intent(in)::x(:,:),mu(:),sigma(:,:)
      integer,intent(out)::info
      real(dp)::ll,ld,q
      integer::i,j,k,p,n,m
      integer,allocatable::idx(:)
      real(dp),allocatable::ss(:,:),si(:,:),d(:)
      n=size(x,1)
      p=size(x,2)
      ll=0
      info=0
      allocate(idx(p))
      do i=1,n
         m=0
         do j=1,p
         if(.not.ieee_is_nan(x(i,j))) then
         m=m+1
         idx(m)=j
         end if
         end do
         if(m==0) cycle
         allocate(ss(m,m),d(m))
         do j=1,m
         d(j)=x(i,idx(j))-mu(idx(j))
         do k=1,m
         ss(j,k)=sigma(idx(j),idx(k))
         end do
         end do
         call inverse_spd(ss,si,info)
         if(info/=0) then
         ll=-huge(1.0_dp)
         return
         end if
         ld=logdet_spd(ss,info)
         q=dot_product(d,matmul(si,d))
         ll=ll-0.5_dp*(real(m,dp)*log(2*acos(-1.0_dp))+ld+q)
         deallocate(ss,d,si)
      end do
   end function mvn_loglik_missing
end module lavaan_objectives
