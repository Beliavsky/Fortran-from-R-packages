module lavaan_twolevel_ml
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_get_free, ram_set_free, ram_sigma, ram_mu
   use lavaan_linalg, only : inverse_spd, inverse_general, logdet_spd
   use lavaan_optimizer, only : bfgs_minimize
   use numderiv, only : hessian, nd_success
   implicit none
   private

   type, public :: twolevel_ml_result
      real(dp), allocatable :: par(:), se(:), vcov(:, :), sigma_within(:, :), sigma_between(:, :), mu(:), icc(:)
      real(dp) :: loglik=-huge(1.0_dp), aic=huge(1.0_dp), bic=huge(1.0_dp)
      real(dp) :: h1_loglik=-huge(1.0_dp), chisq=huge(1.0_dp), df=0.0_dp
      integer :: ncluster=0, iterations=0, status=0
      logical :: converged=.false., h1_converged=.false.
   end type twolevel_ml_result

   public :: fit_ram_twolevel_ml

contains

   subroutine fit_ram_twolevel_ml(within_template,within_map,between_template,between_map,data,cluster,result)
      type(ram_model),intent(in)::within_template,between_template
      type(ram_free_map),intent(in)::within_map,between_map
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::cluster(:)
      type(twolevel_ml_result),intent(out)::result
      type(ram_model)::wm,bm
      real(dp),allocatable::xw(:),xb(:),x(:),hess(:,:),hi(:,:),sw(:,:),sb(:,:),mu(:)
      integer,allocatable::labels(:)
      real(dp)::fval
      integer::kw,kb,k,n,p,info,status,i
      n=size(data,1)
      p=size(data,2)
      if(size(cluster)/=n .or. n<2 .or. p<1) then
      result%status=-1
      return
      end if
      if(size(within_template%observed)/=p .or. size(between_template%observed)/=p) then
      result%status=-2
      return
      end if
      call unique_labels(cluster,labels)
      result%ncluster=size(labels)
      if(result%ncluster<2) then
      result%status=-3
      return
      end if
      xw=ram_get_free(within_template,within_map)
      xb=ram_get_free(between_template,between_map)
      kw=size(xw)
      kb=size(xb)
      k=kw+kb
      allocate(x(k))
      x=[xw,xb]
      call bfgs_minimize(nll,x,fval,result%converged,result%iterations,maxiter=1800,tol=1.0e-7_dp)
      result%par=x
      result%loglik=-fval
      result%aic=2.0_dp*fval+2.0_dp*real(k,dp)
      result%bic=2.0_dp*fval+log(real(n,dp))*real(k,dp)
      call set_models(x,wm,bm)
      call ram_sigma(wm,sw,info)
      call ram_sigma(bm,sb,info)
      call ram_mu(bm,mu,info)
      if(info/=0) then
      result%status=info
      return
      end if
      result%sigma_within=sw
      result%sigma_between=sb
      result%mu=mu
      allocate(result%icc(p))
      do i=1,p
         result%icc(i)=sb(i,i)/max(sb(i,i)+sw(i,i),tiny(1.0_dp))
         result%icc(i)=max(0.0_dp,min(1.0_dp,result%icc(i)))
      end do
      call hessian(nll,x,hess,status=status)
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(status==nd_success) then
         call inverse_general(hess,hi,info)
         if(info==0) then
            result%vcov=hi
            do i=1,k
            if(hi(i,i)>=0.0_dp) result%se(i)=sqrt(hi(i,i))
            end do
         end if
      end if
      result%status=0
   contains
      subroutine set_models(z,wmod,bmod)
         real(dp),intent(in)::z(:)
         type(ram_model),intent(out)::wmod,bmod
         wmod=within_template
         bmod=between_template
         if(kw>0) call ram_set_free(wmod,within_map,z(1:kw))
         if(kb>0) call ram_set_free(bmod,between_map,z(kw+1:kw+kb))
      end subroutine set_models

      function nll(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::ssw(:,:),ssb(:,:),muu(:),wi(:,:),cmat(:,:),ci(:,:),yg(:,:),bar(:),d(:)
         real(dp)::ldw,ldc,q,pi2
         integer::g,m,r,idx,istat
         type(ram_model)::wmod,bmod
         call set_models(z,wmod,bmod)
         call ram_sigma(wmod,ssw,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         call ram_sigma(bmod,ssb,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         call ram_mu(bmod,muu,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         call inverse_spd(ssw,wi,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         ldw=logdet_spd(ssw,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         pi2=log(2.0_dp*acos(-1.0_dp))
         v=0.0_dp
         do g=1,size(labels)
            m=count(cluster==labels(g))
            allocate(yg(m,p))
            idx=0
            do r=1,n
            if(cluster(r)==labels(g)) then
            idx=idx+1
            yg(idx,:)=data(r,:)
            end if
            end do
            bar=sum(yg,dim=1)/real(m,dp)
            q=0.0_dp
            do r=1,m
            d=yg(r,:)-bar
            q=q+dot_product(d,matmul(wi,d))
            end do
            cmat=ssw+real(m,dp)*ssb
            call inverse_spd(cmat,ci,istat)
            if(istat/=0) then
            v=huge(1.0_dp)/100.0_dp
            return
            end if
            ldc=logdet_spd(cmat,istat)
            d=bar-muu
            q=q+real(m,dp)*dot_product(d,matmul(ci,d))
            v=v+0.5_dp*(real(m*p,dp)*pi2+real(m-1,dp)*ldw+ldc+q)
            deallocate(yg)
         end do
      end function nll
   end subroutine fit_ram_twolevel_ml

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
end module lavaan_twolevel_ml
