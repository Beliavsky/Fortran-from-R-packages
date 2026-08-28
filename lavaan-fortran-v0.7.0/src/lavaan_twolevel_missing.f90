module lavaan_twolevel_missing
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_get_free, ram_set_free, ram_sigma, ram_mu
   use lavaan_twolevel_ml, only : twolevel_ml_result
   use lavaan_linalg, only : inverse_spd, inverse_general, logdet_spd, chol_lower
   use lavaan_optimizer, only : bfgs_minimize
   use numderiv, only : hessian, nd_success
   implicit none
   private
   public :: fit_ram_twolevel_fiml
contains
   subroutine fit_ram_twolevel_fiml(within_template,within_map,between_template,between_map,data,cluster,result,fit_h1)
      type(ram_model),intent(in)::within_template,between_template
      type(ram_free_map),intent(in)::within_map,between_map
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::cluster(:)
      type(twolevel_ml_result),intent(out)::result
      logical,intent(in),optional::fit_h1
      type(ram_model)::wm,bm
      real(dp),allocatable::xw(:),xb(:),x(:),sw(:,:),sb(:,:),mu(:),hess(:,:),hi(:,:),xh1(:)
      real(dp)::fval,h1f,h1start
      integer::kw,kb,k,n,p,info,status,i,kh1
      logical::do_h1,h1conv
      n=size(data,1)
      p=size(data,2)
      do_h1=.true.
      if(present(fit_h1)) do_h1=fit_h1
      if(size(cluster)/=n .or. n<2 .or. p<1) then
      result%status=-1
      return
      end if
      if(size(within_template%observed)/=p .or. size(between_template%observed)/=p) then
      result%status=-2
      return
      end if
      result%ncluster=count_unique(cluster)
      if(result%ncluster<2) then
      result%status=-3
      return
      end if
      xw=ram_get_free(within_template,within_map)
      xb=ram_get_free(between_template,between_map)
      kw=size(xw)
      kb=size(xb)
      k=kw+kb
      x=[xw,xb]
      call bfgs_minimize(nll_h0,x,fval,result%converged,result%iterations,maxiter=1800,tol=1.0e-7_dp)
      result%par=x
      result%loglik=-fval
      result%aic=2.0_dp*fval+2.0_dp*real(k,dp)
      result%bic=2.0_dp*fval+log(real(n,dp))*real(k,dp)
      call set_models(x,wm,bm)
      call ram_sigma(wm,sw,info)
      if(info/=0) then
      result%status=info
      return
      end if
      call ram_sigma(bm,sb,info)
      if(info/=0) then
      result%status=info
      return
      end if
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
      call hessian(nll_h0,x,hess,status=status)
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
      if(do_h1) then
         call covariance_start(sw,sb,mu,xh1,info)
         if(info==0) then
            h1start=nll_h1(xh1)
            call bfgs_minimize(nll_h1,xh1,h1f,h1conv,i,maxiter=2200,tol=1.5e-7_dp)
            if(.not.(h1f<huge(1.0_dp)/10.0_dp) .or. h1f>h1start) then
               h1f=h1start
               h1conv=.false.
            end if
            result%h1_loglik=-h1f
            kh1=p*(p+1)+p
            result%chisq=max(0.0_dp,2.0_dp*(result%h1_loglik-result%loglik))
            result%df=real(max(0,kh1-k),dp)
            result%h1_converged=h1conv
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

      function nll_h0(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::ssw(:,:),ssb(:,:),muu(:)
         type(ram_model)::wmod,bmod
         integer::istat
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
         v=cluster_nll(data,cluster,ssw,ssb,muu,istat)
         if(istat/=0) v=huge(1.0_dp)/100.0_dp
      end function nll_h0

      function nll_h1(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::ssw(:,:),ssb(:,:),muu(:)
         integer::istat
         call unpack_covariance(z,p,ssw,ssb,muu,istat)
         if(istat/=0) then
         v=huge(1.0_dp)/100.0_dp
         return
         end if
         v=cluster_nll(data,cluster,ssw,ssb,muu,istat)
         if(istat/=0) v=huge(1.0_dp)/100.0_dp
      end function nll_h1
   end subroutine fit_ram_twolevel_fiml

   function cluster_nll(data,cluster,sw,sb,mu,info) result(v)
      real(dp),intent(in)::data(:,:),sw(:,:),sb(:,:),mu(:)
      integer,intent(in)::cluster(:)
      integer,intent(out)::info
      real(dp)::v
      integer,allocatable::labels(:),rows(:),obs(:)
      real(dp),allocatable::cfull(:,:),mfull(:),yfull(:),cobs(:,:),mobs(:),yobs(:),ci(:,:),d(:)
      integer::g,m,p,nv,r,s,i,j,a,b,idx,istat
      real(dp)::ld,pi2
      p=size(data,2)
      call unique_labels(cluster,labels)
      v=0.0_dp
      info=0
      pi2=log(2.0_dp*acos(-1.0_dp))
      do g=1,size(labels)
         rows=pack([(r,r=1,size(data,1))],cluster==labels(g))
         m=size(rows)
         allocate(cfull(m*p,m*p),mfull(m*p),yfull(m*p))
         cfull=0.0_dp
         do r=1,m
            do i=1,p
               a=(r-1)*p+i
               mfull(a)=mu(i)
               yfull(a)=data(rows(r),i)
               do s=1,m
                  do j=1,p
                     b=(s-1)*p+j
                     cfull(a,b)=sb(i,j)
                     if(r==s) cfull(a,b)=cfull(a,b)+sw(i,j)
                  end do
               end do
            end do
         end do
         obs=pack([(idx,idx=1,m*p)],.not.ieee_is_nan(yfull))
         nv=size(obs)
         if(nv>0) then
            allocate(cobs(nv,nv),mobs(nv),yobs(nv))
            do j=1,nv
            do i=1,nv
            cobs(i,j)=cfull(obs(i),obs(j))
            end do
            end do
            mobs=mfull(obs)
            yobs=yfull(obs)
            call inverse_spd(cobs,ci,istat)
            if(istat/=0) then
            info=istat
            return
            end if
            ld=logdet_spd(cobs,istat)
            if(istat/=0) then
            info=istat
            return
            end if
            d=yobs-mobs
            v=v+0.5_dp*(real(nv,dp)*pi2+ld+dot_product(d,matmul(ci,d)))
            deallocate(cobs,mobs,yobs)
         end if
         deallocate(cfull,mfull,yfull)
      end do
   end function cluster_nll

   subroutine covariance_start(sw,sb,mu,z,info)
      real(dp),intent(in)::sw(:,:),sb(:,:),mu(:)
      real(dp),allocatable,intent(out)::z(:)
      integer,intent(out)::info
      real(dp),allocatable::lw(:,:),lb(:,:),aw(:,:),ab(:,:)
      integer::p,k,i,j,iw,ib
      p=size(mu)
      aw=sw
      ab=sb
      do i=1,p
      aw(i,i)=aw(i,i)+1.0e-8_dp
      ab(i,i)=ab(i,i)+1.0e-8_dp
      end do
      call chol_lower(aw,lw,info)
      if(info/=0) return
      call chol_lower(ab,lb,info)
      if(info/=0) return
      k=p*(p+1)/2
      allocate(z(2*k+p))
      iw=0
      do i=1,p
      do j=1,i
      iw=iw+1
      if(i==j) then
      z(iw)=log(max(lw(i,j),1.0e-8_dp))
      else
      z(iw)=lw(i,j)
      end if
      end do
      end do
      ib=k
      do i=1,p
      do j=1,i
      ib=ib+1
      if(i==j) then
      z(ib)=log(max(lb(i,j),1.0e-8_dp))
      else
      z(ib)=lb(i,j)
      end if
      end do
      end do
      z(2*k+1:)=mu
      info=0
   end subroutine covariance_start

   subroutine unpack_covariance(z,p,sw,sb,mu,info)
      real(dp),intent(in)::z(:)
      integer,intent(in)::p
      real(dp),allocatable,intent(out)::sw(:,:),sb(:,:),mu(:)
      integer,intent(out)::info
      real(dp),allocatable::lw(:,:),lb(:,:)
      integer::k,i,j,pos
      k=p*(p+1)/2
      if(size(z)/=2*k+p) then
      info=-1
      return
      end if
      allocate(lw(p,p),lb(p,p))
      lw=0.0_dp
      lb=0.0_dp
      pos=0
      do i=1,p
      do j=1,i
      pos=pos+1
      if(i==j) then
      lw(i,j)=exp(z(pos))
      else
      lw(i,j)=z(pos)
      end if
      end do
      end do
      do i=1,p
      do j=1,i
      pos=pos+1
      if(i==j) then
      lb(i,j)=exp(z(pos))
      else
      lb(i,j)=z(pos)
      end if
      end do
      end do
      sw=matmul(lw,transpose(lw))
      sb=matmul(lb,transpose(lb))
      mu=z(2*k+1:)
      info=0
   end subroutine unpack_covariance

   integer function count_unique(x) result(nu)
      integer,intent(in)::x(:)
      integer,allocatable::u(:)
      call unique_labels(x,u)
      nu=size(u)
   end function count_unique

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
end module lavaan_twolevel_missing
