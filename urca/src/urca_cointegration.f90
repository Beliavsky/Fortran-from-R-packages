module urca_cointegration
   use urca_kinds, only : dp
   use urca_types, only : johansen_result, po_result, vecm_result, lm_result
   use urca_linalg, only : invert_matrix, invert_spd, chol_lower, symmetric_eigen, trace_matrix
   use urca_regression, only : lm_fit, lm_fit_multi
   implicit none
   private
   integer, parameter, public :: JO_NONE=0, JO_CONST=1, JO_TREND=2
   integer, parameter, public :: JO_LONGRUN=1, JO_TRANSITORY=2
   integer, parameter, public :: JO_EIGEN=1, JO_TRACE=2
   integer, parameter, public :: PO_NONE=0, PO_CONST=1, PO_TREND=2
   integer, parameter, public :: PO_PU=1, PO_PZ=2
   public :: johansen_test, phillips_ouliaris, cajools_fit, alphaols_fit, cajorls_fit
contains
   subroutine seasonal_dummies(n,season,d)
      integer,intent(in)::n,season
      real(dp),allocatable,intent(out)::d(:,:)
      integer::i,j,k
      if(season<=1)then
      allocate(d(n,0))
      return
      end if
      allocate(d(n,season-1))
      d=0.0_dp
      do i=1,n
         k=mod(i-1,season)+1
         do j=1,season-1
            d(i,j)=-1.0_dp/real(season,dp)
            if(k==j)d(i,j)=d(i,j)+1.0_dp
         end do
      end do
   end subroutine seasonal_dummies

   subroutine build_johansen_matrices(x,k,ecdet,spec,z0,z1,zk,season,dumvar,info)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::k,ecdet,spec
      real(dp),allocatable,intent(out)::z0(:,:),z1(:,:),zk(:,:)
      integer,intent(in),optional::season
      real(dp),intent(in),optional::dumvar(:,:)
      integer,intent(out)::info
      real(dp),allocatable::dx(:,:),sd(:,:)
      integer::nt,p,nobs,ns,nd,nbase,nz1,i,j,c,ss
      nt=size(x,1)
      p=size(x,2)
      nobs=nt-k
      info=0
      if(k<2 .or. nobs<=p)then
      info=-1
      allocate(z0(0,0),z1(0,0),zk(0,0))
      return
      end if
      allocate(dx(nt-1,p))
      dx=x(2:nt,:)-x(1:nt-1,:)
      allocate(z0(nobs,p))
      z0=dx(k:nt-1,:)
      ss=0
      if(present(season))ss=season
      call seasonal_dummies(nt,ss,sd)
      ns=size(sd,2)
      nd=0
      if(present(dumvar))then
      if(size(dumvar,1)/=nt)then
      info=-2
      return
      end if
      nd=size(dumvar,2)
      end if
      nbase=(k-1)*p
      if(ecdet==JO_CONST)then
      nz1=ns+nd+nbase
      else
      nz1=1+ns+nd+nbase
      end if
      allocate(z1(nobs,nz1))
      c=0
      if(ecdet/=JO_CONST)then
      c=c+1
      z1(:,c)=1.0_dp
      end if
      if(ns>0)then
      z1(:,c+1:c+ns)=sd(k+1:nt,:)
      c=c+ns
      end if
      if(nd>0)then
      z1(:,c+1:c+nd)=dumvar(k+1:nt,:)
      c=c+nd
      end if
      do j=1,k-1
         z1(:,c+(j-1)*p+1:c+j*p)=dx(k-j:nt-1-j,:)
      end do
      if(ecdet==JO_CONST .or. ecdet==JO_TREND)then
      allocate(zk(nobs,p+1))
      else
      allocate(zk(nobs,p))
      end if
      if(spec==JO_LONGRUN)then
         zk(:,1:p)=x(1:nobs,:)
         if(size(zk,2)>p)then
            if(ecdet==JO_CONST)then
            zk(:,p+1)=1.0_dp
            else
            do i=1,nobs
            zk(i,p+1)=real(i,dp)
            end do
            end if
         end if
      else
         zk(:,1:p)=x(k:nt-1,:)
         if(size(zk,2)>p)then
            if(ecdet==JO_CONST)then
            zk(:,p+1)=1.0_dp
            else
            do i=1,nobs
            zk(i,p+1)=real(k+i-1,dp)
            end do
            end if
         end if
      end if
   end subroutine build_johansen_matrices

   function johansen_test(x,type,ecdet,k,spec,season,dumvar) result(out)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::type,ecdet,k,spec
      integer,intent(in),optional::season
      real(dp),intent(in),optional::dumvar(:,:)
      type(johansen_result)::out
      real(dp),allocatable::z0(:,:),z1(:,:),zk(:,:),m11(:,:),m11i(:,:),m01(:,:),mk1(:,:),m10(:,:),m1k(:,:),r0(:,:), &
         & rk(:,:)
      real(dp),allocatable::s00(:,:),s0k(:,:),sk0(:,:),skk(:,:),l(:,:),linv(:,:),s00i(:,:),a(:,:),eval(:),evec(:,:), &
         & v(:,:),vn(:,:),tmp(:,:),vv(:,:),vvi(:,:),skki(:,:)
      real(dp),allocatable::pi(:,:),delta(:,:),gamma(:,:),m01m(:,:),mk1m(:,:)
      integer::info,n,p,q,i,j,r,nobs
      real(dp)::lam
      call build_johansen_matrices(x,k,ecdet,spec,z0,z1,zk,season,dumvar,info)
      if(info/=0)then
      out%info=info
      return
      end if
      nobs=size(z0,1)
      p=size(z0,2)
      q=size(zk,2)
      m11=matmul(transpose(z1),z1)/real(nobs,dp)
      call invert_spd(m11,m11i,info)
      if(info/=0)then
      out%info=info
      return
      end if
      m01=matmul(transpose(z0),z1)/real(nobs,dp)
      m10=transpose(m01)
      mk1=matmul(transpose(zk),z1)/real(nobs,dp)
      m1k=transpose(mk1)
      r0=z0-matmul(z1,matmul(m11i,m10))
      rk=zk-matmul(z1,matmul(m11i,m1k))
      s00=matmul(transpose(r0),r0)/real(nobs,dp)
      s0k=matmul(transpose(r0),rk)/real(nobs,dp)
      sk0=transpose(s0k)
      skk=matmul(transpose(rk),rk)/real(nobs,dp)
      call chol_lower(skk,l,info)
      if(info/=0)then
      out%info=100+info
      return
      end if
      call invert_matrix(l,linv,info)
      if(info/=0)then
      out%info=200+info
      return
      end if
      call invert_spd(s00,s00i,info)
      if(info/=0)then
      out%info=300+info
      return
      end if
      a=matmul(linv,matmul(sk0,matmul(s00i,matmul(s0k,transpose(linv)))))
      call symmetric_eigen(a,eval,evec,info,.true.)
      if(info/=0)then
      out%info=400+info
      return
      end if
      v=matmul(transpose(linv),evec)
      allocate(vn(q,q))
      vn=v
      do j=1,q
      if(abs(vn(1,j))>1e-14_dp)vn(:,j)=vn(:,j)/vn(1,j)
      end do
      vv=matmul(transpose(vn),matmul(skk,vn))
      call invert_matrix(vv,vvi,info)
      if(info/=0)then
      out%info=500+info
      return
      end if
      out%w=matmul(s0k,matmul(vn,vvi))
      call invert_spd(skk,skki,info)
      if(info/=0)call invert_matrix(skk,skki,info)
      pi=matmul(s0k,skki)
      delta=s00-matmul(s0k,matmul(vn,matmul(vvi,matmul(transpose(vn),sk0))))
      gamma=matmul(m01,m11i)-matmul(pi,matmul(mk1,m11i))
      allocate(out%teststat(p))
      out%teststat=0.0_dp
      do i=1,p
         r=p-i
         if(type==JO_TRACE)then
            do j=r+1,q
            lam=max(0.0_dp,min(1.0_dp-1e-14_dp,eval(j)))
            out%teststat(i)=out%teststat(i)-real(nobs,dp)*log(1.0_dp-lam)
            end do
         else
            lam=max(0.0_dp,min(1.0_dp-1e-14_dp,eval(r+1)))
            out%teststat(i)=-real(nobs,dp)*log(1.0_dp-lam)
         end if
      end do
      call johansen_critical_values(p,type,ecdet,out%critical_values)
      out%x=x
      out%z0=z0
      out%z1=z1
      out%zk=zk
      out%lambda=eval
      out%vorg=v
      out%v=vn
      out%pi=pi
      out%delta=delta
      out%gamma=gamma
      out%r0=r0
      out%rk=rk
      out%p=p
      out%lag=k
      out%ecdet=ecdet
      out%spec=spec
      out%info=0
   end function johansen_test

   subroutine johansen_critical_values(p,type,ecdet,cv)
      integer,intent(in)::p,type,ecdet
      real(dp),allocatable,intent(out)::cv(:,:)
      real(dp),parameter::none_e(11,3)=reshape([6.5_dp,12.91_dp,18.9_dp,24.78_dp,30.84_dp,36.25_dp,42.06_dp, &
         & 48.43_dp,54.01_dp,59._dp,65.07_dp,8.18_dp,14.90_dp,21.07_dp,27.14_dp,33.32_dp,39.43_dp,44.91_dp,51.07_dp, &
         & 57.00_dp,62.42_dp,68.27_dp,11.65_dp,19.19_dp,25.75_dp,32.14_dp,38.78_dp,44.59_dp,51.30_dp,57.07_dp, &
         & 63.37_dp,68.61_dp,74.36_dp],[11,3])
      real(dp),parameter::const_e(11,3)=reshape([7.52_dp,13.75_dp,19.77_dp,25.56_dp,31.66_dp,37.45_dp,43.25_dp, &
         & 48.91_dp,54.35_dp,60.25_dp,66.02_dp,9.24_dp,15.67_dp,22.00_dp,28.14_dp,34.40_dp,40.30_dp,46.45_dp, &
         & 52.00_dp,57.42_dp,63.57_dp,69.74_dp,12.97_dp,20.20_dp,26.81_dp,33.24_dp,39.79_dp,46.82_dp,51.91_dp, &
         & 57.95_dp,63.71_dp,69.94_dp,76.63_dp],[11,3])
      real(dp),parameter::trend_e(11,3)=reshape([10.49_dp,16.85_dp,23.11_dp,29.12_dp,34.75_dp,40.91_dp,46.32_dp, &
         & 52.16_dp,57.87_dp,63.18_dp,69.26_dp,12.25_dp,18.96_dp,25.54_dp,31.46_dp,37.52_dp,43.97_dp,49.42_dp, &
         & 55.50_dp,61.29_dp,66.23_dp,72.72_dp,16.26_dp,23.65_dp,30.34_dp,36.65_dp,42.36_dp,49.51_dp,54.71_dp, &
         & 62.46_dp,67.88_dp,73.73_dp,79.23_dp],[11,3])
      real(dp),parameter::none_t(11,3)=reshape([6.50_dp,15.66_dp,28.71_dp,45.23_dp,66.49_dp,85.18_dp,118.99_dp, &
         & 151.38_dp,186.54_dp,226.34_dp,269.53_dp,8.18_dp,17.95_dp,31.52_dp,48.28_dp,70.6_dp,90.39_dp,124.25_dp, &
         & 157.11_dp,192.84_dp,232.49_dp,277.39_dp,11.65_dp,23.52_dp,37.22_dp,55.43_dp,78.87_dp,104.20_dp,136.06_dp, &
         & 168.92_dp,204.79_dp,246.27_dp,292.65_dp],[11,3])
      real(dp),parameter::const_t(11,3)=reshape([7.52_dp,17.85_dp,32.00_dp,49.65_dp,71.86_dp,97.18_dp,126.58_dp, &
         & 159.48_dp,196.37_dp,236.54_dp,282.45_dp,9.24_dp,19.96_dp,34.91_dp,53.12_dp,76.07_dp,102.14_dp,131.70_dp, &
         & 165.58_dp,202.92_dp,244.15_dp,291.40_dp,12.97_dp,24.60_dp,41.07_dp,60.16_dp,84.45_dp,111.01_dp,143.09_dp, &
         & 177.20_dp,215.74_dp,257.68_dp,307.64_dp],[11,3])
      real(dp),parameter::trend_t(11,3)=reshape([10.49_dp,22.76_dp,39.06_dp,59.14_dp,83.20_dp,110.42_dp,141.01_dp, &
         & 176.67_dp,215.17_dp,256.72_dp,303.13_dp,12.25_dp,25.32_dp,42.44_dp,62.99_dp,87.31_dp,114.90_dp,146.76_dp, &
         & 182.82_dp,222.21_dp,263.42_dp,310.81_dp,16.26_dp,30.45_dp,48.45_dp,70.05_dp,96.58_dp,124.75_dp,158.49_dp, &
         & 196.08_dp,234.41_dp,279.07_dp,327.45_dp],[11,3])
      integer::i
      if(p>11)then
      allocate(cv(0,0))
      return
      end if
      allocate(cv(p,3))
      do i=1,p
         if(type==JO_EIGEN)then
            select case(ecdet)
            case(JO_NONE)
            cv(i,:)=none_e(p-i+1,:)
            case(JO_CONST)
            cv(i,:)=const_e(p-i+1,:)
            case default
            cv(i,:)=trend_e(p-i+1,:)
            end select
         else
            select case(ecdet)
            case(JO_NONE)
            cv(i,:)=none_t(p-i+1,:)
            case(JO_CONST)
            cv(i,:)=const_t(p-i+1,:)
            case default
            cv(i,:)=trend_t(p-i+1,:)
            end select
         end if
      end do
   end subroutine johansen_critical_values

   function phillips_ouliaris(z,demean,lag_mode,type) result(out)
      real(dp),intent(in)::z(:,:)
      integer,intent(in)::demean,lag_mode,type
      type(po_result)::out
      integer::nt,m,n,lmax,k,p,info
      real(dp),allocatable::zl(:,:),zr(:,:),x(:,:),beta(:,:),res(:,:),sigma(:,:),omega(:,:),mzz(:,:),mzzi(:,:), &
         & resu(:),xu(:,:),w21(:),o22(:,:),o22i(:,:)
      real(dp)::wsl,w112,ssqr
      type(lm_result)::ufit
      nt=size(z,1)
      m=size(z,2)
      n=nt-1
      if(m<2 .or. m>6 .or. n<10)then
      out%info=-1
      return
      end if
      lmax=int(merge(12.0_dp,4.0_dp,lag_mode==2)*(real(n,dp)/100.0_dp)**0.25_dp)
      allocate(zl(n,m),zr(n,m))
      zl=z(2:,:)
      zr=z(:nt-1,:)
      select case(demean)
      case(PO_NONE)
      p=m
      allocate(x(n,p))
      x=zr
      case(PO_CONST)
      p=m+1
      allocate(x(n,p))
      x(:,1)=1
      x(:,2:)=zr
      case default
      p=m+2
      allocate(x(n,p))
      x(:,1)=1
      x(:,2:m+1)=zr
      do k=1,n
      x(k,p)=real(k,dp)
      end do
      end select
      call lm_fit_multi(x,zl,beta,res,sigma,info)
      if(info/=0)then
      out%info=info
      return
      end if
      omega=matmul(transpose(res),res)/real(n,dp)
      do k=1,lmax
         wsl=1.0_dp-real(k,dp)/real(lmax+1,dp)
         omega=omega+wsl*(matmul(transpose(res(k+1:n,:)),res(1:n-k,:))+matmul(transpose(res(1:n-k,:)),res(k+1:n, &
            & :)))/real(n,dp)
      end do
      if(type==PO_PZ)then
         mzz=matmul(transpose(zl),zl)/real(n,dp)
         call invert_matrix(mzz,mzzi,info)
         out%statistic=real(n,dp)*trace_matrix(matmul(omega,mzzi))
         out%residuals=res
      else
         select case(demean)
         case(PO_NONE)
         allocate(xu(nt,m-1))
         xu=z(:,2:)
         case(PO_CONST)
         allocate(xu(nt,m))
         xu(:,1)=1
         xu(:,2:)=z(:,2:)
         case default
         allocate(xu(nt,m+1))
         xu(:,1)=1
         xu(:,2:m)=z(:,2:)
         do k=1,nt
         xu(k,m+1)=real(k,dp)
         end do
         end select
         ufit=lm_fit(xu,z(:,1))
         allocate(w21(m-1))
         w21=omega(2:m,1)
         o22=omega(2:m,2:m)
         call invert_matrix(o22,o22i,info)
         w112=omega(1,1)-dot_product(w21,matmul(o22i,w21))
         ssqr=sum(ufit%residuals**2)/real(n,dp)
         out%statistic=real(n,dp)*w112/ssqr
         allocate(out%residuals(size(ufit%residuals),1))
         out%residuals(:,1)=ufit%residuals
      end if
      call po_critical(m,demean,type,out%critical_values)
      out%lags=lmax
   end function phillips_ouliaris

   subroutine po_critical(m,demean,type,cv)
      integer,intent(in)::m,demean,type
      real(dp),intent(out)::cv(3)
      real(dp),parameter::pu(5,3,3)=reshape([20.3933_dp,26.7022_dp,33.5359_dp,39.2826_dp,44.3725_dp,25.9711_dp, &
         & 32.9392_dp,40.1220_dp,46.2691_dp,51.8614_dp,38.3413_dp,46.4097_dp,55.7341_dp,63.2149_dp,69.4939_dp, &
         & 27.8536_dp,33.6955_dp,39.6949_dp,45.3308_dp,50.3537_dp,33.713_dp,40.5252_dp,46.7281_dp,53.2502_dp, &
         & 57.7855_dp,48.0021_dp,53.8731_dp,63.4128_dp,71.5214_dp,76.7705_dp,41.2488_dp,46.1061_dp,52.0015_dp, &
         & 57.3667_dp,61.6155_dp,48.8439_dp,53.8300_dp,60.2384_dp,65.8706_dp,70.7416_dp,65.1714_dp,69.2629_dp, &
         & 78.3470_dp,84.5480_dp,91.0392_dp],[5,3,3])
      real(dp),parameter::pz(5,3,3)=reshape([33.9267_dp,62.1436_dp,99.2664_dp,143.0775_dp,195.6202_dp,40.8217_dp, &
         & 71.2751_dp,109.7426_dp,155.8019_dp,210.2910_dp,55.1911_dp,89.6679_dp,131.5716_dp,180.4845_dp,237.7723_dp, &
         & 47.5877_dp,80.2034_dp,120.3035_dp,168.8572_dp,225.2303_dp,55.2202_dp,89.7619_dp,132.2207_dp,182.0749_dp, &
         & 241.3316_dp,71.9273_dp,109.4525_dp,153.4504_dp,209.8054_dp,270.5018_dp,71.9586_dp,113.4929_dp, &
         & 163.1050_dp,219.5098_dp,284.0100_dp,81.3812_dp,124.3933_dp,175.9902_dp,234.2865_dp,301.0949_dp, &
         & 102.0167_dp,145.8644_dp,201.0905_dp,264.4988_dp,335.9054_dp],[5,3,3])
      integer::d
      d=demean+1
      if(type==PO_PU)then
      cv=pu(m-1,:,d)
      else
      cv=pz(m-1,:,d)
      end if
   end subroutine po_critical

   function cajools_fit(z) result(out)
      type(johansen_result),intent(in)::z
      type(vecm_result)::out
      real(dp),allocatable::x(:,:),b(:,:),r(:,:),s(:,:)
      integer::info,p1,p2
      p1=size(z%z1,2)
      p2=size(z%zk,2)
      allocate(x(size(z%z0,1),p1+p2))
      x(:,1:p1)=z%z1
      x(:,p1+1:)=z%zk
      call lm_fit_multi(x,z%z0,b,r,s,info)
      out%coefficients=b
      out%residuals=r
      out%sigma=s
      out%info=info
   end function cajools_fit

   function alphaols_fit(z,rank) result(out)
      type(johansen_result),intent(in)::z
      integer,intent(in),optional::rank
      type(vecm_result)::out
      integer::rnk,info
      real(dp),allocatable::x(:,:),b(:,:),r(:,:),s(:,:)
      rnk=size(z%v,2)
      if(present(rank))rnk=min(rank,rnk)
      x=matmul(z%rk,z%v(:,1:rnk))
      call lm_fit_multi(x,z%r0,b,r,s,info)
      out%coefficients=b
      out%residuals=r
      out%sigma=s
      out%rank=rnk
      out%info=info
   end function alphaols_fit

   function cajorls_fit(z,rank) result(out)
      type(johansen_result),intent(in)::z
      integer,intent(in)::rank
      type(vecm_result)::out
      real(dp),allocatable::beta(:,:),top(:,:),topi(:,:),bn(:,:),ect(:,:),x(:,:),b(:,:),r(:,:),s(:,:)
      integer::info,n
      if(rank<1 .or. rank>=z%p)then
      out%info=-1
      return
      end if
      beta=z%v(:,1:rank)
      top=beta(1:rank,:)
      call invert_matrix(top,topi,info)
      if(info/=0)then
      out%info=info
      return
      end if
      bn=matmul(beta,topi)
      ect=matmul(z%zk,bn)
      allocate(x(size(ect,1),rank+size(z%z1,2)))
      x(:,1:rank)=ect
      x(:,rank+1:)=z%z1
      call lm_fit_multi(x,z%z0,b,r,s,info)
      out%beta=bn
      out%coefficients=b
      out%residuals=r
      out%sigma=s
      out%rank=rank
      out%info=info
   end function cajorls_fit
end module urca_cointegration
