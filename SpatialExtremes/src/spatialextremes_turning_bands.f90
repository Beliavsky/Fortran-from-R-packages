module spatialextremes_turning_bands
   use spatialextremes_base, only: dp,pi,exp_rand,chisq_rand
   use r_compat, only: runif1,rnorm1,rgamma,rbeta
   implicit none
   private
   integer,parameter,public :: TBM_MATERN=1,TBM_CAUCHY=2,TBM_POWEREXP=3
   integer,parameter,public :: TBM_BESSEL=4,TBM_GAUSSIAN=5,TBM_FBM=6
   public :: van_der_corput_lines,rotate_lines,simulate_gaussian_process_tbm
   public :: simulate_schlather_tbm,simulate_geomgauss_tbm,simulate_extremalt_tbm
contains
   function van_der_corput_lines(n) result(lines)
      ! Upstream randomlines.c Van der Corput construction in R^3.
      integer,intent(in)::n
      real(dp)::lines(n,3),base,u,v
      integer::i,k,r
      do i=0,n-1
         k=i
         u=0.0_dp
         base=2.0_dp
         do while(k>0)
         r=mod(k,2)
         u=u+real(r,dp)/base
         base=2.0_dp*base
         k=k/2
         end do
         k=i
         v=0.0_dp
         base=3.0_dp
         do while(k>0)
         r=mod(k,3)
         v=v+real(r,dp)/base
         base=3.0_dp*base
         k=k/3
         end do
         lines(i+1,1)=cos(2.0_dp*pi*u)*sqrt(max(0.0_dp,1.0_dp-v*v))
         lines(i+1,2)=sin(2.0_dp*pi*u)*sqrt(max(0.0_dp,1.0_dp-v*v))
         lines(i+1,3)=v
      end do
   end function van_der_corput_lines

   subroutine rotate_lines(lines,axis,angle)
      real(dp),intent(inout)::lines(:,:)
      real(dp),intent(in)::axis(3),angle
      real(dp)::p(3),parallel(3),b(3),c(3),r,ca,sa
      integer::i
      ca=cos(angle)
      sa=sin(angle)
      do i=1,size(lines,1)
         p=lines(i,:)
         parallel=dot_product(p,axis)*axis
         b=p-parallel
         r=sqrt(sum(b*b))
         if(r<=tiny(1.0_dp))cycle
         b=b/r
         c=[axis(2)*b(3)-axis(3)*b(2),axis(3)*b(1)-axis(1)*b(3),axis(1)*b(2)-axis(2)*b(1)]
         lines(i,:)=parallel+r*(ca*b+sa*c)
      end do
   end subroutine rotate_lines

   function simulate_gaussian_process_tbm(n,coord,covmod,nugget,sill,range,smooth,nlines,grid,info) result(ans)
      ! Turning-band Gaussian-process simulator translated from upstream
      ! turningbands.c.  For grid=.true., coord(:,d) supplies the axis
      ! coordinates and the Cartesian product is returned flattened.
      integer,intent(in)::n,covmod
      real(dp),intent(in)::coord(:,:),nugget,sill,range,smooth
      integer,intent(in),optional::nlines
      logical,intent(in),optional::grid
      integer,intent(out),optional::info
      real(dp),allocatable::ans(:,:)
      real(dp),allocatable::points(:,:),lines(:,:)
      real(dp)::axis(3),angle,nrm,eprod,freq,phase,g,u1,u2,r,theta,norm_const
      integer::nl,model,dim,neff,i,j,k,istat
      logical::isgrid

      istat=0
      dim=size(coord,2)
      nl=1000
      if(present(nlines))nl=nlines
      isgrid=.false.
      if(present(grid))isgrid=grid
      if(n<0.or.nl<=0.or.(dim/=2.and.dim/=3).or.range<=0.0_dp.or.sill<0.0_dp.or.nugget<0.0_dp)then
         allocate(ans(max(0,n),0))
         istat=1
         if(present(info))info=istat
         return
      end if
      model=covmod
      if(model==TBM_POWEREXP.and.abs(smooth-2.0_dp)<=10*epsilon(1.0_dp))model=TBM_GAUSSIAN
      if(model<TBM_MATERN.or.model>TBM_FBM)then
         allocate(ans(n,0))
         istat=2
         if(present(info))info=istat
         return
      end if
      if((model==TBM_MATERN.or.model==TBM_CAUCHY).and.smooth<=0.0_dp)then
         allocate(ans(n,0))
         istat=3
         if(present(info))info=istat
         return
      end if
      if(model==TBM_POWEREXP.and.(smooth<=0.0_dp.or.smooth>2.0_dp))then
         allocate(ans(n,0))
         istat=3
         if(present(info))info=istat
         return
      end if
      if(model==TBM_BESSEL.and.smooth<=0.5_dp)then
         allocate(ans(n,0))
         istat=3
         if(present(info))info=istat
         return
      end if
      if(model==TBM_FBM.and.(smooth<=0.0_dp.or.smooth>=2.0_dp))then
         allocate(ans(n,0))
         istat=4
         if(present(info))info=istat
         return
      end if

      if(isgrid)then
         call make_grid_points(coord,points)
      else
         points=coord
      end if
      neff=size(points,1)
      allocate(ans(n,neff))
      ans=0.0_dp
      points=points/range
      lines=van_der_corput_lines(nl)

      do i=1,n
         do
            axis=[runif1()-0.5_dp,runif1()-0.5_dp,runif1()-0.5_dp]
            nrm=sqrt(sum(axis*axis))
            if(nrm>sqrt(tiny(1.0_dp)))exit
         end do
         axis=axis/nrm
         angle=2.0_dp*pi*runif1()
         call rotate_lines(lines,axis,angle)
         do j=1,nl
            select case(model)
            case(TBM_MATERN)
               freq=sqrt(0.5_dp*chisq_rand(3.0_dp)/gamma_draw(smooth))
               theta=1.0_dp
            case(TBM_CAUCHY)
               freq=sqrt(2.0_dp*chisq_rand(3.0_dp)*gamma_draw(smooth))
               theta=1.0_dp
            case(TBM_POWEREXP)
               u1=exp_rand()
               u2=pi*(runif1()-0.5_dp)
               g=abs(sin(0.5_dp*smooth*(u2-0.5_dp*pi))*cos(u2)**(-2.0_dp/smooth)* &
                  (cos(u2-0.5_dp*smooth*(u2-0.5_dp*pi))/u1)**((2.0_dp-smooth)/smooth))
               freq=sqrt(2.0_dp*chisq_rand(3.0_dp)*g)
               theta=1.0_dp
            case(TBM_BESSEL)
               freq=sqrt(beta_draw(1.5_dp,smooth-0.5_dp))
               theta=1.0_dp
            case(TBM_GAUSSIAN)
               freq=sqrt(2.0_dp*chisq_rand(3.0_dp))
               theta=1.0_dp
            case(TBM_FBM)
               r=gamma_draw(1.0_dp-0.5_dp*smooth)/gamma_draw(0.5_dp*smooth)
               theta=sqrt((1.0_dp+r)/r**(0.5_dp*smooth+1.0_dp))
               freq=2.0_dp*pi*r
            end select
            phase=2.0_dp*pi*runif1()
            do k=1,neff
               if(dim==2)then
               eprod=points(k,1)*lines(j,2)+points(k,2)*lines(j,3)
               else
               eprod=dot_product(points(k,1:3),lines(j,1:3))
               end if
               ans(i,k)=ans(i,k)+theta*cos(freq*eprod+phase)
               if(model==TBM_FBM.and.r<1.0e-3_dp)ans(i,k)=ans(i,k)-theta*cos(phase)
            end do
         end do
      end do

      norm_const=sqrt(sill/real(nl,dp))
      if(model/=TBM_FBM)then
         norm_const=norm_const*sqrt(2.0_dp)
      else
         norm_const=norm_const*sqrt(4.0_dp*gamma(0.5_dp*smooth+1.0_dp)* &
            gamma(0.5_dp*(real(dim,dp)+smooth))/(pi**smooth*gamma(0.5_dp*real(dim,dp))))
      end if
      ans=ans*norm_const
      if(nugget>0.0_dp)then
         do i=1,n
         do k=1,neff
         ans(i,k)=ans(i,k)+sqrt(nugget)*rnorm1()
         end do
         end do
      end if
      if(present(info))info=istat
   end function simulate_gaussian_process_tbm


   function simulate_schlather_tbm(n,coord,covmod,nugget,range,smooth,nlines,grid,u_bound,info) result(ans)
      ! Turning-band Schlather simulation corresponding to upstream
      ! rschlathertbm.  u_bound defaults to the R wrapper value 3.5.
      integer,intent(in)::n,covmod
      real(dp),intent(in)::coord(:,:),nugget,range,smooth
      integer,intent(in),optional::nlines
      logical,intent(in),optional::grid
      real(dp),intent(in),optional::u_bound
      integer,intent(out),optional::info
      real(dp),allocatable::ans(:,:),gp(:,:)
      real(dp)::poisson,ipoisson,thresh,ub
      integer::i,neff,nko,istat,nl
      logical::isgrid

      nl=1000
      if(present(nlines))nl=nlines
      isgrid=.false.
      if(present(grid))isgrid=grid
      ub=3.5_dp
      if(present(u_bound))ub=u_bound
      neff=effective_site_count(coord,isgrid)
      allocate(ans(max(0,n),max(0,neff)))
      ans=0.0_dp
      istat=0
      if(n<0.or.neff<=0.or.ub<=0.0_dp)then
         istat=1
         if(present(info))info=istat
         return
      end if
      do i=1,n
         poisson=0.0_dp
         nko=neff
         do while(nko>0)
            poisson=poisson+exp_rand()
            ipoisson=1.0_dp/poisson
            thresh=ub*ipoisson
            gp=simulate_gaussian_process_tbm(1,coord,covmod,nugget,1.0_dp-nugget,range,smooth,nl,isgrid,istat)
            if(istat/=0)then
            if(present(info))info=istat
            return
            end if
            ans(i,:)=max(ans(i,:),gp(1,:)*ipoisson)
            nko=count(ans(i,:)<thresh)
         end do
      end do
      ans=ans*sqrt(2.0_dp*pi)
      if(present(info))info=istat
   end function simulate_schlather_tbm

   function simulate_geomgauss_tbm(n,coord,covmod,sigma2,nugget,range,smooth,nlines,grid,u_bound,info) result(ans)
      ! Turning-band geometric-Gaussian simulation corresponding to
      ! upstream rgeomtbm; result has unit-Frechet margins.
      integer,intent(in)::n,covmod
      real(dp),intent(in)::coord(:,:),sigma2,nugget,range,smooth
      integer,intent(in),optional::nlines
      logical,intent(in),optional::grid
      real(dp),intent(in),optional::u_bound
      integer,intent(out),optional::info
      real(dp),allocatable::ans(:,:),gp(:,:)
      real(dp)::poisson,ipoisson,thresh,ub,half_sigma2,sigma
      integer::i,neff,nko,istat,nl
      logical::isgrid

      nl=1000
      if(present(nlines))nl=nlines
      isgrid=.false.
      if(present(grid))isgrid=grid
      ub=exp(3.5_dp*sqrt(max(0.0_dp,sigma2))-0.5_dp*sigma2)
      if(present(u_bound))ub=u_bound
      neff=effective_site_count(coord,isgrid)
      allocate(ans(max(0,n),max(0,neff)))
      ans=0.0_dp
      istat=0
      if(n<0.or.neff<=0.or.sigma2<=0.0_dp.or.ub<=0.0_dp)then
         istat=1
         if(present(info))info=istat
         return
      end if
      sigma=sqrt(sigma2)
      half_sigma2=0.5_dp*sigma2
      do i=1,n
         poisson=0.0_dp
         nko=neff
         do while(nko>0)
            poisson=poisson+exp_rand()
            ipoisson=-log(poisson)
            thresh=log(ub)+ipoisson
            gp=simulate_gaussian_process_tbm(1,coord,covmod,nugget,1.0_dp-nugget,range,smooth,nl,isgrid,istat)
            if(istat/=0)then
            if(present(info))info=istat
            return
            end if
            ans(i,:)=max(ans(i,:),sigma*gp(1,:)+ipoisson-half_sigma2)
            nko=count(ans(i,:)<thresh)
         end do
      end do
      ans=exp(ans)
      if(present(info))info=istat
   end function simulate_geomgauss_tbm

   function simulate_extremalt_tbm(n,coord,covmod,nugget,range,smooth,nu,nlines,grid,u_bound,info) result(ans)
      ! Turning-band extremal-t simulation corresponding to upstream
      ! rextremalttbm.  u_bound defaults to 3**nu as in the R wrapper.
      integer,intent(in)::n,covmod
      real(dp),intent(in)::coord(:,:),nugget,range,smooth,nu
      integer,intent(in),optional::nlines
      logical,intent(in),optional::grid
      real(dp),intent(in),optional::u_bound
      integer,intent(out),optional::info
      real(dp),allocatable::ans(:,:),gp(:,:)
      real(dp)::poisson,ipoisson,thresh,ub,cnu
      integer::i,neff,nko,istat,nl
      logical::isgrid

      nl=1000
      if(present(nlines))nl=nlines
      isgrid=.false.
      if(present(grid))isgrid=grid
      ub=3.0_dp**nu
      if(present(u_bound))ub=u_bound
      neff=effective_site_count(coord,isgrid)
      allocate(ans(max(0,n),max(0,neff)))
      ans=0.0_dp
      istat=0
      if(n<0.or.neff<=0.or.nu<=0.0_dp.or.ub<=0.0_dp)then
         istat=1
         if(present(info))info=istat
         return
      end if
      cnu=sqrt(pi)*2.0_dp**(1.0_dp-0.5_dp*nu)/gamma(0.5_dp*(nu+1.0_dp))
      do i=1,n
         poisson=0.0_dp
         nko=neff
         do while(nko>0)
            poisson=poisson+exp_rand()
            ipoisson=1.0_dp/poisson
            thresh=ub*ipoisson
            gp=simulate_gaussian_process_tbm(1,coord,covmod,nugget,1.0_dp-nugget,range,smooth,nl,isgrid,istat)
            if(istat/=0)then
            if(present(info))info=istat
            return
            end if
            ans(i,:)=max(ans(i,:),max(gp(1,:),0.0_dp)**nu*ipoisson)
            nko=count(ans(i,:)<thresh)
         end do
      end do
      ans=ans*cnu
      if(present(info))info=istat
   end function simulate_extremalt_tbm

   integer function effective_site_count(coord,isgrid) result(neff)
      real(dp),intent(in)::coord(:,:)
      logical,intent(in)::isgrid
      integer::dim,m
      dim=size(coord,2)
      m=size(coord,1)
      if(.not.isgrid)then
         neff=m
      else if(dim==2)then
         neff=m*m
      else if(dim==3)then
         neff=m*m*m
      else
         neff=0
      end if
   end function effective_site_count

   subroutine make_grid_points(coord,points)
      real(dp),intent(in)::coord(:,:)
      real(dp),allocatable,intent(out)::points(:,:)
      integer::m,i,j,k,c
      m=size(coord,1)
      if(size(coord,2)==2)then
         allocate(points(m*m,2))
         c=0
         do j=1,m
         do i=1,m
         c=c+1
         points(c,:)=[coord(i,1),coord(j,2)]
         end do
         end do
      else
         allocate(points(m*m*m,3))
         c=0
         do k=1,m
         do j=1,m
         do i=1,m
            c=c+1
            points(c,:)=[coord(i,1),coord(j,2),coord(k,3)]
         end do
         end do
         end do
      end if
   end subroutine make_grid_points

   real(dp) function gamma_draw(shape) result(x)
      real(dp),intent(in)::shape
      real(dp),allocatable::tmp(:)
      tmp=rgamma(1,shape,1.0_dp)
      x=tmp(1)
   end function gamma_draw

   real(dp) function beta_draw(a,b) result(x)
      real(dp),intent(in)::a,b
      real(dp),allocatable::tmp(:)
      tmp=rbeta(1,a,b)
      x=tmp(1)
   end function beta_draw
end module spatialextremes_turning_bands
