module spatialextremes_simulation
   use spatialextremes_base, only: dp,pi,chol_upper,solve_spd,exp_rand,chisq_rand,nan_dp
   use spatialextremes_covariance, only: fbm_covariance,covariance_matrix,covariance_values
   use r_compat, only: rnorm1
   implicit none
   private
   public :: simulate_gaussian_process,simulate_schlather_spectral,simulate_geomgauss_spectral
   public :: simulate_extremalt_spectral,simulate_brownresnick_spectral
   public :: simulate_brownresnick_exact,simulate_brownresnick_exact_hitting
   public :: simulate_schlather_exact,simulate_extremalt_exact
   public :: conditional_gaussian_process
contains
   function simulate_gaussian_process(n,cov,mean) result(x)
      integer,intent(in)::n
      real(dp),intent(in)::cov(:,:)
      real(dp),intent(in),optional::mean(:)
      real(dp)::x(n,size(cov,1)),r(size(cov,1),size(cov,2)),z(size(cov,1))
      integer::info,i,j
      call chol_upper(cov,r,info)
      if(info/=0)then
      x=nan_dp()
      return
      end if
      do i=1,n
         do j=1,size(z)
         z(j)=rnorm1()
         end do
         x(i,:)=matmul(transpose(r),z)
         if(present(mean))x(i,:)=x(i,:)+mean
      end do
   end function simulate_gaussian_process

   function conditional_gaussian_process(n,coord,data_coord,data,covmod,sill,range,smooth,mean,smooth2) result(x)
      ! Direct conditional Gaussian simulation.  This is distributionally
      ! equivalent to upstream condrgp's unconditional-draw plus kriging-
      ! residual correction, but works directly with the conditional mean
      ! and Schur-complement covariance.
      integer,intent(in)::n,covmod
      real(dp),intent(in)::coord(:,:),data_coord(:,:),data(:),sill,range,smooth
      real(dp),intent(in),optional::mean,smooth2
      real(dp)::x(n,size(coord,1))
      real(dp)::coo(size(data_coord,1),size(data_coord,1))
      real(dp)::bot(size(data_coord,1),size(coord,1)),w(size(data_coord,1),size(coord,1))
      real(dp)::ctt(size(coord,1),size(coord,1)),cc(size(coord,1),size(coord,1))
      real(dp)::cm(size(coord,1)),res(size(data_coord,1)),r(size(coord,1),size(coord,1))
      real(dp)::z(size(coord,1)),h,s2,mu
      integer::i,j,k,info
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      mu=0.0_dp
      if(present(mean))mu=mean
      if(size(data)/=size(data_coord,1) .or. size(coord,2)/=size(data_coord,2)) then
         x=nan_dp()
         return
      end if
      do i=1,size(data_coord,1)
         do j=1,size(data_coord,1)
            h=sqrt(sum((data_coord(i,:)-data_coord(j,:))**2))
            coo(i,j)=covariance_values(h,covmod,0.0_dp,sill,range,smooth,s2,size(coord,2))
         end do
         do j=1,size(coord,1)
            h=sqrt(sum((data_coord(i,:)-coord(j,:))**2))
            bot(i,j)=covariance_values(h,covmod,0.0_dp,sill,range,smooth,s2,size(coord,2))
         end do
      end do
      do i=1,size(coord,1)
         do j=1,size(coord,1)
            h=sqrt(sum((coord(i,:)-coord(j,:))**2))
            ctt(i,j)=covariance_values(h,covmod,0.0_dp,sill,range,smooth,s2,size(coord,2))
         end do
      end do
      call solve_spd(coo,bot,w,info)
      if(info/=0)then
      x=nan_dp()
      return
      end if
      res=data-mu
      cm=mu+matmul(transpose(w),res)
      cc=ctt-matmul(transpose(bot),w)
      cc=0.5_dp*(cc+transpose(cc))
      do i=1,size(cc,1)
         if(cc(i,i)<0.0_dp .and. abs(cc(i,i))<1.0e-12_dp*max(1.0_dp,sill))cc(i,i)=0.0_dp
      end do
      call chol_upper(cc+1.0e-14_dp*max(1.0_dp,sill)*identity_local(size(cc,1)),r,info)
      if(info/=0)then
      x=nan_dp()
      return
      end if
      do i=1,n
         do k=1,size(z)
         z(k)=rnorm1()
         end do
         x(i,:)=cm+matmul(transpose(r),z)
      end do
   end function conditional_gaussian_process

   function simulate_schlather_spectral(n,cor,n_storms) result(x)
      integer,intent(in)::n,n_storms
      real(dp),intent(in)::cor(:,:)
      real(dp)::x(n,size(cor,1)),r(size(cor,1),size(cor,2)),z(size(cor,1)),y(size(cor,1)),g,xi
      integer::info,i,j,s
      call chol_upper(cor,r,info)
      if(info/=0)then
      x=nan_dp()
      return
      end if
      x=0.0_dp
      do i=1,n
         g=0.0_dp
         do s=1,n_storms
            g=g+exp_rand()
            xi=1.0_dp/g
            do j=1,size(z)
            z(j)=rnorm1()
            end do
            z=matmul(transpose(r),z)
            y=sqrt(2.0_dp*pi)*max(z,0.0_dp)
            x(i,:)=max(x(i,:),xi*y)
         end do
      end do
   end function simulate_schlather_spectral

   function simulate_geomgauss_spectral(n,cor,sigma2,n_storms) result(x)
      integer,intent(in)::n,n_storms
      real(dp),intent(in)::cor(:,:),sigma2
      real(dp)::x(n,size(cor,1)),r(size(cor,1),size(cor,2)),z(size(cor,1)),y(size(cor,1)),g,xi
      integer::info,i,j,s
      call chol_upper(cor,r,info)
      if(info/=0)then
      x=nan_dp()
      return
      end if
      x=0.0_dp
      do i=1,n
         g=0.0_dp
         do s=1,n_storms
            g=g+exp_rand()
            xi=1.0_dp/g
            do j=1,size(z)
            z(j)=rnorm1()
            end do
            z=sqrt(sigma2)*matmul(transpose(r),z)
            y=exp(z-0.5_dp*sigma2)
            x(i,:)=max(x(i,:),xi*y)
         end do
      end do
   end function simulate_geomgauss_spectral

   function simulate_extremalt_spectral(n,cor,nu,n_storms) result(x)
      integer,intent(in)::n,n_storms
      real(dp),intent(in)::cor(:,:),nu
      real(dp)::x(n,size(cor,1)),r(size(cor,1),size(cor,2)),z(size(cor,1)),y(size(cor,1)),g,xi,cnu
      integer::info,i,j,s
      call chol_upper(cor,r,info)
      if(info/=0)then
      x=nan_dp()
      return
      end if
      cnu=sqrt(pi)*2.0_dp**(1.0_dp-0.5_dp*nu)/gamma(0.5_dp*(nu+1.0_dp))
      x=0.0_dp
      do i=1,n
         g=0.0_dp
         do s=1,n_storms
            g=g+exp_rand()
            xi=1.0_dp/g
            do j=1,size(z)
            z(j)=rnorm1()
            end do
            z=matmul(transpose(r),z)
            y=cnu*max(z,0.0_dp)**nu
            x(i,:)=max(x(i,:),xi*y)
         end do
      end do
   end function simulate_extremalt_spectral

   function simulate_brownresnick_spectral(n,coord,range,smooth,n_storms) result(x)
      integer,intent(in)::n,n_storms
      real(dp),intent(in)::coord(:,:),range,smooth
      real(dp)::x(n,size(coord,1)),cov(size(coord,1),size(coord,1)),r(size(coord,1),size(coord,2))
      real(dp)::z(size(coord,1)),y(size(coord,1)),v(size(coord,1)),g,xi
      integer::info,i,j,s
      cov=fbm_covariance(coord,range,smooth,1.0_dp)
      v=[(cov(j,j),j=1,size(coord,1))]
      call chol_upper(cov+1.0e-12_dp*identity(size(cov,1)),r,info)
      if(info/=0)then
      x=nan_dp()
      return
      end if
      x=0.0_dp
      do i=1,n
         g=0.0_dp
         do s=1,n_storms
            g=g+exp_rand()
            xi=1.0_dp/g
            do j=1,size(z)
            z(j)=rnorm1()
            end do
            z=matmul(transpose(r),z)
            y=exp(z-0.5_dp*v)
            x(i,:)=max(x(i,:),xi*y)
         end do
      end do
   contains
      pure function identity(m) result(a)
         integer,intent(in)::m
         real(dp)::a(m,m)
         integer::k
         a=0
         do k=1,m
         a(k,k)=1
         end do
      end function
   end function simulate_brownresnick_spectral

   function simulate_brownresnick_exact(n,coord,range,smooth) result(ans)
      integer,intent(in)::n
      real(dp),intent(in)::coord(:,:),range,smooth
      real(dp)::ans(n,size(coord,1)),cov(size(coord,1),size(coord,1)),r(size(coord,1),size(coord,1))
      real(dp)::gp(size(coord,1)),vario(size(coord,1)),shift(size(coord,1),size(coord,2)),g,ipo,dummy
      integer::info,i,j,l,d
      cov=fbm_covariance(coord,range,smooth,1.0_dp)
      call chol_upper(cov+1.0e-12_dp*identity_local(size(cov,1)),r,info)
      if(info/=0)then
      ans=nan_dp()
      return
      end if
      ans=-1.0e10_dp
      do j=1,size(coord,1)
         do l=1,size(coord,1)
         shift(l,:)=coord(l,:)-coord(j,:)
         end do
         do l=1,size(coord,1)
         vario(l)=(sqrt(sum(shift(l,:)**2))/range)**smooth
         end do
         do i=1,n
            g=exp_rand()
            ipo=-log(g)
            do while(ans(i,j)<ipo)
               do l=1,size(gp)
               gp(l)=rnorm1()
               end do
               gp=matmul(transpose(r),gp)
               dummy=gp(j)
               gp=gp-dummy-vario
               if(j==1 .or. all(ipo+gp(1:j-1)<=ans(i,1:j-1))) then
                  ans(i,j:)=max(ans(i,j:),ipo+gp(j:))
               end if
               g=g+exp_rand()
               ipo=-log(g)
            end do
         end do
      end do
      ans=exp(ans)
   end function simulate_brownresnick_exact


   subroutine simulate_brownresnick_exact_hitting(nobs,coord,range,smooth,ans,hitting,info)
      ! Exact Brown-Resnick simulation with the associated hitting scenario.
      ! This is the computational content of upstream rhitscenbrown/newCode.c.
      integer,intent(in)::nobs
      real(dp),intent(in)::coord(:,:),range,smooth
      real(dp),allocatable,intent(out)::ans(:,:)
      integer,allocatable,intent(out)::hitting(:,:)
      integer,intent(out)::info
      real(dp),allocatable::cov(:,:),r(:,:),gp(:),vario(:),shift(:,:)
      real(dp)::g,ipo,dummy
      integer::i,j,l,storm
      if(nobs<1 .or. range<=0.0_dp .or. smooth<=0.0_dp .or. smooth>2.0_dp)then
         allocate(ans(0,0),hitting(0,0))
         info=1
         return
      end if
      allocate(ans(nobs,size(coord,1)),hitting(nobs,size(coord,1)))
      allocate(cov(size(coord,1),size(coord,1)),r(size(coord,1),size(coord,1)))
      allocate(gp(size(coord,1)),vario(size(coord,1)),shift(size(coord,1),size(coord,2)))
      cov=fbm_covariance(coord,range,smooth,1.0_dp)
      call chol_upper(cov+1.0e-12_dp*identity_local(size(coord,1)),r,info)
      if(info/=0)return
      ans=-1.0e10_dp
      hitting=0
      do i=1,nobs
         storm=0
         do j=1,size(coord,1)
            do l=1,size(coord,1)
            shift(l,:)=coord(l,:)-coord(j,:)
            end do
            do l=1,size(coord,1)
            vario(l)=(sqrt(sum(shift(l,:)**2))/range)**smooth
            end do
            g=exp_rand()
            ipo=-log(g)
            do while(ans(i,j)<ipo)
               do l=1,size(gp)
               gp(l)=rnorm1()
               end do
               gp=matmul(transpose(r),gp)
               dummy=gp(j)
               gp=gp-dummy-vario
               if(j==1 .or. all(ipo+gp(1:j-1)<=ans(i,1:j-1)))then
                  storm=storm+1
                  do l=j,size(coord,1)
                     if(ipo+gp(l)>ans(i,l))then
                        ans(i,l)=ipo+gp(l)
                        hitting(i,l)=storm
                     end if
                  end do
               end if
               g=g+exp_rand()
               ipo=-log(g)
            end do
         end do
      end do
      ans=exp(ans)
      info=0
   end subroutine simulate_brownresnick_exact_hitting

   function simulate_schlather_exact(n,coord,covmod,nugget,range,smooth) result(ans)
      integer,intent(in)::n,covmod
      real(dp),intent(in)::coord(:,:),nugget,range,smooth
      real(dp)::ans(n,size(coord,1)),cov(size(coord,1),size(coord,1)),smat(size(coord,1),size(coord,1))
      real(dp)::r(size(coord,1),size(coord,1)),gp(size(coord,1)),mu(size(coord,1)),g,scale
      integer::info,i,j,l,m
      cov=covariance_matrix(coord,covmod,nugget,1.0_dp-nugget,range,smooth)
      ans=0.0_dp
      do j=1,size(coord,1)
         mu=cov(:,j)
         do l=1,size(coord,1)
            do m=1,size(coord,1)
            smat(l,m)=(cov(l,m)-cov(j,l)*cov(j,m))/2.0_dp
            end do
         end do
         smat(j,j)=1.0e-12_dp
         call chol_upper(smat,r,info)
         if(info/=0)then
         ans=nan_dp()
         return
         end if
         r(j,j)=0.0_dp
         do i=1,n
            g=exp_rand()
            do while(g*ans(i,j)<1.0_dp)
               do l=1,size(gp)
               gp(l)=rnorm1()
               end do
               gp=matmul(transpose(r),gp)
               scale=sqrt(2.0_dp/chisq_rand(2.0_dp))
               gp=mu+gp*scale
               if(j==1 .or. all(gp(1:j-1)<=g*ans(i,1:j-1)))ans(i,j:)=max(ans(i,j:),gp(j:)/g)
               g=g+exp_rand()
            end do
         end do
      end do
   end function simulate_schlather_exact

   function simulate_extremalt_exact(n,coord,covmod,nugget,range,smooth,nu) result(ans)
      integer,intent(in)::n,covmod
      real(dp),intent(in)::coord(:,:),nugget,range,smooth,nu
      real(dp)::ans(n,size(coord,1)),cov(size(coord,1),size(coord,1)),smat(size(coord,1),size(coord,1))
      real(dp)::r(size(coord,1),size(coord,1)),gp(size(coord,1)),mu(size(coord,1)),g,point,scale
      integer::info,i,j,l,m
      cov=covariance_matrix(coord,covmod,nugget,1.0_dp-nugget,range,smooth)
      ans=0.0_dp
      do j=1,size(coord,1)
         mu=cov(:,j)
         do l=1,size(coord,1)
            do m=1,size(coord,1)
            smat(l,m)=(cov(l,m)-cov(j,l)*cov(j,m))/(1.0_dp+nu)
            end do
         end do
         smat(j,j)=1.0e-12_dp
         call chol_upper(smat,r,info)
         if(info/=0)then
         ans=nan_dp()
         return
         end if
         r(j,j)=0.0_dp
         do i=1,n
            g=exp_rand()
            point=g**(-1.0_dp/nu)
            do while(ans(i,j)<point)
               do l=1,size(gp)
               gp(l)=rnorm1()
               end do
               gp=matmul(transpose(r),gp)
               scale=sqrt((1.0_dp+nu)/chisq_rand(1.0_dp+nu))
               gp=mu+gp*scale
               if(j==1 .or. all(point*gp(1:j-1)<=ans(i,1:j-1)))ans(i,j:)=max(ans(i,j:),point*gp(j:))
               g=g+exp_rand()
               point=g**(-1.0_dp/nu)
            end do
         end do
      end do
      ans=max(ans,0.0_dp)**nu
   end function simulate_extremalt_exact

   pure function identity_local(m) result(a)
      integer,intent(in)::m
      real(dp)::a(m,m)
      integer::k
      a=0.0_dp
      do k=1,m
      a(k,k)=1.0_dp
      end do
   end function identity_local
end module spatialextremes_simulation
