module spatialextremes_circulant
   use spatialextremes_base, only: dp,pi,exp_rand,nan_dp
   use spatialextremes_covariance, only: covariance_values
   use r_compat, only: rnorm1
   implicit none
   private
   type :: circ_embedding_t
      integer :: m=0,ngrid=0
      real(dp) :: nugget=0.0_dp
      complex(dp),allocatable :: sqrtlam(:,:)
   end type circ_embedding_t
   public :: simulate_gaussian_grid_circulant, simulate_schlather_circulant
   public :: simulate_geomgauss_circulant, simulate_extremalt_circulant
contains
   function simulate_gaussian_grid_circulant(n,ngrid,steps,covmod,nugget,sill,range,smooth,info) result(ans)
      integer,intent(in)::n,ngrid,covmod
      real(dp),intent(in)::steps(2),nugget,sill,range,smooth
      integer,intent(out),optional::info
      real(dp)::ans(n,ngrid*ngrid)
      type(circ_embedding_t)::emb
      integer::k,istat
      call prepare_embedding(ngrid,steps,covmod,nugget,sill,range,smooth,emb,istat)
      if(istat/=0)then
      ans=nan_dp()
      if(present(info))info=istat
      return
      end if
      do k=1,n
      call draw_embedding(emb,ans(k,:))
      end do
      if(present(info))info=0
   end function simulate_gaussian_grid_circulant

   function simulate_schlather_circulant(n,ngrid,steps,covmod,nugget,range,smooth,ubound,info) result(ans)
      integer,intent(in)::n,ngrid,covmod
      real(dp),intent(in)::steps(2),nugget,range,smooth
      real(dp),intent(in),optional::ubound
      integer,intent(out),optional::info
      real(dp)::ans(n,ngrid*ngrid),gp(1,ngrid*ngrid),poisson,ip,thr,ub
      integer::i,istat
      type(circ_embedding_t)::emb
      ub=sqrt(2.0_dp*pi)
      if(present(ubound))ub=ubound
      call prepare_embedding(ngrid,steps,covmod,nugget,1.0_dp-nugget,range,smooth,emb,istat)
      if(istat/=0)then
      ans=nan_dp()
      if(present(info))info=istat
      return
      end if
      ans=0.0_dp
      do i=1,n
         poisson=0.0_dp
         do
            poisson=poisson+exp_rand()
            ip=1.0_dp/poisson
            thr=ub*ip
            call draw_embedding(emb,gp(1,:))
            ans(i,:)=max(ans(i,:),gp(1,:)*ip)
            if(all(ans(i,:)>=thr))exit
         end do
      end do
      ans=ans*sqrt(2.0_dp*pi)
      if(present(info))info=0
   end function simulate_schlather_circulant

   function simulate_geomgauss_circulant(n,ngrid,steps,covmod,nugget,range,smooth,sigma2,ubound,info) result(ans)
      integer,intent(in)::n,ngrid,covmod
      real(dp),intent(in)::steps(2),nugget,range,smooth,sigma2
      real(dp),intent(in),optional::ubound
      integer,intent(out),optional::info
      real(dp)::ans(n,ngrid*ngrid),gp(1,ngrid*ngrid),poisson,lp,thr,ub
      integer::i,istat
      type(circ_embedding_t)::emb
      ub=exp(0.5_dp*sigma2+3.0_dp*sqrt(sigma2))
      if(present(ubound))ub=ubound
      call prepare_embedding(ngrid,steps,covmod,nugget,1.0_dp-nugget,range,smooth,emb,istat)
      if(istat/=0)then
      ans=nan_dp()
      if(present(info))info=istat
      return
      end if
      ans=-huge(1.0_dp)
      do i=1,n
         poisson=0.0_dp
         do
            poisson=poisson+exp_rand()
            lp=-log(poisson)
            thr=log(ub)+lp
            call draw_embedding(emb,gp(1,:))
            ans(i,:)=max(ans(i,:),sqrt(sigma2)*gp(1,:)+lp-0.5_dp*sigma2)
            if(all(ans(i,:)>=thr))exit
         end do
      end do
      ans=exp(ans)
      if(present(info))info=0
   end function simulate_geomgauss_circulant

   function simulate_extremalt_circulant(n,ngrid,steps,covmod,nugget,range,smooth,nu,ubound,info) result(ans)
      integer,intent(in)::n,ngrid,covmod
      real(dp),intent(in)::steps(2),nugget,range,smooth,nu
      real(dp),intent(in),optional::ubound
      integer,intent(out),optional::info
      real(dp)::ans(n,ngrid*ngrid),gp(1,ngrid*ngrid),poisson,ip,thr,ub,cnu
      integer::i,istat
      type(circ_embedding_t)::emb
      ub=5.0_dp**nu
      if(present(ubound))ub=ubound
      call prepare_embedding(ngrid,steps,covmod,nugget,1.0_dp-nugget,range,smooth,emb,istat)
      if(istat/=0)then
      ans=nan_dp()
      if(present(info))info=istat
      return
      end if
      ans=0.0_dp
      do i=1,n
         poisson=0.0_dp
         do
            poisson=poisson+exp_rand()
            ip=1.0_dp/poisson
            thr=ub*ip
            call draw_embedding(emb,gp(1,:))
            ans(i,:)=max(ans(i,:),max(gp(1,:),0.0_dp)**nu*ip)
            if(all(ans(i,:)>=thr))exit
         end do
      end do
      cnu=sqrt(pi)*2.0_dp**(1.0_dp-0.5_dp*nu)/gamma(0.5_dp*(nu+1.0_dp))
      ans=ans*cnu
      if(present(info))info=0
   end function simulate_extremalt_circulant

   subroutine prepare_embedding(ngrid,steps,covmod,nugget,sill,range,smooth,emb,info)
      integer,intent(in)::ngrid,covmod
      real(dp),intent(in)::steps(2),nugget,sill,range,smooth
      type(circ_embedding_t),intent(out)::emb
      integer,intent(out)::info
      complex(dp),allocatable::lam(:,:)
      real(dp)::tol,h
      integer::m,i,j,ii,jj
      m=1
      do while(m<2*max(1,ngrid-1))
      m=2*m
      end do
      info=1
      do while(m<=4096 .and. info/=0)
         allocate(lam(m,m))
         do j=1,m
            jj=j-1
            if(jj>m/2)jj=jj-m
            do i=1,m
               ii=i-1
               if(ii>m/2)ii=ii-m
               h=hypot(steps(1)*real(ii,dp),steps(2)*real(jj,dp))
               lam(i,j)=cmplx(covariance_values(h,covmod,0.0_dp,sill,range,smooth,dim=2),0.0_dp,dp)
            end do
         end do
         call fft2(lam,.false.)
         tol=1.0e-10_dp*max(1.0_dp,maxval(abs(real(lam,dp))))
         if(maxval(abs(aimag(lam)))<=100.0_dp*tol .and. minval(real(lam,dp))>=-tol)then
            info=0
            emb%m=m
            emb%ngrid=ngrid
            emb%nugget=nugget
            allocate(emb%sqrtlam(m,m))
            emb%sqrtlam=cmplx(sqrt(max(real(lam,dp),0.0_dp)),0.0_dp,dp)
         else
            m=2*m
         end if
         deallocate(lam)
      end do
   end subroutine prepare_embedding

   subroutine draw_embedding(emb,x)
      type(circ_embedding_t),intent(in)::emb
      real(dp),intent(out)::x(:)
      complex(dp),allocatable::work(:,:)
      integer::i,j,idx
      allocate(work(emb%m,emb%m))
      do j=1,emb%m
      do i=1,emb%m
      work(i,j)=cmplx(rnorm1(),0.0_dp,dp)
      end do
      end do
      call fft2(work,.false.)
      work=work*emb%sqrtlam
      call fft2(work,.true.)
      idx=0
      do j=1,emb%ngrid
         do i=1,emb%ngrid
            idx=idx+1
            x(idx)=real(work(i,j),dp)+sqrt(max(0.0_dp,emb%nugget))*rnorm1()
         end do
      end do
   end subroutine draw_embedding

   subroutine fft2(a,inverse)
      complex(dp),intent(inout)::a(:,:)
      logical,intent(in)::inverse
      complex(dp),allocatable::v(:)
      integer::i,j,n1,n2
      n1=size(a,1)
      n2=size(a,2)
      allocate(v(max(n1,n2)))
      do j=1,n2
         v(1:n1)=a(:,j)
         call fft1(v(1:n1),inverse)
         a(:,j)=v(1:n1)
      end do
      do i=1,n1
         v(1:n2)=a(i,:)
         call fft1(v(1:n2),inverse)
         a(i,:)=v(1:n2)
      end do
   end subroutine fft2

   subroutine fft1(a,inverse)
      complex(dp),intent(inout)::a(:)
      logical,intent(in)::inverse
      integer::n,i,j,k,m,half
      complex(dp)::tmp,w,wm,u,t
      real(dp)::sgn
      n=size(a)
      j=1
      do i=1,n
         if(i<j)then
         tmp=a(i)
         a(i)=a(j)
         a(j)=tmp
         end if
         k=n/2
         do while(k>=1 .and. j>k)
         j=j-k
         k=k/2
         end do
         j=j+k
      end do
      m=2
      sgn=-1.0_dp
      if(inverse)sgn=1.0_dp
      do while(m<=n)
         half=m/2
         wm=exp(cmplx(0.0_dp,sgn*2.0_dp*pi/real(m,dp),dp))
         do k=1,n,m
            w=cmplx(1.0_dp,0.0_dp,dp)
            do j=0,half-1
               u=a(k+j)
               t=w*a(k+j+half)
               a(k+j)=u+t
               a(k+j+half)=u-t
               w=w*wm
            end do
         end do
         m=2*m
      end do
      if(inverse)a=a/real(n,dp)
   end subroutine fft1
end module spatialextremes_circulant
