module spatialextremes_fit
   use spatialextremes_base, only: dp,distance_vectors,euclidean_distances,eye
   use spatialextremes_covariance
   use spatialextremes_pairwise
   implicit none
   private
   public :: maxstab_fit_t, univ_fit_t, gevmle, gpdmle
   public :: fit_smith_frechet, fit_schlather_frechet, fit_brownresnick_frechet
   public :: fit_extremalt_frechet
   public :: nelder_mead

   type :: univ_fit_t
      real(dp) :: loc=0.0_dp, scale=1.0_dp, shape=0.0_dp
      real(dp) :: loglik=-huge(1.0_dp)
      integer :: iterations=0, convergence=1
   end type

   type :: maxstab_fit_t
      real(dp), allocatable :: par(:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: convergence = 1
   end type

   abstract interface
      function scalar_objective(x) result(v)
         import dp
         real(dp),intent(in)::x(:)
         real(dp)::v
      end function
   end interface
contains
   subroutine nelder_mead(fn,x0,xbest,fbest,iterations,convergence,maxit,tol,step)
      procedure(scalar_objective)::fn
      real(dp),intent(in)::x0(:)
      real(dp),intent(out)::xbest(size(x0)),fbest
      integer,intent(out)::iterations,convergence
      integer,intent(in),optional::maxit
      real(dp),intent(in),optional::tol,step
      integer::n,mx,i,j,ilo,ihi,isecond,it
      real(dp)::eps,st,fr,fe,fc
      real(dp),allocatable::x(:,:),f(:),cent(:),xr(:),xe(:),xc(:)
      n=size(x0)
      mx=1000
      if(present(maxit))mx=maxit
      eps=1.0e-8_dp
      if(present(tol))eps=tol
      st=0.1_dp
      if(present(step))st=step
      allocate(x(n,n+1),f(n+1),cent(n),xr(n),xe(n),xc(n))
      x(:,1)=x0
      do j=1,n
         x(:,j+1)=x0
         x(j,j+1)=x(j,j+1)+st*max(1.0_dp,abs(x0(j)))
      end do
      do j=1,n+1
      f(j)=fn(x(:,j))
      end do
      convergence=1
      do it=1,mx
         ilo=minloc(f,1)
         ihi=maxloc(f,1)
         isecond=ilo
         do j=1,n+1
            if(j/=ihi .and. (isecond==ilo .or. f(j)>f(isecond)))isecond=j
         end do
         if(maxval(abs(f-f(ilo)))<=eps*(1.0_dp+abs(f(ilo))))then
         convergence=0
         exit
         end if
         cent=0
         do j=1,n+1
         if(j/=ihi)cent=cent+x(:,j)
         end do
         cent=cent/n
         xr=cent+(cent-x(:,ihi))
         fr=fn(xr)
         if(fr<f(ilo))then
            xe=cent+2.0_dp*(xr-cent)
            fe=fn(xe)
            if(fe<fr)then
            x(:,ihi)=xe
            f(ihi)=fe
            else
            x(:,ihi)=xr
            f(ihi)=fr
            end if
         else if(fr<f(isecond))then
            x(:,ihi)=xr
            f(ihi)=fr
         else
            if(fr<f(ihi))then
            xc=cent+0.5_dp*(xr-cent)
            else
            xc=cent+0.5_dp*(x(:,ihi)-cent)
            end if
            fc=fn(xc)
            if(fc<min(fr,f(ihi)))then
               x(:,ihi)=xc
               f(ihi)=fc
            else
               do j=1,n+1
                  if(j/=ilo)then
                  x(:,j)=x(:,ilo)+0.5_dp*(x(:,j)-x(:,ilo))
                  f(j)=fn(x(:,j))
                  end if
               end do
            end if
         end if
      end do
      iterations=min(it,mx)
      ilo=minloc(f,1)
      xbest=x(:,ilo)
      fbest=f(ilo)
   end subroutine nelder_mead


   function gevmle(data,start) result(fit)
      use spatialextremes_univariate, only: gev_loglik
      real(dp),intent(in)::data(:)
      real(dp),intent(in),optional::start(:)
      type(univ_fit_t)::fit
      real(dp)::x0(3),xb(3),fb,m,s
      m=sum(data)/real(size(data),dp)
      s=sqrt(sum((data-m)**2)/max(1.0_dp,real(size(data)-1,dp)))
      x0=[m,log(max(s,1.0e-3_dp)),0.0_dp]
      if(present(start))x0=[start(1),log(start(2)),start(3)]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1200,tol=1.0e-9_dp)
      fit%loc=xb(1)
      fit%scale=exp(xb(2))
      fit%shape=xb(3)
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v
         if(x(3)<-0.999_dp .or. abs(x(3))>5.0_dp)then
         v=huge(1.0_dp)/100.0_dp
         else
         v=-gev_loglik(data,x(1),exp(x(2)),x(3))
         end if
      end function
   end function gevmle

   function gpdmle(data,threshold,start) result(fit)
      use spatialextremes_univariate, only: gpd_loglik
      real(dp),intent(in)::data(:),threshold
      real(dp),intent(in),optional::start(:)
      type(univ_fit_t)::fit
      real(dp)::x0(2),xb(2),fb,s
      real(dp),allocatable::exceed(:)
      integer::nexc
      nexc=count(data>threshold)
      if(nexc<=0) then
         fit%loc=threshold
         fit%scale=0.0_dp
         fit%shape=0.0_dp
         fit%loglik=-huge(1.0_dp)
         fit%convergence=2
         return
      end if
      allocate(exceed(nexc))
      exceed=pack(data,data>threshold)
      s=sum(exceed-threshold)/real(nexc,dp)
      x0=[log(max(s,1.0e-3_dp)),0.0_dp]
      if(present(start))x0=[log(start(1)),start(2)]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1000,tol=1.0e-9_dp)
      fit%loc=threshold
      fit%scale=exp(xb(1))
      fit%shape=xb(2)
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v
         if(x(2)<-0.999_dp .or. abs(x(2))>5.0_dp)then
         v=huge(1.0_dp)/100.0_dp
         else
         v=-gpd_loglik(exceed,threshold,exp(x(1)),x(2))
         end if
      end function
   end function gpdmle

   function fit_brownresnick_frechet(data,coord,start) result(fit)
      real(dp),intent(in)::data(:,:),coord(:,:)
      real(dp),intent(in),optional::start(:)
      type(maxstab_fit_t)::fit
      real(dp)::x0(2),xb(2),fb,jac(size(data,1),size(data,2))
      real(dp),allocatable::d(:),a(:)
      integer::i
      jac=0
      d=euclidean_distances(coord)
      allocate(a(size(d)))
      x0=[log(max(1.0_dp,maxval(d)/2)),0.0_dp]
      if(present(start))x0=[log(start(1)),log(start(2)/(2-start(2)))]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=800,tol=1e-9_dp)
      allocate(fit%par(2))
      fit%par=[exp(xb(1)),2.0_dp/(1.0_dp+exp(-xb(2)))]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,ran,sm
         ran=exp(x(1))
         sm=2.0_dp/(1.0_dp+exp(-x(2)))
         do i=1,size(d)
         a(i)=brown_resnick_a(d(i),ran,sm)
         end do
         v=-lplik_smith(data,a,jac)
      end function
   end function fit_brownresnick_frechet

   function fit_schlather_frechet(data,coord,covmod,start,smooth2) result(fit)
      real(dp),intent(in)::data(:,:),coord(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in),optional::start(:),smooth2
      type(maxstab_fit_t)::fit
      real(dp)::x0(3),xb(3),fb,jac(size(data,1),size(data,2)),s2
      real(dp),allocatable::d(:),rho(:)
      integer::i
      jac=0
      d=euclidean_distances(coord)
      allocate(rho(size(d)))
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      x0=[-2.0_dp,log(max(1.0_dp,maxval(d)/2)),0.0_dp]
      if(present(start))x0=[log(start(1)/(1-start(1))),log(start(2)),log(start(3)/(2-start(3)))]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1000,tol=1e-9_dp)
      allocate(fit%par(3))
      fit%par=[1/(1+exp(-xb(1))),exp(xb(2)),2/(1+exp(-xb(3)))]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,nug,ran,sm
         nug=1/(1+exp(-x(1)))
         ran=exp(x(2))
         sm=2/(1+exp(-x(3)))
         do i=1,size(d)
         rho(i)=covariance_values(d(i),covmod,nug,1-nug,ran,sm,s2,size(coord,2))
         end do
         if(any(rho>0.99999999_dp))then
         v=huge(1.0_dp)/100
         else
         v=-lplik_schlather(data,rho,jac)
         end if
      end function
   end function fit_schlather_frechet

   function fit_extremalt_frechet(data,coord,covmod,start,smooth2) result(fit)
      real(dp),intent(in)::data(:,:),coord(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in),optional::start(:),smooth2
      type(maxstab_fit_t)::fit
      real(dp)::x0(4),xb(4),fb,jac(size(data,1),size(data,2)),s2
      real(dp),allocatable::d(:),rho(:)
      integer::i
      jac=0
      d=euclidean_distances(coord)
      allocate(rho(size(d)))
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      x0=[-2.0_dp,log(max(1.0_dp,maxval(d)/2)),0.0_dp,log(3.0_dp)]
      if(present(start))x0=[log(start(1)/(1-start(1))),log(start(2)),log(start(3)/(2-start(3))),log(start(4))]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1200,tol=1e-8_dp)
      allocate(fit%par(4))
      fit%par=[1/(1+exp(-xb(1))),exp(xb(2)),2/(1+exp(-xb(3))),exp(xb(4))]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,nug,ran,sm,nu
         nug=1/(1+exp(-x(1)))
         ran=exp(x(2))
         sm=2/(1+exp(-x(3)))
         nu=exp(x(4))
         do i=1,size(d)
         rho(i)=covariance_values(d(i),covmod,nug,1-nug,ran,sm,s2,size(coord,2))
         end do
         v=-lplik_extremalt(data,rho,nu,jac)
      end function
   end function fit_extremalt_frechet

   function fit_smith_frechet(data,coord,start) result(fit)
      real(dp),intent(in)::data(:,:),coord(:,:)
      real(dp),intent(in),optional::start(:)
      type(maxstab_fit_t)::fit
      real(dp)::x0(3),xb(3),fb,jac(size(data,1),size(data,2)),cov(2,2),l11,l21,l22
      real(dp),allocatable::dv(:,:),a(:)
      jac=0
      dv=distance_vectors(coord)
      allocate(a(size(dv,1)))
      x0=[0.0_dp,0.0_dp,0.0_dp]
      if(present(start))x0=start(1:3)
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1000,tol=1e-9_dp)
      l11=exp(xb(1))
      l21=xb(2)
      l22=exp(xb(3))
      allocate(fit%par(3))
      fit%par=[l11*l11,l11*l21,l21*l21+l22*l22]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v
         l11=exp(x(1))
         l21=x(2)
         l22=exp(x(3))
         cov=reshape([l11*l11,l11*l21,l11*l21,l21*l21+l22*l22],[2,2])
         a=mahalanobis_distances_2d(dv,cov)
         if(any(a<=0))then
         v=huge(1.0_dp)/100
         else
         v=-lplik_smith(data,a,jac)
         end if
      end function
   end function fit_smith_frechet
end module spatialextremes_fit
