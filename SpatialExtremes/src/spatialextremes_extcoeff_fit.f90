module spatialextremes_extcoeff_fit
   use spatialextremes_base, only: dp,pair_count,pair_indices,euclidean_distances,distance_vectors
   use spatialextremes_dependence, only: smith_extcoeff
   use spatialextremes_covariance
   use spatialextremes_models
   use spatialextremes_fit, only: maxstab_fit_t,nelder_mead,gevmle
   use spatialextremes_univariate, only: gev_to_frechet
   implicit none
   private
   public :: smith_extcoeff_frechet, st_extcoeff_frechet, empirical_frechet_margins
   public :: fitextcoeff_smith_empirical, fitextcoeff_st_empirical
   public :: lsfit_schlather,lsfit_schlather_ind,lsfit_brownresnick
   public :: lsfit_geomgauss,lsfit_extremalt,lsfit_smith_2d
contains
   function smith_extcoeff_frechet(frech) result(theta)
      ! High-level Smith extremal-coefficient estimator for unit-Frechet
      ! observations.  Upstream native extCoeffSmith receives 1/Z, i.e.
      ! standard-exponential margins, so invert here before calling the
      ! exact low-level recurrence.
      real(dp),intent(in)::frech(:,:)
      real(dp)::theta(pair_count(size(frech,2)))
      theta=smith_extcoeff(1.0_dp/frech)
   end function smith_extcoeff_frechet

   function st_extcoeff_frechet(frech,prob) result(theta)
      ! Schlather-Tawn estimator.  This minimizes the same one-parameter
      ! censored objective as upstream extCoeffST, analytically and then
      ! applies its [1,2] optimization interval.
      real(dp),intent(in)::frech(:,:)
      real(dp),intent(in),optional::prob
      real(dp)::theta(pair_count(size(frech,2)))
      real(dp)::xbar(size(frech,2)),p,z,scaled,s
      integer::k,i,j,t,m
      p=0.0_dp
      if(present(prob))p=prob
      if(p>0.0_dp)then
      z=-1.0_dp/log(p)
      else
      z=0.0_dp
      end if
      do j=1,size(frech,2)
      xbar(j)=sum(1.0_dp/frech(:,j))/real(size(frech,1),dp)
      end do
      do k=1,size(theta)
         call pair_indices(k,size(frech,2),i,j)
         m=0
         s=0.0_dp
         do t=1,size(frech,1)
            scaled=max(frech(t,i)*xbar(i),frech(t,j)*xbar(j))
            if(scaled>z)then
               m=m+1
               s=s+1.0_dp/scaled
            else if(z>0.0_dp)then
               s=s+1.0_dp/z
            end if
         end do
         if(s<=0.0_dp .or. m<=0)then
         theta(k)=2.0_dp
         else
         theta(k)=min(2.0_dp,max(1.0_dp,real(m,dp)/s))
         end if
      end do
   end function st_extcoeff_frechet

   function empirical_frechet_margins(data) result(frech)
      real(dp),intent(in)::data(:,:)
      real(dp)::frech(size(data,1),size(data,2)),u
      integer::i,j,k,r,n
      n=size(data,1)
      do j=1,size(data,2)
         do i=1,n
            r=1
            do k=1,n
               if(data(k,j)<data(i,j))r=r+1
            end do
            ! R's rank() averages ties.  For exact ties compute the average
            ! of the strict-lower and lower-or-equal ranks.
            k=count(data(:,j)==data(i,j))
            if(k>1)then
               u=(real(r,dp)+0.5_dp*real(k-1,dp))/real(n+1,dp)
            else
               u=real(r,dp)/real(n+1,dp)
            end if
            frech(i,j)=-1.0_dp/log(u)
         end do
      end do
   end function empirical_frechet_margins

   function fitextcoeff_smith_empirical(data) result(theta)
      real(dp),intent(in)::data(:,:)
      real(dp)::theta(pair_count(size(data,2)))
      theta=smith_extcoeff_frechet(empirical_frechet_margins(data))
   end function fitextcoeff_smith_empirical

   function fitextcoeff_st_empirical(data,prob) result(theta)
      real(dp),intent(in)::data(:,:)
      real(dp),intent(in),optional::prob
      real(dp)::theta(pair_count(size(data,2)))
      if(present(prob))then
         theta=st_extcoeff_frechet(empirical_frechet_margins(data),prob)
      else
         theta=st_extcoeff_frechet(empirical_frechet_margins(data))
      end if
   end function fitextcoeff_st_empirical

   function lsfit_schlather(extcoeff,coord,covmod,weights,start,smooth2) result(fit)
      real(dp),intent(in)::extcoeff(:),coord(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in),optional::weights(:),start(:),smooth2
      type(maxstab_fit_t)::fit
      real(dp)::x0(3),xb(3),fb,s2
      real(dp),allocatable::d(:),w(:),rho(:)
      integer::i
      d=euclidean_distances(coord)
      allocate(w(size(d)),rho(size(d)))
      w=1.0_dp
      if(present(weights))w=weights
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      x0=[-2.0_dp,log(max(1.0_dp,maxval(d)/2.0_dp)),0.0_dp]
      if(present(start))x0=[log(start(1)/(1-start(1))),log(start(2)),log(start(3)/(2-start(3)))]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1500,tol=1.0e-10_dp)
      allocate(fit%par(3))
      fit%par=[logistic(xb(1)),exp(xb(2)),2.0_dp*logistic(xb(3))]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,nug,ran,sm,res
         nug=logistic(x(1))
         ran=exp(x(2))
         sm=2.0_dp*logistic(x(3))
         v=0.0_dp
         do i=1,size(d)
            rho(i)=covariance_values(d(i),covmod,nug,1.0_dp-nug,ran,sm,s2,size(coord,2))
            res=extremal_coefficient_schlather(rho(i))-extcoeff(i)
            v=v+(res/w(i))**2
         end do
      end function
   end function lsfit_schlather

   function lsfit_schlather_ind(extcoeff,coord,covmod,weights,start,smooth2) result(fit)
      real(dp),intent(in)::extcoeff(:),coord(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in),optional::weights(:),start(:),smooth2
      type(maxstab_fit_t)::fit
      real(dp)::x0(4),xb(4),fb,s2
      real(dp),allocatable::d(:),w(:),rho(:)
      integer::i
      d=euclidean_distances(coord)
      allocate(w(size(d)),rho(size(d)))
      w=1.0_dp
      if(present(weights))w=weights
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      x0=[-2.0_dp,-2.0_dp,log(max(1.0_dp,maxval(d)/2.0_dp)),0.0_dp]
      if(present(start))x0=[logit01(start(1)),logit01(start(2)),log(start(3)),logit2(start(4))]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1800,tol=1.0e-10_dp)
      allocate(fit%par(4))
      fit%par=[logistic(xb(1)),logistic(xb(2)),exp(xb(3)),2.0_dp*logistic(xb(4))]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,alpha,nug,ran,sm,res
         alpha=logistic(x(1))
         nug=logistic(x(2))
         ran=exp(x(3))
         sm=2.0_dp*logistic(x(4))
         v=0.0_dp
         do i=1,size(d)
            rho(i)=covariance_values(d(i),covmod,nug,1.0_dp-nug,ran,sm,s2,size(coord,2))
            res=extremal_coefficient_schlather_ind(alpha,rho(i))-extcoeff(i)
            v=v+(res/w(i))**2
         end do
      end function
   end function lsfit_schlather_ind

   function lsfit_brownresnick(extcoeff,coord,weights,start) result(fit)
      real(dp),intent(in)::extcoeff(:),coord(:,:)
      real(dp),intent(in),optional::weights(:),start(:)
      type(maxstab_fit_t)::fit
      real(dp)::x0(2),xb(2),fb
      real(dp),allocatable::d(:),w(:)
      integer::i
      d=euclidean_distances(coord)
      allocate(w(size(d)))
      w=1.0_dp
      if(present(weights))w=weights
      x0=[log(max(1.0_dp,maxval(d)/2.0_dp)),0.0_dp]
      if(present(start))x0=[log(start(1)),logit2(start(2))]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1200,tol=1.0e-10_dp)
      allocate(fit%par(2))
      fit%par=[exp(xb(1)),2.0_dp*logistic(xb(2))]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,ran,sm,res,a
         ran=exp(x(1))
         sm=2.0_dp*logistic(x(2))
         v=0.0_dp
         do i=1,size(d)
            a=brown_resnick_a(d(i),ran,sm)
            res=extremal_coefficient_smith(a)-extcoeff(i)
            v=v+(res/w(i))**2
         end do
      end function
   end function lsfit_brownresnick

   function lsfit_geomgauss(extcoeff,coord,covmod,weights,start,smooth2) result(fit)
      real(dp),intent(in)::extcoeff(:),coord(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in),optional::weights(:),start(:),smooth2
      type(maxstab_fit_t)::fit
      real(dp)::x0(4),xb(4),fb,s2
      real(dp),allocatable::d(:),w(:)
      integer::i
      d=euclidean_distances(coord)
      allocate(w(size(d)))
      w=1.0_dp
      if(present(weights))w=weights
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      x0=[log(1.0_dp),-2.0_dp,log(max(1.0_dp,maxval(d)/2.0_dp)),0.0_dp]
      if(present(start))x0=[log(start(1)),logit01(start(2)),log(start(3)),logit2(start(4))]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1800,tol=1.0e-10_dp)
      allocate(fit%par(4))
      fit%par=[exp(xb(1)),logistic(xb(2)),exp(xb(3)),2.0_dp*logistic(xb(4))]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,sig,nug,ran,sm,res,a
         sig=exp(x(1))
         nug=logistic(x(2))
         ran=exp(x(3))
         sm=2.0_dp*logistic(x(4))
         v=0.0_dp
         do i=1,size(d)
            a=geom_gauss_a(d(i),covmod,sig,nug,ran,sm,s2,size(coord,2))
            res=extremal_coefficient_smith(a)-extcoeff(i)
            v=v+(res/w(i))**2
         end do
      end function
   end function lsfit_geomgauss

   function lsfit_extremalt(extcoeff,coord,covmod,weights,start,smooth2) result(fit)
      real(dp),intent(in)::extcoeff(:),coord(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in),optional::weights(:),start(:),smooth2
      type(maxstab_fit_t)::fit
      real(dp)::x0(4),xb(4),fb,s2
      real(dp),allocatable::d(:),w(:)
      integer::i
      d=euclidean_distances(coord)
      allocate(w(size(d)))
      w=1.0_dp
      if(present(weights))w=weights
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      x0=[-2.0_dp,log(max(1.0_dp,maxval(d)/2.0_dp)),0.0_dp,log(3.0_dp)]
      if(present(start))x0=[logit01(start(1)),log(start(2)),logit2(start(3)),log(start(4))]
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1800,tol=1.0e-10_dp)
      allocate(fit%par(4))
      fit%par=[logistic(xb(1)),exp(xb(2)),2.0_dp*logistic(xb(3)),exp(xb(4))]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,nug,ran,sm,nu,res,rho
         nug=logistic(x(1))
         ran=exp(x(2))
         sm=2.0_dp*logistic(x(3))
         nu=exp(x(4))
         v=0.0_dp
         do i=1,size(d)
            rho=covariance_values(d(i),covmod,nug,1.0_dp-nug,ran,sm,s2,size(coord,2))
            res=extremal_coefficient_extremalt(rho,nu)-extcoeff(i)
            v=v+(res/w(i))**2
         end do
      end function
   end function lsfit_extremalt

   function lsfit_smith_2d(extcoeff,coord,weights,start) result(fit)
      real(dp),intent(in)::extcoeff(:),coord(:,:)
      real(dp),intent(in),optional::weights(:),start(:)
      type(maxstab_fit_t)::fit
      real(dp)::x0(3),xb(3),fb,cov(2,2),l11,l21,l22
      real(dp),allocatable::dv(:,:),w(:),a(:)
      integer::i
      dv=distance_vectors(coord)
      allocate(w(size(dv,1)),a(size(dv,1)))
      w=1.0_dp
      if(present(weights))w=weights
      x0=0.0_dp
      if(present(start))x0=start(1:3)
      call nelder_mead(obj,x0,xb,fb,fit%iterations,fit%convergence,maxit=1500,tol=1.0e-10_dp)
      l11=exp(xb(1))
      l21=xb(2)
      l22=exp(xb(3))
      allocate(fit%par(3))
      fit%par=[l11*l11,l11*l21,l21*l21+l22*l22]
      fit%loglik=-fb
   contains
      function obj(x) result(v)
         real(dp),intent(in)::x(:)
         real(dp)::v,res
         l11=exp(x(1))
         l21=x(2)
         l22=exp(x(3))
         cov=reshape([l11*l11,l11*l21,l11*l21,l21*l21+l22*l22],[2,2])
         a=mahalanobis_distances_2d(dv,cov)
         v=0.0_dp
         do i=1,size(a)
         res=extremal_coefficient_smith(a(i))-extcoeff(i)
         v=v+(res/w(i))**2
         end do
      end function
   end function lsfit_smith_2d

   pure elemental real(dp) function logistic(x) result(y)
      real(dp),intent(in)::x
      if(x>=0.0_dp)then
      y=1.0_dp/(1.0_dp+exp(-x))
      else
      y=exp(x)/(1.0_dp+exp(x))
      end if
   end function
   pure real(dp) function logit01(x) result(y)
      real(dp),intent(in)::x
      y=log(x/(1.0_dp-x))
   end function
   pure real(dp) function logit2(x) result(y)
      real(dp),intent(in)::x
      y=log(x/(2.0_dp-x))
   end function
end module spatialextremes_extcoeff_fit
