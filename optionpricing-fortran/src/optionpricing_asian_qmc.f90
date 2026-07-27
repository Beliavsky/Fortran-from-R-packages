! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_asian_qmc
   use optionpricing_kinds, only : dp
   use optionpricing_math, only : normal_quantile, mean_value, sample_sd, normal_cdf, bisection_root
   use optionpricing_random, only : seed_rng
   use optionpricing_linalg, only : symmetric_eigen, identity_matrix, lower_ones, orthonormal_complete, least_squares
   use optionpricing_types, only : greeks_result, conditional_result, moments_result
   use optionpricing_asian_analytic, only : eval_ecv, eval_lb, eval_eqcv
   use optionpricing_asian_mc, only : asian_call_naive_greeks_z, conditional_estimates_z, &
      asian_call_naive_mc, asian_call_best_mc
   implicit none
   private
   public :: korobov_lattice, randomized_korobov_normals
   public :: naive_pca_matrix, conditional_generation_matrix
   public :: asian_call_naive_qmc, asian_call_best_qmc, asian_call
contains

   function korobov_lattice(n,a,d) result(u)
      integer, intent(in) :: n,a,d
      real(dp), allocatable :: u(:,:)
      integer, allocatable :: z(:)
      integer :: i,j
      allocate(u(d,n),z(d)); u=0.0_dp; z=1
      do j=2,d
         z(j)=modulo(a*z(j-1),n)
      end do
      do i=2,n
         u(:,i)=modulo(u(:,i-1)+real(z,dp)/real(n,dp),1.0_dp)
      end do
   end function korobov_lattice

   subroutine randomized_korobov_normals(lattice,z,baker)
      real(dp), intent(in) :: lattice(:,:)
      real(dp), intent(out) :: z(:,:)
      logical, intent(in), optional :: baker
      real(dp), allocatable :: shift(:),u(:,:)
      logical :: do_baker
      integer :: i,j,d,n
      d=size(lattice,1); n=size(lattice,2)
      if(any(shape(z)/=[d,n])) error stop 'randomized_korobov_normals: shape mismatch'
      allocate(shift(d),u(d,n)); call random_number(shift)
      u=modulo(lattice+spread(shift,2,n),1.0_dp)
      do_baker=.true.; if(present(baker))do_baker=baker
      if(do_baker) then
         where(u<=0.5_dp)
            u=2.0_dp*u
         elsewhere
            u=2.0_dp*(1.0_dp-u)
         end where
      end if
      do j=1,n
         do i=1,d
            z(i,j)=normal_quantile(min(max(u(i,j),1.0_dp/real(n*n,dp)),1.0_dp-1.0_dp/real(n*n,dp)))
         end do
      end do
   end subroutine randomized_korobov_normals

   function naive_pca_matrix(d,status) result(qmat)
      integer, intent(in) :: d
      integer, intent(out), optional :: status
      real(dp), allocatable :: qmat(:,:)
      real(dp), allocatable :: l(:,:),values(:),vectors(:,:)
      integer :: istat,j
      allocate(qmat(d,d),l(d,d),values(d),vectors(d,d))
      l=lower_ones(d)
      call symmetric_eigen(matmul(l,transpose(l)),values,vectors,istat)
      do j=1,d
         qmat(:,j)=vectors(:,j)*sqrt(max(0.0_dp,values(j)))
      end do
      if(present(status))status=istat
   end function naive_pca_matrix

   function conditional_generation_matrix(t,d,k,r,sigma,s0,mode,dirnum,status) result(qmat)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      character(len=*), intent(in) :: mode
      integer, intent(in), optional :: dirnum
      integer, intent(out), optional :: status
      real(dp), allocatable :: qmat(:,:)
      real(dp), allocatable :: l(:,:),v(:),bmat(:,:),values(:),vectors(:,:),first(:,:),amat(:,:)
      integer :: istat,nd
      nd=1; if(present(dirnum))nd=max(1,min(dirnum,d))
      allocate(l(d,d),v(d),bmat(d,d),values(d),vectors(d,d))
      l=lower_ones(d); call conditional_vector(d,v)
      bmat=matmul(l,identity_matrix(d)-outer(v,v))
      select case(trim(adjustl(mode)))
      case('pca')
         call symmetric_eigen(matmul(bmat,transpose(bmat)),values,vectors,istat)
         allocate(qmat(d,max(1,d-1)))
         qmat=0.0_dp
         if(d>1) call scaled_eigen_columns(vectors,values,qmat,d-1)
      case('pcamain')
         call symmetric_eigen(matmul(bmat,transpose(bmat)),values,vectors,istat)
         allocate(first(d,nd),amat(d,d),qmat(d,d))
         first=matmul(transpose(bmat),vectors(:,1:nd))
         call orthonormal_complete(first,amat,istat)
         qmat=matmul(bmat,amat)
      case('lt')
         allocate(amat(d,d),qmat(d,d))
         call lt_rotation(bmat,t,d,k,r,sigma,s0,nd,amat,istat)
         qmat=matmul(bmat,amat)
      case('ltpca')
         call symmetric_eigen(matmul(bmat,transpose(bmat)),values,vectors,istat)
         if(d<=1) then
            allocate(qmat(d,1)); qmat=0.0_dp
         else
            allocate(first(d,d-1))
            call scaled_eigen_columns(vectors,values,first,d-1)
            allocate(amat(d-1,d-1),qmat(d,d-1))
            call lt_rotation(first,t,d,k,r,sigma,s0,d-1,amat,istat)
            qmat=matmul(first,amat)
         end if
      case default
         istat=1; allocate(qmat(d,d)); qmat=0.0_dp
      end select
      if(present(status))status=istat
   end function conditional_generation_matrix

   function asian_call_naive_qmc(nout,n,a_gen,t,d,k,r,sigma,s0,genmethod,baker,seed) result(res)
      integer, intent(in) :: nout,n,a_gen,d
      real(dp), intent(in) :: t,k,r,sigma,s0
      character(len=*), intent(in), optional :: genmethod
      logical, intent(in), optional :: baker
      integer, intent(in), optional :: seed
      type(greeks_result) :: res
      real(dp), allocatable :: lattice(:,:),z(:,:),qmat(:,:),rep(:,:)
      real(dp) :: estimate(3)
      character(len=16) :: mode
      logical :: do_baker
      integer :: i,status
      if(nout<2 .or. n<2 .or. d<1) then
         res%status=1; res%message='invalid QMC dimensions'; return
      end if
      mode='pca'; if(present(genmethod))mode=adjustl(genmethod)
      if(trim(mode)/='pca' .and. trim(mode)/='std') then
         res%status=2; res%message='naive QMC supports pca or std'; return
      end if
      do_baker=.true.; if(present(baker))do_baker=baker
      if(present(seed))call seed_rng(seed)
      lattice=korobov_lattice(n,a_gen,d); allocate(z(d,n),rep(nout,3))
      if(trim(mode)=='pca')qmat=naive_pca_matrix(d,status)
      do i=1,nout
         call randomized_korobov_normals(lattice,z,do_baker)
         if(trim(mode)=='pca') then
            call asian_call_naive_greeks_z(z,t,k,r,sigma,s0,estimate,qmat,'pca',status)
         else
            call asian_call_naive_greeks_z(z,t,k,r,sigma,s0,estimate,mode='std',status=status)
         end if
         if(status/=0) then
            res%status=status; res%message='QMC path transform failed'; return
         end if
         rep(i,:)=estimate
      end do
      do i=1,3
         res%estimate(i)=mean_value(rep(:,i))
         res%error95(i)=1.96_dp*sample_sd(rep(:,i))/sqrt(real(nout,dp))
      end do
   end function asian_call_naive_qmc

   function asian_call_best_qmc(nout,n,a_gen,t,d,k,r,sigma,s0,genmethod,dirnum,baker,cvmethod,seed,noutp,maxiter,tol) result(res)
      integer, intent(in) :: nout,n,a_gen,d
      real(dp), intent(in) :: t,k,r,sigma,s0
      character(len=*), intent(in), optional :: genmethod,cvmethod
      integer, intent(in), optional :: dirnum,seed,noutp,maxiter
      logical, intent(in), optional :: baker
      real(dp), intent(in), optional :: tol
      type(greeks_result) :: res
      real(dp), allocatable :: lattice(:,:),z(:,:),qmat(:,:),y(:,:),x(:,:),adjusted(:,:),beta(:,:)
      real(dp), allocatable :: xtrain(:,:),ytrain(:,:)
      real(dp) :: expected(6),ewvec(3),eps
      type(conditional_result) :: one
      type(moments_result) :: lb,eq,ew
      character(len=16) :: mode,cv
      logical :: do_baker
      integer :: total,pilot,i,j,status,nd,iters,ii
      logical, allocatable :: keep(:)
      mode='pca'; if(present(genmethod))mode=adjustl(genmethod)
      cv='splitting'; if(present(cvmethod))cv=adjustl(cvmethod)
      nd=1; if(present(dirnum))nd=dirnum
      do_baker=.true.; if(present(baker))do_baker=baker
      eps=1.0e-14_dp; if(present(tol))eps=tol
      iters=100; if(present(maxiter))iters=maxiter
      pilot=0
      if(trim(cv)=='pilotrun') then
         pilot=max(8,nout); if(present(noutp))pilot=noutp
      end if
      total=nout+pilot
      if(nout<2 .or. n<2 .or. total<2) then
         res%status=1; res%message='invalid QMC dimensions'; return
      end if
      if((trim(cv)=='direct' .or. trim(cv)=='splitting') .and. nout<8) then
         res%status=2; res%message='at least eight QMC replications are required for six controls'; return
      end if
      if(present(seed))call seed_rng(seed)
      lattice=korobov_lattice(n,a_gen,d); allocate(z(d,n),y(total,3),x(total,6))
      if(trim(mode)/='std') then
         qmat=conditional_generation_matrix(t,d,k,r,sigma,s0,mode,nd,status)
         if(status/=0) then
            res%status=3; res%message='generation matrix construction failed'; return
         end if
      end if
      do i=1,total
         call randomized_korobov_normals(lattice,z,do_baker)
         if(trim(mode)=='std') then
            one=conditional_estimates_z(z,t,d,k,r,sigma,s0,'std',maxiter=iters,tol=eps)
         else
            one=conditional_estimates_z(z,t,d,k,r,sigma,s0,mode,qmat,iters,eps)
         end if
         if(one%status/=0) then
            res%status=10+one%status; res%message=one%message; return
         end if
         y(i,:)=one%y; x(i,:)=one%controls
      end do
      lb=eval_lb(t,d,k,r,sigma,s0,.false.); eq=eval_eqcv(t,d,k,r,sigma,s0); ew=eval_ecv(t,d,k,r,sigma,s0)
      expected=[lb%price,eq%price,exp(r*t)*lb%delta,exp(r*t)*lb%gamma,eq%delta,eq%gamma]
      ewvec=[ew%price,ew%delta,ew%gamma]
      allocate(adjusted(nout,3),beta(7,3))
      select case(trim(cv))
      case('direct')
         call least_squares(x(1:nout,:),y(1:nout,:),beta,status,.true.)
         if(status==1) then
            res%status=20; res%message='direct control regression failed'; return
         end if
         adjusted=y(1:nout,:)-matmul(x(1:nout,:)-spread(expected,1,nout),beta(2:7,:))
      case('splitting')
         allocate(keep(nout),xtrain(nout-1,6),ytrain(nout-1,3))
         do i=1,nout
            keep=.true.; keep(i)=.false.; ii=0
            do j=1,nout
               if(keep(j)) then
                  ii=ii+1; xtrain(ii,:)=x(j,:); ytrain(ii,:)=y(j,:)
               end if
            end do
            call least_squares(xtrain,ytrain,beta,status,.true.)
            if(status==1) then
               res%status=21; res%message='split control regression failed'; return
            end if
            adjusted(i,:)=y(i,:)-matmul(x(i,:)-expected,beta(2:7,:))
         end do
      case('pilotrun')
         if(pilot<8) then
            res%status=22; res%message='pilotrun requires at least eight pilot replications'; return
         end if
         call least_squares(x(1:pilot,:),y(1:pilot,:),beta,status,.true.)
         if(status==1) then
            res%status=23; res%message='pilot control regression failed'; return
         end if
         adjusted=y(pilot+1:total,:)-matmul(x(pilot+1:total,:)-spread(expected,1,nout),beta(2:7,:))
      case default
         res%status=24; res%message='unknown control-variate method'; return
      end select
      adjusted=adjusted+spread(ewvec,1,nout)
      do j=1,3
         res%estimate(j)=mean_value(adjusted(:,j))
         res%error95(j)=1.96_dp*sample_sd(adjusted(:,j))/sqrt(real(nout,dp))
      end do
   end function asian_call_best_qmc

   function asian_call(t,d,k,r,sigma,s0,method,sampling,n,nout,a_gen,np,seed,baker, &
      genmethod,dirnum,cvmethod,maxiter,tol) result(res)
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in) :: d
      character(len=*), intent(in) :: method,sampling
      integer, intent(in), optional :: n,nout,a_gen,np,seed,dirnum,maxiter
      logical, intent(in), optional :: baker
      character(len=*), intent(in), optional :: genmethod,cvmethod
      real(dp), intent(in), optional :: tol
      type(greeks_result) :: res
      integer :: nn,nr,ag,pilot,sd,nd,iters
      real(dp) :: eps
      logical :: bk
      character(len=16) :: gm,cv
      nn=10000; if(present(n))nn=n
      nr=50; if(present(nout))nr=nout
      ag=1487; if(present(a_gen))ag=a_gen
      pilot=100; if(present(np))pilot=np
      sd=4711; if(present(seed))sd=seed
      nd=1; if(present(dirnum))nd=dirnum
      iters=100; if(present(maxiter))iters=maxiter
      eps=1.0e-14_dp; if(present(tol))eps=tol
      bk=.true.; if(present(baker))bk=baker
      gm='pca'; if(present(genmethod))gm=adjustl(genmethod)
      cv='splitting'; if(present(cvmethod))cv=adjustl(cvmethod)
      if(trim(adjustl(method))=='best' .and. trim(adjustl(sampling))=='QMC') then
         res=asian_call_best_qmc(nr,nn,ag,t,d,k,r,sigma,s0,gm,nd,bk,cv,sd,maxiter=iters,tol=eps)
      else if(trim(adjustl(method))=='best' .and. trim(adjustl(sampling))=='MC') then
         res=asian_call_best_mc(nn,t,d,k,r,sigma,s0,pilot,iters,eps,sd)
      else if(trim(adjustl(method))=='naive' .and. trim(adjustl(sampling))=='QMC') then
         res=asian_call_naive_qmc(nr,nn,ag,t,d,k,r,sigma,s0,gm,bk,sd)
      else if(trim(adjustl(method))=='naive' .and. trim(adjustl(sampling))=='MC') then
         res=asian_call_naive_mc(nn,t,d,k,r,sigma,s0,sd)
      else
         res%status=1; res%message='method must be best/naive and sampling MC/QMC'
      end if
   end function asian_call

   subroutine conditional_vector(d,v)
      integer, intent(in) :: d
      real(dp), intent(out) :: v(d)
      real(dp) :: varx
      integer :: i
      varx=real(d*(d+1)*(2*d+1),dp)/6.0_dp
      do i=1,d
         v(i)=real(d-i+1,dp)/sqrt(varx)
      end do
   end subroutine conditional_vector

   pure function outer(x,y) result(a)
      real(dp), intent(in) :: x(:),y(:)
      real(dp) :: a(size(x),size(y))
      a=spread(x,2,size(y))*spread(y,1,size(x))
   end function outer

   subroutine scaled_eigen_columns(vectors,values,qmat,ncols)
      real(dp), intent(in) :: vectors(:,:),values(:)
      real(dp), intent(out) :: qmat(:,:)
      integer, intent(in) :: ncols
      integer :: j
      qmat=0.0_dp
      do j=1,ncols
         qmat(:,j)=vectors(:,j)*sqrt(max(0.0_dp,values(j)))
      end do
   end subroutine scaled_eigen_columns

   subroutine lt_rotation(base,t,d,k,r,sigma,s0,ndir,amat,status)
      real(dp), intent(in) :: base(:,:),t,k,r,sigma,s0
      integer, intent(in) :: d,ndir
      real(dp), intent(out) :: amat(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: eps(:),z(:),s(:),avec(:),m(:),gradient(:),prefix(:,:),completed(:,:)
      real(dp) :: dt,varx,vol,zcut,mus,sigmas,broot
      integer :: mbase,i,j,istat,nuse
      mbase=size(base,2); status=0; amat=identity_matrix(mbase)
      allocate(eps(mbase),z(d),s(d),avec(d),m(d),gradient(mbase))
      dt=t/real(d,dp); varx=real(d*(d+1)*(2*d+1),dp)/6.0_dp; vol=sigma*sqrt(dt)
      mus=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(d+1,dp)/2.0_dp
      sigmas=sigma/real(d,dp)*sqrt(dt*varx); zcut=(log(k)-mus)/sigmas
      do i=1,d
         avec(i)=vol*sum([(real(d-j+1,dp)/sqrt(varx),j=1,i)])
      end do
      nuse=min(ndir,mbase)
      amat=1.0_dp
      do i=1,nuse
         eps=0.0_dp
         if(i>1)eps(1:i-1)=1.0_dp
         z=matmul(base,matmul(amat,eps))
         do j=1,d
            s(j)=s0*exp((r-0.5_dp*sigma*sigma)*real(j,dp)*dt+vol*z(j))
         end do
         broot=bisection_root(rootfun,-20.0_dp,max(20.0_dp,zcut+5.0_dp),1.0e-12_dp,200,istat)
         if(istat/=0) then
            status=1; return
         end if
         m=exp(0.5_dp*avec*avec)*(normal_cdf(zcut-avec)-normal_cdf(broot-avec))
         do j=1,mbase
            gradient(j)=vol/real(d,dp)*sum(m*s*base(:,j))
         end do
         amat(:,i)=gradient
         allocate(prefix(mbase,i),completed(mbase,mbase))
         prefix=amat(:,1:i)
         call orthonormal_complete(prefix,completed,istat)
         amat(:,1:i)=completed(:,1:i)
         deallocate(prefix,completed)
      end do
      if(nuse<mbase) then
         allocate(prefix(mbase,nuse),completed(mbase,mbase))
         prefix=amat(:,1:nuse)
         call orthonormal_complete(prefix,completed,istat)
         amat=completed
      end if
   contains
      function rootfun(x) result(value)
         real(dp), intent(in) :: x
         real(dp) :: value
         value=sum(exp(avec*x)*s)/real(d,dp)-k
      end function rootfun
   end subroutine lt_rotation
end module optionpricing_asian_qmc
