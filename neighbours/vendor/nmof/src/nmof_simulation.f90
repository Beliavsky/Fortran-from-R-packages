! SPDX-License-Identifier: GPL-3.0-only
module nmof_simulation
   use nmof_kinds, only: dp, i8
   use nmof_rng, only: rng_state, rng_seed, rng_normal
   use nmof_math, only: normal_cdf, sort_real
   use nmof_linalg, only: covariance_matrix, cholesky_lower, solve_linear_matrix, eigen_symmetric, column_means, column_sds
   use nmof_types, only: cppi_result, nmof_ok, nmof_invalid_input, nmof_numerical_failure
   implicit none
   private
   public :: geometric_brownian_motion, geometric_brownian_bridge
   public :: random_returns, resample_correlated, cppi
contains
   function geometric_brownian_motion(npaths,timesteps,r,v,tau,s0,exponentiate,antithetic,seed,status) result(paths)
      integer,intent(in)::npaths,timesteps
      real(dp),intent(in)::r,v,tau
      real(dp),intent(in),optional::s0
      logical,intent(in),optional::exponentiate,antithetic
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::status
      real(dp),allocatable::paths(:,:)
      type(rng_state)::rng
      real(dp)::dt,initial,inc
      integer::i,j,half
      logical::expa,anti
      integer(i8)::seedv
      expa=.true.; if(present(exponentiate)) expa=exponentiate
      anti=.false.; if(present(antithetic)) anti=antithetic
      initial=0.0_dp; if(present(s0)) initial=log(s0)
      if(npaths<1.or.timesteps<1.or.v<0.0_dp.or.tau<=0.0_dp.or.(present(s0).and.s0<=0.0_dp)) then
         allocate(paths(0,0)); if(present(status)) status=nmof_invalid_input; return
      end if
      if(anti.and.mod(npaths,2)/=0) then; allocate(paths(0,0)); if(present(status)) status=nmof_invalid_input; return; end if
      seedv=7919_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      dt=tau/real(timesteps,dp); allocate(paths(timesteps+1,npaths)); paths(1,:)=initial
      if(.not.anti) then
         do j=1,npaths; do i=2,timesteps+1
            inc=(r-0.5_dp*v)*dt+sqrt(v*dt)*rng_normal(rng); paths(i,j)=paths(i-1,j)+inc
         end do; end do
      else
         half=npaths/2
         do j=1,half; do i=2,timesteps+1
            inc=sqrt(v*dt)*rng_normal(rng)
            paths(i,j)=paths(i-1,j)+(r-0.5_dp*v)*dt+inc
            paths(i,j+half)=paths(i-1,j+half)+(r-0.5_dp*v)*dt-inc
         end do; end do
      end if
      if(expa) paths=exp(paths)
      if(present(status)) status=nmof_ok
   end function geometric_brownian_motion

   function geometric_brownian_bridge(npaths,timesteps,s0,st,v,tau,log_input,exponentiate,seed,status) result(paths)
      integer,intent(in)::npaths,timesteps
      real(dp),intent(in)::s0,st,v,tau
      logical,intent(in),optional::log_input,exponentiate
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::status
      real(dp),allocatable::paths(:,:),base(:,:)
      real(dp)::a0,aT,w
      logical::logi,expa
      integer::i,j,info
      logi=.false.; if(present(log_input)) logi=log_input
      expa=.true.; if(present(exponentiate)) expa=exponentiate
      if(logi) then; a0=log(s0); aT=log(st); else; a0=s0; aT=st; end if
      allocate(base(timesteps+1,npaths))
      if(logi) then
         base=geometric_brownian_motion(npaths,timesteps,0.0_dp,v,tau,s0,.false.,.false.,seed,info)
      else
         base=geometric_brownian_motion(npaths,timesteps,0.0_dp,v,tau,exp(a0),.true.,.false.,seed,info)
      end if
      if(info/=nmof_ok) then; allocate(paths(0,0)); if(present(status)) status=info; return; end if
      allocate(paths(timesteps+1,npaths)); paths=base
      do j=1,npaths
         do i=2,timesteps+1
            w=real(i-1,dp)/real(timesteps,dp)
            paths(i,j)=base(i,j)+w*(aT-base(timesteps+1,j))
         end do
         paths(1,j)=a0
      end do
      if(logi.and.expa) paths=exp(paths)
      if(present(status)) status=nmof_ok
   end function geometric_brownian_bridge

   function random_returns(n_assets,n_samples,sd,mean,rho,corr,exact,seed,status) result(ret)
      integer,intent(in)::n_assets,n_samples
      real(dp),intent(in)::sd(:)
      real(dp),intent(in),optional::mean(:),rho,corr(:,:)
      logical,intent(in),optional::exact
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::status
      real(dp),allocatable::ret(:,:)
      real(dp),allocatable::z(:,:),c(:,:),l(:,:),white(:,:),cov(:,:)
      real(dp)::rh
      type(rng_state)::rng
      integer::i,j,info
      integer(i8)::seedv
      logical::ex
      ex=.false.; if(present(exact)) ex=exact
      if(n_assets<1.or.n_samples<2.or.(size(sd)/=1.and.size(sd)/=n_assets)) then
         allocate(ret(0,0)); if(present(status)) status=nmof_invalid_input; return
      end if
      seedv=104729_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      allocate(z(n_assets,n_samples)); do j=1,n_samples; do i=1,n_assets; z(i,j)=rng_normal(rng); end do; end do
      if(ex) then
         cov=covariance_matrix(transpose(z)); allocate(l(n_assets,n_assets),white(n_assets,n_samples))
         call cholesky_lower(cov,l,info)
         if(info/=0) then; allocate(ret(0,0)); if(present(status)) status=nmof_numerical_failure; return; end if
         call solve_linear_matrix(l,z,white,info); if(info/=0) then; allocate(ret(0,0)); if(present(status)) status=nmof_numerical_failure; return; end if
         z=white
         do i=1,n_assets
            z(i,:)=z(i,:)-sum(z(i,:))/real(n_samples,dp)
            z(i,:)=z(i,:)/sqrt(sum(z(i,:)**2)/real(n_samples-1,dp))
         end do
      end if
      if(present(corr)) then
         if(any(shape(corr)/=[n_assets,n_assets])) then; allocate(ret(0,0)); if(present(status)) status=nmof_invalid_input; return; end if
         c=corr
      else
         rh=0.0_dp; if(present(rho)) rh=rho
         allocate(c(n_assets,n_assets)); c=rh; do i=1,n_assets; c(i,i)=1.0_dp; end do
      end if
      if(n_assets>1.and.maxval(abs(c-identity(n_assets)))>epsilon(1.0_dp)) then
         if(allocated(l)) deallocate(l)
         allocate(l(n_assets,n_assets)); call cholesky_lower(c,l,info)
         if(info/=0) then; allocate(ret(0,0)); if(present(status)) status=nmof_numerical_failure; return; end if
         z=matmul(l,z)
      end if
      allocate(ret(n_samples,n_assets)); ret=transpose(z)
      do i=1,n_assets
         if(size(sd)==1) then; ret(:,i)=ret(:,i)*sd(1); else; ret(:,i)=ret(:,i)*sd(i); end if
         if(present(mean)) then
            if(size(mean)==1) then; ret(:,i)=ret(:,i)+mean(1); else if(size(mean)==n_assets) then; ret(:,i)=ret(:,i)+mean(i); end if
         end if
      end do
      if(present(status)) status=nmof_ok
   contains
      pure function identity(n) result(a)
         integer,intent(in)::n; real(dp)::a(n,n); integer::ii
         a=0.0_dp; do ii=1,n; a(ii,ii)=1.0_dp; end do
      end function identity
   end function random_returns

   function resample_correlated(data,sample_size,cormat,lengths,seed,status) result(sample)
      real(dp),intent(in)::data(:,:),cormat(:,:)
      integer,intent(in)::sample_size
      integer,intent(in),optional::lengths(:)
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::status
      real(dp),allocatable::sample(:,:)
      real(dp),allocatable::sorted(:,:),x(:,:),pearson(:,:),eval(:),evec(:,:),transform(:,:),u(:,:)
      integer,allocatable::lens(:)
      type(rng_state)::rng
      integer::n,i,j,info,idx,maxlen
      integer(i8)::seedv
      n=size(data,2); maxlen=size(data,1)
      if(sample_size<1.or.any(shape(cormat)/=[n,n])) then; allocate(sample(0,0)); if(present(status)) status=nmof_invalid_input; return; end if
      allocate(lens(n)); if(present(lengths)) then; lens=lengths; else; lens=maxlen; end if
      if(any(lens<1).or.any(lens>maxlen)) then; allocate(sample(0,0)); if(present(status)) status=nmof_invalid_input; return; end if
      allocate(sorted(maxlen,n)); sorted=data
      do j=1,n; call sort_real(sorted(1:lens(j),j)); end do
      seedv=1299709_i8; if(present(seed)) seedv=seed; call rng_seed(rng,seedv)
      allocate(x(sample_size,n)); do j=1,n; do i=1,sample_size; x(i,j)=rng_normal(rng); end do; end do
      pearson=2.0_dp*sin(cormat*acos(-1.0_dp)/6.0_dp)
      allocate(eval(n),evec(n,n)); call eigen_symmetric(pearson,eval,evec,info)
      if(info/=0) then; allocate(sample(0,0)); if(present(status)) status=nmof_numerical_failure; return; end if
      eval=max(eval,0.0_dp); allocate(transform(n,n)); transform=0.0_dp
      do i=1,n; transform(:,i)=evec(:,i)*sqrt(eval(i)); end do
      u=normal_cdf(matmul(x,transpose(transform))); allocate(sample(sample_size,n))
      do j=1,n; do i=1,sample_size
         idx=ceiling(u(i,j)*real(lens(j),dp)); idx=max(1,min(lens(j),idx)); sample(i,j)=sorted(idx,j)
      end do; end do
      if(present(status)) status=nmof_ok
   end function resample_correlated

   function cppi(spot,multiplier,floor_value,r,tau,gap) result(ans)
      real(dp),intent(in)::spot(:),multiplier,floor_value,r
      real(dp),intent(in),optional::tau
      integer,intent(in),optional::gap
      type(cppi_result)::ans
      real(dp)::tt,dt
      integer::t,n,g
      tt=1.0_dp; if(present(tau)) tt=tau; g=1; if(present(gap)) g=gap
      n=size(spot)-1
      if(n<1.or.g<1.or.any(spot<=0.0_dp)) then; ans%status=nmof_invalid_input; return; end if
      allocate(ans%value(n+1),ans%cushion(n+1),ans%bond(n+1),ans%floor(n+1), &
               ans%exposure(n+1),ans%units(n+1),ans%spot(n+1)); ans%spot=spot; dt=tt/real(n,dp)
      do t=0,n; ans%floor(t+1)=floor_value*exp(-r*tt*(1.0_dp-real(t,dp)/real(n,dp))); end do
      ans%value(1)=1.0_dp; ans%cushion(1)=ans%value(1)-ans%floor(1); ans%exposure(1)=multiplier*ans%cushion(1)
      ans%units(1)=ans%exposure(1)/spot(1); ans%bond(1)=ans%value(1)-ans%exposure(1)
      do t=2,n+1
         ans%bond(t)=ans%bond(t-1)*exp(r*dt); ans%value(t)=ans%units(t-1)*spot(t)+ans%bond(t)
         ans%cushion(t)=max(0.0_dp,ans%value(t)-ans%floor(t))
         if(mod(t-1,g)==0) then
            ans%exposure(t)=min(multiplier*ans%cushion(t),ans%value(t)); ans%units(t)=ans%exposure(t)/spot(t)
            ans%bond(t)=ans%value(t)-ans%exposure(t)
         else
            ans%exposure(t)=ans%value(t)-ans%bond(t); ans%units(t)=ans%units(t-1)
         end if
      end do
      ans%status=nmof_ok
   end function cppi
end module nmof_simulation
