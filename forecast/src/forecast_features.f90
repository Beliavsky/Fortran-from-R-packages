module forecast_features
   use forecast_kinds, only : dp, pi
   use forecast_stats, only : acf_values, pacf_values, mean_value, quantile_type8
   use forecast_linalg, only : symmetric_eigen, cholesky_lower, solve_linear
   implicit none
   private
   public :: fourier_terms, seasonal_dummy, moving_average, difference_series, seasonal_difference
   public :: acf_values, pacf_values, findfrequency, fourier_terms_multi, fourier_forecast_terms
   public :: seasonal_dummy_forecast, seasonal_index_forecast, tapered_acf, tapered_pacf
   public :: linear_process_bootstrap, tapered_correlation_ci
contains
   function fourier_terms_multi(n,periods,k,start_index) result(x)
      integer,intent(in)::n,k(:)
      real(dp),intent(in)::periods(:)
      integer,intent(in),optional::start_index
      real(dp),allocatable::x(:,:),freq(:),uniq(:)
      integer::j,h,nf,u,col,i,s
      real(dp)::f,ang,tol
      if(size(periods)/=size(k))error stop 'fourier_terms_multi: periods/orders mismatch'
      if(any(k<0) .or. any(2.0_dp*real(k,dp)>periods+1.0e-12_dp))error stop 'fourier_terms_multi: invalid order'
      allocate(freq(sum(k)))
      nf=0
      do j=1,size(periods)
         do h=1,k(j)
            nf=nf+1
            freq(nf)=real(h,dp)/periods(j)
         end do
      end do
      allocate(uniq(nf))
      u=0
      tol=100.0_dp*epsilon(1.0_dp)
      do j=1,nf
         if(u==0 .or. all(abs(uniq(1:u)-freq(j))>tol*max(1.0_dp,abs(freq(j)))))then
            u=u+1
            uniq(u)=freq(j)
         end if
      end do
      col=0
      do j=1,u
         if(abs(2.0_dp*uniq(j)-anint(2.0_dp*uniq(j)))>epsilon(1.0_dp))col=col+1
         col=col+1
      end do
      allocate(x(n,col))
      s=1
      if(present(start_index))s=start_index
      col=0
      do j=1,u
         f=uniq(j)
         if(abs(2.0_dp*f-anint(2.0_dp*f))>epsilon(1.0_dp))then
            col=col+1
            do i=1,n
               ang=2.0_dp*pi*f*real(s+i-1,dp)
               x(i,col)=sin(ang)
            end do
         end if
         col=col+1
         do i=1,n
            ang=2.0_dp*pi*f*real(s+i-1,dp)
            x(i,col)=cos(ang)
         end do
      end do
   end function fourier_terms_multi

   function fourier_forecast_terms(n_history,h,period,k,start_index) result(x)
      integer,intent(in)::n_history,h,period,k
      integer,intent(in),optional::start_index
      real(dp),allocatable::x(:,:)
      integer::s
      s=1
      if(present(start_index))s=start_index
      x=fourier_terms(h,period,k,s+n_history)
   end function fourier_forecast_terms

   function seasonal_dummy_forecast(n_history,h,m,start_index) result(x)
      integer,intent(in)::n_history,h,m
      integer,intent(in),optional::start_index
      real(dp),allocatable::x(:,:)
      integer::s
      s=1
      if(present(start_index))s=start_index
      x=seasonal_dummy(h,m,s+n_history)
   end function seasonal_dummy_forecast

   function seasonal_index_forecast(seasonal,period,h) result(out)
      real(dp),intent(in)::seasonal(:)
      integer,intent(in)::period,h
      real(dp),allocatable::out(:)
      integer::i,n
      if(period<1 .or. size(seasonal)<period)error stop 'seasonal_index_forecast: invalid period'
      n=size(seasonal)
      allocate(out(h))
      do i=1,h
         out(i)=seasonal(n-period+1+mod(i-1,period))
      end do
   end function seasonal_index_forecast

   function tapered_acf(x,lag_max) result(gamma)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::lag_max
      real(dp),allocatable::gamma(:),raw(:),gmat(:,:),eval(:),evec(:,:),g2(:,:)
      real(dp)::upper,kap,md
      integer::lag,n,l,k,s,i,j,info
      n=size(x)
      lag=max(0,n-1)
      if(present(lag_max))lag=min(max(0,lag_max),n-1)
      raw=acf_values(x,lag,.true.)
      s=lag+1
      allocate(gamma(s))
      gamma=raw
      upper=2.0_dp*sqrt(log10(real(max(2,n),dp))/real(max(1,n),dp))
      l=0
      if(s>=5)then
         do k=1,s-4
            if(all(abs(gamma(k:k+4))<upper))then
               l=k
               exit
            end if
         end do
      end if
      if(l==0)l=s
      do i=1,s
         if(real(i,dp)/real(l,dp)<=1.0_dp)then
            kap=1.0_dp
         else if(real(i,dp)/real(l,dp)<=2.0_dp)then
            kap=2.0_dp-real(i,dp)/real(l,dp)
         else
            kap=0.0_dp
         end if
         gamma(i)=gamma(i)*kap
      end do
      allocate(gmat(s,s))
      do i=1,s
         do j=1,s
            gmat(i,j)=gamma(abs(i-j)+1)
         end do
      end do
      call symmetric_eigen(gmat,eval,evec,info)
      if(info/=0)return
      eval=max(eval,20.0_dp/real(n,dp))
      md=sum(eval)/real(s,dp)
      g2=matmul(evec*spread(eval,1,s),transpose(evec))/max(md,tiny(1.0_dp))
      gamma(1)=sum([(g2(i,i),i=1,s)])/real(s,dp)
      do k=1,lag
         gamma(k+1)=0.0_dp
         do i=1,s-k
            gamma(k+1)=gamma(k+1)+g2(i+k,i)
         end do
         gamma(k+1)=gamma(k+1)/real(s-k,dp)
      end do
   end function tapered_acf

   function tapered_pacf(x,lag_max) result(pacf)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::lag_max
      real(dp),allocatable::pacf(:),a(:),oldphi(:),phi(:)
      real(dp)::pev
      integer::lag,k
      lag=max(0,size(x)-1)
      if(present(lag_max))lag=min(max(0,lag_max),size(x)-1)
      allocate(pacf(lag+1))
      pacf=0.0_dp
      pacf(1)=1.0_dp
      if(lag==0)return
      a=tapered_acf(x,lag)
      allocate(oldphi(1))
      oldphi(1)=a(2)/max(a(1),tiny(1.0_dp))
      pacf(2)=oldphi(1)
      pev=a(1)*(1.0_dp-oldphi(1)**2)
      do k=2,lag
         allocate(phi(k))
         phi=0.0_dp
         phi(k)=(a(k+1)-dot_product(oldphi,a(k:2:-1)))/max(pev,tiny(1.0_dp))
         phi(1:k-1)=oldphi-phi(k)*oldphi(k-1:1:-1)
         pacf(k+1)=phi(k)
         pev=pev*(1.0_dp-phi(k)**2)
         call move_alloc(phi,oldphi)
      end do
   end function tapered_pacf

   function linear_process_bootstrap(x,nsim) result(out)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::nsim
      real(dp),allocatable::out(:,:),gamma(:),gmat(:,:),l(:,:),centered(:),w(:),sampled(:,:)
      real(dp)::mx,u
      integer::n,nb,i,j,k,info
      n=size(x)
      nb=100
      if(present(nsim))nb=max(1,nsim)
      if(n<2)then
         allocate(out(n,nb))
         out=spread(x,2,nb)
         return
      end if
      mx=mean_value(x)
      centered=x-mx
      gamma=tapered_acf(centered,n-1)
      allocate(gmat(n,n))
      do i=1,n
         do j=1,n
            gmat(i,j)=gamma(abs(i-j)+1)
         end do
      end do
      call cholesky_lower(gmat,l,info)
      if(info/=0)error stop 'linear_process_bootstrap: covariance Cholesky failed'
      call solve_linear(l,centered,w,info)
      if(info/=0)error stop 'linear_process_bootstrap: whitening solve failed'
      allocate(sampled(n,nb))
      do j=1,nb
         do i=1,n
            call random_number(u)
            k=1+min(n-1,int(u*real(n,dp)))
            sampled(i,j)=w(k)
         end do
      end do
      out=matmul(l,sampled)+mx
   end function linear_process_bootstrap

   subroutine tapered_correlation_ci(x,lag_max,partial,level,nsim,estimate,lower,upper)
      real(dp),intent(in)::x(:)
      integer,intent(in)::lag_max
      logical,intent(in),optional::partial
      real(dp),intent(in),optional::level
      integer,intent(in),optional::nsim
      real(dp),allocatable,intent(out)::estimate(:),lower(:),upper(:)
      real(dp),allocatable::boots(:,:),vals(:,:),tmp(:)
      real(dp)::lev,prob
      integer::lag,nb,j,k
      logical::part
      lag=min(max(0,lag_max),size(x)-1)
      nb=100
      if(present(nsim))nb=max(2,nsim)
      lev=95.0_dp
      if(present(level))lev=max(0.0_dp,min(100.0_dp,level))
      part=.false.
      if(present(partial))part=partial
      allocate(estimate(lag),lower(lag),upper(lag))
      if(lag==0)return
      if(part)then
         tmp=tapered_pacf(x,lag)
         estimate=tmp(2:lag+1)
      else
         tmp=tapered_acf(x,lag)
         estimate=tmp(2:lag+1)
      end if
      boots=linear_process_bootstrap(x,nb)
      allocate(vals(lag,nb))
      do j=1,nb
         if(part)then
            tmp=tapered_pacf(boots(:,j),lag)
            vals(:,j)=tmp(2:lag+1)
         else
            tmp=tapered_acf(boots(:,j),lag)
            vals(:,j)=tmp(2:lag+1)
         end if
      end do
      prob=(100.0_dp-lev)/200.0_dp
      do k=1,lag
         lower(k)=quantile_type8(vals(k,:),prob)
         upper(k)=quantile_type8(vals(k,:),1.0_dp-prob)
      end do
   end subroutine tapered_correlation_ci

   function fourier_terms(n,period,k,start_index) result(x)
      integer,intent(in)::n,period,k
      integer,intent(in),optional::start_index
      real(dp),allocatable::x(:,:)
      integer::i,j,s,col
      real(dp)::ang
      s=1
      if(present(start_index))s=start_index
      allocate(x(n,2*k))
      x=0.0_dp
      do i=1,n
         do j=1,k
            ang=2.0_dp*pi*real(j*(s+i-1),dp)/real(period,dp)
            col=2*j-1
            x(i,col)=sin(ang)
            x(i,col+1)=cos(ang)
         end do
      end do
   end function
   function seasonal_dummy(n,m,start_index) result(x)
      integer,intent(in)::n,m
      integer,intent(in),optional::start_index
      real(dp),allocatable::x(:,:)
      integer::i,s,season
      s=1
      if(present(start_index))s=start_index
      if(m<=1) then
      allocate(x(n,0))
      return
      end if
      allocate(x(n,m-1))
      x=0.0_dp
      do i=1,n
         season=mod(s+i-2,m)+1
         if(season<m)x(i,season)=1.0_dp
      end do
   end function
   function moving_average(x,order,centered) result(y)
      real(dp),intent(in)::x(:)
      integer,intent(in)::order
      logical,intent(in),optional::centered
      real(dp),allocatable::y(:)
      logical::c
      integer::i,j1,j2,n
      n=size(x)
      allocate(y(n))
      y=0.0_dp
      c=.true.
      if(present(centered))c=centered
      do i=1,n
         if(c)then
         j1=max(1,i-order/2)
         j2=min(n,i+(order-1)/2)
         else
         j1=max(1,i-order+1)
         j2=i
         end if
         y(i)=sum(x(j1:j2))/real(j2-j1+1,dp)
      end do
   end function
   function difference_series(x,d) result(y)
      real(dp),intent(in)::x(:)
      integer,intent(in)::d
      real(dp),allocatable::y(:),z(:)
      integer::j
      z=x
      do j=1,d
      y=z(2:)-z(:size(z)-1)
      z=y
      end do
      y=z
   end function
   function seasonal_difference(x,m,D) result(y)
      real(dp),intent(in)::x(:)
      integer,intent(in)::m,D
      real(dp),allocatable::y(:),z(:)
      integer::j,n
      z=x
      do j=1,D
      n=size(z)
      if(n<=m)then
      allocate(y(0))
      return
      end if
      y=z(1+m:n)-z(1:n-m)
      z=y
      end do
      y=z
   end function
   integer function findfrequency(x,max_period) result(period)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::max_period
      real(dp),allocatable::a(:)
      real(dp)::best
      integer::k,mx
      mx=min(size(x)/2,350)
      if(present(max_period))mx=min(mx,max_period)
      if(mx<2)then
      period=1
      return
      end if
      a=acf_values(x,mx,.true.)
      best=0.0_dp
      period=1
      do k=2,mx
         if(a(k)>best .and. a(k)>0.3_dp)then
            if(k==mx .or. (a(k)>=a(k-1).and.a(k)>=a(k+1)))then
            best=a(k)
            period=k
            end if
         end if
      end do
   end function
end module forecast_features
