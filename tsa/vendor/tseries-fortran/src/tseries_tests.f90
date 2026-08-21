! Experimental modern Fortran translation of computational routines from
! the R package tseries 0.10-62. Original authors include Adrian Trapletti
! and Kurt Hornik; Blake LeBaron contributed the original BDS code.
! Licensed under GPL-2.0-only OR GPL-3.0-only. See LICENSE and NOTICE.

module tseries_tests
   use tseries_kinds, only : dp
   use tseries_types, only : test_result, bds_result
   use tseries_stats, only : mean_value, variance_value, standard_deviation, long_run_variance, &
      cumulative_sum, linear_interpolate
   use tseries_special, only : normal_cdf, chi_square_cdf, f_cdf
   use tseries_linalg, only : least_squares, standardize_columns, covariance_matrix, jacobi_eigen
   use tseries_random, only : seed_random
   implicit none
   private

   public :: runs_test
   public :: jarque_bera_test
   public :: adf_test
   public :: pp_test
   public :: kpss_test
   public :: po_test
   public :: bds_test
   public :: terasvirta_test
   public :: white_test

contains

   function runs_test(x,alternative) result(res)
      integer, intent(in) :: x(:)
      character(len=*), intent(in), optional :: alternative
      type(test_result) :: res
      character(len=:), allocatable :: alt
      integer :: n,n1,n2,r,i
      real(dp) :: mu,sd,z
      alt='two.sided'; if(present(alternative)) alt=trim(alternative)
      n=size(x)
      if(n<2 .or. any((x/=0).and.(x/=1))) then
         res%status=1; res%message='runs_test requires at least two binary observations'; return
      end if
      n1=count(x==0); n2=count(x==1)
      if(n1==0 .or. n2==0) then
         res%status=1; res%message='both binary levels must occur'; return
      end if
      r=1
      do i=2,n
         if(x(i)/=x(i-1)) r=r+1
      end do
      mu=1.0_dp+2.0_dp*real(n1*n2,dp)/real(n1+n2,dp)
      sd=sqrt(2.0_dp*real(n1*n2,dp)*(2.0_dp*real(n1*n2,dp)-real(n1+n2,dp))/ &
         (real(n1+n2,dp)**2*real(n1+n2-1,dp)))
      z=(real(r,dp)-mu)/sd
      res%statistic=z
      select case(alt)
      case('less')
         res%p_value=normal_cdf(z)
      case('greater')
         res%p_value=1.0_dp-normal_cdf(z)
      case default
         res%p_value=2.0_dp*normal_cdf(-abs(z))
      end select
      res%method='Runs Test'; res%status=0; res%message='ok'
   end function runs_test

   function jarque_bera_test(x) result(res)
      real(dp), intent(in) :: x(:)
      type(test_result) :: res
      real(dp) :: mu,m2,m3,m4,b1,b2
      integer :: n
      n=size(x)
      if(n<3) then
         res%status=1; res%message='at least three observations are required'; return
      end if
      mu=mean_value(x)
      m2=sum((x-mu)**2)/real(n,dp)
      if(m2<=tiny(1.0_dp)) then
         res%status=1; res%message='zero variance'; return
      end if
      m3=sum((x-mu)**3)/real(n,dp)
      m4=sum((x-mu)**4)/real(n,dp)
      b1=(m3/m2**1.5_dp)**2
      b2=m4/m2**2
      res%statistic=real(n,dp)*b1/6.0_dp+real(n,dp)*(b2-3.0_dp)**2/24.0_dp
      res%p_value=1.0_dp-chi_square_cdf(res%statistic,2.0_dp)
      allocate(res%parameters(1)); res%parameters=2.0_dp
      res%method='Jarque Bera Test'; res%status=0; res%message='ok'
   end function jarque_bera_test

   function adf_test(x,lags,explosive) result(res)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: lags
      logical, intent(in), optional :: explosive
      type(test_result) :: res
      real(dp), allocatable :: dy(:),design(:,:),y(:),beta(:),resid(:),cov(:,:),crit(:)
      real(dp), parameter :: table_t(6)=[25.0_dp,50.0_dp,100.0_dp,250.0_dp,500.0_dp,100000.0_dp]
      real(dp), parameter :: table_p(8)=[0.01_dp,0.025_dp,0.05_dp,0.10_dp,0.90_dp,0.95_dp,0.975_dp,0.99_dp]
      real(dp), parameter :: table(6,8)=reshape(-[ &
         4.38_dp,4.15_dp,4.04_dp,3.99_dp,3.98_dp,3.96_dp, &
         3.95_dp,3.80_dp,3.73_dp,3.69_dp,3.68_dp,3.66_dp, &
         3.60_dp,3.50_dp,3.45_dp,3.43_dp,3.42_dp,3.41_dp, &
         3.24_dp,3.18_dp,3.15_dp,3.13_dp,3.13_dp,3.12_dp, &
         1.14_dp,1.19_dp,1.22_dp,1.23_dp,1.24_dp,1.25_dp, &
         0.80_dp,0.87_dp,0.90_dp,0.92_dp,0.93_dp,0.94_dp, &
         0.50_dp,0.58_dp,0.62_dp,0.64_dp,0.65_dp,0.66_dp, &
         0.15_dp,0.24_dp,0.28_dp,0.31_dp,0.32_dp,0.33_dp ],[6,8])
      integer :: k,n0,nobs,i,j,status
      logical :: use_explosive

      n0=size(x); k=int(real(max(0,n0-1),dp)**(1.0_dp/3.0_dp)); if(present(lags)) k=lags
      use_explosive=.false.; if(present(explosive)) use_explosive=explosive
      if(k<0 .or. n0<=k+4) then
         res%status=1; res%message='invalid ADF lag or sample size'; return
      end if
      allocate(dy(n0-1)); dy=x(2:)-x(:n0-1)
      nobs=n0-1-k
      allocate(design(nobs,3+k),y(nobs),beta(3+k),resid(nobs),cov(3+k,3+k),crit(8))
      do i=1,nobs
         y(i)=dy(k+i)
         design(i,1)=1.0_dp
         design(i,2)=x(k+i)
         design(i,3)=real(k+i,dp)
         do j=1,k
            design(i,3+j)=dy(k+i-j)
         end do
      end do
      call least_squares(design,y,beta,resid,cov,status)
      if(status/=0 .or. cov(2,2)<=0.0_dp) then
         res%status=2; res%message='singular ADF regression'; return
      end if
      res%statistic=beta(2)/sqrt(cov(2,2))
      do j=1,8
         crit(j)=linear_interpolate(table_t,table(:,j),real(nobs,dp))
      end do
      res%p_value=linear_interpolate(crit,table_p,res%statistic)
      if(use_explosive) res%p_value=1.0_dp-res%p_value
      allocate(res%parameters(1)); res%parameters=real(k,dp)
      res%method='Augmented Dickey-Fuller Test'; res%status=0; res%message='ok'
   end function adf_test

   function pp_test(x,use_t_statistic,long_lag,explosive) result(res)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: use_t_statistic,long_lag,explosive
      type(test_result) :: res
      real(dp), allocatable :: design(:,:),yt(:),beta(:),u(:),cov(:,:),crit(:),yt1(:)
      real(dp), parameter :: table_t(6)=[25.0_dp,50.0_dp,100.0_dp,250.0_dp,500.0_dp,100000.0_dp]
      real(dp), parameter :: table_p(8)=[0.01_dp,0.025_dp,0.05_dp,0.10_dp,0.90_dp,0.95_dp,0.975_dp,0.99_dp]
      real(dp), parameter :: table_zt(6,8)=reshape(-[ &
         4.38_dp,4.15_dp,4.04_dp,3.99_dp,3.98_dp,3.96_dp, &
         3.95_dp,3.80_dp,3.73_dp,3.69_dp,3.68_dp,3.66_dp, &
         3.60_dp,3.50_dp,3.45_dp,3.43_dp,3.42_dp,3.41_dp, &
         3.24_dp,3.18_dp,3.15_dp,3.13_dp,3.13_dp,3.12_dp, &
         1.14_dp,1.19_dp,1.22_dp,1.23_dp,1.24_dp,1.25_dp, &
         0.80_dp,0.87_dp,0.90_dp,0.92_dp,0.93_dp,0.94_dp, &
         0.50_dp,0.58_dp,0.62_dp,0.64_dp,0.65_dp,0.66_dp, &
         0.15_dp,0.24_dp,0.28_dp,0.31_dp,0.32_dp,0.33_dp ],[6,8])
      real(dp), parameter :: table_za(6,8)=reshape(-[ &
         22.5_dp,25.7_dp,27.4_dp,28.4_dp,28.9_dp,29.5_dp, &
         19.9_dp,22.4_dp,23.6_dp,24.4_dp,24.8_dp,25.1_dp, &
         17.9_dp,19.8_dp,20.7_dp,21.3_dp,21.5_dp,21.8_dp, &
         15.6_dp,16.8_dp,17.5_dp,18.0_dp,18.1_dp,18.3_dp, &
         3.66_dp,3.71_dp,3.74_dp,3.75_dp,3.76_dp,3.77_dp, &
         2.51_dp,2.60_dp,2.62_dp,2.64_dp,2.65_dp,2.66_dp, &
         1.53_dp,1.66_dp,1.73_dp,1.78_dp,1.78_dp,1.79_dp, &
         0.43_dp,0.65_dp,0.75_dp,0.82_dp,0.84_dp,0.87_dp ],[6,8])
      integer :: n,n0,i,j,l,status
      real(dp) :: ssqru,ssqrtl,n2,trm1,trm2,trm3,trm4,dx,alpha,tstat
      logical :: ttype,llag,expalt

      n0=size(x); n=n0-1; ttype=.false.; llag=.false.; expalt=.false.
      if(present(use_t_statistic)) ttype=use_t_statistic
      if(present(long_lag)) llag=long_lag
      if(present(explosive)) expalt=explosive
      if(n<6) then
         res%status=1; res%message='sample too short for PP test'; return
      end if
      allocate(design(n,3),yt(n),yt1(n),beta(3),u(n),cov(3,3),crit(8))
      yt=x(2:); yt1=x(:n0-1)
      do i=1,n
         design(i,:)=[1.0_dp,real(i,dp)-real(n,dp)/2.0_dp,yt1(i)]
      end do
      call least_squares(design,yt,beta,u,cov,status)
      if(status/=0) then
         res%status=2; res%message='singular PP regression'; return
      end if
      ssqru=sum(u*u)/real(n,dp)
      if(llag) then
         l=int(12.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      else
         l=int(4.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      end if
      ssqrtl=long_run_variance(u,l)
      n2=real(n,dp)**2
      trm1=n2*(n2-1.0_dp)*sum(yt1**2)/12.0_dp
      trm2=real(n,dp)*sum(yt1*[(real(i,dp),i=1,n)])**2
      trm3=real(n*(n+1),dp)*sum(yt1*[(real(i,dp),i=1,n)])*sum(yt1)
      trm4=real(n*(n+1)*(2*n+1),dp)*sum(yt1)**2/6.0_dp
      dx=trm1-trm2+trm3-trm4
      alpha=beta(3)
      if(ttype) then
         tstat=(alpha-1.0_dp)/sqrt(max(cov(3,3),tiny(1.0_dp)))
         res%statistic=sqrt(ssqru/ssqrtl)*tstat-real(n,dp)**3*(ssqrtl-ssqru)/ &
            (4.0_dp*sqrt(3.0_dp*dx*ssqrtl))
         do j=1,8
            crit(j)=linear_interpolate(table_t,table_zt(:,j),real(n,dp))
         end do
      else
         res%statistic=real(n,dp)*(alpha-1.0_dp)-real(n,dp)**6*(ssqrtl-ssqru)/(24.0_dp*dx)
         do j=1,8
            crit(j)=linear_interpolate(table_t,table_za(:,j),real(n,dp))
         end do
      end if
      res%p_value=linear_interpolate(crit,table_p,res%statistic)
      if(expalt) res%p_value=1.0_dp-res%p_value
      allocate(res%parameters(1)); res%parameters=real(l,dp)
      res%method='Phillips-Perron Unit Root Test'; res%status=0; res%message='ok'
   end function pp_test

   function kpss_test(x,trend,long_lag) result(res)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: trend,long_lag
      type(test_result) :: res
      real(dp), allocatable :: design(:,:),beta(:),e(:),s(:)
      real(dp), parameter :: ptab(4)=[0.01_dp,0.025_dp,0.05_dp,0.10_dp]
      real(dp), parameter :: level_tab(4)=[0.739_dp,0.574_dp,0.463_dp,0.347_dp]
      real(dp), parameter :: trend_tab(4)=[0.216_dp,0.176_dp,0.146_dp,0.119_dp]
      integer :: n,i,l,status,p
      real(dp) :: eta,s2
      logical :: use_trend,llag

      n=size(x); use_trend=.false.; llag=.false.
      if(present(trend)) use_trend=trend
      if(present(long_lag)) llag=long_lag
      p=merge(2,1,use_trend)
      allocate(design(n,p),beta(p),e(n),s(n))
      design(:,1)=1.0_dp
      if(use_trend) design(:,2)=[(real(i,dp),i=1,n)]
      call least_squares(design,x,beta,e,status=status)
      if(status/=0) then
         res%status=2; res%message='singular KPSS regression'; return
      end if
      call cumulative_sum(e,s)
      eta=sum(s*s)/real(n*n,dp)
      if(llag) then
         l=int(12.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      else
         l=int(4.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      end if
      s2=long_run_variance(e,l)
      res%statistic=eta/s2
      if(use_trend) then
         res%p_value=linear_interpolate(trend_tab,ptab,res%statistic)
         res%method='KPSS Test for Trend Stationarity'
      else
         res%p_value=linear_interpolate(level_tab,ptab,res%statistic)
         res%method='KPSS Test for Level Stationarity'
      end if
      allocate(res%parameters(1)); res%parameters=real(l,dp)
      res%status=0; res%message='ok'
   end function kpss_test

   function po_test(x,demean,long_lag) result(res)
      real(dp), intent(in) :: x(:, :)
      logical, intent(in), optional :: demean,long_lag
      type(test_result) :: res
      real(dp), allocatable :: design(:,:),beta(:),u(:),ut(:),ut1(:),arx(:,:),arb(:),k(:)
      real(dp), parameter :: ptab(7)=[0.01_dp,0.025_dp,0.05_dp,0.075_dp,0.10_dp,0.125_dp,0.15_dp]
      real(dp), parameter :: demean_tab(5,7)=reshape(-[ &
         28.32_dp,34.17_dp,41.13_dp,47.51_dp,52.17_dp, &
         23.81_dp,29.74_dp,35.71_dp,41.64_dp,46.53_dp, &
         20.49_dp,26.09_dp,32.06_dp,37.15_dp,41.94_dp, &
         18.48_dp,23.87_dp,29.51_dp,34.71_dp,39.11_dp, &
         17.04_dp,22.19_dp,27.58_dp,32.74_dp,37.01_dp, &
         15.93_dp,21.04_dp,26.23_dp,31.15_dp,35.48_dp, &
         14.91_dp,19.95_dp,25.05_dp,29.88_dp,34.20_dp ],[5,7])
      real(dp), parameter :: nodemean_tab(5,7)=reshape(-[ &
         22.83_dp,29.27_dp,36.16_dp,42.87_dp,48.52_dp, &
         18.89_dp,25.21_dp,31.54_dp,37.48_dp,42.55_dp, &
         15.64_dp,21.48_dp,27.85_dp,33.48_dp,38.09_dp, &
         13.81_dp,19.61_dp,25.52_dp,30.93_dp,35.51_dp, &
         12.54_dp,18.18_dp,23.92_dp,28.85_dp,33.80_dp, &
         11.57_dp,17.01_dp,22.62_dp,27.40_dp,32.27_dp, &
         10.74_dp,16.02_dp,21.53_dp,26.17_dp,30.90_dp ],[5,7])
      integer :: n,d,p,l,status
      real(dp) :: ssqrk,ssqrtl,alpha
      logical :: dm,llag

      n=size(x,1); d=size(x,2); dm=.true.; llag=.false.
      if(present(demean)) dm=demean
      if(present(long_lag)) llag=long_lag
      if(d<2 .or. d>6 .or. n<5) then
         res%status=1; res%message='PO test requires 2 to 6 series'; return
      end if
      p=d-1+merge(1,0,dm)
      allocate(design(n,p),beta(p),u(n))
      if(dm) then
         design(:,1)=1.0_dp; design(:,2:)=x(:,2:)
      else
         design=x(:,2:)
      end if
      call least_squares(design,x(:,1),beta,u,status=status)
      if(status/=0) then
         res%status=2; res%message='singular cointegrating regression'; return
      end if
      allocate(ut(n-1),ut1(n-1),arx(n-1,1),arb(1),k(n-1))
      ut=u(2:); ut1=u(:n-1); arx(:,1)=ut1
      call least_squares(arx,ut,arb,k,status=status)
      if(status/=0) then
         res%status=2; res%message='singular residual AR regression'; return
      end if
      ssqrk=sum(k*k)/real(n-1,dp)
      if(llag) then
         l=int(real(n-1,dp)/30.0_dp)
      else
         l=int(real(n-1,dp)/100.0_dp)
      end if
      ssqrtl=long_run_variance(k,l)
      alpha=arb(1)
      res%statistic=real(n-1,dp)*(alpha-1.0_dp)-0.5_dp*real(n-1,dp)**2*(ssqrtl-ssqrk)/sum(ut1**2)
      if(dm) then
         res%p_value=linear_interpolate(-demean_tab(d-1,:),ptab,res%statistic)
      else
         res%p_value=linear_interpolate(-nodemean_tab(d-1,:),ptab,res%statistic)
      end if
      allocate(res%parameters(1)); res%parameters=real(l,dp)
      res%method='Phillips-Ouliaris Cointegration Test'; res%status=0; res%message='ok'
   end function po_test

   function bds_test(x,max_embedding,eps) result(res)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: max_embedding
      real(dp), intent(in), optional :: eps(:)
      type(bds_result) :: res
      real(dp), allocatable :: epsilon_values(:)
      real(dp) :: c1,kvalue,cm,sigma,std
      integer :: n,mmax,ne,i,m,j,nobs,count_close,tcount,phi_count,a,b,d

      n=size(x); mmax=3; if(present(max_embedding)) mmax=max_embedding
      if(present(eps)) then
         allocate(epsilon_values(size(eps))); epsilon_values=eps
      else
         allocate(epsilon_values(4))
         epsilon_values=[0.5_dp,1.0_dp,1.5_dp,2.0_dp]*standard_deviation(x)
      end if
      ne=size(epsilon_values)
      if(n<=mmax+2 .or. mmax<2 .or. any(epsilon_values<=0.0_dp)) then
         res%status=1; res%message='invalid BDS dimensions or epsilon'; return
      end if
      allocate(res%statistic(mmax-1,ne),res%p_value(mmax-1,ne),res%eps(ne))
      res%eps=epsilon_values; res%max_embedding=mmax
      do i=1,ne
         nobs=n-mmax+1
         count_close=0; phi_count=0
         do a=1,nobs
            tcount=0
            do b=1,nobs
               if(abs(x(a)-x(b))<=epsilon_values(i)) tcount=tcount+1
            end do
            phi_count=phi_count+tcount*tcount
            count_close=count_close+tcount
         end do
         count_close=count_close-nobs
         phi_count=phi_count-nobs-3*count_close
         c1=real(count_close,dp)/real(nobs*(nobs-1),dp)
         kvalue=real(phi_count,dp)/real(nobs*(nobs-1)*(nobs-2),dp)
         do m=2,mmax
            count_close=0
            do a=1,nobs-1
               do b=a+1,nobs
                  d=0
                  do j=0,m-1
                     if(abs(x(a+j)-x(b+j))>epsilon_values(i)) then
                        d=1; exit
                     end if
                  end do
                  if(d==0) count_close=count_close+1
               end do
            end do
            cm=2.0_dp*real(count_close,dp)/real(nobs*(nobs-1),dp)
            sigma=0.0_dp
            do j=1,m-1
               sigma=sigma+2.0_dp*kvalue**(m-j)*c1**(2*j)
            end do
            sigma=4.0_dp*(sigma+kvalue**m+real((m-1)*(m-1),dp)*c1**(2*m)-real(m*m,dp)*kvalue*c1**(2*m-2))
            std=sqrt(max(sigma/real(nobs,dp),tiny(1.0_dp)))
            res%statistic(m-1,i)=(cm-c1**m)/std
            res%p_value(m-1,i)=2.0_dp*normal_cdf(-abs(res%statistic(m-1,i)))
         end do
      end do
      res%status=0; res%message='ok'
   end function bds_test

   function terasvirta_test(x,lag,use_f,scale_data) result(res)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: lag
      logical, intent(in), optional :: use_f,scale_data
      type(test_result) :: res
      real(dp), allocatable :: lags(:,:),y(:),zlags(:,:),zy(:,:),base(:,:),b(:),u(:),aug(:,:),b2(:),v(:)
      integer :: n,l,nobs,i,j,k,m,col,status,nin
      real(dp) :: ssr0,ssr
      logical :: ftype,scale
      l=1; if(present(lag)) l=lag
      ftype=.false.; if(present(use_f)) ftype=use_f
      scale=.true.; if(present(scale_data)) scale=scale_data
      n=size(x); nobs=n-l; nin=l
      if(l<1 .or. nobs<=l+4) then
         res%status=1; res%message='invalid Terasvirta lag/sample'; return
      end if
      allocate(lags(nobs,l),y(nobs))
      do i=1,nobs
         y(i)=x(l+i)
         do j=1,l
            lags(i,j)=x(l+i-j)
         end do
      end do
      if(scale) then
         allocate(zlags(nobs,l),zy(nobs,1))
         call standardize_columns(lags,zlags)
         call standardize_columns(reshape(y,[nobs,1]),zy)
         lags=zlags; y=zy(:,1)
      end if
      allocate(base(nobs,1+l),b(1+l),u(nobs))
      base(:,1)=1.0_dp; base(:,2:)=lags
      call least_squares(base,y,b,u,status=status)
      if(status/=0) then
         res%status=2; res%message='singular base regression'; return
      end if
      ssr0=sum(u*u)
      m=l*(l+1)/2+l*(l+1)*(l+2)/6
      allocate(aug(nobs,1+l+m),b2(1+l+m),v(nobs))
      aug(:,1:1+l)=base; col=1+l
      do i=1,l
         do j=i,l
            col=col+1; aug(:,col)=lags(:,i)*lags(:,j)
         end do
      end do
      do i=1,l
         do j=i,l
            do k=j,l
               col=col+1; aug(:,col)=lags(:,i)*lags(:,j)*lags(:,k)
            end do
         end do
      end do
      call least_squares(aug,u,b2,v,status=status)
      if(status/=0) then
         res%status=2; res%message='singular auxiliary regression'; return
      end if
      ssr=sum(v*v)
      if(ftype) then
         res%statistic=((ssr0-ssr)/real(m,dp))/(ssr/real(nobs-l-m,dp))
         res%p_value=1.0_dp-f_cdf(res%statistic,real(m,dp),real(nobs-l-m,dp))
         allocate(res%parameters(2)); res%parameters=[real(m,dp),real(nobs-l-m,dp)]
      else
         res%statistic=real(n,dp)*log(ssr0/ssr)
         res%p_value=1.0_dp-chi_square_cdf(res%statistic,real(m,dp))
         allocate(res%parameters(1)); res%parameters=real(m,dp)
      end if
      res%method='Terasvirta Neural Network Test'; res%status=0; res%message='ok'
   end function terasvirta_test

   function white_test(x,lag,qstar,q,weight_range,seed,use_f,scale_data) result(res)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: lag,qstar,q,seed
      real(dp), intent(in), optional :: weight_range
      logical, intent(in), optional :: use_f,scale_data
      type(test_result) :: res
      real(dp), allocatable :: lags(:,:),y(:),zlags(:,:),zy(:,:),base(:,:),b(:),u(:),gamma(:,:),phantom(:,:), &
         zphantom(:,:),covp(:,:),eval(:),evec(:,:),scores(:,:),aug(:,:),b2(:),v(:)
      real(dp) :: ssr0,ssr,wr,u01
      integer :: n,l,qs,nh,nobs,i,j,status,seedv
      logical :: ftype,scale
      l=1; qs=2; nh=10; wr=4.0_dp; seedv=12345
      if(present(lag)) l=lag
      if(present(qstar)) qs=qstar
      if(present(q)) nh=q
      if(present(weight_range)) wr=weight_range
      if(present(seed)) seedv=seed
      ftype=.false.; if(present(use_f)) ftype=use_f
      scale=.true.; if(present(scale_data)) scale=scale_data
      n=size(x); nobs=n-l
      if(l<1 .or. qs<1 .or. nh<=qs .or. nobs<=l+qs+3) then
         res%status=1; res%message='invalid White test dimensions'; return
      end if
      allocate(lags(nobs,l),y(nobs))
      do i=1,nobs
         y(i)=x(l+i)
         do j=1,l
            lags(i,j)=x(l+i-j)
         end do
      end do
      if(scale) then
         allocate(zlags(nobs,l),zy(nobs,1))
         call standardize_columns(lags,zlags)
         call standardize_columns(reshape(y,[nobs,1]),zy)
         lags=zlags; y=zy(:,1)
      end if
      allocate(base(nobs,1+l),b(1+l),u(nobs))
      base(:,1)=1.0_dp; base(:,2:)=lags
      call least_squares(base,y,b,u,status=status)
      if(status/=0) then
         res%status=2; res%message='singular base regression'; return
      end if
      ssr0=sum(u*u)
      allocate(gamma(l+1,nh),phantom(nobs,nh),zphantom(nobs,nh),covp(nh,nh),eval(nh),evec(nh,nh),scores(nobs,qs))
      call seed_random(seedv)
      do j=1,nh
         do i=1,l+1
            call random_number(u01); gamma(i,j)=wr*(u01-0.5_dp)
         end do
      end do
      phantom=1.0_dp/(1.0_dp+exp(-matmul(base,gamma)))
      call standardize_columns(phantom,zphantom)
      call covariance_matrix(zphantom,covp)
      call jacobi_eigen(covp,eval,evec,status)
      if(status/=0) then
         res%status=2; res%message='PCA did not converge'; return
      end if
      scores=matmul(zphantom,evec(:,2:qs+1))
      allocate(aug(nobs,1+l+qs),b2(1+l+qs),v(nobs))
      aug(:,1:1+l)=base; aug(:,2+l:)=scores
      call least_squares(aug,u,b2,v,status=status)
      if(status/=0) then
         res%status=2; res%message='singular auxiliary regression'; return
      end if
      ssr=sum(v*v)
      if(ftype) then
         res%statistic=((ssr0-ssr)/real(qs,dp))/(ssr/real(nobs-qs-l,dp))
         res%p_value=1.0_dp-f_cdf(res%statistic,real(qs,dp),real(nobs-qs-l,dp))
         allocate(res%parameters(2)); res%parameters=[real(qs,dp),real(nobs-qs-l,dp)]
      else
         res%statistic=real(n,dp)*log(ssr0/ssr)
         res%p_value=1.0_dp-chi_square_cdf(res%statistic,real(qs,dp))
         allocate(res%parameters(1)); res%parameters=real(qs,dp)
      end if
      res%method='White Neural Network Test'; res%status=0; res%message='ok'
   end function white_test

end module tseries_tests
