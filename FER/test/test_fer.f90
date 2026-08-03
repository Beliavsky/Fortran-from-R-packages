program test_fer
   use fer
   implicit none
   integer :: failures
   failures = 0
   call test_vanilla(failures)
   call test_cev(failures)
   call test_sabr(failures)
   call test_nsvh(failures)
   call test_spreads(failures)
   if (failures /= 0) error stop 'FER tests failed'
   print '(a)', 'All FER tests passed.'
contains
   subroutine check_close(name,x,y,tol,failures)
      character(*), intent(in) :: name
      real(dp), intent(in) :: x,y,tol
      integer, intent(inout) :: failures
      if (abs(x-y) > tol*max(1.0_dp,abs(x),abs(y))) then
         failures=failures+1
         print '(a,2es24.14)', 'FAIL '//trim(name)//': ',x,y
      end if
   end subroutine
   subroutine test_vanilla(f)
      integer,intent(inout)::f
      integer :: i, cp
      real(dp)::k,s,t,r,q,fw,df,p,iv,p2
      do i=1,200
         k=50.0_dp+100.0_dp*real(mod(37*i,199),dp)/198.0_dp
         s=0.001_dp+1.999_dp*real(mod(53*i,197),dp)/196.0_dp
         t=1.0e-4_dp+9.9999_dp*real(mod(71*i,193),dp)/192.0_dp
         r=0.1_dp*real(mod(29*i,191),dp)/190.0_dp
         q=0.1_dp*real(mod(43*i,189),dp)/188.0_dp
         cp=merge(1,-1,mod(i,2)==0); df=exp(-r*t); fw=100.0_dp*exp((r-q)*t)
         p=black_scholes_price(k,fw,t,s,cp,df); iv=black_scholes_impvol(p,k,fw,t,cp,df)
         p2=black_scholes_price(k,fw,t,iv,cp,df); call check_close('BS roundtrip',p,p2,1e-10_dp,f)
         s=1.0_dp+19.0_dp*real(mod(47*i,181),dp)/180.0_dp
         p=bachelier_price(k,fw,t,s,cp,df); iv=bachelier_impvol(p,k,fw,t,cp,df)
         p2=bachelier_price(k,fw,t,iv,cp,df); call check_close('Bachelier roundtrip',p,p2,1e-10_dp,f)
      end do
   end subroutine
   subroutine test_cev(f)
      integer,intent(inout)::f
      real(dp),parameter::ref(5)=[0.04608_dp,0.04229_dp,0.03868_dp,0.03525_dp,0.03203_dp]
      integer::i
      real(dp)::v
      do i=1,5
         v=cev_price(0.4_dp*real(i,dp)*0.05_dp,0.05_dp,1.0_dp,0.4_dp,0.3_dp,1,1.0_dp)
         call check_close('CEV reference',anint(v*1e5_dp)/1e5_dp,ref(i),1e-6_dp,f)
      end do
      do i=1,10
         v=1000.0_dp*cev_price(0.001_dp,1.0_dp,real(i,dp),0.5_dp,0.2_dp,-1,1.0_dp)
         call check_close('CEV mass zero',v,cev_mass_zero(1.0_dp,real(i,dp),0.5_dp,0.2_dp),1e-5_dp,f)
      end do
   end subroutine
   subroutine test_sabr(f)
      integer,intent(inout)::f
      real(dp), parameter :: r1(20) = [ &
         .7176_dp, .5725_dp, .4886_dp, .4293_dp, .3835_dp, &
         .3462_dp, .3148_dp, .2876_dp, .2638_dp, .2427_dp, &
         .2238_dp, .2068_dp, .1916_dp, .1781_dp, .1663_dp, &
         .1562_dp, .1478_dp, .1412_dp, .1360_dp, .1322_dp ]
      integer::i; real(dp)::v
      do i=1,20
         v=sabr_hagan_2002(0.1_dp*i,1.0_dp,10.0_dp,.25_dp,.3_dp,-.8_dp,.3_dp)
         call check_close('SABR table',anint(v*1e4_dp)/1e4_dp,r1(i),1e-5_dp,f)
      end do
   end subroutine
   subroutine test_nsvh(f)
      integer,intent(inout)::f
      integer::i; real(dp)::k,fw,df,c,p,target
      df=exp(-.1_dp*2.3_dp); fw=120.0_dp*exp((.1_dp-.05_dp)*2.3_dp)
      do i=-3,3
         k=120.0_dp+10.0_dp*i
         c=nsvh1_choi_2019(k,fw,2.3_dp,20.0_dp,.2_dp,-.5_dp,1,df)
         p=nsvh1_choi_2019(k,fw,2.3_dp,20.0_dp,.2_dp,-.5_dp,-1,df)
         target=exp(-.05_dp*2.3_dp)*120.0_dp-df*k
         call check_close('NSVh parity',c-p,target,1e-12_dp,f)
      end do
   end subroutine
   subroutine test_spreads(f)
      integer,intent(inout)::f
      integer::i; real(dp)::k,f1,f2,df,c,p,target
      df=exp(-.1_dp*1.3_dp); f1=110.0_dp*exp((.1_dp-.05_dp)*1.3_dp); f2=100.0_dp*exp((.1_dp-.07_dp)*1.3_dp)
      do i=-3,3
         k=10.0_dp+5.0_dp*i; target=df*(f1-f2-k)
         c=spread_kirk(k,f1,f2,1.3_dp,.2_dp,.3_dp,-.5_dp,1,df); p=spread_kirk(k,f1,f2,1.3_dp,.2_dp,.3_dp,-.5_dp,-1,df)
         call check_close('Kirk parity',c-p,target,1e-11_dp,f)
         c = spread_bjerksund_2014(k,f1,f2,1.3_dp,.2_dp,.3_dp,-.5_dp,1,df)
         p = spread_bjerksund_2014(k,f1,f2,1.3_dp,.2_dp,.3_dp,-.5_dp,-1,df)
         call check_close('Bjerksund parity',c-p,target,1e-10_dp,f)
         c = spread_bachelier(k,f1,f2,1.3_dp,22.0_dp,30.0_dp,-.5_dp,1,df)
         p = spread_bachelier(k,f1,f2,1.3_dp,22.0_dp,30.0_dp,-.5_dp,-1,df)
         call check_close('Spread Bachelier parity',c-p,target,1e-11_dp,f)
      end do
   end subroutine
end program test_fer
