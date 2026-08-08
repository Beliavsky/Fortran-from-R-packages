! SPDX-License-Identifier: GPL-3.0-only
program test_core
   use nmof
   implicit none
   integer :: failures,status,i
   real(dp) :: ma(4),x(4),d(5),cov(3,3),sample_cov(3,3),w(3),ret(200,3),meanv(3),sdev(3),repaired(3,3),mrc1(3),mrc2(3)
   real(dp),allocatable :: paths(:,:),sample(:,:),intervals(:,:),resampled(:,:)
   real(dp) :: nodes4(4),weights4(4),integral
   type(quadrature_rule) :: rule
   type(cppi_result) :: cp
   type(drawdown_summary) :: di
   type(pbo_result) :: pbo
   real(dp) :: m(12,3)

   failures=0
   x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   ma=moving_average(x,2)
   call check_array('moving average',ma,[0.5_dp,1.5_dp,2.5_dp,3.5_dp],1.0e-14_dp)
   call check_close('partial moment',partial_moment([-2.0_dp,-1.0_dp,1.0_dp,2.0_dp],xp=1.0_dp),0.75_dp,1.0e-14_dp)
   d=drawdown_series([100.0_dp,110.0_dp,90.0_dp,95.0_dp,120.0_dp])
   call check_close('drawdown',maxval(d),20.0_dp/110.0_dp,1.0e-14_dp)
   di=drawdown_info([100.0_dp,110.0_dp,90.0_dp,95.0_dp,120.0_dp])
   call check_true('drawdown positions',di%high_position==2.and.di%low_position==3)

   call check_close('Ackley zero',test_ackley([0.0_dp,0.0_dp]),0.0_dp,1.0e-14_dp)
   call check_close('Griewank zero',test_griewank([0.0_dp,0.0_dp]),0.0_dp,1.0e-14_dp)
   call check_close('Rosenbrock one',test_rosenbrock([1.0_dp,1.0_dp,1.0_dp]),0.0_dp,1.0e-14_dp)
   call check_close('Trefethen',test_trefethen([-0.0244_dp,0.2106_dp]),-3.306868_dp,2.0e-5_dp)

   rule=xw_gauss(4,'legendre'); nodes4=rule%nodes; weights4=rule%weights
   integral=sum(weights4*nodes4**6)
   call check_close('Gauss-Legendre',integral,2.0_dp/7.0_dp,2.0e-14_dp)
   intervals=bracketing(cubic_root,-2.0_dp,2.0_dp,41)
   call check_true('bracketing count',size(intervals,1)>=1)

   repaired=repair_matrix(reshape([1.0_dp,1.2_dp,0.2_dp,1.2_dp,1.0_dp,0.3_dp,0.2_dp,0.3_dp,1.0_dp],[3,3]),eps=1.0e-8_dp,status=status)
   call check_true('repair status',status==nmof_ok)
   call check_array('repair diagonal',[repaired(1,1),repaired(2,2),repaired(3,3)],[1.0_dp,1.0_dp,1.0_dp],2.0e-12_dp)

   cov=reshape([1.0_dp,0.2_dp,0.1_dp,0.2_dp,2.0_dp,0.3_dp,0.1_dp,0.3_dp,1.5_dp],[3,3])
   call minimum_variance(cov,w,status=status)
   call check_true('minvar status',status==nmof_ok)
   call check_close('minvar budget',sum(w),1.0_dp,2.0e-9_dp)
   call check_true('minvar bounds',all(w>=-1.0e-10_dp).and.all(w<=1.0_dp+1.0e-10_dp))

   call equal_risk_contribution(cov,w,status=status)
   call check_true('ERC status',status==nmof_ok)
   call check_close('ERC budget',sum(w),1.0_dp,2.0e-10_dp)
   call check_true('ERC equal RC',maxval(abs(w*matmul(cov,w)/dot_product(w,matmul(cov,w))-1.0_dp/3.0_dp))<2.0e-6_dp)

   ret=random_returns(3,200,[0.01_dp,0.02_dp,0.03_dp],mean=[0.001_dp,0.002_dp,0.003_dp],rho=0.25_dp,exact=.true.,seed=123_i8,status=status)
   meanv=column_means(ret); sdev=column_sds(ret); sample_cov=covariance_matrix(ret)
   call check_array('random means',meanv,[0.001_dp,0.002_dp,0.003_dp],2.0e-12_dp)
   call check_array('random sds',sdev,[0.01_dp,0.02_dp,0.03_dp],2.0e-12_dp)
   call check_close('random corr12',sample_cov(1,2)/(sdev(1)*sdev(2)),0.25_dp,2.0e-10_dp)
   w=[0.2_dp,0.3_dp,0.5_dp]
   mrc1=marginal_risk_contributions(w,sample_cov)
   mrc2=marginal_risk_contributions_fd(w,ret,sd_metric)
   call check_array('MRC finite difference',mrc2,mrc1/sqrt(dot_product(w,matmul(sample_cov,w))),3.0e-5_dp)

   resampled=resample_correlated(reshape([(real(i,dp),i=1,30)],[10,3]),40, &
      reshape([1.0_dp,0.4_dp,0.2_dp,0.4_dp,1.0_dp,0.3_dp,0.2_dp,0.3_dp,1.0_dp],[3,3]),seed=321_i8,status=status)
   call check_true('resample status',status==nmof_ok.and.all(shape(resampled)==[40,3]))
   call check_true('resample support',all(resampled(:,1)>=1.0_dp).and.all(resampled(:,1)<=10.0_dp))

   paths=geometric_brownian_motion(4,8,0.02_dp,0.04_dp,1.0_dp,s0=100.0_dp,antithetic=.true.,seed=77_i8,status=status)
   call check_true('gbm shape',all(shape(paths)==[9,4]).and.status==nmof_ok)
   sample=geometric_brownian_bridge(3,10,100.0_dp,120.0_dp,0.04_dp,1.0_dp,log_input=.true.,exponentiate=.true.,seed=88_i8,status=status)
   call check_array('bridge start',sample(1,:),[100.0_dp,100.0_dp,100.0_dp],1.0e-12_dp)
   call check_array('bridge end',sample(11,:),[120.0_dp,120.0_dp,120.0_dp],1.0e-10_dp)

   cp=cppi([100.0_dp,105.0_dp,95.0_dp,110.0_dp],3.0_dp,0.8_dp,0.01_dp)
   call check_true('CPPI status',cp%status==nmof_ok)
   call check_true('CPPI nonnegative',all(cp%value>=0.0_dp).and.all(cp%floor>0.0_dp))

   do i=1,12
      m(i,1)=0.01_dp*real(i,dp)
      m(i,2)=0.02_dp*real(13-i,dp)
      m(i,3)=0.005_dp*(-1.0_dp)**i
   end do
   pbo=probability_backtest_overfitting(m,s=6)
   call check_true('PBO status',pbo%status==nmof_ok)
   call check_true('PBO range',pbo%pbo>=0.0_dp.and.pbo%pbo<=1.0_dp)

   if(failures>0) then
      write(*,'(a,i0)') 'test_core failures: ',failures
      error stop 1
   end if
   write(*,'(a)') 'test_core: PASS'
contains
   function cubic_root(z,context) result(f)
      real(dp),intent(in)::z
      class(*),intent(in),optional::context
      real(dp)::f
      f=z**3-z-0.25_dp
   end function cubic_root
   function sd_metric(z,context) result(value)
      real(dp),intent(in)::z(:)
      class(*),intent(in),optional::context
      real(dp)::value
      value=standard_deviation(z)
   end function sd_metric
   subroutine check_array(name,actual,expected,tol)
      character(len=*),intent(in)::name
      real(dp),intent(in)::actual(:),expected(:),tol
      if(size(actual)/=size(expected).or.maxval(abs(actual-expected))>tol) then
         failures=failures+1; write(*,'(a)') trim(name)//' failed'
      end if
   end subroutine check_array
   subroutine check_close(name,actual,expected,tol)
      character(len=*),intent(in)::name
      real(dp),intent(in)::actual,expected,tol
      if(abs(actual-expected)>tol) then
         failures=failures+1; write(*,'(a,2(1x,es16.8))') trim(name)//' failed:',actual,expected
      end if
   end subroutine check_close
   subroutine check_true(name,condition)
      character(len=*),intent(in)::name
      logical,intent(in)::condition
      if(.not.condition) then; failures=failures+1; write(*,'(a)') trim(name)//' failed'; end if
   end subroutine check_true
end program test_core
