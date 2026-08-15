program test_v09
   use gamlss
   use gamlss_kinds, only : dp
   use gamlss_fit, only : GAMLSS_NO,GAMLSS_GA
   implicit none
   integer :: failures
   failures=0
   call test_identity_link(failures)
   call test_log_link_oracle(failures)
   call test_sigma_parameter(failures)
   call test_object_adapter(failures)
   call test_fitted_random_scale(failures)
   if(failures/=0)then
      write(*,'(a,i0)') 'test_v09: FAIL ',failures
      error stop 1
   end if
   write(*,'(a)') 'test_v09: PASS'
contains
   subroutine test_identity_link(failures)
      integer,intent(inout) :: failures
      real(dp) :: eta(3),sb
      type(marginal_prediction_result_t) :: r1,r2
      eta=[-0.7_dp,0.2_dp,1.4_dp];sb=0.65_dp
      call marginal_predict_eta(eta,sb,GAMLSS_NO,1,r1,MARGINAL_INTEGRATE)
      call marginal_predict_eta(eta,sb,GAMLSS_NO,1,r2,MARGINAL_QFUNCTION)
      if(r1%status/=0.or.maxval(abs(r1%fitted-eta))>2.0e-9_dp)failures=failures+1
      if(r2%status/=0.or.maxval(abs(r2%fitted-eta))>2.0e-12_dp)failures=failures+1
   end subroutine test_identity_link

   subroutine test_log_link_oracle(failures)
      integer,intent(inout) :: failures
      real(dp) :: eta(3),sb,truth(3)
      type(marginal_prediction_result_t) :: ri,rq,rr
      eta=[-0.5_dp,0.1_dp,0.8_dp];sb=0.55_dp
      truth=exp(eta+0.5_dp*sb*sb)
      call marginal_predict_eta(eta,sb,GAMLSS_GA,1,ri,MARGINAL_INTEGRATE)
      call marginal_predict_eta(eta,sb,GAMLSS_GA,1,rq,MARGINAL_QFUNCTION)
      call seed_rng(909)
      call marginal_predict_eta(eta,sb,GAMLSS_GA,1,rr,MARGINAL_RANDOM,n_random=30000)
      if(ri%status/=0.or.maxval(abs(ri%fitted-truth)/truth)>2.0e-8_dp)failures=failures+1
      if(rq%status/=0.or.maxval(abs(rq%fitted-truth)/truth)>4.0e-3_dp)failures=failures+1
      if(rr%status/=0.or.maxval(abs(rr%fitted-truth)/truth)>2.5e-2_dp)failures=failures+1
      if(rq%integration_points/=999.or.rr%integration_points/=30000)failures=failures+1
   end subroutine test_log_link_oracle

   subroutine test_sigma_parameter(failures)
      integer,intent(inout) :: failures
      real(dp) :: eta(2),sb,truth(2)
      type(marginal_prediction_result_t) :: r
      eta=[-1.1_dp,-0.3_dp];sb=0.4_dp;truth=exp(eta+0.5_dp*sb*sb)
      call marginal_predict_eta(eta,sb,GAMLSS_NO,2,r,MARGINAL_INTEGRATE)
      if(r%status/=0.or.maxval(abs(r%fitted-truth)/truth)>2.0e-8_dp)failures=failures+1
   end subroutine test_sigma_parameter

   subroutine test_object_adapter(failures)
      integer,intent(inout) :: failures
      type(random_intercept_result_t) :: fit
      type(marginal_prediction_result_t) :: r
      integer :: group(4)
      real(dp) :: base(4)
      group=[10,10,20,20];base=[0.2_dp,0.4_dp,-0.1_dp,0.7_dp]
      fit%status=0;fit%parameter=1;fit%model%family=GAMLSS_NO;fit%sigma_b=0.3_dp
      allocate(fit%levels(2),fit%effects(2),fit%model%mu%eta(4))
      fit%levels=[10,20];fit%effects=[0.25_dp,-0.15_dp]
      fit%model%mu%eta=base+[0.25_dp,0.25_dp,-0.15_dp,-0.15_dp]
      call get_marginal_random_intercept(fit,group,r,MARGINAL_NONE)
      if(r%status/=0.or.maxval(abs(r%fitted-base))>1.0e-12_dp)failures=failures+1
      if(abs(r%sigma_b-0.3_dp)>1.0e-12_dp)failures=failures+1
   end subroutine test_object_adapter

   subroutine test_fitted_random_scale(failures)
      integer,intent(inout) :: failures
      integer,parameter :: ng=4,m=5,n=ng*m
      real(dp) :: y(n),x(n,1)
      integer :: group(n),g,i,k
      type(random_intercept_result_t) :: fit
      type(gamlss_control_t) :: ctl
      k=0;x=1.0_dp
      do g=1,ng
         do i=1,m
            k=k+1;group(k)=g
            y(k)=1.0_dp+0.18_dp*real(g-2,dp)+0.04_dp*sin(real(k,dp))
         end do
      end do
      ctl=gamlss_control_t();ctl%n_cyc=10;ctl%inner_cyc=20
      call fit_gamlss_random_intercept(y,x,group,GAMLSS_NO,fit,control=ctl,use_nlme_start=.false.)
      if(fit%status/=0.or.fit%random_edf<=0.0_dp.or.fit%sigma_b<0.0_dp)failures=failures+1
      if(.not.(fit%sigma_b<huge(1.0_dp)))failures=failures+1
   end subroutine test_fitted_random_scale

   subroutine seed_rng(base)
      integer,intent(in) :: base
      integer,allocatable :: s(:)
      integer :: n,i
      call random_seed(size=n);allocate(s(n))
      do i=1,n;s(i)=base+97*i;end do
      call random_seed(put=s)
   end subroutine seed_rng
end program test_v09
