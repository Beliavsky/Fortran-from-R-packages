program v05_extended
   use gamlss
   use gamlss_continuous, only : rGA
   use nlme_types, only : correlation_spec,variance_spec,nlme_control,COR_AR1,VAR_CONSTANT
   implicit none
   integer,parameter::n=72
   real(dp)::y(n),x(n,2),time(n),eta,e,rho
   integer::group(n),fold(n),i
   type(correlation_spec)::cor
   type(variance_spec)::var
   type(nlme_control)::nctl
   type(gamlss_control_t)::ctl
   type(correlated_rs_result_t)::fit
   type(gamlss_cv_result_t)::cv

   call seed_rng(505)
   rho=0.45_dp;e=0.15_dp*randn()
   do i=1,n
      x(i,1)=1.0_dp;x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
      time(i)=real(i,dp);group(i)=1;fold(i)=1+mod(i-1,3)
      if(i>1)e=rho*e+0.15_dp*sqrt(1.0_dp-rho*rho)*randn()
      eta=0.35_dp+0.65_dp*x(i,2)+e;y(i)=rGA(exp(eta),0.28_dp)
   end do
   cor%kind=COR_AR1;cor%fixed=.true.;allocate(cor%par(1));cor%par=rho
   var%kind=VAR_CONSTANT;var%fixed=.true.;allocate(var%par(0))
   nctl=nlme_control();nctl%optimize_covariance=.false.;nctl%reml=.false.
   ctl=gamlss_control_t();ctl%n_cyc=6;ctl%inner_cyc=30;ctl%c_crit=1.0e-4_dp
   call fit_gamlss_correlated_rs(y,x,GAMLSS_GA,fit,correlation=cor,variance=var,time=time,group=group, &
      control=ctl,nlme_control_in=nctl,max_outer=6)
   call cross_validate_gamlss(y,x,fold,GAMLSS_GA,cv,control=ctl)
   write(*,'(a,2f10.4)')'Correlated Gamma mu coefficients: ',fit%model%mu%coefficients
   write(*,'(a,f10.4)')'Fixed AR(1) correlation: ',fit%correlation_parameters(1)
   write(*,'(a,f10.4)')'Three-fold mean log score: ',cv%mean_log_score
contains
   subroutine seed_rng(base)
      integer,intent(in)::base
      integer,allocatable::s(:)
      integer::m,k
      call random_seed(size=m);allocate(s(m));do k=1,m;s(k)=base+31*k;end do;call random_seed(put=s)
   end subroutine seed_rng
   real(dp) function randn() result(z)
      real(dp)::u1,u2
      call random_number(u1);call random_number(u2);u1=max(u1,1.0e-12_dp)
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function randn
end program v05_extended
