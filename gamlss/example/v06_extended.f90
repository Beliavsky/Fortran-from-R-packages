program v06_extended
   use gamlss
   use gamlss_continuous, only : qGA
   use gamlss_special, only : normal_cdf
   use nlme_types, only : correlation_spec,COR_AR1
   implicit none
   integer,parameter :: n=72,ng=8,m=10,nr=ng*m
   real(dp) :: y(n),x(n,2),xs(n,1),time(n),z,u,mu,rho
   real(dp) :: yr(nr),xm(nr,2),xsr(nr,1),zr(nr,1,2),umu(ng),usig(ng),xx
   integer :: group(n),gr(nr),i,g,k
   logical :: act(4)
   type(correlation_spec) :: cor
   type(gaussian_copula_result_t) :: cf
   type(joint_random_effects_result_t) :: rf
   type(gamlss_control_t) :: ctl
   call seed_rng(2606)
   rho=0.55_dp;z=randn()
   do i=1,n
      x(i,:)=[1.0_dp,-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)];xs(i,1)=1.0_dp
      time(i)=real(i,dp);group(i)=1
      if(i>1)z=rho*z+sqrt(1.0_dp-rho*rho)*randn()
      u=min(1.0_dp-1.0e-9_dp,max(1.0e-9_dp,normal_cdf(z)))
      mu=exp(0.40_dp+0.60_dp*x(i,2));y(i)=qGA(u,mu,0.35_dp)
   end do
   cor%kind=COR_AR1;cor%fixed=.false.;allocate(cor%par(1));cor%par=0.20_dp
   ctl=gamlss_control_t();ctl%n_cyc=15;ctl%inner_cyc=30;ctl%c_crit=1.0e-4_dp
   call fit_gamlss_gaussian_copula(y,x,GAMLSS_GA,cf,correlation=cor,x_sigma=xs,time=time,group=group, &
      control=ctl,max_iter=100,tolerance=5.0e-6_dp)
   print '(a,2f10.4)','Copula Gamma mu coefficients: ',cf%model%mu%coefficients
   print '(a,f10.4)','Gaussian-copula AR(1) rho: ',cf%correlation_parameters(1)
   print '(a,f10.4)','Copula log-likelihood contribution: ',cf%copula_log_likelihood

   do g=1,ng
      umu(g)=0.38_dp*sin(0.72_dp*real(g,dp))
      usig(g)=0.50_dp*umu(g)+0.07_dp*cos(0.53_dp*real(g,dp))
   end do
   k=0
   do g=1,ng;do i=1,m
      k=k+1;gr(k)=g;xx=-1.0_dp+2.0_dp*real(i-1,dp)/real(m-1,dp)
      xm(k,:)=[1.0_dp,xx];xsr(k,1)=1.0_dp;zr(k,1,1)=1.0_dp;zr(k,1,2)=1.0_dp
      yr(k)=1.10_dp+0.70_dp*xx+umu(g)+exp(-1.0_dp+usig(g))*randn()
   end do;end do
   act=[.true.,.true.,.false.,.false.]
   call fit_gamlss_joint_random_effects(yr,xm,zr,gr,GAMLSS_NO,rf,active_parameters=act,x_sigma=xsr, &
      control=ctl,max_outer=5,max_inner=70,tol_cov=3.0e-4_dp,tolerance=5.0e-6_dp)
   print '(a,2f10.4)','Joint-RE fixed mu coefficients: ',rf%model%mu%coefficients
   print '(a,f10.5)','Estimated Cov(b_mu,b_sigma): ',rf%joint_covariance(1,2)
contains
   subroutine seed_rng(base)
      integer,intent(in) :: base
      integer,allocatable :: s(:)
      integer :: ns,j
      call random_seed(size=ns);allocate(s(ns));do j=1,ns;s(j)=base+61*j;end do;call random_seed(put=s)
   end subroutine seed_rng
   real(dp) function randn() result(v)
      real(dp) :: u1,u2
      call random_number(u1);call random_number(u2);u1=max(u1,1.0e-12_dp)
      v=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function randn
end program v06_extended
