program v04_extended
   use gamlss
   use nlme_types, only : correlation_spec, variance_spec, nlme_control, COR_AR1, VAR_CONSTANT
   implicit none
   integer, parameter :: n = 120, ng = 8, m = 15
   real(dp) :: y(n), x(n,2), time(n), e, rho, sd, xx
   real(dp) :: xm(n,2), xs(n,1), u_mu(ng), u_sig(ng)
   integer :: group(n), i, j, k
   logical :: active(4)
   type(correlation_spec) :: cor
   type(variance_spec) :: var
   type(nlme_control) :: nctl
   type(gamlss_control_t) :: gctl
   type(correlated_no_result_t) :: glsfit
   type(multi_random_intercept_result_t) :: refit

   call seed_rng(404)
   rho = 0.60_dp
   sd = 0.40_dp
   e = sd*randn()
   do i = 1, n
      xx = -1.0_dp + 2.0_dp*real(i-1,dp)/real(n-1,dp)
      x(i,:) = [1.0_dp, xx]
      time(i) = real(i,dp)
      group(i) = 1
      if (i > 1) e = rho*e + sd*sqrt(1.0_dp-rho*rho)*randn()
      y(i) = 0.7_dp + 1.15_dp*xx + e
   end do
   cor%kind = COR_AR1
   allocate(cor%par(1)); cor%par = 0.25_dp
   var%kind = VAR_CONSTANT
   allocate(var%par(0))
   nctl = nlme_control(); nctl%reml = .false.; nctl%max_iter = 160
   call fit_gamlss_no_gls(y,x,glsfit,correlation=cor,variance=var,time=time,group=group,control=nctl)
   write(*,'(a,2f10.4)') 'Correlated-NO coefficients: ', glsfit%gls%beta
   write(*,'(a,f10.4)') 'Fitted AR(1) correlation:  ', glsfit%gls%correlation_parameters(1)

   do j = 1, ng
      u_mu(j) = 0.28_dp*sin(0.8_dp*real(j,dp))
      u_sig(j) = 0.18_dp*cos(0.7_dp*real(j,dp))
   end do
   k = 0
   do j = 1, ng
      do i = 1, m
         k = k + 1
         xx = -1.0_dp + 2.0_dp*real(i-1,dp)/real(m-1,dp)
         group(k) = j
         xm(k,:) = [1.0_dp, xx]
         xs(k,1) = 1.0_dp
         sd = exp(-0.85_dp + u_sig(j))
         y(k) = 1.1_dp + 0.65_dp*xx + u_mu(j) + sd*randn()
      end do
   end do
   active = [.true.,.true.,.false.,.false.]
   gctl = gamlss_control_t(); gctl%n_cyc = 35; gctl%inner_cyc = 70
   call fit_gamlss_multi_random_intercept(y,xm,group,GAMLSS_NO,refit, &
      active_parameters=active,x_sigma=xs,control=gctl)
   write(*,'(a,2f10.4)') 'mu/sigma random-effect correlations: ', &
      correlation(refit%effects(:,1),u_mu), correlation(refit%effects(:,2),u_sig)
contains
   subroutine seed_rng(base)
      integer, intent(in) :: base
      integer, allocatable :: seed(:)
      integer :: ns, ii
      call random_seed(size=ns)
      allocate(seed(ns))
      do ii=1,ns
         seed(ii)=base+37*ii
      end do
      call random_seed(put=seed)
   end subroutine seed_rng

   real(dp) function randn() result(z)
      real(dp) :: u1,u2
      call random_number(u1); call random_number(u2)
      u1=max(u1,1.0e-12_dp)
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function randn

   real(dp) function correlation(a,b) result(c)
      real(dp),intent(in) :: a(:),b(:)
      real(dp) :: am,bm,da,db
      am=sum(a)/real(size(a),dp); bm=sum(b)/real(size(b),dp)
      da=sum((a-am)**2); db=sum((b-bm)**2)
      if(da<=tiny(1.0_dp).or.db<=tiny(1.0_dp))then
         c=0.0_dp
      else
         c=sum((a-am)*(b-bm))/sqrt(da*db)
      end if
   end function correlation
end program v04_extended
