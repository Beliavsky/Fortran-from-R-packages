! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran computational translation of tvgarch 2.4.3.
module tvgarch_multivariate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchx_kinds, only : dp
   use garchx_math, only : random_normal_vector
   use garchx_linalg, only : invert_matrix, cholesky_lower
   use garchx_optimize, only : bounded_nelder_mead, numerical_hessian
   use tvgarch_transition, only : tv_component
   use tvgarch_model, only : tvgarch_spec, tvgarch_fit, tvgarch_simulation, &
                             fit_tvgarch, tvgarch_simulate, tvgarch_forecast
   implicit none
   private
   real(dp), parameter :: log_two_pi = log(2.0_dp*acos(-1.0_dp))

   type, public :: dcc_fit
      real(dp) :: alpha = 0.0_dp
      real(dp) :: beta = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp), allocatable :: correlations(:, :, :)
      real(dp), allocatable :: vcov(:, :), se(:)
      integer :: status = 1
      integer :: iterations = 0
   end type dcc_fit

   type, public :: mtvgarch_fit
      type(tvgarch_fit), allocatable :: margins(:)
      real(dp), allocatable :: y(:, :), sigma2(:, :), residuals(:, :), g(:, :), h(:, :)
      real(dp), allocatable :: correlation(:, :), correlations(:, :, :)
      type(dcc_fit) :: dcc
      real(dp) :: loglik = -huge(1.0_dp)
      logical :: dynamic_correlation = .false.
      integer :: status = 1
   end type mtvgarch_fit

   type, public :: mtvgarch_simulation
      real(dp), allocatable :: y(:, :), sigma2(:, :), innovations(:, :), g(:, :), h(:, :)
      real(dp), allocatable :: correlations(:, :, :)
      integer :: status = 1
   end type mtvgarch_simulation

   public :: correlation_matrix, dcc_objective, dcc_filter, fit_dcc
   public :: mtvgarch_simulate, fit_mtvgarch, fit_mtvgarch_spillover
   public :: mtvgarch_forecast, lower_triangle_series
contains
   subroutine correlation_matrix(x, r, status)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: r(:, :)
      integer, intent(out) :: status
      integer :: n, m, i, j
      real(dp), allocatable :: centered(:, :), sd(:)
      real(dp) :: xbar
      n = size(x,1); m = size(x,2)
      if (n < 2 .or. m < 1 .or. .not. all(ieee_is_finite(x))) then
         status=1; allocate(r(0,0)); return
      end if
      allocate(centered(n,m), sd(m), r(m,m))
      do j=1,m
         xbar=sum(x(:,j))/real(n,dp)
         centered(:,j)=x(:,j)-xbar
         sd(j)=sqrt(sum(centered(:,j)**2)/real(n-1,dp))
         if (sd(j) <= epsilon(1.0_dp)) then
            status=2; r=0.0_dp; return
         end if
         centered(:,j)=centered(:,j)/sd(j)
      end do
      do i=1,m
         do j=1,m
            r(i,j)=dot_product(centered(:,i),centered(:,j))/real(n-1,dp)
         end do
      end do
      do i=1,m; r(i,i)=1.0_dp; end do
      status=0
   end subroutine correlation_matrix

   subroutine normalize_covariance(q, r, status)
      real(dp), intent(in) :: q(:, :)
      real(dp), intent(out) :: r(:, :)
      integer, intent(out) :: status
      integer :: i,j,m
      m=size(q,1)
      if (size(q,2)/=m .or. any([(q(i,i),i=1,m)]<=0.0_dp)) then
         status=1; r=0.0_dp; return
      end if
      do i=1,m
         do j=1,m
            r(i,j)=q(i,j)/sqrt(q(i,i)*q(j,j))
         end do
      end do
      status=0
   end subroutine normalize_covariance

   subroutine logdet_inverse(a, logdet, ainv, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: logdet
      real(dp), allocatable, intent(out) :: ainv(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: l(:, :)
      integer :: i
      call cholesky_lower(a,l,status)
      if (status/=0) then
         logdet=huge(1.0_dp); allocate(ainv(0,0)); return
      end if
      logdet=0.0_dp
      do i=1,size(a,1); logdet=logdet+2.0_dp*log(l(i,i)); end do
      call invert_matrix(a,ainv,status)
   end subroutine logdet_inverse

   subroutine dcc_filter(par, z, correlations, loglik_terms, status, sigma2, initial_window)
      real(dp), intent(in) :: par(2), z(:, :)
      real(dp), allocatable, intent(out) :: correlations(:, :, :), loglik_terms(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: sigma2(:, :)
      integer, intent(in), optional :: initial_window
      integer :: n,m,t,w,st
      real(dp), allocatable :: rbar(:,:), qbar(:,:), qt(:,:), zt2(:,:), rt(:,:), rinv(:,:)
      real(dp) :: ldet, quad, marginal_term
      n=size(z,1); m=size(z,2)
      if (n<2 .or. m<2 .or. par(1)<0.0_dp .or. par(2)<0.0_dp .or. sum(par)>=1.0_dp) then
         status=1; allocate(correlations(0,0,0),loglik_terms(0)); return
      end if
      if (present(sigma2)) then
         if (size(sigma2,1)/=n .or. size(sigma2,2)/=m .or. any(sigma2<=0.0_dp)) then
            status=2; allocate(correlations(0,0,0),loglik_terms(0)); return
         end if
      end if
      call correlation_matrix(z,rbar,st)
      if (st/=0) then; status=3; allocate(correlations(0,0,0),loglik_terms(0)); return; end if
      allocate(qbar(m,m),qt(m,m),zt2(m,m),rt(m,m),correlations(m,m,n),loglik_terms(n))
      qbar=(1.0_dp-par(1)-par(2))*rbar
      w=min(n,100); if(present(initial_window)) w=min(n,max(2,initial_window))
      call correlation_matrix(z(1:w,:),qt,st)
      if(st/=0) qt=rbar
      zt2=qt
      do t=1,n
         qt=qbar+par(1)*zt2+par(2)*qt
         call normalize_covariance(qt,rt,st)
         if(st/=0) then; status=4; return; end if
         correlations(:,:,t)=rt
         call logdet_inverse(rt,ldet,rinv,st)
         if(st/=0) then; status=5; return; end if
         quad=dot_product(z(t,:),matmul(rinv,z(t,:)))
         marginal_term=0.0_dp
         if(present(sigma2)) marginal_term=sum(log(sigma2(t,:)))
         loglik_terms(t)=-0.5_dp*(real(m,dp)*log_two_pi+marginal_term+ldet+quad)
         zt2=spread(z(t,:),2,m)*spread(z(t,:),1,m)
      end do
      status=0
   end subroutine dcc_filter

   subroutine dcc_objective(par, z, value, status, sigma2, per_observation)
      real(dp), intent(in) :: par(2), z(:, :)
      real(dp), intent(out) :: value
      integer, intent(out) :: status
      real(dp), intent(in), optional :: sigma2(:, :)
      real(dp), allocatable, intent(out), optional :: per_observation(:)
      real(dp), allocatable :: correlations(:,:,:), terms(:)
      call dcc_filter(par,z,correlations,terms,status,sigma2)
      if(status/=0) then
         value=huge(1.0_dp)*0.01_dp
         if(present(per_observation)) allocate(per_observation(0))
         return
      end if
      value=-sum(terms)
      if(present(per_observation)) then
         allocate(per_observation(size(terms))); per_observation=-terms
      end if
   end subroutine dcc_objective

   subroutine fit_dcc(z, fit, sigma2, initial, turbo)
      real(dp), intent(in) :: z(:, :)
      type(dcc_fit), intent(out) :: fit
      real(dp), intent(in), optional :: sigma2(:, :), initial(2)
      logical, intent(in), optional :: turbo
      real(dp) :: x0(2),lo(2),hi(2),fbest,step
      real(dp), allocatable :: best(:),terms(:),hess(:,:),hinv(:,:),scores(:,:),meat(:,:),pp(:),pm(:)
      integer :: st,i,invst,n
      logical :: fast
      x0=[0.05_dp,0.90_dp]; if(present(initial)) x0=initial
      lo=[5.0e-5_dp,5.0e-5_dp]; hi=[0.9998_dp,0.9998_dp]
      call bounded_nelder_mead(obj,x0,lo,hi,best,fbest,fit%status,fit%iterations,3000,1.0e-8_dp)
      fit%alpha=best(1); fit%beta=best(2)
      call dcc_filter(best,z,fit%correlations,terms,st,sigma2)
      if(st/=0) then; fit%status=st; return; end if
      fit%loglik=sum(terms)
      allocate(fit%vcov(2,2),fit%se(2)); fit%vcov=0.0_dp; fit%se=0.0_dp
      fast=.false.; if(present(turbo)) fast=turbo
      if(.not.fast) then
         call numerical_hessian(obj,best,hess)
         call invert_matrix(hess,hinv,invst)
         if(invst==0) then
            n=size(z,1); allocate(scores(n,2),meat(2,2))
            do i=1,2
               step=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(best(i)))
               pp=best; pm=best; pp(i)=pp(i)+step; pm(i)=pm(i)-step
               call perobs(pp,terms)
               scores(:,i)=terms
               call perobs(pm,terms)
               scores(:,i)=(scores(:,i)-terms)/(2.0_dp*step)
            end do
            meat=matmul(transpose(scores),scores)
            fit%vcov=matmul(hinv,matmul(meat,hinv))
            fit%se=sqrt(max([fit%vcov(1,1),fit%vcov(2,2)],0.0_dp))
         end if
      end if
      fit%status=0
   contains
      function obj(p) result(v)
         real(dp),intent(in)::p(:); real(dp)::v; integer::s
         if(size(p)/=2 .or. p(1)+p(2)>=0.99995_dp) then; v=huge(1.0_dp)*0.01_dp; return; end if
         call dcc_objective(p,z,v,s,sigma2)
         if(s/=0) v=huge(1.0_dp)*0.01_dp
      end function obj
      subroutine perobs(p,v)
         real(dp),intent(in)::p(:); real(dp),allocatable,intent(out)::v(:)
         real(dp)::total; integer::s
         call dcc_objective(p,z,total,s,sigma2,v)
         if(s/=0) then; if(allocated(v))deallocate(v); allocate(v(size(z,1))); v=huge(1.0_dp)*0.01_dp; end if
      end subroutine perobs
   end subroutine fit_dcc

   subroutine mtvgarch_simulate(n, specs, par_h, correlation, result, xtv, xreg, dcc_par, innovations)
      integer,intent(in)::n
      type(tvgarch_spec),intent(in)::specs(:)
      real(dp),intent(in)::par_h(:,:),correlation(:,:)
      type(mtvgarch_simulation),intent(out)::result
      real(dp),intent(in),optional::xtv(:),xreg(:,:,:),dcc_par(2),innovations(:,:)
      integer::m,t,j,st
      real(dp),allocatable::base(:,:),z(:,:),qbar(:,:),qt(:,:),zt2(:,:),rt(:,:),l(:,:),draw(:)
      type(tvgarch_simulation)::one
      m=size(specs)
      if(m<2 .or. size(par_h,2)/=m .or. size(correlation,1)/=m .or. size(correlation,2)/=m)then
         result%status=1;return
      end if
      allocate(base(n,m),z(n,m),result%correlations(m,m,n))
      if(present(innovations))then
         if(size(innovations,1)/=n .or. size(innovations,2)/=m)then;result%status=2;return;end if
         base=innovations
      else
         do j=1,m; call random_normal_vector(base(:,j)); end do
      end if
      allocate(draw(m),rt(m,m),l(m,m))
      if(present(dcc_par))then
         if(sum(dcc_par)>=1.0_dp .or. any(dcc_par<0.0_dp))then;result%status=3;return;end if
         allocate(qbar(m,m),qt(m,m),zt2(m,m))
         qbar=(1.0_dp-sum(dcc_par))*correlation;qt=correlation;zt2=correlation
         do t=1,n
            qt=qbar+dcc_par(1)*zt2+dcc_par(2)*qt
            call normalize_covariance(qt,rt,st); if(st/=0)then;result%status=4;return;end if
            call cholesky_lower(rt,l,st);if(st/=0)then;result%status=5;return;end if
            draw=matmul(l,base(t,:));z(t,:)=draw;zt2=spread(draw,2,m)*spread(draw,1,m)
            result%correlations(:,:,t)=rt
         end do
      else
         call cholesky_lower(correlation,l,st);if(st/=0)then;result%status=6;return;end if
         do t=1,n;z(t,:)=matmul(l,base(t,:));result%correlations(:,:,t)=correlation;end do
      end if
      allocate(result%y(n,m),result%sigma2(n,m),result%innovations(n,m),result%g(n,m),result%h(n,m))
      do j=1,m
         if(present(xreg))then
            call tvgarch_simulate(n,specs(j),par_h(:,j),one,xtv=xtv, &
                 xreg=xreg(:,1:specs(j)%garch%xreg_count,j),innovations=z(:,j))
         else
            call tvgarch_simulate(n,specs(j),par_h(:,j),one,xtv=xtv,innovations=z(:,j))
         end if
         if(one%status/=0)then;result%status=10+j;return;end if
         result%y(:,j)=one%y;result%sigma2(:,j)=one%sigma2;result%innovations(:,j)=one%innovations
         result%g(:,j)=one%g;result%h(:,j)=one%h
      end do
      result%status=0
   end subroutine mtvgarch_simulate

   subroutine fit_mtvgarch(y,specs,fit,xtv,xreg,dynamic_correlation,turbo)
      real(dp),intent(in)::y(:,:)
      type(tvgarch_spec),intent(in)::specs(:)
      type(mtvgarch_fit),intent(out)::fit
      real(dp),intent(in),optional::xtv(:),xreg(:,:,:)
      logical,intent(in),optional::dynamic_correlation,turbo
      integer::n,m,j,st
      logical::dcc,fast
      n=size(y,1);m=size(y,2)
      if(size(specs)/=m)then;fit%status=1;return;end if
      dcc=.false.;if(present(dynamic_correlation))dcc=dynamic_correlation
      fast=.false.;if(present(turbo))fast=turbo
      allocate(fit%margins(m),fit%y(n,m),fit%sigma2(n,m),fit%residuals(n,m),fit%g(n,m),fit%h(n,m))
      fit%y=y
      do j=1,m
         if(present(xreg))then
            call fit_tvgarch(y(:,j),specs(j),fit%margins(j),xtv=xtv, &
                 xreg=xreg(:,1:specs(j)%garch%xreg_count,j),turbo=fast)
         else
            call fit_tvgarch(y(:,j),specs(j),fit%margins(j),xtv=xtv,turbo=fast)
         end if
         if(fit%margins(j)%status/=0 .and. .not.allocated(fit%margins(j)%sigma2))then;fit%status=10+j;return;end if
         fit%sigma2(:,j)=fit%margins(j)%sigma2;fit%residuals(:,j)=fit%margins(j)%residuals
         fit%g(:,j)=fit%margins(j)%g;fit%h(:,j)=fit%margins(j)%h
      end do
      call correlation_matrix(fit%residuals,fit%correlation,st)
      if(st/=0)then;fit%status=30;return;end if
      fit%dynamic_correlation=dcc
      if(dcc)then
         call fit_dcc(fit%residuals,fit%dcc,sigma2=fit%sigma2,turbo=fast)
         if(fit%dcc%status/=0)then;fit%status=31;return;end if
         fit%correlations=fit%dcc%correlations;fit%loglik=fit%dcc%loglik
      else
         allocate(fit%correlations(m,m,n))
         do j=1,n;fit%correlations(:,:,j)=fit%correlation;end do
         call multivariate_loglik(fit%residuals,fit%sigma2,fit%correlations,fit%loglik,st)
      end if
      fit%status=0
   end subroutine fit_mtvgarch

   subroutine fit_mtvgarch_spillover(y,specs,order_x,fit,xtv,dynamic_correlation,max_outer,turbo)
      real(dp),intent(in)::y(:,:)
      type(tvgarch_spec),intent(in)::specs(:)
      logical,intent(in)::order_x(:,:)
      type(mtvgarch_fit),intent(out)::fit
      real(dp),intent(in),optional::xtv(:)
      logical,intent(in),optional::dynamic_correlation,turbo
      integer,intent(in),optional::max_outer
      integer::n,m,i,j,k,it,maxit,nx
      real(dp),allocatable::phi2(:,:),xall(:,:,:),old_sigma(:,:)
      type(tvgarch_spec),allocatable::local_specs(:)
      n=size(y,1);m=size(y,2)
      if(size(order_x,1)/=m .or. size(order_x,2)/=m .or. size(specs)/=m)then;fit%status=1;return;end if
      nx=maxval([(count(order_x(i,:)),i=1,m)])
      if(nx==0)then;call fit_mtvgarch(y,specs,fit,xtv=xtv,dynamic_correlation=dynamic_correlation,turbo=turbo);return;end if
      allocate(phi2(n,m),xall(n,nx,m),local_specs(m),old_sigma(n,m));phi2=y*y
      local_specs=specs
      do i=1,m;local_specs(i)%garch%xreg_count=count(order_x(i,:));end do
      maxit=10;if(present(max_outer))maxit=max_outer
      do it=1,maxit
         xall=0.0_dp
         do i=1,m
            k=0
            do j=1,m
               if(order_x(i,j))then
                  k=k+1;xall(1,k,i)=sum(phi2(:,j))/real(n,dp);xall(2:,k,i)=phi2(1:n-1,j)
               end if
            end do
         end do
         old_sigma=0.0_dp;if(allocated(fit%sigma2))old_sigma=fit%sigma2
         call fit_mtvgarch(y,local_specs,fit,xtv=xtv,xreg=xall,dynamic_correlation=dynamic_correlation,turbo=turbo)
         if(fit%status/=0)return
         phi2=y*y/max(fit%g,epsilon(1.0_dp))
         if(it>1 .and. maxval(abs(fit%sigma2-old_sigma))/max(1.0_dp,maxval(abs(old_sigma)))<1.0e-3_dp)exit
      end do
   end subroutine fit_mtvgarch_spillover

   subroutine multivariate_loglik(z,sigma2,rpath,value,status)
      real(dp),intent(in)::z(:,:),sigma2(:,:),rpath(:,:,:)
      real(dp),intent(out)::value
      integer,intent(out)::status
      integer::t,n,m,st
      real(dp)::ldet
      real(dp),allocatable::rinv(:,:)
      n=size(z,1);m=size(z,2);value=0.0_dp
      do t=1,n
         call logdet_inverse(rpath(:,:,t),ldet,rinv,st)
         if(st/=0)then;status=1;value=-huge(1.0_dp);return;end if
         value=value-0.5_dp*(real(m,dp)*log_two_pi+sum(log(sigma2(t,:)))+ldet+ &
               dot_product(z(t,:),matmul(rinv,z(t,:))))
      end do
      status=0
   end subroutine multivariate_loglik

   subroutine mtvgarch_forecast(fit,n_ahead,forecast,status,new_xtv,new_xreg,n_sim)
      type(mtvgarch_fit),intent(in)::fit
      integer,intent(in)::n_ahead
      real(dp),allocatable,intent(out)::forecast(:,:)
      integer,intent(out)::status
      real(dp),intent(in),optional::new_xtv(:),new_xreg(:,:,:)
      integer,intent(in),optional::n_sim
      integer::m,j,st
      real(dp),allocatable::one(:)
      m=size(fit%margins);allocate(forecast(n_ahead,m))
      do j=1,m
         if(present(new_xreg))then
            call tvgarch_forecast(fit%margins(j),n_ahead,one,st,new_xtv, &
                 new_xreg(:,1:fit%margins(j)%spec%garch%xreg_count,j),n_sim)
         else
            call tvgarch_forecast(fit%margins(j),n_ahead,one,st,new_xtv,n_sim=n_sim)
         end if
         if(st/=0)then;status=st;return;end if
         forecast(:,j)=one
      end do
      status=0
   end subroutine mtvgarch_forecast

   subroutine lower_triangle_series(rpath,values,status)
      real(dp),intent(in)::rpath(:,:,:)
      real(dp),allocatable,intent(out)::values(:,:)
      integer,intent(out)::status
      integer::m,n,t,i,j,k
      m=size(rpath,1);n=size(rpath,3)
      if(size(rpath,2)/=m)then;status=1;allocate(values(0,0));return;end if
      allocate(values(n,m*(m-1)/2))
      do t=1,n
         k=0
         do j=1,m-1
            do i=j+1,m;k=k+1;values(t,k)=rpath(i,j,t);end do
         end do
      end do
      status=0
   end subroutine lower_triangle_series
end module tvgarch_multivariate
