! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_fitting
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ghyp_kinds, only : dp
   use ghyp_model, only : ghyp_model_type, make_ghyp, ghyp_uv, gaussian_uv, &
      gaussian_mv, model_ghyp, model_hyp, model_nig, model_student, model_vg, &
      model_gaussian
   use ghyp_distribution, only : log_dghyp
   use ghyp_linalg, only : cholesky_lower, solve_linear, symmetrize
   use ghyp_optimize, only : nelder_mead, numerical_hessian
   implicit none
   private

   type, public :: fit_result
      type(ghyp_model_type) :: model
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: standard_error(:)
      integer :: iterations = 0
      integer :: n_parameters = 0
      logical :: converged = .false.
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type fit_result

   type, public :: likelihood_ratio_result
      real(dp) :: statistic = 0.0_dp
      integer :: degrees_freedom = 0
      real(dp) :: p_value = 1.0_dp
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type likelihood_ratio_result

   type, public :: model_selection_result
      type(fit_result) :: best
      character(len=16) :: family = ''
      real(dp), allocatable :: aic(:)
      character(len=16), allocatable :: names(:)
      logical :: ok = .false.
   end type model_selection_result

   public :: fit_ghyp_uv, fit_ghyp_mv, fit_gaussian_uv, fit_gaussian_mv
   public :: likelihood_ratio_test, step_aic_ghyp

   real(dp), allocatable, save :: active_uv_data(:)
   real(dp), allocatable, save :: active_mv_data(:,:)
   integer, save :: active_family = model_ghyp
   integer, save :: active_dimension = 1

contains

   pure function sample_mean(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sum(x)/real(size(x),dp)
   end function sample_mean

   pure function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, m
      m = sample_mean(x)
      value = sum((x-m)**2)/real(max(1,size(x)-1),dp)
   end function sample_variance

   function fit_gaussian_uv(data) result(result)
      real(dp), intent(in) :: data(:)
      type(fit_result) :: result
      real(dp) :: m, s2
      integer :: n
      n = size(data)
      if (n < 2 .or. .not. all(ieee_is_finite(data))) then
         result%message = 'at least two finite observations are required'
         return
      end if
      m = sample_mean(data)
      s2 = sum((data-m)**2)/real(n,dp)
      result%model = gaussian_uv(m,sqrt(max(s2,tiny(1.0_dp))))
      result%n_parameters = 2
      result%log_likelihood = -0.5_dp*real(n,dp)*(log(2.0_dp*acos(-1.0_dp)*s2)+1.0_dp)
      result%aic = -2.0_dp*result%log_likelihood+4.0_dp
      result%bic = -2.0_dp*result%log_likelihood+2.0_dp*log(real(n,dp))
      allocate(result%theta(2),result%standard_error(2),result%covariance(2,2))
      result%theta = [m,log(sqrt(s2))]
      result%covariance = 0.0_dp
      result%covariance(1,1) = s2/real(n,dp)
      result%covariance(2,2) = 1.0_dp/(2.0_dp*real(n,dp))
      result%standard_error = sqrt(max(0.0_dp,[result%covariance(1,1),result%covariance(2,2)]))
      result%converged = .true.; result%ok = .true.
   end function fit_gaussian_uv

   function decode_uv(theta, family) result(model)
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: family
      type(ghyp_model_type) :: model
      real(dp) :: lambda, chi, psi, mu, sigma, gamma, nu
      select case(family)
      case(model_ghyp)
         lambda=theta(1); chi=exp(theta(2)); psi=exp(theta(3)); mu=theta(4)
         sigma=exp(theta(5)); gamma=theta(6)
      case(model_hyp)
         lambda=1.0_dp; chi=exp(theta(1)); psi=exp(theta(2)); mu=theta(3)
         sigma=exp(theta(4)); gamma=theta(5)
      case(model_nig)
         lambda=-0.5_dp; chi=exp(theta(1)); psi=exp(theta(2)); mu=theta(3)
         sigma=exp(theta(4)); gamma=theta(5)
      case(model_student)
         nu=2.0_dp+exp(theta(1)); lambda=-0.5_dp*nu; chi=exp(theta(2)); psi=0.0_dp
         mu=theta(3); sigma=exp(theta(4)); gamma=theta(5)
      case(model_vg)
         lambda=exp(theta(1)); chi=0.0_dp; psi=exp(theta(2)); mu=theta(3)
         sigma=exp(theta(4)); gamma=theta(5)
      case default
         model%message='unknown family'; return
      end select
      model = ghyp_uv(lambda,chi,psi,mu,sigma,gamma)
   end function decode_uv

   function uv_objective(theta) result(value)
      real(dp), intent(in) :: theta(:)
      real(dp) :: value, lv
      type(ghyp_model_type) :: model
      integer :: i
      model = decode_uv(theta,active_family)
      if (.not. model%ok) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      value = 0.0_dp
      do i = 1, size(active_uv_data)
         lv = log_dghyp(active_uv_data(i),model)
         if (.not. ieee_is_finite(lv)) then
            value = huge(1.0_dp)/100.0_dp
            return
         end if
         value = value-lv
      end do
   end function uv_objective

   function fit_ghyp_uv(data, family, max_iter) result(result)
      real(dp), intent(in) :: data(:)
      character(len=*), intent(in), optional :: family
      integer, intent(in), optional :: max_iter
      type(fit_result) :: result
      real(dp), allocatable :: start(:), optimum(:), hessian(:,:), column(:), solution(:)
      real(dp) :: objective_value, m, sd
      integer :: fam, p, i, maxit
      logical :: converged, inv_ok
      character(len=:), allocatable :: key

      if (size(data) < 6 .or. .not. all(ieee_is_finite(data))) then
         result%message = 'at least six finite observations are required'
         return
      end if
      key = 'ghyp'
      if (present(family)) key = lower_string(trim(family))
      select case(key)
      case('gaussian','gauss','normal')
         result = fit_gaussian_uv(data); return
      case('hyp','hyperbolic'); fam=model_hyp; p=5
      case('nig'); fam=model_nig; p=5
      case('student','t','student-t'); fam=model_student; p=5
      case('vg','variance-gamma'); fam=model_vg; p=5
      case default; fam=model_ghyp; p=6
      end select
      m = sample_mean(data)
      sd = sqrt(max(sample_variance(data),1.0e-8_dp))
      allocate(start(p),active_uv_data(size(data)))
      active_uv_data = data; active_family = fam
      select case(fam)
      case(model_ghyp); start=[0.5_dp,log(1.0_dp),log(2.0_dp),m,log(sd),0.0_dp]
      case(model_hyp); start=[log(1.0_dp),log(2.0_dp),m,log(sd),0.0_dp]
      case(model_nig); start=[log(2.0_dp),log(2.0_dp),m,log(sd),0.0_dp]
      case(model_student); start=[log(6.0_dp),log(6.0_dp),m,log(sd),0.0_dp]
      case(model_vg); start=[log(1.0_dp),log(2.0_dp),m,log(sd),0.0_dp]
      end select
      maxit = 1000; if (present(max_iter)) maxit=max_iter
      call nelder_mead(uv_objective,start,optimum,objective_value,converged, &
         result%iterations,maxit,1.0e-7_dp)
      result%model = decode_uv(optimum,fam)
      result%theta = optimum
      result%n_parameters = p
      result%log_likelihood = -objective_value
      result%aic = 2.0_dp*real(p,dp)+2.0_dp*objective_value
      result%bic = log(real(size(data),dp))*real(p,dp)+2.0_dp*objective_value
      result%converged = converged
      result%ok = result%model%ok .and. ieee_is_finite(objective_value)
      if (.not. result%ok) then
         result%message = 'optimization failed'
         return
      end if
      call numerical_hessian(uv_objective,optimum,hessian)
      allocate(result%covariance(p,p),result%standard_error(p),column(p))
      result%covariance = 0.0_dp; inv_ok=.true.
      do i=1,p
         column=0.0_dp; column(i)=1.0_dp
         call solve_linear(hessian,column,solution,inv_ok)
         if (.not. inv_ok) exit
         result%covariance(:,i)=solution
      end do
      if (inv_ok) then
         call symmetrize(result%covariance)
         do i=1,p
            result%standard_error(i)=sqrt(max(0.0_dp,result%covariance(i,i)))
         end do
      else
         result%covariance=0.0_dp; result%standard_error=0.0_dp
      end if
   end function fit_ghyp_uv

   function fit_gaussian_mv(data) result(result)
      real(dp), intent(in) :: data(:,:)
      type(fit_result) :: result
      real(dp), allocatable :: mean(:), covariance(:,:)
      integer :: i, n, d
      n=size(data,1); d=size(data,2)
      if (n < 2 .or. d < 1 .or. .not. all(ieee_is_finite(data))) then
         result%message='invalid multivariate data'; return
      end if
      allocate(mean(d),covariance(d,d))
      mean=sum(data,dim=1)/real(n,dp)
      covariance=0.0_dp
      do i=1,n
         covariance=covariance+spread(data(i,:)-mean,2,d)*spread(data(i,:)-mean,1,d)
      end do
      covariance=covariance/real(n,dp)
      result%model=gaussian_mv(mean,covariance)
      result%n_parameters=d+d*(d+1)/2
      result%log_likelihood=sum([(log_dghyp(data(i,:),result%model),i=1,n)])
      result%aic=-2.0_dp*result%log_likelihood+2.0_dp*real(result%n_parameters,dp)
      result%bic=-2.0_dp*result%log_likelihood+log(real(n,dp))*real(result%n_parameters,dp)
      result%converged=.true.; result%ok=result%model%ok
   end function fit_gaussian_mv

   function decode_mv(theta, family, d) result(model)
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: family, d
      type(ghyp_model_type) :: model
      real(dp), allocatable :: mu(:), gamma(:), l(:,:), scatter(:,:)
      real(dp) :: lambda,chi,psi,nu
      integer :: k,i,j,pos
      select case(family)
      case(model_ghyp); k=3; lambda=theta(1);chi=exp(theta(2));psi=exp(theta(3))
      case(model_hyp); k=2;lambda=0.5_dp*real(d+1,dp);chi=exp(theta(1));psi=exp(theta(2))
      case(model_nig); k=2;lambda=-0.5_dp;chi=exp(theta(1));psi=exp(theta(2))
      case(model_student);k=2;nu=2.0_dp+exp(theta(1));lambda=-0.5_dp*nu;chi=exp(theta(2));psi=0.0_dp
      case(model_vg);k=2;lambda=exp(theta(1));chi=0.0_dp;psi=exp(theta(2))
      case default;model%message='unknown family';return
      end select
      allocate(mu(d),gamma(d),l(d,d),scatter(d,d))
      mu=theta(k+1:k+d); pos=k+d
      l=0.0_dp
      do i=1,d
         do j=1,i
            pos=pos+1
            if (i==j) then
               l(i,j)=exp(theta(pos))
            else
               l(i,j)=theta(pos)
            end if
         end do
      end do
      gamma=theta(pos+1:pos+d)
      scatter=matmul(l,transpose(l))
      model=make_ghyp(lambda,chi,psi,mu,scatter,gamma)
   end function decode_mv

   function mv_objective(theta) result(value)
      real(dp), intent(in) :: theta(:)
      real(dp) :: value,lv
      type(ghyp_model_type)::model
      integer::i
      model=decode_mv(theta,active_family,active_dimension)
      if(.not.model%ok)then;value=huge(1.0_dp)/100.0_dp;return;end if
      value=0.0_dp
      do i=1,size(active_mv_data,1)
         lv=log_dghyp(active_mv_data(i,:),model)
         if(.not.ieee_is_finite(lv))then;value=huge(1.0_dp)/100.0_dp;return;end if
         value=value-lv
      end do
   end function mv_objective

   function fit_ghyp_mv(data,family,max_iter) result(result)
      real(dp),intent(in)::data(:,:)
      character(len=*),intent(in),optional::family
      integer,intent(in),optional::max_iter
      type(fit_result)::result
      type(fit_result)::gauss_fit
      real(dp),allocatable::start(:),optimum(:),l(:,:)
      real(dp)::obj
      integer::d,n,fam,k,p,i,j,pos,maxit
      logical::conv,chol_ok
      character(len=:),allocatable::key
      n=size(data,1);d=size(data,2)
      if(n<max(8,d+2).or.d<1.or..not.all(ieee_is_finite(data)))then
         result%message='insufficient finite multivariate data';return
      end if
      key='ghyp';if(present(family))key=lower_string(trim(family))
      if(key=='gaussian'.or.key=='gauss'.or.key=='normal')then
         result=fit_gaussian_mv(data);return
      end if
      select case(key)
      case('hyp','hyperbolic');fam=model_hyp;k=2
      case('nig');fam=model_nig;k=2
      case('student','t','student-t');fam=model_student;k=2
      case('vg','variance-gamma');fam=model_vg;k=2
      case default;fam=model_ghyp;k=3
      end select
      p=k+d+d*(d+1)/2+d
      allocate(start(p),active_mv_data(n,d));active_mv_data=data
      active_family=fam;active_dimension=d
      gauss_fit=fit_gaussian_mv(data)
      select case(fam)
      case(model_ghyp);start(1:3)=[0.5_dp,log(1.0_dp),log(2.0_dp)]
      case(model_hyp);start(1:2)=[log(1.0_dp),log(2.0_dp)]
      case(model_nig);start(1:2)=[log(2.0_dp),log(2.0_dp)]
      case(model_student);start(1:2)=[log(6.0_dp),log(6.0_dp)]
      case(model_vg);start(1:2)=[log(1.0_dp),log(2.0_dp)]
      end select
      start(k+1:k+d)=gauss_fit%model%mu
      call cholesky_lower(gauss_fit%model%scatter,l,chol_ok)
      pos=k+d
      do i=1,d
         do j=1,i
            pos=pos+1
            if(i==j)then;start(pos)=log(l(i,j));else;start(pos)=l(i,j);end if
         end do
      end do
      start(pos+1:pos+d)=0.0_dp
      maxit=1500;if(present(max_iter))maxit=max_iter
      call nelder_mead(mv_objective,start,optimum,obj,conv,result%iterations,maxit,1.0e-6_dp)
      result%model=decode_mv(optimum,fam,d);result%theta=optimum
      result%n_parameters=p;result%log_likelihood=-obj
      result%aic=2.0_dp*real(p,dp)+2.0_dp*obj
      result%bic=log(real(n,dp))*real(p,dp)+2.0_dp*obj
      result%converged=conv;result%ok=result%model%ok.and.ieee_is_finite(obj)
      if(.not.result%ok)result%message='optimization failed'
   end function fit_ghyp_mv

   function likelihood_ratio_test(full_fit,restricted_fit) result(result)
      type(fit_result),intent(in)::full_fit,restricted_fit
      type(likelihood_ratio_result)::result
      real(dp)::x
      if(.not.full_fit%ok.or..not.restricted_fit%ok)then
         result%message='both fits must be valid';return
      end if
      result%degrees_freedom=full_fit%n_parameters-restricted_fit%n_parameters
      if(result%degrees_freedom<1)then;result%message='nonpositive degrees of freedom';return;end if
      result%statistic=max(0.0_dp,2.0_dp*(full_fit%log_likelihood-restricted_fit%log_likelihood))
      x=0.5_dp*result%statistic
      result%p_value=gamma_upper_regularized(0.5_dp*real(result%degrees_freedom,dp),x)
      result%ok=.true.
   end function likelihood_ratio_test

   function step_aic_ghyp(data) result(result)
      real(dp),intent(in)::data(:)
      type(model_selection_result)::result
      character(len=16),parameter::families(6)=[character(len=16):: &
         'ghyp','hyp','nig','student','vg','gaussian']
      type(fit_result)::fit
      integer::i,best_i
      allocate(result%aic(6),result%names(6));result%names=families
      result%aic=huge(1.0_dp);best_i=1
      do i=1,6
         fit=fit_ghyp_uv(data,trim(families(i)),max_iter=500)
         if(fit%ok)result%aic(i)=fit%aic
         if(result%aic(i)<result%aic(best_i))best_i=i
      end do
      result%best=fit_ghyp_uv(data,trim(families(best_i)),max_iter=1000)
      result%family=families(best_i);result%ok=result%best%ok
   end function step_aic_ghyp

   pure function lower_string(s) result(out)
      character(len=*),intent(in)::s
      character(len=len(s))::out
      integer::i,c
      do i=1,len(s)
         c=iachar(s(i:i));if(c>=iachar('A').and.c<=iachar('Z'))then
            out(i:i)=achar(c+32)
         else;out(i:i)=s(i:i);end if
      end do
   end function lower_string

   pure function gamma_upper_regularized(a,x) result(q)
      real(dp),intent(in)::a,x
      real(dp)::q,ap,del,sumv,b,c,d,h,an
      integer::n
      if(x<=0.0_dp)then;q=1.0_dp;return;end if
      if(x<a+1.0_dp)then
         ap=a;del=1.0_dp/a;sumv=del
         do n=1,300
            ap=ap+1.0_dp;del=del*x/ap;sumv=sumv+del
            if(abs(del)<abs(sumv)*1.0e-14_dp)exit
         end do
         q=1.0_dp-sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b=x+1.0_dp-a;c=1.0e300_dp;d=1.0_dp/b;h=d
         do n=1,300
            an=-real(n,dp)*(real(n,dp)-a);b=b+2.0_dp
            d=an*d+b;if(abs(d)<1.0e-300_dp)d=1.0e-300_dp
            c=b+an/c;if(abs(c)<1.0e-300_dp)c=1.0e-300_dp
            d=1.0_dp/d;del=d*c;h=h*del
            if(abs(del-1.0_dp)<1.0e-14_dp)exit
         end do
         q=exp(-x+a*log(x)-log_gamma(a))*h
      end if
      q=min(1.0_dp,max(0.0_dp,q))
   end function gamma_upper_regularized

end module ghyp_fitting
