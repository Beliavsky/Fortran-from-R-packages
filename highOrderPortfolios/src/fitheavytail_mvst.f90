! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_mvst
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use fitheavytail_kinds, only: dp
   use fitheavytail_status, only: ht_success, ht_invalid_argument, &
      ht_too_few_observations, ht_singular_matrix, ht_no_convergence, &
      ht_numerical_error
   use fitheavytail_types, only: heavy_tail_fit, clear_fit
   use fitheavytail_linalg, only: clean_complete_rows, column_mean, &
      sample_covariance, inverse_matrix, quadratic_forms, &
      weighted_covariance, outer_product, symmetrize, logdet_spd
   use fitheavytail_special, only: digamma_dp, bessel_k_ratio, &
      bessel_order_derivative, log_bessel_k
   use fitheavytail_tail, only: default_nu_min, default_nu_max, cap_nu
   implicit none
   private
   public :: fit_mvst, mvst_log_likelihood, sample_skewness

contains

   subroutine fit_mvst(x,result,fixed_nu,fixed_gamma,initial_nu, &
      initial_mu,initial_gamma,initial_scatter,max_iter,ptol,ftol, &
      pxem,nu_min,nu_max)
      real(dp),intent(in)::x(:,:)
      type(heavy_tail_fit),intent(out)::result
      real(dp),intent(in),optional::fixed_nu,fixed_gamma(:),initial_nu
      real(dp),intent(in),optional::initial_mu(:),initial_gamma(:)
      real(dp),intent(in),optional::initial_scatter(:,:)
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::ptol,ftol,nu_min,nu_max
      logical,intent(in),optional::pxem
      real(dp),allocatable::data(:,:),mu(:),old_mu(:),gamma_u(:)
      real(dp),allocatable::old_gamma(:),scatter_u(:,:),old_scatter(:,:)
      real(dp),allocatable::scatter(:,:),covariance(:,:),mean_vec(:)
      real(dp),allocatable::etau(:),einvtau(:),elogtau(:),xc(:,:)
      real(dp)::nu,old_nu,alpha,old_alpha,lo,hi,p_tol,f_tol
      real(dp)::ll,old_ll,start_time,end_time,nan_value
      real(dp)::sum_etau,sum_einvtau
      real(dp)::sum_xc(size(x,2))
      integer::t,n,niter,iter,status
      logical::optimize_nu,optimize_gamma,do_px,params_ok,fun_ok
      logical::converged_flag

      call clear_fit(result)
      call clean_complete_rows(x,data)
      t=size(data,1)
      n=size(data,2)
      if(n<1.or.t<=n) then
         call set_error(result,ht_too_few_observations, &
            'fit_mvst requires T > N after dropping NaN rows')
         return
      end if
      lo=default_nu_min
      hi=default_nu_max
      if(present(nu_min)) lo=nu_min
      if(present(nu_max)) hi=nu_max
      niter=500
      if(present(max_iter)) niter=max_iter
      p_tol=1.0e-3_dp
      if(present(ptol)) p_tol=ptol
      f_tol=huge(1.0_dp)
      if(present(ftol)) f_tol=ftol
      do_px=.true.
      if(present(pxem)) do_px=pxem
      if(niter<1.or.lo<=2.0_dp.or.hi<=lo) then
         call set_error(result,ht_invalid_argument,'invalid fitting options')
         return
      end if
      optimize_nu=.not.present(fixed_nu)
      if(present(fixed_nu)) then
         if(fixed_nu<=2.0_dp) then
            call set_error(result,ht_invalid_argument,'fixed_nu must exceed two')
            return
         end if
         nu=fixed_nu
         if(nu>=huge(1.0_dp)/100.0_dp) nu=1.0e5_dp
      else
         nu=4.0_dp
         if(present(initial_nu)) nu=initial_nu
         nu=cap_nu(nu,lo,hi)
      end if
      optimize_gamma=.not.present(fixed_gamma)
      allocate(mu(n),old_mu(n),gamma_u(n),old_gamma(n), &
         scatter_u(n,n),old_scatter(n,n),scatter(n,n), &
         covariance(n,n),mean_vec(n),etau(t),einvtau(t), &
         elogtau(t),xc(t,n))
      alpha=1.0_dp
      if(present(fixed_gamma)) then
         if(size(fixed_gamma)/=n) then
            call set_error(result,ht_invalid_argument, &
               'fixed_gamma has the wrong size')
            return
         end if
         gamma_u=fixed_gamma
      else if(present(initial_gamma)) then
         if(size(initial_gamma)/=n) then
            call set_error(result,ht_invalid_argument, &
               'initial_gamma has the wrong size')
            return
         end if
         gamma_u=initial_gamma
      else
         gamma_u=0.0_dp
      end if
      if(present(initial_mu)) then
         if(size(initial_mu)/=n) then
            call set_error(result,ht_invalid_argument, &
               'initial_mu has the wrong size')
            return
         end if
         mu=initial_mu
      else
         mu=column_mean(data)-max(nu,2.1_dp)/ &
            (max(nu,2.1_dp)-2.0_dp)*gamma_u/alpha
      end if
      if(present(initial_scatter)) then
         if(any(shape(initial_scatter)/=[n,n])) then
            call set_error(result,ht_invalid_argument, &
               'initial_scatter has the wrong shape')
            return
         end if
         scatter_u=initial_scatter
      else
         scatter_u=(max(nu,2.1_dp)-2.0_dp)/max(nu,2.1_dp) * &
            sample_covariance(data)
      end if
      ll=0.0_dp
      if(f_tol<huge(1.0_dp)/2.0_dp) then
         ll=mvst_log_likelihood(data,nu,gamma_u/alpha,mu, &
            scatter_u/alpha,status)
         if(status/=ht_success) then
            call set_error(result,status,'initial skew-t likelihood failed')
            return
         end if
      end if
      converged_flag = .false.
      call cpu_time(start_time)
      status=ht_success
      do iter=1,niter
         old_mu=mu
         old_gamma=gamma_u/alpha
         old_scatter=scatter_u/alpha
         old_nu=nu
         old_alpha=alpha
         old_ll=ll
         call estep_mvst(data,nu,gamma_u,mu,scatter_u,alpha, &
            etau,einvtau,elogtau,status)
         if(status/=ht_success) exit
         if(optimize_nu) then
            nu=optimize_mvst_nu(etau,elogtau,alpha,lo,hi)
         end if
         sum_etau=sum(etau)
         sum_einvtau=sum(einvtau)
         mu=(matmul(etau,data)-real(t,dp)*gamma_u)/sum_etau
         if(optimize_gamma) then
            gamma_u=(sum(data,dim=1)-real(t,dp)*mu)/sum_einvtau
         end if
         call center_rows(data,mu,xc)
         sum_xc=sum(xc,dim=1)
         scatter_u=weighted_covariance(xc,etau,1.0_dp) - &
            2.0_dp*outer_product(sum_xc,gamma_u) + &
            sum_einvtau*outer_product(gamma_u,gamma_u)
         scatter_u=symmetrize(scatter_u/real(t,dp))
         if(do_px) alpha=sum_etau/real(t,dp)
         if(.not.all_finite(scatter_u).or.alpha<=0.0_dp) then
            status=ht_numerical_error
            exit
         end if
         if(f_tol<huge(1.0_dp)/2.0_dp) then
            ll=mvst_log_likelihood(data,nu,gamma_u/alpha,mu, &
               scatter_u/alpha,status)
            if(status/=ht_success) exit
            fun_ok=relative_scalar_converged(ll,old_ll,f_tol)
         else
            fun_ok=.true.
         end if
         params_ok=relative_scalar_converged(fnu(nu),fnu(old_nu),p_tol) .and. &
            relative_vector_converged(mu,old_mu,p_tol) .and. &
            relative_vector_converged(gamma_u/alpha,old_gamma,p_tol) .and. &
            relative_matrix_converged(scatter_u/alpha,old_scatter,p_tol)
         if (params_ok .and. fun_ok) then
            converged_flag = .true.
            exit
         end if
         if(abs(alpha-old_alpha)>huge(1.0_dp)) exit
      end do
      call cpu_time(end_time)
      if(status/=ht_success) then
         call set_error(result,status,'numerical failure in fit_mvst')
         return
      end if
      scatter=scatter_u/alpha
      result%mu=mu
      result%gamma=gamma_u/alpha
      result%scatter=scatter
      result%nu=nu
      nan_value=ieee_value(0.0_dp,ieee_quiet_nan)
      if(nu>2.0_dp) then
         mean_vec=mu+nu/(nu-2.0_dp)*result%gamma
      else
         mean_vec=nan_value
      end if
      if(nu>4.0_dp) then
         covariance=nu/(nu-2.0_dp)*scatter + &
            2.0_dp*nu*nu/((nu-2.0_dp)**2*(nu-4.0_dp)) * &
            outer_product(result%gamma,result%gamma)
      else
         covariance=nan_value
      end if
      result%mean=mean_vec
      result%covariance=covariance
      result%latent_weights=etau
      result%num_iterations = min(iter,niter)
      result%converged = converged_flag
      result%cpu_time=end_time-start_time
      result%log_likelihood=ll
      if(result%converged) then
         result%status=ht_success
         result%message='success'
      else
         result%status=ht_no_convergence
         result%message='maximum iterations reached'
      end if
   end subroutine fit_mvst

   subroutine estep_mvst(x,nu,gamma_u,mu,scatter_u,alpha, &
      etau,einvtau,elogtau,status)
      real(dp),intent(in)::x(:,:),nu,gamma_u(:),mu(:)
      real(dp),intent(in)::scatter_u(:,:),alpha
      real(dp),intent(out)::etau(:),einvtau(:),elogtau(:)
      integer,intent(out)::status
      real(dp),allocatable::scatter(:,:),invs(:,:),xc(:,:),q(:)
      real(dp)::gamma(size(gamma_u)),delta,lambda,kappa,ratio,arg
      integer::i,n,t
      n=size(x,2)
      t=size(x,1)
      status=ht_success
      allocate(scatter(n,n),invs(n,n),xc(t,n),q(t))
      gamma=gamma_u/alpha
      scatter=scatter_u/alpha
      call inverse_matrix(scatter,invs,status)
      if(status/=ht_success) return
      delta=sqrt(max(0.0_dp,dot_product(gamma,matmul(invs,gamma))))
      call center_rows(x,mu,xc)
      do i=1,t
         q(i)=dot_product(xc(i,:),matmul(invs,xc(i,:)))
      end do
      lambda=0.5_dp*(nu+real(n,dp))
      do i=1,t
         kappa=sqrt(max(nu+q(i),tiny(1.0_dp)))
         if(delta<=1.0e-10_dp) then
            etau(i)=(nu+real(n,dp))/(kappa*kappa)
            einvtau(i)=kappa*kappa/(nu+real(n,dp)-2.0_dp)
            elogtau(i)=digamma_dp(lambda)-log(0.5_dp*kappa*kappa)
         else
            arg=delta*kappa
            ratio=bessel_k_ratio(arg,lambda)
            etau(i)=delta/kappa*ratio
            einvtau(i)=kappa/delta * &
               (ratio-2.0_dp*lambda/max(arg,tiny(1.0_dp)))
            elogtau(i)=log(delta/kappa) + &
               bessel_order_derivative(arg,lambda)
            if(einvtau(i)<=0.0_dp) then
               einvtau(i)=kappa/delta / &
                  max(bessel_k_ratio(arg,lambda-1.0_dp),tiny(1.0_dp))
            end if
         end if
      end do
      etau=etau*alpha
      einvtau=einvtau/alpha
      elogtau=elogtau+log(alpha)
      if(.not.all_finite_vector(etau).or. &
         .not.all_finite_vector(einvtau).or. &
         .not.all_finite_vector(elogtau)) status=ht_numerical_error
   end subroutine estep_mvst

   function optimize_mvst_nu(etau,elogtau,alpha,lo,hi) result(root)
      real(dp),intent(in)::etau(:),elogtau(:),alpha,lo,hi
      real(dp)::root,a,b,c,d,fc,fd,summary
      real(dp),parameter::gr=0.6180339887498948482_dp
      integer::iter,t
      t=size(etau)
      summary=sum(elogtau-log(alpha)-etau/alpha)
      a=lo
      b=hi
      c=b-gr*(b-a)
      d=a+gr*(b-a)
      fc=-mvst_qnu(c,summary,t)
      fd=-mvst_qnu(d,summary,t)
      do iter=1,120
         if(abs(b-a)<=1.0e-8_dp*(1.0_dp+abs(a)+abs(b))) exit
         if(fc<fd) then
            b=d
            d=c
            fd=fc
            c=b-gr*(b-a)
            fc=-mvst_qnu(c,summary,t)
         else
            a=c
            c=d
            fc=fd
            d=a+gr*(b-a)
            fd=-mvst_qnu(d,summary,t)
         end if
      end do
      root=0.5_dp*(a+b)
   end function optimize_mvst_nu

   pure function mvst_qnu(nu,summary,t) result(value)
      real(dp),intent(in)::nu,summary
      integer,intent(in)::t
      real(dp)::value
      value=0.5_dp*nu*summary + real(t,dp) * &
         (0.5_dp*nu*log(0.5_dp*nu)-log_gamma(0.5_dp*nu))
   end function mvst_qnu

   function mvst_log_likelihood(x,nu,gamma,mu,scatter,status) result(value)
      real(dp),intent(in)::x(:,:),nu,gamma(:),mu(:),scatter(:,:)
      integer,intent(out),optional::status
      real(dp)::value
      real(dp),allocatable::invs(:,:),xc(:,:),q(:)
      real(dp)::delta,kappa,first,second,third,logdet,pi
      integer::i,n,t,istat
      n=size(x,2)
      t=size(x,1)
      pi=acos(-1.0_dp)
      allocate(invs(n,n),xc(t,n),q(t))
      call inverse_matrix(scatter,invs,istat)
      if(istat/=ht_success) then
         value=-huge(1.0_dp)
         if(present(status)) status=istat
         return
      end if
      logdet=logdet_spd(scatter,istat)
      call center_rows(x,mu,xc)
      do i=1,t
      q(i)=dot_product(xc(i,:),matmul(invs,xc(i,:)))
      end do
      if(nu>=50000.0_dp) then
         value=-0.5_dp*real(t,dp)*(real(n,dp)*log(2.0_dp*pi)+logdet) - &
            0.5_dp*sum(q)
      else
         delta=sqrt(max(0.0_dp,dot_product(gamma,matmul(invs,gamma))))
         if(delta<=1.0e-12_dp) then
            first=log_gamma(0.5_dp*(nu+real(n,dp))) - &
               log_gamma(0.5_dp*nu) - 0.5_dp*real(n,dp)*log(nu*pi) - &
               0.5_dp*logdet
            value=real(t,dp)*first - &
               0.5_dp*(nu+real(n,dp))*sum(log(1.0_dp+q/nu))
         else
            second=log(2.0_dp)+0.5_dp*nu*log(0.5_dp*nu) - &
               log_gamma(0.5_dp*nu)
            value=0.0_dp
            do i=1,t
               kappa=sqrt(nu+q(i))
               first=dot_product(xc(i,:),matmul(invs,gamma)) - &
                  0.5_dp*real(n,dp)*log(2.0_dp*pi) - 0.5_dp*logdet
               third=-0.5_dp*(nu+real(n,dp))*log(kappa/delta) + &
                  log_bessel_k(delta*kappa,-0.5_dp*(nu+real(n,dp)))
               value=value+first+second+third
            end do
         end if
      end if
      if(present(status)) status=istat
   end function mvst_log_likelihood

   function sample_skewness(x) result(skew)
      real(dp),intent(in)::x(:,:)
      real(dp)::skew(size(x,2)),mu(size(x,2)),centered(size(x,1),size(x,2))
      integer::i
      mu=column_mean(x)
      do i=1,size(x,1)
      centered(i,:)=x(i,:)-mu
      end do
      skew=(sum(centered**3,dim=1)/real(size(x,1),dp)) / &
         max(sum(centered**2,dim=1)/real(size(x,1),dp),tiny(1.0_dp))**1.5_dp
   end function sample_skewness

   subroutine center_rows(x,mu,xc)
      real(dp),intent(in)::x(:,:),mu(:)
      real(dp),intent(out)::xc(:,:)
      integer::i
      do i=1,size(x,1)
      xc(i,:)=x(i,:)-mu
      end do
   end subroutine center_rows

   function all_finite(a) result(ok)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      real(dp),intent(in)::a(:,:)
      logical::ok
      integer::i,j
      ok=.true.
      do j=1,size(a,2)
      do i=1,size(a,1)
         if(.not.ieee_is_finite(a(i,j))) then
         ok=.false.
         return
         end if
      end do
      end do
   end function all_finite

   function all_finite_vector(a) result(ok)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      real(dp),intent(in)::a(:)
      logical::ok
      integer::i
      ok=.true.
      do i=1,size(a)
         if(.not.ieee_is_finite(a(i))) then
         ok=.false.
         return
         end if
      end do
   end function all_finite_vector

   pure function fnu(nu) result(value)
      real(dp),intent(in)::nu
      real(dp)::value
      value=nu/(nu-2.0_dp)
   end function fnu

   pure function relative_vector_converged(a,b,tol) result(ok)
      real(dp),intent(in)::a(:),b(:),tol
      logical::ok
      ok=all(abs(a-b)<=0.5_dp*tol*(abs(a)+abs(b)+tiny(1.0_dp)))
   end function relative_vector_converged

   pure function relative_matrix_converged(a,b,tol) result(ok)
      real(dp),intent(in)::a(:,:),b(:,:),tol
      logical::ok
      ok=all(abs(a-b)<=0.5_dp*tol*(abs(a)+abs(b)+tiny(1.0_dp)))
   end function relative_matrix_converged

   pure function relative_scalar_converged(a,b,tol) result(ok)
      real(dp),intent(in)::a,b,tol
      logical::ok
      ok=abs(a-b)<=0.5_dp*tol*(abs(a)+abs(b)+tiny(1.0_dp))
   end function relative_scalar_converged

   subroutine set_error(result,status,message)
      type(heavy_tail_fit),intent(inout)::result
      integer,intent(in)::status
      character(len=*),intent(in)::message
      result%status=status
      result%message=message
      result%converged=.false.
   end subroutine set_error

end module fitheavytail_mvst
