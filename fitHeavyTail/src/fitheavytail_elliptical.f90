! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_elliptical
   use fitheavytail_kinds, only: dp
   use fitheavytail_status, only: ht_success, ht_invalid_argument, &
      ht_too_few_observations, ht_singular_matrix, ht_no_convergence
   use fitheavytail_types, only: heavy_tail_fit, clear_fit
   use fitheavytail_linalg, only: clean_complete_rows, column_mean, sample_covariance, &
      spatial_median, trace_matrix, quadratic_forms, weighted_covariance, logdet_spd
   use fitheavytail_tail, only: cap_nu, default_nu_min, default_nu_max
   implicit none
   private
   public :: fit_tyler, fit_cauchy, recover_scaled_scatter_and_nu

contains

   subroutine fit_tyler(x, result, estimate_mu, initial_mu, initial_covariance, &
      max_iter, ptol, ftol, nu_min, nu_max)
      real(dp), intent(in) :: x(:,:)
      type(heavy_tail_fit), intent(out) :: result
      logical, intent(in), optional :: estimate_mu
      real(dp), intent(in), optional :: initial_mu(:), initial_covariance(:,:)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: ptol, ftol, nu_min, nu_max
      real(dp), allocatable :: data(:,:), xc(:,:), sigma(:,:), old_sigma(:,:), q(:), weights(:)
      real(dp), allocatable :: mu(:), covariance(:,:)
      real(dp) :: p_tol, f_tol, ll, old_ll, start_time, finish_time, tr, lo, hi
      integer :: niter, iter, t, n, status
      logical :: do_mu, param_ok, fun_ok, converged_flag

      call clear_fit(result)
      call clean_complete_rows(x,data)
      t=size(data,1)
      n=size(data,2)
      if(n<1 .or. t<=n) then
         call set_error(result,ht_too_few_observations,'fit_tyler requires T > N after dropping NaN rows')
         return
      end if
      niter=200
      if(present(max_iter)) niter=max_iter
      p_tol=1.0e-3_dp
      if(present(ptol)) p_tol=ptol
      f_tol=huge(1.0_dp)
      if(present(ftol)) f_tol=ftol
      lo=default_nu_min
      hi=default_nu_max
      if(present(nu_min)) lo=nu_min
      if(present(nu_max)) hi=nu_max
      do_mu=.true.
      if(present(estimate_mu)) do_mu=estimate_mu
      if(niter<1) then
         call set_error(result,ht_invalid_argument,'max_iter must be positive')
         return
      end if
      allocate(mu(n),sigma(n,n),old_sigma(n,n),xc(t,n),q(t),weights(t),covariance(n,n))
      if(do_mu) then
         if(present(initial_mu)) then
            if(size(initial_mu)/=n) then
               call set_error(result,ht_invalid_argument,'initial_mu has the wrong size')
               return
            end if
            mu=initial_mu
         else
            call spatial_median(data,mu)
         end if
      else
         mu=0.0_dp
      end if
      if(present(initial_covariance)) then
         if(any(shape(initial_covariance)/=[n,n])) then
            call set_error(result,ht_invalid_argument,'initial_covariance has the wrong shape')
            return
         end if
         sigma=initial_covariance
      else
         sigma=sample_covariance(data)
      end if
      tr=trace_matrix(sigma)
      if(tr<=0.0_dp) then
         call set_error(result,ht_singular_matrix,'initial covariance is singular')
         return
      end if
      sigma=sigma/tr
      call center_rows(data,mu,xc)
      call quadratic_forms(xc,sigma,q,status)
      if(status/=ht_success .or. any(q<=0.0_dp)) then
         call set_error(result,ht_singular_matrix,'failed to invert initial covariance')
         return
      end if
      weights=1.0_dp/q
      if(f_tol<huge(1.0_dp)/2.0_dp) then
         ll=0.5_dp*real(n,dp)*sum(log(weights))-0.5_dp*real(t,dp)*logdet_spd(sigma,status)
      else
         ll=0.0_dp
      end if
      converged_flag = .false.
      call cpu_time(start_time)
      do iter=1,niter
         old_sigma=sigma
         old_ll=ll
         sigma=real(n,dp)*weighted_covariance(xc,weights,real(t,dp))
         tr=trace_matrix(sigma)
         if(tr<=0.0_dp) exit
         sigma=sigma/tr
         call quadratic_forms(xc,sigma,q,status)
         if(status/=ht_success .or. any(q<=0.0_dp)) exit
         weights=1.0_dp/q
         param_ok=relative_matrix_converged(sigma,old_sigma,p_tol)
         if(f_tol<huge(1.0_dp)/2.0_dp) then
            ll=0.5_dp*real(n,dp)*sum(log(weights))-0.5_dp*real(t,dp)*logdet_spd(sigma,status)
            fun_ok=relative_scalar_converged(ll,old_ll,f_tol)
         else
            fun_ok=.true.
         end if
         if (param_ok .and. fun_ok) then
            converged_flag = .true.
            exit
         end if
      end do
      call cpu_time(finish_time)
      if(status/=ht_success) then
         call set_error(result,ht_singular_matrix,'numerical failure during Tyler iterations')
         return
      end if
      call recover_scaled_scatter_and_nu(sigma,xc,result%nu,lo,hi,status)
      if(status/=ht_success) then
         call set_error(result,status,'could not recover Tyler scatter scale')
         return
      end if
      covariance=result%nu/(result%nu-2.0_dp)*sigma
      result%mu=mu
      result%scatter=sigma
      result%covariance=covariance
      result%mean=mu
      result%num_iterations = min(iter,niter)
      result%converged = converged_flag
      result%cpu_time=finish_time-start_time
      result%log_likelihood=ll
      if(result%converged) then
         result%status=ht_success
         result%message='success'
      else
         result%status=ht_no_convergence
         result%message='maximum iterations reached'
      end if
   end subroutine fit_tyler

   subroutine fit_cauchy(x, result, initial_mu, initial_covariance, max_iter, &
      ptol, ftol, nu_min, nu_max)
      real(dp), intent(in) :: x(:,:)
      type(heavy_tail_fit), intent(out) :: result
      real(dp), intent(in), optional :: initial_mu(:), initial_covariance(:,:)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: ptol, ftol, nu_min, nu_max
      real(dp), allocatable :: data(:,:), xc(:,:), sigma(:,:), old_sigma(:,:), q(:), weights(:)
      real(dp), allocatable :: mu(:), old_mu(:), covariance(:,:)
      real(dp) :: p_tol,f_tol,ll,old_ll,start_time,finish_time,beta,lo,hi
      integer :: niter,iter,t,n,status
      logical :: param_ok,fun_ok,converged_flag

      call clear_fit(result)
      call clean_complete_rows(x,data)
      t=size(data,1)
      n=size(data,2)
      if(n<1.or.t<=n) then
         call set_error(result,ht_too_few_observations,'fit_cauchy requires T > N after dropping NaN rows')
         return
      end if
      niter=200
      if(present(max_iter)) niter=max_iter
      p_tol=1.0e-3_dp
      if(present(ptol)) p_tol=ptol
      f_tol=huge(1.0_dp)
      if(present(ftol)) f_tol=ftol
      lo=default_nu_min
      hi=default_nu_max
      if(present(nu_min)) lo=nu_min
      if(present(nu_max)) hi=nu_max
      allocate(mu(n),old_mu(n),sigma(n,n),old_sigma(n,n),xc(t,n),q(t),weights(t),covariance(n,n))
      if(present(initial_mu)) then
         if(size(initial_mu)/=n) then
         call set_error(result,ht_invalid_argument,'initial_mu has the wrong size')
         return
         end if
         mu=initial_mu
      else
         mu=column_mean(data)
      end if
      if(present(initial_covariance)) then
         if(any(shape(initial_covariance)/=[n,n])) then
         call set_error(result,ht_invalid_argument,'initial covariance shape mismatch')
         return
         end if
         sigma=initial_covariance
      else
         sigma=sample_covariance(data)
      end if
      call center_rows(data,mu,xc)
      call quadratic_forms(xc,sigma,q,status)
      if(status/=ht_success) then
      call set_error(result,status,'initial covariance is singular')
      return
      end if
      weights=1.0_dp/(1.0_dp+q)
      if(f_tol<huge(1.0_dp)/2.0_dp) then
         ll=0.5_dp*real(n+1,dp)*sum(log(weights))-0.5_dp*real(t,dp)*logdet_spd(sigma,status)
      else
         ll=0.0_dp
      end if
      converged_flag = .false.
      call cpu_time(start_time)
      do iter=1,niter
         old_mu=mu
         old_sigma=sigma
         old_ll=ll
         mu=matmul(weights,data)/sum(weights)
         call center_rows(data,mu,xc)
         beta=real(t,dp)/(real(n+1,dp)*sum(weights))
         sigma=beta*real(n+1,dp)*weighted_covariance(xc,weights,real(t,dp))
         call quadratic_forms(xc,sigma,q,status)
         if(status/=ht_success) exit
         weights=1.0_dp/(1.0_dp+q)
         param_ok=relative_vector_converged(mu,old_mu,p_tol).and.relative_matrix_converged(sigma,old_sigma,p_tol)
         if(f_tol<huge(1.0_dp)/2.0_dp) then
            ll=0.5_dp*real(n+1,dp)*sum(log(weights))-0.5_dp*real(t,dp)*logdet_spd(sigma,status)
            fun_ok=relative_scalar_converged(ll,old_ll,f_tol)
         else
            fun_ok=.true.
         end if
         if (param_ok .and. fun_ok) then
            converged_flag = .true.
            exit
         end if
      end do
      call cpu_time(finish_time)
      if(status/=ht_success) then
      call set_error(result,status,'numerical failure during Cauchy iterations')
      return
      end if
      call recover_scaled_scatter_and_nu(sigma,xc,result%nu,lo,hi,status)
      covariance=result%nu/(result%nu-2.0_dp)*sigma
      result%mu=mu
      result%scatter=sigma
      result%covariance=covariance
      result%mean=mu
      result%num_iterations = min(iter,niter)
      result%converged = converged_flag
      result%cpu_time=finish_time-start_time
      result%log_likelihood=ll
      if(result%converged) then
      result%status=ht_success
      result%message='success'
      else
      result%status=ht_no_convergence
      result%message='maximum iterations reached'
      end if
   end subroutine fit_cauchy

   subroutine recover_scaled_scatter_and_nu(sigma,xc,nu,nu_min,nu_max,status)
      real(dp),intent(inout)::sigma(:,:)
      real(dp),intent(in)::xc(:,:)
      real(dp),intent(out)::nu
      real(dp),intent(in),optional::nu_min,nu_max
      integer,intent(out),optional::status
      real(dp),allocatable::q(:)
      real(dp)::eta_scatter,eta_cov,theta,tr
      integer::n,t,istat
      n=size(xc,2)
      t=size(xc,1)
      allocate(q(t))
      tr=trace_matrix(sigma)
      if(tr<=0.0_dp) then
      istat=ht_singular_matrix
      nu=default_nu_max
      else
         sigma=sigma*real(n,dp)/tr
         call quadratic_forms(xc,sigma,q,istat)
         if(istat==ht_success) then
            eta_scatter=1.0_dp/(sum(real(n,dp)/max(q,tiny(1.0_dp)))/real(t,dp))
            sigma=eta_scatter*sigma
            eta_cov=sum(sum(xc*xc,dim=1)/real(t-1,dp))/real(n,dp)
            theta=eta_cov/max(eta_scatter,tiny(1.0_dp))
            if(abs(theta-1.0_dp)<1.0e-12_dp) then
               nu=default_nu_max
            else
               nu=cap_nu(2.0_dp*theta/(theta-1.0_dp),nu_min,nu_max)
            end if
         else
            nu=default_nu_max
         end if
      end if
      if(present(status)) status=istat
   end subroutine recover_scaled_scatter_and_nu

   subroutine center_rows(x,mu,xc)
      real(dp),intent(in)::x(:,:),mu(:)
      real(dp),intent(out)::xc(:,:)
      integer::i
      do i=1,size(x,1)
      xc(i,:)=x(i,:)-mu
      end do
   end subroutine center_rows

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

end module fitheavytail_elliptical
