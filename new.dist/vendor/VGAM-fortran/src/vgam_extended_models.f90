! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_extended_models
   use vgam_kinds, only : dp
   use vgam_links, only : link_inverse, link_logit
   use vgam_linalg, only : invert_matrix
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   implicit none
   private

   type, public :: ordinal_result_t
      real(dp), allocatable :: cutpoints(:)
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_probabilities(:, :)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: link = link_logit
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict_proba => predict_ordinal
   end type ordinal_result_t

   type, public :: beta_regression_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp) :: precision = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_beta_regression
   end type beta_regression_result_t

   type, public :: negative_binomial_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp) :: size = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_negative_binomial
   end type negative_binomial_result_t

   type, public :: zip_result_t
      real(dp), allocatable :: count_coefficients(:)
      real(dp), allocatable :: zero_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_mean(:)
      real(dp), allocatable :: zero_probability(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_zip
   end type zip_result_t

   public :: fit_ordinal, ordinal_probabilities
   public :: fit_beta_regression, fit_negative_binomial
   public :: fit_zero_inflated_poisson

contains

   subroutine fit_ordinal(y, x, ncat, result, link_id, weights, max_iter, tol)
      integer, intent(in) :: y(:), ncat
      real(dp), intent(in) :: x(:, :)
      type(ordinal_result_t), intent(out) :: result
      integer, intent(in), optional :: link_id, max_iter
      real(dp), intent(in), optional :: weights(:), tol
      real(dp), allocatable :: w(:), par(:), hess(:,:), covu(:,:), jac(:,:), covn(:,:)
      real(dp), allocatable :: cuts(:), beta(:), prob(:,:)
      real(dp) :: fval, tolerance
      integer :: n, p, nc, i, j, stat, stat2, niter

      n = size(y)
      p = size(x,2)
      nc = ncat - 1
      if (n <= 0 .or. p <= 0 .or. size(x,1) /= n .or. ncat < 2 .or. &
          any(y < 1) .or. any(y > ncat)) then
         result%status = 1
         return
      end if
      result%link = link_logit
      if (present(link_id)) result%link = link_id
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         w = weights
      else
         allocate(w(n))
         w = 1.0_dp
      end if
      niter = 250
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      allocate(par(nc+p))
      call ordinal_initial(y, ncat, result%link, par(1:nc))
      par(nc+1:) = 0.0_dp
      call bfgs_minimize(objective, par, fval, stat, niter, tolerance)
      result%status = stat
      result%converged = stat == 0
      call unpack_ordinal(par, nc, cuts, beta)
      call ordinal_probabilities(x, cuts, beta, result%link, prob)
      result%cutpoints = cuts
      result%coefficients = beta
      result%fitted_probabilities = prob
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(size(par),dp)

      allocate(hess(size(par),size(par)), covu(size(par),size(par)))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, covu, stat2)
      if (stat2 == 0) then
         allocate(jac(size(par),size(par)), covn(size(par),size(par)))
         jac = 0.0_dp
         do i = 1, nc
            jac(i,1) = 1.0_dp
            do j = 2, i
               jac(i,j) = exp(par(j))
            end do
         end do
         do i = 1, p
            jac(nc+i,nc+i) = 1.0_dp
         end do
         covn = matmul(jac, matmul(covu, transpose(jac)))
         result%covariance = covn
      else
         allocate(result%covariance(size(par),size(par)))
         result%covariance = 0.0_dp
         if (result%status == 0) result%status = 20 + stat2
      end if

   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp), allocatable :: cc(:), bb(:), pp(:,:)
         integer :: ii
         call unpack_ordinal(theta, nc, cc, bb)
         call ordinal_probabilities(x, cc, bb, result%link, pp)
         nll = 0.0_dp
         do ii = 1, n
            nll = nll - w(ii)*log(max(pp(ii,y(ii)), tiny(1.0_dp)))
         end do
      end function objective
   end subroutine fit_ordinal

   subroutine ordinal_initial(y, ncat, link_id, u)
      integer, intent(in) :: y(:), ncat, link_id
      real(dp), intent(out) :: u(:)
      real(dp), allocatable :: cut(:)
      real(dp) :: pr
      integer :: j
      allocate(cut(ncat-1))
      do j = 1, ncat-1
         pr = real(count(y <= j),dp)/real(size(y),dp)
         pr = min(1.0_dp-1.0e-5_dp,max(1.0e-5_dp,pr))
         cut(j) = inverse_cdf_link(pr, link_id)
      end do
      u(1) = cut(1)
      do j = 2, ncat-1
         u(j) = log(max(cut(j)-cut(j-1),1.0e-3_dp))
      end do
   end subroutine ordinal_initial

   real(dp) function inverse_cdf_link(p, link_id) result(eta)
      use vgam_links, only : link_value
      real(dp), intent(in) :: p
      integer, intent(in) :: link_id
      eta = link_value(p, link_id)
   end function inverse_cdf_link

   subroutine unpack_ordinal(par, nc, cuts, beta)
      real(dp), intent(in) :: par(:)
      integer, intent(in) :: nc
      real(dp), allocatable, intent(out) :: cuts(:), beta(:)
      integer :: j
      allocate(cuts(nc), beta(size(par)-nc))
      cuts(1) = par(1)
      do j = 2, nc
         cuts(j) = cuts(j-1) + exp(min(par(j), 50.0_dp))
      end do
      beta = par(nc+1:)
   end subroutine unpack_ordinal

   subroutine ordinal_probabilities(x, cutpoints, beta, link_id, prob)
      real(dp), intent(in) :: x(:,:), cutpoints(:), beta(:)
      integer, intent(in) :: link_id
      real(dp), allocatable, intent(out) :: prob(:,:)
      real(dp), allocatable :: cum(:)
      real(dp) :: eta
      integer :: n, k, i, j
      n = size(x,1)
      k = size(cutpoints)+1
      if (size(x,2) /= size(beta)) error stop "ordinal_probabilities: shape mismatch"
      allocate(prob(n,k), cum(k-1))
      do i = 1, n
         eta = dot_product(x(i,:), beta)
         do j = 1, k-1
            cum(j) = link_inverse(cutpoints(j)-eta, link_id)
            cum(j) = min(1.0_dp,max(0.0_dp,cum(j)))
         end do
         prob(i,1) = cum(1)
         do j = 2, k-1
            prob(i,j) = max(0.0_dp,cum(j)-cum(j-1))
         end do
         prob(i,k) = max(0.0_dp,1.0_dp-cum(k-1))
         prob(i,:) = prob(i,:)/max(sum(prob(i,:)),tiny(1.0_dp))
      end do
   end subroutine ordinal_probabilities

   subroutine predict_ordinal(self, x, prob)
      class(ordinal_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: prob(:,:)
      call ordinal_probabilities(x, self%cutpoints, self%coefficients, self%link, prob)
   end subroutine predict_ordinal

   subroutine fit_beta_regression(y, x, result, weights, max_iter, tol)
      real(dp), intent(in) :: y(:), x(:,:)
      type(beta_regression_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:,:), covu(:,:), jac(:,:), covn(:,:)
      real(dp) :: fval, tolerance, phi, mu0
      integer :: n, p, stat, stat2, i, niter

      n = size(y)
      p = size(x,2)
      if (n <= 0 .or. p <= 0 .or. size(x,1) /= n .or. any(y <= 0.0_dp) .or. any(y >= 1.0_dp)) then
         result%status = 1
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         w = weights
      else
         allocate(w(n))
         w = 1.0_dp
      end if
      niter = 250
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      allocate(par(p+1))
      par = 0.0_dp
      mu0 = min(0.999_dp,max(0.001_dp,sum(w*y)/max(sum(w),tiny(1.0_dp))))
      if (p >= 1 .and. all(abs(x(:,1)-1.0_dp) < 1.0e-12_dp)) par(1) = log(mu0/(1.0_dp-mu0))
      phi = max(2.0_dp, mu0*(1.0_dp-mu0)/max(weighted_variance(y,w),1.0e-4_dp)-1.0_dp)
      par(p+1) = log(phi)
      call bfgs_minimize(objective, par, fval, stat, niter, tolerance)
      result%status = stat
      result%converged = stat == 0
      result%coefficients = par(1:p)
      result%precision = exp(min(par(p+1), 50.0_dp))
      allocate(result%fitted(n))
      result%fitted = logistic_vec(matmul(x, result%coefficients))
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(p+1,dp)

      allocate(hess(p+1,p+1), covu(p+1,p+1))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, covu, stat2)
      if (stat2 == 0) then
         allocate(jac(p+1,p+1), covn(p+1,p+1))
         jac = 0.0_dp
         do i = 1, p
            jac(i,i) = 1.0_dp
         end do
         jac(p+1,p+1) = result%precision
         covn = matmul(jac, matmul(covu, transpose(jac)))
         result%covariance = covn
      else
         allocate(result%covariance(p+1,p+1))
         result%covariance = 0.0_dp
         if (result%status == 0) result%status = 20 + stat2
      end if

   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp), allocatable :: mu(:)
         real(dp) :: ph, aa, bb
         integer :: ii
         ph = exp(min(theta(p+1),50.0_dp))
         mu = logistic_vec(matmul(x,theta(1:p)))
         nll = 0.0_dp
         do ii = 1, n
            aa = max(mu(ii)*ph,tiny(1.0_dp))
            bb = max((1.0_dp-mu(ii))*ph,tiny(1.0_dp))
            nll = nll - w(ii)*(log_gamma(ph)-log_gamma(aa)-log_gamma(bb) + &
               (aa-1.0_dp)*log(y(ii))+(bb-1.0_dp)*log(1.0_dp-y(ii)))
         end do
      end function objective
   end subroutine fit_beta_regression

   subroutine predict_beta_regression(self, x, fitted)
      class(beta_regression_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: fitted(:)
      if (size(x,2) /= size(self%coefficients)) error stop "beta regression predict: shape mismatch"
      fitted = logistic_vec(matmul(x,self%coefficients))
   end subroutine predict_beta_regression

   subroutine fit_negative_binomial(y, x, result, weights, offset, max_iter, tol)
      real(dp), intent(in) :: y(:), x(:,:)
      type(negative_binomial_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), off(:), par(:), hess(:,:), covu(:,:), jac(:,:), covn(:,:)
      real(dp) :: fval, tolerance, mean_y, var_y
      integer :: n, p, stat, stat2, i, niter

      n = size(y)
      p = size(x,2)
      if (n <= 0 .or. p <= 0 .or. size(x,1) /= n .or. any(y < 0.0_dp)) then
         result%status = 1
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         w = weights
      else
         allocate(w(n)); w = 1.0_dp
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = 3
            return
         end if
         off = offset
      else
         allocate(off(n)); off = 0.0_dp
      end if
      niter = 300
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      allocate(par(p+1)); par = 0.0_dp
      mean_y = sum(w*y)/max(sum(w),tiny(1.0_dp))
      var_y = weighted_variance(y,w)
      if (p >= 1 .and. all(abs(x(:,1)-1.0_dp) < 1.0e-12_dp)) &
         par(1) = log(max(mean_y,1.0e-3_dp))
      par(p+1) = log(max(0.1_dp, mean_y*mean_y/max(var_y-mean_y,0.1_dp)))
      call bfgs_minimize(objective, par, fval, stat, niter, tolerance)
      result%status = stat
      result%converged = stat == 0
      result%coefficients = par(1:p)
      result%size = exp(min(par(p+1),50.0_dp))
      result%fitted = exp(min(matmul(x,result%coefficients)+off, 700.0_dp))
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(p+1,dp)

      allocate(hess(p+1,p+1), covu(p+1,p+1))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess,covu,stat2)
      if (stat2 == 0) then
         allocate(jac(p+1,p+1),covn(p+1,p+1)); jac = 0.0_dp
         do i=1,p; jac(i,i)=1.0_dp; end do
         jac(p+1,p+1)=result%size
         covn=matmul(jac,matmul(covu,transpose(jac)))
         result%covariance=covn
      else
         allocate(result%covariance(p+1,p+1)); result%covariance=0.0_dp
         if (result%status==0) result%status=20+stat2
      end if

   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp), allocatable :: mu(:)
         real(dp) :: sz
         integer :: ii
         sz=exp(min(theta(p+1),50.0_dp))
         mu=exp(min(matmul(x,theta(1:p))+off,700.0_dp))
         nll=0.0_dp
         do ii=1,n
            nll=nll-w(ii)*(log_gamma(y(ii)+sz)-log_gamma(sz)-log_gamma(y(ii)+1.0_dp) + &
               sz*(log(sz)-log(sz+mu(ii))) + y(ii)*(log(max(mu(ii),tiny(1.0_dp)))-log(sz+mu(ii))))
         end do
      end function objective
   end subroutine fit_negative_binomial

   subroutine predict_negative_binomial(self, x, fitted, offset)
      class(negative_binomial_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: fitted(:)
      real(dp), intent(in), optional :: offset(:)
      real(dp), allocatable :: eta(:)
      eta=matmul(x,self%coefficients)
      if(present(offset)) eta=eta+offset
      fitted=exp(min(eta,700.0_dp))
   end subroutine predict_negative_binomial

   subroutine fit_zero_inflated_poisson(y, x_count, x_zero, result, weights, offset, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_count(:,:), x_zero(:,:)
      type(zip_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), off(:), par(:), hess(:,:), cov(:,:)
      real(dp) :: fval, tolerance, mean_y
      integer :: n, pc, pz, stat, stat2, niter

      n=size(y); pc=size(x_count,2); pz=size(x_zero,2)
      if(n<=0.or.pc<=0.or.pz<=0.or.size(x_count,1)/=n.or.size(x_zero,1)/=n.or.any(y<0))then
         result%status=1; return
      end if
      if(present(weights))then
         if(size(weights)/=n.or.any(weights<0.0_dp))then
            result%status=2; return
         end if
         w=weights
      else
         allocate(w(n)); w=1.0_dp
      end if
      if(present(offset))then
         if(size(offset)/=n)then
            result%status=3; return
         end if
         off=offset
      else
         allocate(off(n)); off=0.0_dp
      end if
      niter=300; if(present(max_iter))niter=max_iter
      tolerance=1.0e-7_dp; if(present(tol))tolerance=tol
      allocate(par(pc+pz)); par=0.0_dp
      mean_y=sum(w*real(y,dp))/max(sum(w),tiny(1.0_dp))
      if(all(abs(x_count(:,1)-1.0_dp)<1.0e-12_dp))par(1)=log(max(mean_y,1.0e-3_dp))
      if(all(abs(x_zero(:,1)-1.0_dp)<1.0e-12_dp))par(pc+1)=-2.0_dp
      call bfgs_minimize(objective,par,fval,stat,niter,tolerance)
      result%status=stat; result%converged=stat==0
      result%count_coefficients=par(1:pc)
      result%zero_coefficients=par(pc+1:)
      call zip_fitted(x_count,x_zero,result%count_coefficients,result%zero_coefficients,off, &
         result%fitted_mean,result%zero_probability)
      result%loglik=-fval
      result%aic=2.0_dp*fval+2.0_dp*real(size(par),dp)
      allocate(hess(size(par),size(par)),cov(size(par),size(par)))
      call numerical_hessian(objective,par,hess)
      call invert_matrix(hess,cov,stat2)
      if(stat2==0)then
         result%covariance=cov
      else
         allocate(result%covariance(size(par),size(par))); result%covariance=0.0_dp
         if(result%status==0)result%status=20+stat2
      end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp), allocatable :: mu(:),zi(:)
         real(dp)::pr
         integer::ii
         mu=exp(min(matmul(x_count,theta(1:pc))+off,700.0_dp))
         zi=logistic_vec(matmul(x_zero,theta(pc+1:)))
         nll=0.0_dp
         do ii=1,n
            if(y(ii)==0)then
               pr=zi(ii)+(1.0_dp-zi(ii))*exp(-mu(ii))
               nll=nll-w(ii)*log(max(pr,tiny(1.0_dp)))
            else
               nll=nll-w(ii)*(log(max(1.0_dp-zi(ii),tiny(1.0_dp)))-mu(ii) + &
                  real(y(ii),dp)*log(max(mu(ii),tiny(1.0_dp)))-log_gamma(real(y(ii)+1,dp)))
            end if
         end do
      end function objective
   end subroutine fit_zero_inflated_poisson

   subroutine zip_fitted(xc,xz,bc,bz,off,mean,pi0)
      real(dp),intent(in)::xc(:,:),xz(:,:),bc(:),bz(:),off(:)
      real(dp),allocatable,intent(out)::mean(:),pi0(:)
      real(dp),allocatable::mu(:)
      mu=exp(min(matmul(xc,bc)+off,700.0_dp))
      pi0=logistic_vec(matmul(xz,bz))
      mean=(1.0_dp-pi0)*mu
   end subroutine zip_fitted

   subroutine predict_zip(self, x_count, x_zero, fitted_mean, zero_probability, offset)
      class(zip_result_t), intent(in) :: self
      real(dp), intent(in) :: x_count(:,:), x_zero(:,:)
      real(dp), allocatable, intent(out) :: fitted_mean(:), zero_probability(:)
      real(dp), intent(in), optional :: offset(:)
      real(dp), allocatable :: off(:)
      allocate(off(size(x_count,1))); off=0.0_dp
      if(present(offset))off=offset
      call zip_fitted(x_count,x_zero,self%count_coefficients,self%zero_coefficients,off, &
         fitted_mean,zero_probability)
   end subroutine predict_zip

   pure function logistic_vec(eta) result(p)
      real(dp), intent(in) :: eta(:)
      real(dp) :: p(size(eta))
      integer :: i
      do i=1,size(eta)
         if(eta(i)>=0.0_dp)then
            p(i)=1.0_dp/(1.0_dp+exp(-min(eta(i),700.0_dp)))
         else
            p(i)=exp(max(eta(i),-700.0_dp))/(1.0_dp+exp(max(eta(i),-700.0_dp)))
         end if
      end do
   end function logistic_vec

   pure real(dp) function weighted_variance(x,w) result(v)
      real(dp),intent(in)::x(:),w(:)
      real(dp)::m,sw
      sw=sum(w)
      if(sw<=0.0_dp)then
         v=0.0_dp
      else
         m=sum(w*x)/sw
         v=sum(w*(x-m)**2)/sw
      end if
   end function weighted_variance

end module vgam_extended_models
