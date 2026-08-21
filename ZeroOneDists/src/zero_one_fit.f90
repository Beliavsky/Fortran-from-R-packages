! SPDX-License-Identifier: MIT
module zero_one_fit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use zero_one_kinds, only : dp
   use zero_one_families
   implicit none
   private
   type, public :: zero_one_fit_result_t
      real(dp), allocatable :: beta_mu(:), beta_sigma(:), beta_nu(:)
      real(dp), allocatable :: fitted(:,:)
      real(dp), allocatable :: covariance(:,:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
   end type zero_one_fit_result_t
   type :: fit_context_t
      real(dp), allocatable :: y(:), xm(:,:), xs(:,:), xn(:,:), w(:)
      integer :: family=0, pm=0, ps=0, pn=0, np=0
      integer :: links(3)=[link_identity,link_identity,link_identity]
   end type fit_context_t
   type(fit_context_t), save :: ctx
   public :: fit_zero_one
contains
   subroutine fit_zero_one(y, x_mu, family, result, x_sigma, x_nu, weights, links, start, max_iter, tol)
      real(dp), intent(in) :: y(:), x_mu(:,:)
      integer, intent(in) :: family
      type(zero_one_fit_result_t), intent(out) :: result
      real(dp), intent(in), optional :: x_sigma(:,:), x_nu(:,:), weights(:), start(:)
      integer, intent(in), optional :: links(3), max_iter
      real(dp), intent(in), optional :: tol
      integer :: n, npar, pm, ps, pn, np, maxit, status
      integer :: lnk(3)
      real(dp) :: tolerance, fval
      real(dp), allocatable :: theta(:), hinv(:,:)
      n = size(y)
      npar = family_npar(family)
      if (npar == 0 .or. size(x_mu,1) /= n .or. .not.family_valid_y(family,y)) then
         result%status = 10
         return
      end if
      pm = size(x_mu,2)
      ps = 0
      pn = 0
      if (npar >= 2) then
         if (.not.present(x_sigma)) then
            result%status = 11
            return
         end if
         if (size(x_sigma,1) /= n) then
            result%status = 12
            return
         end if
         ps = size(x_sigma,2)
      end if
      if (npar >= 3) then
         if (.not.present(x_nu)) then
            result%status = 13
            return
         end if
         if (size(x_nu,1) /= n) then
            result%status = 14
            return
         end if
         pn = size(x_nu,2)
      end if
      np = pm+ps+pn
      if (np <= 0) then
         result%status = 15
         return
      end if
      call family_default_links(family,lnk)
      if (present(links)) lnk = links
      allocate(theta(np))
      theta = 0.0_dp
      if (present(start)) then
         if (size(start) /= np) then
            result%status = 16
            return
         end if
         theta = start
      else
         call initialize_coefficients(y,x_mu,family,lnk,pm,ps,pn,theta)
      end if
      call set_context(y,x_mu,family,lnk,pm,ps,pn,x_sigma,x_nu,weights,status)
      if (status /= 0) then
         result%status = status
         call clear_context()
         return
      end if
      maxit = 250
      if (present(max_iter)) maxit = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(theta,fval,hinv,result%iterations,result%converged,maxit,tolerance)
      result%status = merge(0,1,result%converged)
      result%loglik = -fval
      result%aic = 2.0_dp*real(np,dp)-2.0_dp*result%loglik
      result%covariance = hinv
      call unpack_result(theta,pm,ps,pn,result)
      call fill_fitted(theta,result)
      call clear_context()
   end subroutine fit_zero_one

   subroutine initialize_coefficients(y, xm, family, links, pm, ps, pn, theta)
      real(dp), intent(in) :: y(:), xm(:,:)
      integer, intent(in) :: family, links(3), pm, ps, pn
      real(dp), intent(out) :: theta(:)
      real(dp) :: par(3)
      integer :: o2, o3
      theta = 0.0_dp
      call family_initial(family,y,par)
      o2 = pm
      o3 = pm+ps
      if (pm > 0 .and. is_intercept(xm(:,1))) theta(1) = linkfun(par(1),links(1))
      if (ps > 0) theta(o2+1) = linkfun(par(2),links(2))
      if (pn > 0) theta(o3+1) = linkfun(par(3),links(3))
   end subroutine initialize_coefficients

   pure logical function is_intercept(x) result(ok)
      real(dp), intent(in) :: x(:)
      ok = maxval(abs(x-1.0_dp)) <= 1.0e-12_dp
   end function is_intercept

   subroutine set_context(y,xm,family,links,pm,ps,pn,xs,xn,weights,status)
      real(dp), intent(in) :: y(:), xm(:,:)
      integer, intent(in) :: family, links(3), pm, ps, pn
      real(dp), intent(in), optional :: xs(:,:), xn(:,:), weights(:)
      integer, intent(out) :: status
      integer :: n
      call clear_context()
      n = size(y)
      ctx%y = y
      ctx%xm = xm
      ctx%family = family
      ctx%links = links
      ctx%pm = pm
      ctx%ps = ps
      ctx%pn = pn
      ctx%np = pm+ps+pn
      allocate(ctx%w(n))
      ctx%w = 1.0_dp
      status = 0
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            status = 20
            return
         end if
         ctx%w = weights
      end if
      if (ps > 0) then
         ctx%xs = xs
      else
         allocate(ctx%xs(n,0))
      end if
      if (pn > 0) then
         ctx%xn = xn
      else
         allocate(ctx%xn(n,0))
      end if
   end subroutine set_context

   subroutine clear_context()
      if (allocated(ctx%y)) deallocate(ctx%y)
      if (allocated(ctx%xm)) deallocate(ctx%xm)
      if (allocated(ctx%xs)) deallocate(ctx%xs)
      if (allocated(ctx%xn)) deallocate(ctx%xn)
      if (allocated(ctx%w)) deallocate(ctx%w)
      ctx%family = 0
      ctx%pm = 0
      ctx%ps = 0
      ctx%pn = 0
      ctx%np = 0
   end subroutine clear_context

   real(dp) function objective(theta) result(nll)
      real(dp), intent(in) :: theta(:)
      real(dp) :: par(3), lp
      integer :: i
      nll = 0.0_dp
      do i = 1, size(ctx%y)
         call observation_parameters(i,theta,par)
         if (.not.family_valid_parameters(ctx%family,par)) then
            nll = 1.0e100_dp
            return
         end if
         lp = family_logpdf(ctx%family,ctx%y(i),par)
         if (.not.ieee_is_finite(lp)) then
            nll = 1.0e100_dp
            return
         end if
         nll = nll-ctx%w(i)*lp
      end do
   end function objective

   subroutine observation_parameters(i,theta,par)
      integer, intent(in) :: i
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: par(3)
      integer :: o2, o3
      par = 1.0_dp
      o2 = ctx%pm
      o3 = ctx%pm+ctx%ps
      par(1) = linkinv(dot_product(ctx%xm(i,:),theta(1:ctx%pm)),ctx%links(1))
      if (ctx%ps > 0) par(2) = linkinv(dot_product(ctx%xs(i,:),theta(o2+1:o3)),ctx%links(2))
      if (ctx%pn > 0) par(3) = linkinv(dot_product(ctx%xn(i,:),theta(o3+1:ctx%np)),ctx%links(3))
   end subroutine observation_parameters

   subroutine numerical_gradient(theta, grad)
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: grad(:)
      real(dp) :: tp(size(theta)), tm(size(theta)), h
      integer :: j
      do j = 1, size(theta)
         h = 1.0e-5_dp*max(1.0_dp,abs(theta(j)))
         tp = theta
         tm = theta
         tp(j) = tp(j)+h
         tm(j) = tm(j)-h
         grad(j) = (objective(tp)-objective(tm))/(2.0_dp*h)
      end do
   end subroutine numerical_gradient

   subroutine bfgs_minimize(theta,fval,hinv,iterations,converged,max_iter,tol)
      real(dp), intent(inout) :: theta(:)
      real(dp), intent(out) :: fval
      real(dp), allocatable, intent(out) :: hinv(:,:)
      integer, intent(out) :: iterations
      logical, intent(out) :: converged
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: tol
      integer :: n, i, j, it
      real(dp) :: fnew, alpha, slope, ys, rho
      real(dp), allocatable :: g(:), gnew(:), p(:), trial(:), s(:), yv(:), eye(:,:), a(:,:), tmp(:,:)
      n = size(theta)
      allocate(g(n),gnew(n),p(n),trial(n),s(n),yv(n),hinv(n,n),eye(n,n),a(n,n),tmp(n,n))
      eye = 0.0_dp
      do i = 1, n
         eye(i,i) = 1.0_dp
      end do
      hinv = eye
      fval = objective(theta)
      call numerical_gradient(theta,g)
      converged = .false.
      iterations = 0
      do it = 1, max_iter
         iterations = it
         if (maxval(abs(g)) <= tol*(1.0_dp+abs(fval))) then
            converged = .true.
            exit
         end if
         p = -matmul(hinv,g)
         slope = dot_product(g,p)
         if (slope >= 0.0_dp) then
            p = -g
            hinv = eye
            slope = -dot_product(g,g)
         end if
         alpha = 1.0_dp
         do j = 1, 50
            trial = theta+alpha*p
            fnew = objective(trial)
            if (ieee_is_finite(fnew) .and. fnew <= fval+1.0e-4_dp*alpha*slope) exit
            alpha = 0.5_dp*alpha
         end do
         if (alpha < 1.0e-12_dp) exit
         s = alpha*p
         theta = trial
         call numerical_gradient(theta,gnew)
         yv = gnew-g
         ys = dot_product(yv,s)
         if (ys > 1.0e-12_dp*sqrt(max(dot_product(yv,yv)*dot_product(s,s),tiny(1.0_dp)))) then
            rho = 1.0_dp/ys
            a = eye
            do i = 1, n
               do j = 1, n
                  a(i,j) = a(i,j)-rho*s(i)*yv(j)
               end do
            end do
            tmp = matmul(a,matmul(hinv,transpose(a)))
            do i = 1, n
               do j = 1, n
                  hinv(i,j) = tmp(i,j)+rho*s(i)*s(j)
               end do
            end do
         end if
         fval = fnew
         g = gnew
      end do
      if (.not.converged) converged = maxval(abs(g)) <= 10.0_dp*tol*(1.0_dp+abs(fval))
   end subroutine bfgs_minimize

   subroutine unpack_result(theta,pm,ps,pn,result)
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: pm, ps, pn
      type(zero_one_fit_result_t), intent(inout) :: result
      integer :: o2, o3
      o2 = pm
      o3 = pm+ps
      result%beta_mu = theta(1:pm)
      if (ps > 0) result%beta_sigma = theta(o2+1:o3)
      if (pn > 0) result%beta_nu = theta(o3+1:o3+pn)
   end subroutine unpack_result

   subroutine fill_fitted(theta,result)
      real(dp), intent(in) :: theta(:)
      type(zero_one_fit_result_t), intent(inout) :: result
      real(dp) :: par(3)
      integer :: i, npar
      npar = family_npar(ctx%family)
      allocate(result%fitted(size(ctx%y),npar))
      do i = 1, size(ctx%y)
         call observation_parameters(i,theta,par)
         result%fitted(i,:) = par(1:npar)
      end do
   end subroutine fill_fitted
end module zero_one_fit
