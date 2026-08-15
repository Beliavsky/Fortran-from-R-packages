! Modern Fortran translation of R package skewunit.
! SPDX-License-Identifier: GPL-2.0-or-later
module skewunit_fit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use skewunit_kinds, only : dp, nan_dp
   use skewunit_distributions, only : family_none, family_asin, family_uquad, &
      n_families, family_has_delta, dskewunit
   use skewunit_optimize, only : brent_minimize, nelder_mead
   implicit none
   private

   type, public :: skewunit_fit_result
      integer :: family1 = family_asin
      integer :: family2 = family_asin
      integer :: npar = 0
      integer :: convergence = 0
      integer :: iterations = 0
      real(dp) :: coefficients(3) = 0.0_dp
      real(dp) :: std_error(3) = 0.0_dp
      logical :: std_error_available = .false.
      real(dp) :: loglik = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
   end type skewunit_fit_result

   type, public :: skewunit_model_row
      integer :: family1 = family_asin
      integer :: family2 = family_none
      integer :: npar = 0
      real(dp) :: criterion = 0.0_dp
   end type skewunit_model_row

   type, public :: skewunit_choice_result
      character(len=3) :: criterion_name = 'AIC'
      type(skewunit_model_row) :: summary(30)
      type(skewunit_fit_result) :: best_fit
   end type skewunit_choice_result

   type :: fit_context
      real(dp), allocatable :: x(:)
      integer :: family1 = family_asin
      integer :: family2 = family_asin
      integer :: case_id = 0
   end type fit_context

   public :: estimate_skewunit, choose_skewunit, coefficient_name

contains

   function coefficient_name(fit, index) result(name)
      type(skewunit_fit_result), intent(in) :: fit
      integer, intent(in) :: index
      character(len=8) :: name
      logical :: h1, h2

      name = ''
      h1 = family_has_delta(fit%family1)
      h2 = family_has_delta(fit%family2)
      if (fit%family2 == family_none) then
         if (h1 .and. index == 1) name = 'delta'
      else
         if (index == 1) name = 'lambda'
         if (h1 .and. h2) then
            if (index == 2) name = 'delta1'
            if (index == 3) name = 'delta2'
         else if ((h1 .or. h2) .and. index == 2) then
            name = 'delta'
         end if
      end if
   end function coefficient_name

   subroutine estimate_skewunit(x, family1, family2, fit, est_var, maxit)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: family1
      integer, intent(in), optional :: family2
      type(skewunit_fit_result), intent(out) :: fit
      logical, intent(in), optional :: est_var
      integer, intent(in), optional :: maxit
      type(fit_context) :: context
      integer :: f2, fit_case, npar, status1, status2, iter1, iter2, limit
      real(dp) :: fmin1, fmin2, tmin1, zero_par(3)
      real(dp), allocatable :: start1(:), start2(:), theta1(:), theta2(:), theta(:)
      logical :: get_var

      f2 = family_asin
      if (present(family2)) f2 = family2
      get_var = .true.
      if (present(est_var)) get_var = est_var
      limit = 10000
      if (present(maxit)) limit = max(10,maxit)

      call initialize_fit(fit,family1,f2)
      if (size(x) == 0 .or. any(x < 0.0_dp) .or. any(x > 1.0_dp) .or. &
          family1 < 1 .or. family1 > n_families .or. &
          f2 < family_none .or. f2 > n_families) then
         fit%convergence = 2
         fit%loglik = -huge(1.0_dp)
         fit%aic = huge(1.0_dp)
         fit%bic = huge(1.0_dp)
         return
      end if

      call determine_case(family1,f2,fit_case,npar)
      fit%npar = npar
      context%x = x
      context%family1 = family1
      context%family2 = f2
      context%case_id = fit_case
      zero_par = 0.0_dp

      if (npar == 0) then
         fmin1 = nll_from_parameters(zero_par,context)
         fit%convergence = 0
         fit%iterations = 0
      else if (fit_case == 1 .or. fit_case == 5) then
         call brent_minimize(objective_scalar_ctx,context,-100.0_dp,100.0_dp, &
            tmin1,fmin1,status1,iter1,tol=1.0e-10_dp,maxit=limit)
         allocate(theta(1))
         theta(1) = tmin1
         fit%convergence = status1
         fit%iterations = iter1
         call transformed_to_natural(theta,fit_case,fit%coefficients)
      else
         allocate(start1(npar),start2(npar),theta1(npar),theta2(npar),theta(npar))
         start1 = -0.2_dp
         start2 = 0.2_dp
         call nelder_mead(objective_vector_ctx,context,start1,theta1,fmin1, &
            status1,iter1,tol=1.0e-10_dp,maxit=limit)
         call nelder_mead(objective_vector_ctx,context,start2,theta2,fmin2, &
            status2,iter2,tol=1.0e-10_dp,maxit=limit)
         if (fmin2 < fmin1) then
            theta = theta2
            fmin1 = fmin2
            fit%convergence = status2
            fit%iterations = iter2
         else
            theta = theta1
            fit%convergence = status1
            fit%iterations = iter1
         end if
         call transformed_to_natural(theta,fit_case,fit%coefficients)
      end if

      fit%loglik = -fmin1
      fit%aic = -2.0_dp*fit%loglik+2.0_dp*real(npar,dp)
      fit%bic = -2.0_dp*fit%loglik+log(real(size(x),dp))*real(npar,dp)

      if (get_var .and. npar > 0 .and. fit%convergence == 0) then
         call estimate_standard_errors(fit%coefficients(1:npar),context,fit)
      end if
   end subroutine estimate_skewunit

   real(dp) function objective_scalar_ctx(t, context_any) result(v)
      real(dp), intent(in) :: t
      class(*), intent(in) :: context_any
      real(dp) :: theta(1)

      theta(1) = t
      select type(context => context_any)
      type is(fit_context)
         v = transformed_nll(theta,context)
      class default
         v = huge(1.0_dp)
      end select
   end function objective_scalar_ctx

   real(dp) function objective_vector_ctx(theta, context_any) result(v)
      real(dp), intent(in) :: theta(:)
      class(*), intent(in) :: context_any
      select type(context => context_any)
      type is(fit_context)
         v = transformed_nll(theta,context)
      class default
         v = huge(1.0_dp)
      end select
   end function objective_vector_ctx

   real(dp) function transformed_nll(theta, context) result(v)
      real(dp), intent(in) :: theta(:)
      type(fit_context), intent(in) :: context
      real(dp) :: par(3)
      call transformed_to_natural(theta,context%case_id,par)
      v = nll_from_parameters(par,context)
   end function transformed_nll

   real(dp) function nll_from_parameters(par, context) result(v)
      real(dp), intent(in) :: par(3)
      type(fit_context), intent(in) :: context
      real(dp) :: lambda, d1, d2, lf
      integer :: i

      lambda = 0.0_dp
      d1 = 1.0_dp
      d2 = 1.0_dp
      select case(context%case_id)
      case(1)
         lambda = par(1)
      case(2)
         lambda = par(1)
         d1 = par(2)
      case(3)
         lambda = par(1)
         d1 = par(2)
         d2 = par(3)
      case(5)
         d1 = par(1)
      end select

      if (abs(lambda) > 1.0_dp .or. d1 <= 0.0_dp .or. d2 <= 0.0_dp) then
         v = huge(1.0_dp)
         return
      end if

      v = 0.0_dp
      do i = 1, size(context%x)
         ! Preserve the upstream exact-0.5 omission for U-quadratic f.
         if (context%family1 == family_uquad .and. context%x(i) == 0.5_dp) cycle
         lf = dskewunit(context%x(i),lambda,d1,d2,context%family1, &
            context%family2,.true.)
         if (.not.ieee_is_finite(lf)) then
            v = huge(1.0_dp)
            return
         end if
         v = v-lf
      end do
   end function nll_from_parameters

   subroutine estimate_standard_errors(par, context, fit)
      real(dp), intent(in) :: par(:)
      type(fit_context), intent(in) :: context
      type(skewunit_fit_result), intent(inout) :: fit
      real(dp) :: h(size(par),size(par)), invh(size(par),size(par))
      logical :: ok
      integer :: i

      call numerical_hessian_ctx(par,context,h)
      call invert_matrix(h,invh,ok)
      if (.not.ok) return
      do i = 1, size(par)
         if (invh(i,i) <= 0.0_dp) return
      end do
      fit%std_error = 0.0_dp
      do i = 1, size(par)
         fit%std_error(i) = sqrt(invh(i,i))
      end do
      fit%std_error_available = .true.
   end subroutine estimate_standard_errors

   subroutine numerical_hessian_ctx(x, context, hess)
      real(dp), intent(in) :: x(:)
      type(fit_context), intent(in) :: context
      real(dp), intent(out) :: hess(size(x),size(x))
      real(dp) :: xp(size(x)), xm(size(x)), xpp(size(x)), xpm(size(x))
      real(dp) :: xmp(size(x)), xmm(size(x)), step(size(x)), f0
      integer :: i, j, n

      n = size(x)
      f0 = nll_natural_vector(x,context)
      do i = 1, n
         step(i) = 1.0e-4_dp*max(1.0_dp,abs(x(i)))
         if (x(i) > 0.0_dp .and. x(i)-step(i) <= 0.0_dp) &
            step(i) = 0.25_dp*x(i)
      end do

      hess = 0.0_dp
      do i = 1, n
         xp = x
         xm = x
         xp(i) = xp(i)+step(i)
         xm(i) = xm(i)-step(i)
         hess(i,i) = (nll_natural_vector(xp,context)-2.0_dp*f0 &
            + nll_natural_vector(xm,context))/(step(i)*step(i))
         do j = i+1, n
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = xpp(i)+step(i)
            xpp(j) = xpp(j)+step(j)
            xpm(i) = xpm(i)+step(i)
            xpm(j) = xpm(j)-step(j)
            xmp(i) = xmp(i)-step(i)
            xmp(j) = xmp(j)+step(j)
            xmm(i) = xmm(i)-step(i)
            xmm(j) = xmm(j)-step(j)
            hess(i,j) = (nll_natural_vector(xpp,context) &
               - nll_natural_vector(xpm,context) &
               - nll_natural_vector(xmp,context) &
               + nll_natural_vector(xmm,context))/(4.0_dp*step(i)*step(j))
            hess(j,i) = hess(i,j)
         end do
      end do
   end subroutine numerical_hessian_ctx

   real(dp) function nll_natural_vector(p, context) result(v)
      real(dp), intent(in) :: p(:)
      type(fit_context), intent(in) :: context
      real(dp) :: full(3)
      full = 0.0_dp
      full(1:min(3,size(p))) = p(1:min(3,size(p)))
      v = nll_from_parameters(full,context)
   end function nll_natural_vector

   subroutine choose_skewunit(x, choice, criteria, est_var, maxit)
      real(dp), intent(in) :: x(:)
      type(skewunit_choice_result), intent(out) :: choice
      character(len=*), intent(in), optional :: criteria
      logical, intent(in), optional :: est_var
      integer, intent(in), optional :: maxit
      character(len=3) :: crit
      type(skewunit_fit_result) :: fit
      logical :: get_var
      integer :: f1, f2, k, i, j

      crit = 'AIC'
      if (present(criteria)) then
         if (trim(adjustl(criteria)) == 'BIC' .or. &
             trim(adjustl(criteria)) == 'bic') crit = 'BIC'
      end if
      get_var = .true.
      if (present(est_var)) get_var = est_var
      choice%criterion_name = crit

      k = 0
      do f1 = 1, n_families
         do f2 = 1, n_families
            k = k+1
            call estimate_for_choice(x,f1,f2,fit,.false.,maxit)
            call set_row(choice%summary(k),fit,crit)
         end do
      end do
      do f1 = 1, n_families
         k = k+1
         call estimate_for_choice(x,f1,family_none,fit,.false.,maxit)
         call set_row(choice%summary(k),fit,crit)
      end do

      do i = 1, 29
         k = i
         do j = i+1, 30
            if (choice%summary(j)%criterion < choice%summary(k)%criterion) k = j
         end do
         if (k /= i) call swap_rows(choice%summary(i),choice%summary(k))
      end do

      call estimate_for_choice(x,choice%summary(1)%family1, &
         choice%summary(1)%family2,choice%best_fit,get_var,maxit)
   end subroutine choose_skewunit

   subroutine estimate_for_choice(x, f1, f2, fit, est_var, maxit)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: f1, f2
      type(skewunit_fit_result), intent(out) :: fit
      logical, intent(in) :: est_var
      integer, intent(in), optional :: maxit
      if (present(maxit)) then
         call estimate_skewunit(x,f1,f2,fit,est_var,maxit)
      else
         call estimate_skewunit(x,f1,f2,fit,est_var)
      end if
   end subroutine estimate_for_choice

   subroutine set_row(row, fit, criterion)
      type(skewunit_model_row), intent(out) :: row
      type(skewunit_fit_result), intent(in) :: fit
      character(len=*), intent(in) :: criterion
      row%family1 = fit%family1
      row%family2 = fit%family2
      row%npar = fit%npar
      if (criterion == 'BIC') then
         row%criterion = fit%bic
      else
         row%criterion = fit%aic
      end if
   end subroutine set_row

   subroutine swap_rows(a,b)
      type(skewunit_model_row), intent(inout) :: a,b
      type(skewunit_model_row) :: tmp
      tmp = a
      a = b
      b = tmp
   end subroutine swap_rows

   subroutine initialize_fit(fit, f1, f2)
      type(skewunit_fit_result), intent(out) :: fit
      integer, intent(in) :: f1, f2
      fit%family1 = f1
      fit%family2 = f2
      fit%npar = 0
      fit%convergence = 0
      fit%iterations = 0
      fit%coefficients = 0.0_dp
      fit%std_error = nan_dp()
      fit%std_error_available = .false.
      fit%loglik = 0.0_dp
      fit%aic = 0.0_dp
      fit%bic = 0.0_dp
   end subroutine initialize_fit

   subroutine determine_case(f1, f2, case_id, npar)
      integer, intent(in) :: f1, f2
      integer, intent(out) :: case_id, npar
      logical :: h1, h2
      h1 = family_has_delta(f1)
      h2 = family_has_delta(f2)
      if (f2 == family_none) then
         if (h1) then
            case_id = 5
            npar = 1
         else
            case_id = 0
            npar = 0
         end if
      else if (.not.h1 .and. .not.h2) then
         case_id = 1
         npar = 1
      else if (h1 .and. h2) then
         case_id = 3
         npar = 3
      else
         case_id = 2
         npar = 2
      end if
   end subroutine determine_case

   subroutine transformed_to_natural(theta, case_id, par)
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: case_id
      real(dp), intent(out) :: par(3)
      par = 0.0_dp
      select case(case_id)
      case(1)
         par(1) = theta(1)/sqrt(1.0_dp+theta(1)*theta(1))
      case(2)
         par(1) = theta(1)/sqrt(1.0_dp+theta(1)*theta(1))
         par(2) = exp(min(theta(2),700.0_dp))
      case(3)
         par(1) = theta(1)/sqrt(1.0_dp+theta(1)*theta(1))
         par(2) = exp(min(theta(2),700.0_dp))
         par(3) = exp(min(theta(3),700.0_dp))
      case(5)
         par(1) = exp(min(theta(1),700.0_dp))
      end select
   end subroutine transformed_to_natural

   subroutine invert_matrix(a, ainv, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(size(a,1),size(a,2))
      logical, intent(out) :: ok
      real(dp), allocatable :: aug(:,:), tmp(:)
      real(dp) :: pivot, factor, maxp
      integer :: n, i, j, k, p

      n = size(a,1)
      if (size(a,2) /= n) then
         ok = .false.
         return
      end if
      allocate(aug(n,2*n),tmp(2*n))
      aug(:,1:n) = a
      aug(:,n+1:2*n) = 0.0_dp
      do i = 1, n
         aug(i,n+i) = 1.0_dp
      end do

      ok = .true.
      do i = 1, n
         p = i
         maxp = abs(aug(i,i))
         do k = i+1, n
            if (abs(aug(k,i)) > maxp) then
               p = k
               maxp = abs(aug(k,i))
            end if
         end do
         if (maxp <= 1.0e-14_dp*max(1.0_dp,maxval(abs(a)))) then
            ok = .false.
            return
         end if
         if (p /= i) then
            tmp = aug(i,:)
            aug(i,:) = aug(p,:)
            aug(p,:) = tmp
         end if
         pivot = aug(i,i)
         aug(i,:) = aug(i,:)/pivot
         do j = 1, n
            if (j == i) cycle
            factor = aug(j,i)
            aug(j,:) = aug(j,:)-factor*aug(i,:)
         end do
      end do
      ainv = aug(:,n+1:2*n)
   end subroutine invert_matrix

end module skewunit_fit
