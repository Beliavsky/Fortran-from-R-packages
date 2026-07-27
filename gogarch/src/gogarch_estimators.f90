! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_estimators
   use gogarch_kinds, only : dp
   use gogarch_linalg, only : identity_matrix, outer_product, symmetric_eigen, inverse_matrix, symmetric_invsqrt
   use gogarch_core, only : initialize_gogarch, cora
   use gogarch_orthogonal, only : uprod_r, umatch, unvech
   use gogarch_optimizer, only : optimizer_result, nelder_mead
   use gogarch_univariate, only : fit_univariate
   use gogarch_ica, only : ica_result, fastica
   use gogarch_model, only : build_gogarch_fit
   use gogarch_types, only : gogarch_fit, garch11_fit, univariate_spec
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)

   type :: nls_context
      real(dp), allocatable :: ssi(:,:,:)
   end type nls_context

   type :: ml_context
      real(dp), allocatable :: whitened(:,:)
      integer :: garch_iterations = 120
      type(univariate_spec) :: factor_spec
   end type ml_context

   public :: fit_gogarch, fit_gogarch_ica, fit_gogarch_mm
   public :: fit_gogarch_nls, fit_gogarch_ml, gonls_objective
   public :: gogarch_from_angles, gogarch_negloglik

contains

   function fit_gogarch(data, method, lag_max, max_outer_iterations, max_garch_iterations, factor_spec) result(fit)
      real(dp), intent(in) :: data(:,:)
      character(len=*), intent(in) :: method
      integer, intent(in), optional :: lag_max, max_outer_iterations, max_garch_iterations
      type(univariate_spec), intent(in), optional :: factor_spec
      type(gogarch_fit) :: fit
      select case (trim(adjustl(method)))
      case ('ica','ICA')
         fit = fit_gogarch_ica(data,max_garch_iterations=max_garch_iterations,factor_spec=factor_spec)
      case ('mm','MM')
         fit = fit_gogarch_mm(data,lag_max=lag_max,max_garch_iterations=max_garch_iterations,factor_spec=factor_spec)
      case ('nls','NLS')
         fit = fit_gogarch_nls(data,max_outer_iterations=max_outer_iterations,max_garch_iterations=max_garch_iterations, &
            factor_spec=factor_spec)
      case ('ml','ML')
         fit = fit_gogarch_ml(data,max_outer_iterations=max_outer_iterations,max_garch_iterations=max_garch_iterations, &
            factor_spec=factor_spec)
      case default
         fit%n = size(data,1)
         fit%m = size(data,2)
         fit%method = 'unknown'
         fit%status = 9
      end select
   end function fit_gogarch

   function gogarch_from_angles(data, angles, max_garch_iterations, factor_spec) result(fit)
      real(dp), intent(in) :: data(:,:), angles(:)
      integer, intent(in), optional :: max_garch_iterations
      type(univariate_spec), intent(in), optional :: factor_spec
      type(gogarch_fit) :: fit
      real(dp) :: rotation(size(data,2),size(data,2))
      logical :: ok
      integer :: git
      git = 350
      if (present(max_garch_iterations)) git = max_garch_iterations
      if (size(angles) /= size(data,2)*(size(data,2)-1)/2) then
         fit%n = size(data,1)
         fit%m = size(data,2)
         fit%method = 'angles'
         fit%status = 8
         return
      end if
      rotation = uprod_r(angles,ok)
      if (.not. ok) then
         fit%n = size(data,1)
         fit%m = size(data,2)
         fit%method = 'angles'
         fit%status = 8
         return
      end if
      call build_gogarch_fit(data,rotation,'angles',fit,max_garch_iterations=git,parameters=angles,factor_spec=factor_spec)
   end function gogarch_from_angles

   function gogarch_negloglik(data, angles, max_garch_iterations, factor_spec) result(value)
      real(dp), intent(in) :: data(:,:), angles(:)
      integer, intent(in), optional :: max_garch_iterations
      type(univariate_spec), intent(in), optional :: factor_spec
      real(dp) :: value
      type(gogarch_fit) :: fit
      fit = gogarch_from_angles(data,angles,max_garch_iterations,factor_spec)
      if (fit%status <= 1) then
         value = -fit%log_likelihood
      else
         value = huge(1.0_dp)/100.0_dp
      end if
   end function gogarch_negloglik

   function fit_gogarch_ica(data, max_ica_iterations, max_garch_iterations, factor_spec) result(fit)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: max_ica_iterations, max_garch_iterations
      type(univariate_spec), intent(in), optional :: factor_spec
      type(gogarch_fit) :: fit
      type(ica_result) :: ica
      integer :: icait, git
      icait = 500
      git = 350
      if (present(max_ica_iterations)) icait = max_ica_iterations
      if (present(max_garch_iterations)) git = max_garch_iterations
      ica = fastica(data,max_iterations=icait)
      call build_gogarch_fit(data,ica%rotation,'ica',fit,max_garch_iterations=git, &
         optimizer_iterations=ica%iterations,factor_spec=factor_spec)
      if (ica%status > 1) fit%status = 5
   end function fit_gogarch_ica

   function fit_gogarch_mm(data, lag_max, max_garch_iterations, factor_spec) result(fit)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: lag_max, max_garch_iterations
      type(univariate_spec), intent(in), optional :: factor_spec
      type(gogarch_fit) :: fit
      real(dp) :: covariance(size(data,2),size(data,2)), eigenvectors(size(data,2),size(data,2))
      real(dp) :: eigenvalues(size(data,2)), covroot(size(data,2),size(data,2)), invroot(size(data,2),size(data,2))
      real(dp), allocatable :: whitened(:,:), ssi(:,:,:), rotations(:,:,:), weights(:), gaps(:)
      real(dp) :: gamma(size(data,2),size(data,2)), values(size(data,2)), vectors(size(data,2),size(data,2))
      real(dp) :: sm(size(data,2),size(data,2)), term(size(data,2),size(data,2))
      real(dp) :: rotation(size(data,2),size(data,2)), inverse(size(data,2),size(data,2))
      real(dp) :: polar(size(data,2),size(data,2)), denominator
      logical :: ok, inv_ok, polar_ok
      integer :: n, m, lmax, lag, git
      n = size(data,1)
      m = size(data,2)
      lmax = 1
      git = 350
      if (present(lag_max)) lmax = max(0,min(abs(lag_max),n-1))
      if (present(max_garch_iterations)) git = max_garch_iterations
      call initialize_gogarch(data,covariance,eigenvectors,eigenvalues,covroot,invroot,ok)
      if (.not. ok) then
         call build_gogarch_fit(data,identity_matrix(m),'mm',fit,max_garch_iterations=git,factor_spec=factor_spec)
         fit%status = 2
         return
      end if
      if (lmax == 0) then
         call build_gogarch_fit(data,identity_matrix(m),'mm',fit,max_garch_iterations=git,weights=[1.0_dp], &
            factor_spec=factor_spec)
         return
      end if
      allocate(whitened(n,m),ssi(m,m,n),rotations(m,m,lmax),weights(lmax),gaps(lmax))
      whitened = matmul(data,invroot)
      do lag = 1, n
         ssi(:,:,lag) = outer_product(whitened(lag,:),whitened(lag,:))-identity_matrix(m)
      end do
      do lag = 1, lmax
         gamma = cora(ssi,lag=lag,standardize=.true.,ok=ok)
         if (.not. ok) then
            vectors = identity_matrix(m)
            values = [(real(m-lag,dp),lag=1,m)]
         else
            call symmetric_eigen(gamma,values,vectors,ok)
            if (.not. ok) vectors = identity_matrix(m)
         end if
         rotations(:,:,lag) = vectors
         gaps(lag) = minimum_eigen_gap(values)
      end do
      rotations(:,:,1) = umatch(identity_matrix(m),rotations(:,:,1),ok)
      do lag = 2, lmax
         rotations(:,:,lag) = umatch(rotations(:,:,1),rotations(:,:,lag),ok)
      end do
      denominator = sum(gaps)
      if (denominator > 1.0e-14_dp) then
         weights = gaps/denominator
      else
         weights = 1.0_dp/real(lmax,dp)
      end if
      sm = 0.0_dp
      do lag = 1, lmax
         inverse = inverse_matrix(identity_matrix(m)+rotations(:,:,lag),inv_ok)
         if (inv_ok) then
            term = matmul(identity_matrix(m)-rotations(:,:,lag),inverse)
            sm = sm+weights(lag)*term
         end if
      end do
      inverse = inverse_matrix(identity_matrix(m)+sm,inv_ok)
      if (inv_ok) then
         rotation = matmul(identity_matrix(m)-sm,inverse)
      else
         rotation = rotations(:,:,1)
      end if
      polar = symmetric_invsqrt(matmul(transpose(rotation),rotation),1.0e-12_dp,polar_ok)
      if (polar_ok) rotation = matmul(rotation,polar)
      call build_gogarch_fit(data,rotation,'mm',fit,max_garch_iterations=git,weights=weights,factor_spec=factor_spec)
   end function fit_gogarch_mm

   function fit_gogarch_nls(data, max_outer_iterations, max_garch_iterations, factor_spec) result(fit)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: max_outer_iterations, max_garch_iterations
      type(univariate_spec), intent(in), optional :: factor_spec
      type(gogarch_fit) :: fit
      type(nls_context) :: context
      type(optimizer_result) :: opt
      real(dp) :: covariance(size(data,2),size(data,2)), eigenvectors(size(data,2),size(data,2))
      real(dp) :: eigenvalues(size(data,2)), covroot(size(data,2),size(data,2)), invroot(size(data,2),size(data,2))
      real(dp), allocatable :: whitened(:,:), initial(:)
      real(dp) :: b(size(data,2),size(data,2)), values(size(data,2)), rotation(size(data,2),size(data,2))
      logical :: ok
      integer :: n, m, t, outer, git
      n = size(data,1)
      m = size(data,2)
      outer = 350
      git = 350
      if (present(max_outer_iterations)) outer = max_outer_iterations
      if (present(max_garch_iterations)) git = max_garch_iterations
      call initialize_gogarch(data,covariance,eigenvectors,eigenvalues,covroot,invroot,ok)
      if (.not. ok) then
         call build_gogarch_fit(data,identity_matrix(m),'nls',fit,max_garch_iterations=git,factor_spec=factor_spec)
         fit%status = 2
         return
      end if
      allocate(whitened(n,m),context%ssi(m,m,n),initial(m*(m+1)/2))
      whitened = matmul(data,invroot)
      do t = 1, n
         context%ssi(:,:,t) = outer_product(whitened(t,:),whitened(t,:))-identity_matrix(m)
      end do
      initial = 0.1_dp
      opt = nelder_mead(nls_objective_vector,initial,context,step=0.12_dp,tolerance=1.0e-7_dp,max_iterations=outer)
      b = unvech(opt%x,ok)
      call symmetric_eigen(b,values,rotation,ok)
      if (.not. ok) rotation = identity_matrix(m)
      rotation = umatch(identity_matrix(m),rotation,ok)
      call build_gogarch_fit(data,rotation,'nls',fit,max_garch_iterations=git,optimizer_iterations=opt%iterations, &
         parameters=opt%x,factor_spec=factor_spec)
   end function fit_gogarch_nls

   function fit_gogarch_ml(data, max_outer_iterations, max_garch_iterations, factor_spec) result(fit)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: max_outer_iterations, max_garch_iterations
      type(univariate_spec), intent(in), optional :: factor_spec
      type(gogarch_fit) :: fit
      type(ml_context) :: context
      type(optimizer_result) :: opt
      real(dp) :: covariance(size(data,2),size(data,2)), eigenvectors(size(data,2),size(data,2))
      real(dp) :: eigenvalues(size(data,2)), covroot(size(data,2),size(data,2)), invroot(size(data,2),size(data,2))
      real(dp), allocatable :: initial(:), angles(:)
      real(dp) :: rotation(size(data,2),size(data,2))
      logical :: ok
      integer :: m, nangles, outer, git, j
      m = size(data,2)
      nangles = m*(m-1)/2
      outer = 160
      git = 160
      if (present(max_outer_iterations)) outer = max_outer_iterations
      if (present(max_garch_iterations)) git = max_garch_iterations
      call initialize_gogarch(data,covariance,eigenvectors,eigenvalues,covroot,invroot,ok)
      if (.not. ok) then
         call build_gogarch_fit(data,identity_matrix(m),'ml',fit,max_garch_iterations=git,factor_spec=factor_spec)
         fit%status = 2
         return
      end if
      if (nangles == 0) then
         call build_gogarch_fit(data,identity_matrix(m),'ml',fit,max_garch_iterations=git, &
            parameters=[real(dp)::],factor_spec=factor_spec)
         return
      end if
      allocate(context%whitened(size(data,1),m),initial(nangles),angles(nangles))
      context%whitened = matmul(data,invroot)
      context%garch_iterations = git
      context%factor_spec = univariate_spec()
      if (present(factor_spec)) context%factor_spec = factor_spec
      initial = 0.0_dp
      opt = nelder_mead(ml_objective,initial,context,step=0.45_dp,tolerance=2.0e-6_dp,max_iterations=outer)
      do j = 1, nangles
         angles(j) = bounded_angle(opt%x(j))
      end do
      rotation = uprod_r(angles,ok)
      if (.not. ok) rotation = identity_matrix(m)
      call build_gogarch_fit(data,rotation,'ml',fit,max_garch_iterations=git,optimizer_iterations=opt%iterations, &
         parameters=angles,factor_spec=factor_spec)
   end function fit_gogarch_ml

   function gonls_objective(parameters, ssi) result(value)
      real(dp), intent(in) :: parameters(:), ssi(:,:,:)
      real(dp) :: value, b(size(ssi,1),size(ssi,2)), difference(size(ssi,1),size(ssi,2))
      logical :: ok
      integer :: t, n
      b = unvech(parameters,ok)
      if (.not. ok .or. size(ssi,3) < 2) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      n = size(ssi,3)-1
      value = 0.0_dp
      do t = 1, n
         difference = ssi(:,:,t+1)-matmul(b,matmul(ssi(:,:,t),b))
         value = value+sum(difference*difference)
      end do
      value = value/real(n,dp)
   end function gonls_objective

   function nls_objective_vector(parameters, generic_context) result(value)
      real(dp), intent(in) :: parameters(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value
      select type (context => generic_context)
      type is (nls_context)
         value = gonls_objective(parameters,context%ssi)
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function nls_objective_vector

   function ml_objective(raw_angles, generic_context) result(value)
      real(dp), intent(in) :: raw_angles(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value
      real(dp), allocatable :: angles(:), factors(:,:)
      real(dp) :: rotation(max(1,nint(0.5_dp+sqrt(0.25_dp+2.0_dp*real(size(raw_angles),dp)))), &
                           max(1,nint(0.5_dp+sqrt(0.25_dp+2.0_dp*real(size(raw_angles),dp)))))
      type(garch11_fit) :: model
      logical :: ok
      integer :: j, m
      select type (context => generic_context)
      type is (ml_context)
         m = size(context%whitened,2)
         allocate(angles(size(raw_angles)),factors(size(context%whitened,1),m))
         do j = 1, size(raw_angles)
            angles(j) = bounded_angle(raw_angles(j))
         end do
         rotation = uprod_r(angles,ok)
         if (.not. ok) then
            value = huge(1.0_dp)/100.0_dp
            return
         end if
         factors = matmul(context%whitened,rotation)
         value = 0.0_dp
         do j = 1, m
            model = fit_univariate(factors(:,j),context%factor_spec,context%garch_iterations)
            if (model%status > 1) then
               value = huge(1.0_dp)/100.0_dp
               return
            end if
            value = value-model%log_likelihood
         end do
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function ml_objective

   pure function bounded_angle(raw) result(theta)
      real(dp), intent(in) :: raw
      real(dp) :: theta, clipped
      clipped = max(-30.0_dp,min(30.0_dp,raw))
      theta = 1.0e-8_dp+(0.5_dp*pi-2.0e-8_dp)/(1.0_dp+exp(-clipped))
   end function bounded_angle

   pure function minimum_eigen_gap(values) result(gap)
      real(dp), intent(in) :: values(:)
      real(dp) :: gap
      integer :: i, j
      if (size(values) < 2) then
         gap = 1.0_dp
         return
      end if
      gap = huge(1.0_dp)
      do i = 1, size(values)-1
         do j = i+1, size(values)
            gap = min(gap,(values(i)-values(j))**2)
         end do
      end do
      if (.not. (gap < huge(1.0_dp))) gap = 0.0_dp
   end function minimum_eigen_gap

end module gogarch_estimators
