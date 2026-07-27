! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from RND 1.2, Copyright (C) 2017 Kam Hamidieh.
module rnd_optimize
   use rnd_kinds, only : dp
   use rnd_types, only : optimizer_result
   implicit none
   private
   public :: objective_function, nelder_mead, numerical_hessian

   abstract interface
      function objective_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface

contains

   function nelder_mead(objective, initial, max_iter, tolerance, initial_step, compute_hessian) result(result)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: initial(:)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tolerance, initial_step
      logical, intent(in), optional :: compute_hessian
      type(optimizer_result) :: result
      integer :: n, limit, iter, i
      real(dp) :: tol, step, f_reflect, f_expand, f_contract
      real(dp), allocatable :: simplex(:, :), values(:), centroid(:), reflected(:), expanded(:), contracted(:)
      logical :: want_hessian
      real(dp), parameter :: alpha = 1.0_dp, gamma = 2.0_dp, rho = 0.5_dp, shrink = 0.5_dp

      n = size(initial)
      limit = 5000
      if (present(max_iter)) limit = max_iter
      tol = 1.0e-9_dp
      if (present(tolerance)) tol = tolerance
      step = 0.05_dp
      if (present(initial_step)) step = initial_step
      want_hessian = .false.
      if (present(compute_hessian)) want_hessian = compute_hessian

      allocate(simplex(n,n+1), values(n+1), centroid(n), reflected(n), expanded(n), contracted(n))
      simplex(:,1) = initial
      do i = 1, n
         simplex(:,i+1) = initial
         simplex(i,i+1) = initial(i)+step*max(1.0_dp, abs(initial(i)))
      end do
      do i = 1, n+1
         values(i) = finite_objective(objective, simplex(:,i))
      end do

      result%convergence = 1
      do iter = 1, limit
         call sort_simplex(simplex, values)
         if (maxval(abs(values-values(1))) <= tol*(1.0_dp+abs(values(1))) .and. &
               maxval(abs(simplex-spread(simplex(:,1),dim=2,ncopies=n+1))) &
               <= sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,1))))) then
            result%convergence = 0
            exit
         end if
         centroid = sum(simplex(:,1:n),dim=2)/real(n,dp)
         reflected = centroid+alpha*(centroid-simplex(:,n+1))
         f_reflect = finite_objective(objective, reflected)
         if (f_reflect < values(1)) then
            expanded = centroid+gamma*(reflected-centroid)
            f_expand = finite_objective(objective, expanded)
            if (f_expand < f_reflect) then
               simplex(:,n+1) = expanded
               values(n+1) = f_expand
            else
               simplex(:,n+1) = reflected
               values(n+1) = f_reflect
            end if
         else if (f_reflect < values(n)) then
            simplex(:,n+1) = reflected
            values(n+1) = f_reflect
         else
            if (f_reflect < values(n+1)) then
               contracted = centroid+rho*(reflected-centroid)
            else
               contracted = centroid-rho*(centroid-simplex(:,n+1))
            end if
            f_contract = finite_objective(objective, contracted)
            if (f_contract < min(values(n+1),f_reflect)) then
               simplex(:,n+1) = contracted
               values(n+1) = f_contract
            else
               do i = 2, n+1
                  simplex(:,i) = simplex(:,1)+shrink*(simplex(:,i)-simplex(:,1))
                  values(i) = finite_objective(objective, simplex(:,i))
               end do
            end if
         end if
      end do
      call sort_simplex(simplex, values)
      allocate(result%par(n))
      result%par = simplex(:,1)
      result%value = values(1)
      result%iterations = min(iter,limit)
      if (want_hessian) then
         result%hessian = numerical_hessian(objective, result%par)
      else
         allocate(result%hessian(0,0))
      end if
   end function nelder_mead

   function numerical_hessian(objective, x) result(hessian)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:), h(:)
      real(dp) :: f0
      integer :: i, j, n
      n = size(x)
      allocate(hessian(n,n), xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n), h(n))
      h = epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x))
      f0 = finite_objective(objective,x)
      hessian = 0.0_dp
      do i = 1, n
         xp = x
         xm = x
         xp(i) = xp(i)+h(i)
         xm(i) = xm(i)-h(i)
         hessian(i,i) = (finite_objective(objective,xp)-2.0_dp*f0 &
            +finite_objective(objective,xm))/(h(i)*h(i))
         do j = i+1, n
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = xpp(i)+h(i)
            xpp(j) = xpp(j)+h(j)
            xpm(i) = xpm(i)+h(i)
            xpm(j) = xpm(j)-h(j)
            xmp(i) = xmp(i)-h(i)
            xmp(j) = xmp(j)+h(j)
            xmm(i) = xmm(i)-h(i)
            xmm(j) = xmm(j)-h(j)
            hessian(i,j) = (finite_objective(objective,xpp)-finite_objective(objective,xpm) &
               -finite_objective(objective,xmp)+finite_objective(objective,xmm))/(4.0_dp*h(i)*h(j))
            hessian(j,i) = hessian(i,j)
         end do
      end do
   end function numerical_hessian

   real(dp) function finite_objective(objective, x) result(value)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      procedure(objective_function) :: objective
      real(dp), intent(in) :: x(:)
      value = objective(x)
      if (.not. ieee_is_finite(value)) value = huge(1.0_dp)/1000.0_dp
   end function finite_objective

   subroutine sort_simplex(simplex, values)
      real(dp), intent(inout) :: simplex(:, :), values(:)
      real(dp), allocatable :: tmp_column(:)
      real(dp) :: tmp_value
      integer :: i, j, minimum
      allocate(tmp_column(size(simplex,1)))
      do i = 1, size(values)-1
         minimum = i
         do j = i+1, size(values)
            if (values(j) < values(minimum)) minimum = j
         end do
         if (minimum /= i) then
            tmp_value = values(i)
            values(i) = values(minimum)
            values(minimum) = tmp_value
            tmp_column = simplex(:,i)
            simplex(:,i) = simplex(:,minimum)
            simplex(:,minimum) = tmp_column
         end if
      end do
   end subroutine sort_simplex

end module rnd_optimize
