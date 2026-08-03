! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_optimization
   use qbc_kinds, only : dp
   use qbc_status, only : qbc_success, qbc_invalid_argument, qbc_no_convergence
   implicit none
   private

   abstract interface
      subroutine qbc_objective(x, value, data)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: value
         class(*), intent(in), optional :: data
      end subroutine qbc_objective
   end interface

   type, public :: qbc_optimizer_result
      real(dp), allocatable :: x(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = qbc_invalid_argument
      character(len=120) :: message = 'not optimized'
   end type qbc_optimizer_result

   public :: qbc_objective, bounded_nelder_mead

contains

   subroutine bounded_nelder_mead(fn, start, lower, upper, result, data, tolerance, max_iterations)
      procedure(qbc_objective) :: fn
      real(dp), intent(in) :: start(:), lower(:), upper(:)
      type(qbc_optimizer_result), intent(out) :: result
      class(*), intent(in), optional :: data
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:)
      real(dp) :: alpha, gamma, rho, sigma, tol, fr, fe, fc, spread_x, spread_f, step
      integer :: n, i, iter, maxit

      n = size(start)
      tol = 1.0e-9_dp
      maxit = 3000
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      if (present(max_iterations)) maxit = max(1, max_iterations)
      alpha = 1.0_dp
      gamma = 2.0_dp
      rho = 0.5_dp
      sigma = 0.5_dp

      if (n == 0 .or. size(lower) /= n .or. size(upper) /= n .or. any(lower > upper)) then
         allocate(result%x(0))
         result%status = qbc_invalid_argument
         result%message = 'invalid bounds or starting vector'
         return
      end if

      allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n),result%x(n))
      simplex(:,1) = min(max(start,lower),upper)
      do i=1,n
         simplex(:,i+1) = simplex(:,1)
         step = 0.05_dp * max(1.0_dp,abs(simplex(i,1)))
         if (upper(i) > lower(i)) step = min(step,0.1_dp*(upper(i)-lower(i)))
         simplex(i,i+1) = min(upper(i),simplex(i,i+1)+step)
         if (abs(simplex(i,i+1)-simplex(i,1)) <= epsilon(1.0_dp)) &
            simplex(i,i+1) = max(lower(i),simplex(i,i+1)-step)
      end do
      do i=1,n+1
         call evaluate(fn,simplex(:,i),f(i),data)
      end do
      result%evaluations = n+1

      do iter=1,maxit
         call sort_simplex(simplex,f)
         spread_x = maxval(abs(simplex(:,2:n+1)-spread(simplex(:,1),2,n)))
         spread_f = maxval(abs(f(2:n+1)-f(1)))
         if (spread_x <= tol*(1.0_dp+maxval(abs(simplex(:,1)))) .and. &
             spread_f <= tol*(1.0_dp+abs(f(1)))) exit

         centroid = sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr = project(centroid + alpha*(centroid-simplex(:,n+1)),lower,upper)
         call evaluate(fn,xr,fr,data); result%evaluations=result%evaluations+1

         if (fr < f(1)) then
            xe = project(centroid + gamma*(xr-centroid),lower,upper)
            call evaluate(fn,xe,fe,data); result%evaluations=result%evaluations+1
            if (fe < fr) then
               simplex(:,n+1)=xe; f(n+1)=fe
            else
               simplex(:,n+1)=xr; f(n+1)=fr
            end if
         else if (fr < f(n)) then
            simplex(:,n+1)=xr; f(n+1)=fr
         else
            if (fr < f(n+1)) then
               xc = project(centroid + rho*(xr-centroid),lower,upper)
            else
               xc = project(centroid + rho*(simplex(:,n+1)-centroid),lower,upper)
            end if
            call evaluate(fn,xc,fc,data); result%evaluations=result%evaluations+1
            if (fc < min(fr,f(n+1))) then
               simplex(:,n+1)=xc; f(n+1)=fc
            else
               do i=2,n+1
                  simplex(:,i)=project(simplex(:,1)+sigma*(simplex(:,i)-simplex(:,1)),lower,upper)
                  call evaluate(fn,simplex(:,i),f(i),data)
               end do
               result%evaluations=result%evaluations+n
            end if
         end if
      end do

      call sort_simplex(simplex,f)
      result%x = simplex(:,1)
      result%value = f(1)
      result%iterations = min(iter,maxit)
      if (iter <= maxit) then
         result%status = qbc_success
         result%message = 'converged'
      else
         result%status = qbc_no_convergence
         result%message = 'maximum iterations reached'
      end if
   end subroutine bounded_nelder_mead

   subroutine evaluate(fn,x,value,data)
      procedure(qbc_objective) :: fn
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      if (present(data)) then
         call fn(x,value,data)
      else
         call fn(x,value)
      end if
      if (.not. (value < huge(1.0_dp))) value = huge(1.0_dp)
   end subroutine evaluate

   pure function project(x,lower,upper) result(y)
      real(dp), intent(in) :: x(:),lower(:),upper(:)
      real(dp) :: y(size(x))
      y=min(max(x,lower),upper)
   end function project

   subroutine sort_simplex(simplex,f)
      real(dp), intent(inout) :: simplex(:,:),f(:)
      real(dp) :: key, col(size(simplex,1))
      integer :: i,j
      do i=2,size(f)
         key=f(i);col=simplex(:,i);j=i-1
         do while(j>=1)
            if(f(j)<=key) exit
            f(j+1)=f(j);simplex(:,j+1)=simplex(:,j);j=j-1
         end do
         f(j+1)=key;simplex(:,j+1)=col
      end do
   end subroutine sort_simplex

end module qbc_optimization
