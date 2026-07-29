! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_optimize
   use ghyp_kinds, only : dp
   implicit none
   private
   public :: nelder_mead, numerical_hessian

   abstract interface
      function objective_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface

contains

   subroutine nelder_mead(objective, start, optimum, value, converged, iterations, &
      max_iter, tolerance)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: start(:)
      real(dp), allocatable, intent(out) :: optimum(:)
      real(dp), intent(out) :: value
      logical, intent(out) :: converged
      integer, intent(out) :: iterations
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
      real(dp) :: fr, fe, fc, tol, scale
      integer :: n, i, j, best, worst, second_worst, maxit

      n = size(start)
      maxit = 1500
      if (present(max_iter)) maxit = max(1,max_iter)
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(tolerance,epsilon(1.0_dp))
      allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n),optimum(n))
      simplex(:,1) = start
      do j = 1, n
         simplex(:,j+1) = start
         scale = 0.05_dp*max(1.0_dp,abs(start(j)))
         simplex(j,j+1) = simplex(j,j+1)+scale
      end do
      do i = 1, n+1
         f(i) = objective(simplex(:,i))
      end do
      converged = .false.
      do iterations = 1, maxit
         best = minloc(f,dim=1)
         worst = maxloc(f,dim=1)
         second_worst = best
         do i = 1, n+1
            if (i /= worst) then
               if (second_worst == best .or. f(i) > f(second_worst)) second_worst = i
            end if
         end do
         if (maxval(abs(f-f(best))) <= tol*(1.0_dp+abs(f(best))) .and. &
             maxval(abs(simplex-spread(simplex(:,best),2,n+1))) <= &
             sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,best))))) then
            converged = .true.
            exit
         end if
         centroid = 0.0_dp
         do i = 1, n+1
            if (i /= worst) centroid = centroid+simplex(:,i)
         end do
         centroid = centroid/real(n,dp)
         xr = centroid+(centroid-simplex(:,worst))
         fr = objective(xr)
         if (fr < f(best)) then
            xe = centroid+2.0_dp*(xr-centroid)
            fe = objective(xe)
            if (fe < fr) then
               simplex(:,worst) = xe; f(worst) = fe
            else
               simplex(:,worst) = xr; f(worst) = fr
            end if
         else if (fr < f(second_worst)) then
            simplex(:,worst) = xr; f(worst) = fr
         else
            if (fr < f(worst)) then
               xc = centroid+0.5_dp*(xr-centroid)
            else
               xc = centroid+0.5_dp*(simplex(:,worst)-centroid)
            end if
            fc = objective(xc)
            if (fc < min(fr,f(worst))) then
               simplex(:,worst) = xc; f(worst) = fc
            else
               do i = 1, n+1
                  if (i /= best) then
                     simplex(:,i) = simplex(:,best)+0.5_dp*(simplex(:,i)-simplex(:,best))
                     f(i) = objective(simplex(:,i))
                  end if
               end do
            end if
         end if
      end do
      best = minloc(f,dim=1)
      optimum = simplex(:,best)
      value = f(best)
   end subroutine nelder_mead

   subroutine numerical_hessian(objective, x, hessian)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: hessian(:,:)
      real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
      real(dp) :: hi, hj, f0
      integer :: i, j, n
      n = size(x)
      allocate(hessian(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n))
      f0 = objective(x)
      do i = 1, n
         hi = epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
         xp = x; xm = x
         xp(i) = xp(i)+hi; xm(i) = xm(i)-hi
         hessian(i,i) = (objective(xp)-2.0_dp*f0+objective(xm))/(hi*hi)
         do j = i+1, n
            hj = epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
            xpp=x; xpm=x; xmp=x; xmm=x
            xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
            hessian(i,j) = (objective(xpp)-objective(xpm)-objective(xmp)+ &
               objective(xmm))/(4.0_dp*hi*hj)
            hessian(j,i) = hessian(i,j)
         end do
      end do
   end subroutine numerical_hessian

end module ghyp_optimize
