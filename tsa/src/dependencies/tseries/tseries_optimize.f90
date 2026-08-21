! Experimental modern Fortran translation of computational routines from
! the R package tseries 0.10-62. Original authors include Adrian Trapletti
! and Kurt Hornik; Blake LeBaron contributed the original BDS code.
! Licensed under GPL-2.0-only OR GPL-3.0-only. See LICENSE and NOTICE.

module tseries_optimize
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use tseries_kinds, only : dp
   implicit none
   private

   public :: objective_function
   public :: nelder_mead
   public :: bfgs
   public :: numerical_hessian
   public :: optim_hessian

   abstract interface
      function objective_function(x) result(value)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface

contains

   subroutine bfgs(fun,x,fvalue,iterations,status,max_iterations,reltol,abstol,parscale,ndeps)
      ! BFGS variable-metric minimizer following the algorithm used by R's
      ! vmmin()/optim(method="BFGS").  Optimization is performed in scaled
      ! coordinates b = x/parscale; numerical derivatives use ndeps in those
      ! scaled coordinates, matching the R optim interface.
      procedure(objective_function) :: fun
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out) :: fvalue
      integer, intent(out) :: iterations,status
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: reltol,abstol,parscale(:),ndeps(:)
      real(dp), allocatable :: b(:), btrial(:), scale(:), deps(:)
      real(dp), allocatable :: g(:), gold(:), direction(:), stepvec(:), by(:)
      real(dp), allocatable :: hessinv(:,:), yvec(:), actual(:)
      real(dp) :: fmin, ftrial, gradproj, steplen, rt, at, d1, d2
      real(dp), parameter :: stepredn = 0.2_dp
      real(dp), parameter :: acctol = 1.0e-4_dp
      real(dp), parameter :: reltest = 10.0_dp
      integer :: n, maxit, i, j, gradcount, ilast, nunchanged
      logical :: accepted, enough, converged

      n = size(x)
      maxit = 100
      if (present(max_iterations)) maxit = max_iterations
      rt = sqrt(epsilon(1.0_dp))
      if (present(reltol)) rt = reltol
      at = -huge(1.0_dp)
      if (present(abstol)) at = abstol

      if (n == 0) then
         fvalue = fun(x)
         iterations = 0
         status = 0
         return
      end if
      if (maxit <= 0 .or. rt < 0.0_dp) then
         fvalue = huge(1.0_dp)
         iterations = 0
         status = 2
         return
      end if

      allocate(b(n),btrial(n),scale(n),deps(n),g(n),gold(n),direction(n), &
         stepvec(n),by(n),hessinv(n,n),yvec(n),actual(n))
      scale = 1.0_dp
      if (present(parscale)) then
         if (size(parscale) /= n .or. any(parscale <= 0.0_dp) .or. &
            any(.not. ieee_is_finite(parscale))) then
            fvalue = huge(1.0_dp)
            iterations = 0
            status = 2
            return
         end if
         scale = parscale
      end if
      deps = 1.0e-3_dp
      if (present(ndeps)) then
         if (size(ndeps) /= n .or. any(ndeps <= 0.0_dp) .or. &
            any(.not. ieee_is_finite(ndeps))) then
            fvalue = huge(1.0_dp)
            iterations = 0
            status = 2
            return
         end if
         deps = ndeps
      end if

      b = x/scale
      call scaled_value(b,fmin)
      if (.not. ieee_is_finite(fmin)) then
         fvalue = fmin
         iterations = 0
         status = 3
         return
      end if
      call scaled_gradient(b,g,status)
      if (status /= 0) then
         fvalue = fmin
         iterations = 0
         return
      end if

      hessinv = 0.0_dp
      do i = 1, n
         hessinv(i,i) = 1.0_dp
      end do
      gradcount = 1
      ilast = gradcount
      iterations = 0
      converged = .false.

      do while (iterations < maxit)
         if (ilast == gradcount) then
            hessinv = 0.0_dp
            do i = 1, n
               hessinv(i,i) = 1.0_dp
            end do
         end if

         gold = g
         direction = -matmul(hessinv,g)
         gradproj = dot_product(direction,g)

         if (gradproj < 0.0_dp .and. ieee_is_finite(gradproj)) then
            steplen = 1.0_dp
            accepted = .false.
            nunchanged = 0
            do
               btrial = b + steplen*direction
               nunchanged = 0
               do i = 1, n
                  if (.not. (reltest+b(i) < reltest+btrial(i) .or. &
                     reltest+b(i) > reltest+btrial(i))) nunchanged = nunchanged+1
               end do
               if (nunchanged == n) exit

               call scaled_value(btrial,ftrial)
               accepted = ieee_is_finite(ftrial) .and. &
                  ftrial <= fmin + gradproj*steplen*acctol
               if (accepted) exit
               steplen = steplen*stepredn
            end do

            if (nunchanged == n) then
               if (ilast == gradcount) then
                  converged = .true.
                  exit
               else
                  ilast = gradcount
                  cycle
               end if
            end if

            enough = (ftrial > at) .and. &
               (abs(ftrial-fmin) > rt*(abs(fmin)+rt))
            b = btrial
            fmin = ftrial
            iterations = iterations + 1
            if (.not. enough) then
               converged = .true.
               exit
            end if

            call scaled_gradient(b,g,status)
            if (status /= 0) exit
            gradcount = gradcount + 1
            stepvec = steplen*direction
            yvec = g-gold
            d1 = dot_product(stepvec,yvec)
            if (d1 > 0.0_dp .and. ieee_is_finite(d1)) then
               by = matmul(hessinv,yvec)
               d2 = 1.0_dp + dot_product(yvec,by)/d1
               do j = 1, n
                  do i = 1, n
                     hessinv(i,j) = hessinv(i,j) + &
                        (d2*stepvec(i)*stepvec(j)-by(i)*stepvec(j)- &
                        stepvec(i)*by(j))/d1
                  end do
               end do
               hessinv = 0.5_dp*(hessinv+transpose(hessinv))
            else
               ilast = gradcount
            end if
            if (gradcount-ilast > 2*n) ilast = gradcount
         else
            if (ilast == gradcount) then
               converged = .true.
               exit
            else
               ilast = gradcount
            end if
         end if
      end do

      x = b*scale
      fvalue = fmin
      if (status /= 0) return
      if (converged) then
         status = 0
      else
         status = 1
      end if

   contains

      subroutine scaled_value(bs,value)
         real(dp), intent(in) :: bs(:)
         real(dp), intent(out) :: value
         actual = bs*scale
         value = fun(actual)
      end subroutine scaled_value

      subroutine scaled_gradient(bs,grad,istat)
         real(dp), intent(in) :: bs(:)
         real(dp), intent(out) :: grad(:)
         integer, intent(out) :: istat
         real(dp) :: fp,fm,old
         integer :: k

         istat = 0
         btrial = bs
         do k = 1, n
            old = btrial(k)
            btrial(k) = old + deps(k)
            call scaled_value(btrial,fp)
            btrial(k) = old - deps(k)
            call scaled_value(btrial,fm)
            btrial(k) = old
            if (.not. ieee_is_finite(fp) .or. .not. ieee_is_finite(fm)) then
               istat = 3
               return
            end if
            grad(k) = (fp-fm)/(2.0_dp*deps(k))
         end do
      end subroutine scaled_gradient

   end subroutine bfgs

   subroutine nelder_mead(fun,x,fvalue,iterations,status,max_iterations,tolerance,step)
      procedure(objective_function) :: fun
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out) :: fvalue
      integer, intent(out) :: iterations,status
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance,step
      real(dp), allocatable :: simplex(:, :), values(:), centroid(:), xr(:),xe(:),xc(:),tmp(:)
      real(dp) :: tol,initial_step,fr,fe,fc,spread
      integer :: n,maxit,j,best,worst,second_worst

      n=size(x); maxit=2000; tol=1.0e-8_dp; initial_step=0.05_dp
      if(present(max_iterations)) maxit=max_iterations
      if(present(tolerance)) tol=tolerance
      if(present(step)) initial_step=step
      allocate(simplex(n,n+1),values(n+1),centroid(n),xr(n),xe(n),xc(n),tmp(n))
      simplex(:,1)=x
      do j=1,n
         simplex(:,j+1)=x
         simplex(j,j+1)=simplex(j,j+1)+initial_step*max(1.0_dp,abs(x(j)))
      end do
      do j=1,n+1
         values(j)=fun(simplex(:,j))
      end do
      status=1
      do iterations=1,maxit
         call order_simplex(simplex,values)
         best=1; worst=n+1; second_worst=n
         spread=maxval(abs(values-values(best)))
         if(spread<=tol*(1.0_dp+abs(values(best))) .and. &
            maxval(abs(simplex-spread_columns(simplex(:,best),n+1)))<=sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,best))))) then
            status=0; exit
         end if
         centroid=sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr=centroid+(centroid-simplex(:,worst))
         fr=fun(xr)
         if(fr<values(best)) then
            xe=centroid+2.0_dp*(xr-centroid)
            fe=fun(xe)
            if(fe<fr) then
               simplex(:,worst)=xe; values(worst)=fe
            else
               simplex(:,worst)=xr; values(worst)=fr
            end if
         else if(fr<values(second_worst)) then
            simplex(:,worst)=xr; values(worst)=fr
         else
            if(fr<values(worst)) then
               xc=centroid+0.5_dp*(xr-centroid)
            else
               xc=centroid+0.5_dp*(simplex(:,worst)-centroid)
            end if
            fc=fun(xc)
            if(fc<min(fr,values(worst))) then
               simplex(:,worst)=xc; values(worst)=fc
            else
               do j=2,n+1
                  simplex(:,j)=simplex(:,best)+0.5_dp*(simplex(:,j)-simplex(:,best))
                  values(j)=fun(simplex(:,j))
               end do
            end if
         end if
      end do
      call order_simplex(simplex,values)
      x=simplex(:,1); fvalue=values(1)
      if(iterations>maxit) iterations=maxit
   end subroutine nelder_mead

   function spread_columns(column,ncol) result(matrix)
      real(dp), intent(in) :: column(:)
      integer, intent(in) :: ncol
      real(dp) :: matrix(size(column),ncol)
      integer :: j
      do j=1,ncol
         matrix(:,j)=column
      end do
   end function spread_columns

   subroutine order_simplex(simplex,values)
      real(dp), intent(inout) :: simplex(:, :),values(:)
      real(dp), allocatable :: col(:)
      real(dp) :: tv
      integer :: i,j,k,m
      m=size(values); allocate(col(size(simplex,1)))
      do i=1,m-1
         k=i
         do j=i+1,m
            if(values(j)<values(k)) k=j
         end do
         if(k/=i) then
            tv=values(i); values(i)=values(k); values(k)=tv
            col=simplex(:,i); simplex(:,i)=simplex(:,k); simplex(:,k)=col
         end if
      end do
   end subroutine order_simplex


   subroutine optim_hessian(fun,x,hessian,status,parscale,ndeps)
      ! Numerical Hessian matching R's C_optimhess path when no analytical
      ! gradient is supplied. The outer perturbation is ndeps in the actual
      ! parameter coordinate; each inner gradient uses ndeps in the scaled
      ! BFGS coordinate, exactly as optimHess/optim(hessian=TRUE) does.
      procedure(objective_function) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hessian(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: parscale(:), ndeps(:)
      real(dp), allocatable :: scale(:), deps(:), xp(:), xm(:), gp(:), gm(:)
      integer :: n, i, j

      n = size(x)
      status = 0
      hessian = 0.0_dp
      if (size(hessian,1) /= n .or. size(hessian,2) /= n) then
         status = 1
         return
      end if
      allocate(scale(n),deps(n),xp(n),xm(n),gp(n),gm(n))
      scale = 1.0_dp
      deps = 1.0e-3_dp
      if (present(parscale)) then
         if (size(parscale) /= n .or. any(parscale <= 0.0_dp) .or. &
            any(.not. ieee_is_finite(parscale))) then
            status = 2
            return
         end if
         scale = parscale
      end if
      if (present(ndeps)) then
         if (size(ndeps) /= n .or. any(ndeps <= 0.0_dp) .or. &
            any(.not. ieee_is_finite(ndeps))) then
            status = 2
            return
         end if
         deps = ndeps
      end if

      do i = 1, n
         xp = x
         xm = x
         xp(i) = xp(i) + deps(i)
         xm(i) = xm(i) - deps(i)
         call scaled_numeric_gradient(xp,gp,status)
         if (status /= 0) return
         call scaled_numeric_gradient(xm,gm,status)
         if (status /= 0) return
         do j = 1, n
            hessian(j,i) = (gp(j)-gm(j))/(2.0_dp*deps(i)*scale(j))
         end do
      end do
      hessian = 0.5_dp*(hessian+transpose(hessian))

   contains

      subroutine scaled_numeric_gradient(point,grad,istat)
         real(dp), intent(in) :: point(:)
         real(dp), intent(out) :: grad(:)
         integer, intent(out) :: istat
         real(dp), allocatable :: pplus(:), pminus(:)
         real(dp) :: fp, fm, h
         integer :: k

         istat = 0
         allocate(pplus(n),pminus(n))
         do k = 1, n
            pplus = point
            pminus = point
            h = deps(k)*scale(k)
            pplus(k) = pplus(k)+h
            pminus(k) = pminus(k)-h
            fp = fun(pplus)
            fm = fun(pminus)
            if (.not. ieee_is_finite(fp) .or. .not. ieee_is_finite(fm)) then
               istat = 3
               return
            end if
            grad(k) = (fp-fm)/(2.0_dp*deps(k))
         end do
      end subroutine scaled_numeric_gradient

   end subroutine optim_hessian

   subroutine numerical_hessian(fun,x,hessian,status)
      procedure(objective_function) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hessian(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:),h(:)
      real(dp) :: f0
      integer :: n,i,j
      n=size(x); status=0
      if(size(hessian,1)/=n .or. size(hessian,2)/=n) then
         status=1; return
      end if
      allocate(xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n),h(n))
      do i=1,n
         h(i)=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
      end do
      f0=fun(x)
      do i=1,n
         xp=x; xm=x
         xp(i)=xp(i)+h(i); xm(i)=xm(i)-h(i)
         hessian(i,i)=(fun(xp)-2.0_dp*f0+fun(xm))/(h(i)*h(i))
         do j=i+1,n
            xpp=x; xpm=x; xmp=x; xmm=x
            xpp(i)=xpp(i)+h(i); xpp(j)=xpp(j)+h(j)
            xpm(i)=xpm(i)+h(i); xpm(j)=xpm(j)-h(j)
            xmp(i)=xmp(i)-h(i); xmp(j)=xmp(j)+h(j)
            xmm(i)=xmm(i)-h(i); xmm(j)=xmm(j)-h(j)
            hessian(i,j)=(fun(xpp)-fun(xpm)-fun(xmp)+fun(xmm))/(4.0_dp*h(i)*h(j))
            hessian(j,i)=hessian(i,j)
         end do
      end do
   end subroutine numerical_hessian

end module tseries_optimize
