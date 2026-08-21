! Experimental modern Fortran translation of computational routines from
! the R package tseries 0.10-62. Original authors include Adrian Trapletti
! and Kurt Hornik; Blake LeBaron contributed the original BDS code.
! Licensed under GPL-2.0-only OR GPL-3.0-only. See LICENSE and NOTICE.

module tseries_optimize
   use tseries_kinds, only : dp
   implicit none
   private

   public :: objective_function
   public :: nelder_mead
   public :: numerical_hessian

   abstract interface
      function objective_function(x) result(value)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface

contains

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
