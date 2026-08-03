! SPDX-License-Identifier: Artistic-2.0
module mts_optimize
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mts_kinds, only : dp
   use mts_types, only : mts_success, mts_no_convergence
   implicit none
   private
   abstract interface
      function objective_function(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_function
   end interface
   public :: objective_function, nelder_mead, bfgs_minimize, numerical_gradient, numerical_hessian
contains

   subroutine numerical_gradient(f,x,g,step)
      procedure(objective_function) :: f
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(size(x))
      real(dp), intent(in), optional :: step
      real(dp) :: xp(size(x)), xm(size(x)), h
      integer :: i
      do i = 1, size(x)
         h = sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x(i)))
         if (present(step)) h = max(step,h)
         xp = x
         xm = x
         xp(i) = xp(i)+h
         xm(i) = xm(i)-h
         g(i) = (f(xp)-f(xm))/(2.0_dp*h)
      end do
   end subroutine numerical_gradient

   subroutine numerical_hessian(f,x,hess,step)
      procedure(objective_function) :: f
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hess(size(x),size(x))
      real(dp), intent(in), optional :: step
      real(dp) :: xpp(size(x)), xpm(size(x)), xmp(size(x)), xmm(size(x))
      real(dp) :: hi, hj, f0
      integer :: i, j
      f0 = f(x)
      hess = 0.0_dp
      do i = 1, size(x)
         hi = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
         if (present(step)) hi = max(step,hi)
         xpp = x
         xpm = x
         xpp(i) = xpp(i)+hi
         xpm(i) = xpm(i)-hi
         hess(i,i) = (f(xpp)-2.0_dp*f0+f(xpm))/(hi*hi)
         do j = i+1, size(x)
            hj = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(j)))
            if (present(step)) hj = max(step,hj)
            xpp = x; xpm = x; xmp = x; xmm = x
            xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
            hess(i,j) = (f(xpp)-f(xpm)-f(xmp)+f(xmm))/(4.0_dp*hi*hj)
            hess(j,i) = hess(i,j)
         end do
      end do
   end subroutine numerical_hessian

   subroutine nelder_mead(f,x,value,status,iterations,max_iterations,tolerance,lower,upper)
      procedure(objective_function) :: f
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out) :: value
      integer, intent(out), optional :: status, iterations
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance, lower(:), upper(:)
      real(dp), allocatable :: simplex(:,:), fv(:), centroid(:), xr(:), xe(:), xc(:)
      real(dp), allocatable :: tempcol(:)
      real(dp) :: alpha, gamma, rho, shrink, tol, fspread, step
      integer :: n, i, j, iter, maxit, best, worst, second_worst, istat

      n = size(x)
      allocate(simplex(n,n+1),fv(n+1),centroid(n),xr(n),xe(n),xc(n),tempcol(n))
      alpha=1.0_dp; gamma=2.0_dp; rho=0.5_dp; shrink=0.5_dp
      maxit = 1000
      if (present(max_iterations)) maxit = max_iterations
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      simplex(:,1) = projected(x,lower,upper)
      do j = 2, n+1
         simplex(:,j) = simplex(:,1)
         step = 0.05_dp*max(1.0_dp,abs(x(j-1)))
         simplex(j-1,j) = simplex(j-1,j)+step
         simplex(:,j) = projected(simplex(:,j),lower,upper)
      end do
      do j = 1, n+1
         fv(j) = safe_value(f,simplex(:,j))
      end do
      istat = mts_no_convergence
      do iter = 1, maxit
         call sort_simplex(simplex,fv)
         best = 1; worst = n+1; second_worst = n
         fspread = maxval(abs(fv-fv(best)))
         if (fspread <= tol*(1.0_dp+abs(fv(best))) .and. &
             maxval(abs(simplex-spread(simplex(:,best),dim=2,ncopies=n+1))) <= sqrt(tol)*(1.0_dp+maxval(abs(x)))) then
            istat = mts_success
            exit
         end if
         centroid = sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr = projected(centroid+alpha*(centroid-simplex(:,worst)),lower,upper)
         if (safe_value(f,xr) < fv(best)) then
            xe = projected(centroid+gamma*(xr-centroid),lower,upper)
            if (safe_value(f,xe) < safe_value(f,xr)) then
               simplex(:,worst)=xe; fv(worst)=safe_value(f,xe)
            else
               simplex(:,worst)=xr; fv(worst)=safe_value(f,xr)
            end if
         else if (safe_value(f,xr) < fv(second_worst)) then
            simplex(:,worst)=xr; fv(worst)=safe_value(f,xr)
         else
            if (safe_value(f,xr) < fv(worst)) then
               xc = projected(centroid+rho*(xr-centroid),lower,upper)
            else
               xc = projected(centroid-rho*(centroid-simplex(:,worst)),lower,upper)
            end if
            if (safe_value(f,xc) < min(fv(worst),safe_value(f,xr))) then
               simplex(:,worst)=xc; fv(worst)=safe_value(f,xc)
            else
               do j = 2, n+1
                  simplex(:,j)=projected(simplex(:,best)+shrink*(simplex(:,j)-simplex(:,best)),lower,upper)
                  fv(j)=safe_value(f,simplex(:,j))
               end do
            end if
         end if
      end do
      call sort_simplex(simplex,fv)
      x = simplex(:,1)
      value = fv(1)
      if (present(status)) status = istat
      if (present(iterations)) iterations = min(iter,maxit)
   end subroutine nelder_mead

   subroutine bfgs_minimize(f,x,value,status,iterations,max_iterations,tolerance,lower,upper)
      procedure(objective_function) :: f
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out) :: value
      integer, intent(out), optional :: status, iterations
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance, lower(:), upper(:)
      real(dp), allocatable :: h(:,:), g(:), gnew(:), p(:), xnew(:), s(:), y(:)
      real(dp) :: alpha, fnew, ys, tol
      integer :: n, i, iter, maxit, istat
      n = size(x)
      allocate(h(n,n),g(n),gnew(n),p(n),xnew(n),s(n),y(n))
      h = 0.0_dp
      do i = 1, n
         h(i,i)=1.0_dp
      end do
      maxit=500
      if (present(max_iterations)) maxit=max_iterations
      tol=1.0e-7_dp
      if (present(tolerance)) tol=tolerance
      x=projected(x,lower,upper)
      value=safe_value(f,x)
      call numerical_gradient(f,x,g)
      istat=mts_no_convergence
      do iter=1,maxit
         if (maxval(abs(projected_gradient(x,g,lower,upper))) <= tol) then
            istat=mts_success
            exit
         end if
         p=-matmul(h,g)
         if (dot_product(p,g) >= -epsilon(1.0_dp)) p=-g
         alpha=1.0_dp
         do
            xnew=projected(x+alpha*p,lower,upper)
            fnew=safe_value(f,xnew)
            if (fnew <= value+1.0e-4_dp*dot_product(g,xnew-x)) exit
            alpha=0.5_dp*alpha
            if (alpha < 1.0e-12_dp) exit
         end do
         if (alpha < 1.0e-12_dp) exit
         call numerical_gradient(f,xnew,gnew)
         s=xnew-x
         y=gnew-g
         ys=dot_product(y,s)
         if (ys > 1.0e-12_dp*sqrt(max(dot_product(y,y)*dot_product(s,s),tiny(1.0_dp)))) then
            h = bfgs_update(h,s,y,ys)
         else
            h=0.0_dp
            do i=1,n
               h(i,i)=1.0_dp
            end do
         end if
         x=xnew; value=fnew; g=gnew
         if (maxval(abs(s)) <= tol*(1.0_dp+maxval(abs(x)))) then
            istat=mts_success
            exit
         end if
      end do
      if (present(status)) status=istat
      if (present(iterations)) iterations=min(iter,maxit)
   end subroutine bfgs_minimize

   function bfgs_update(h,s,y,ys) result(hnew)
      real(dp), intent(in) :: h(:,:), s(:), y(:), ys
      real(dp) :: hnew(size(h,1),size(h,2)), hy(size(y)), yhy, rho
      integer :: i,j
      hy=matmul(h,y)
      yhy=dot_product(y,hy)
      rho=1.0_dp/ys
      hnew=h
      do i=1,size(s)
         do j=1,size(s)
            hnew(i,j)=h(i,j)+(1.0_dp+yhy*rho)*rho*s(i)*s(j)-rho*(s(i)*hy(j)+hy(i)*s(j))
         end do
      end do
   end function bfgs_update

   function projected(x,lower,upper) result(y)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: lower(:), upper(:)
      real(dp) :: y(size(x))
      y=x
      if (present(lower)) y=max(y,lower)
      if (present(upper)) y=min(y,upper)
   end function projected

   function projected_gradient(x,g,lower,upper) result(pg)
      real(dp), intent(in) :: x(:),g(:)
      real(dp), intent(in), optional :: lower(:),upper(:)
      real(dp) :: pg(size(x))
      integer :: i
      pg=g
      do i=1,size(x)
         if (present(lower)) then
            if (x(i) <= lower(i)+1.0e-12_dp .and. g(i)>0.0_dp) pg(i)=0.0_dp
         end if
         if (present(upper)) then
            if (x(i) >= upper(i)-1.0e-12_dp .and. g(i)<0.0_dp) pg(i)=0.0_dp
         end if
      end do
   end function projected_gradient

   function safe_value(f,x) result(v)
      procedure(objective_function) :: f
      real(dp), intent(in) :: x(:)
      real(dp) :: v
      v=f(x)
      if (.not. ieee_is_finite(v)) v=huge(1.0_dp)/100.0_dp
   end function safe_value

   subroutine sort_simplex(simplex,fv)
      real(dp), intent(inout) :: simplex(:,:),fv(:)
      real(dp), allocatable :: col(:)
      real(dp) :: tv
      integer :: i,j,imin
      allocate(col(size(simplex,1)))
      do i=1,size(fv)-1
         imin=i-1+minloc(fv(i:),dim=1)
         if (imin/=i) then
            tv=fv(i);fv(i)=fv(imin);fv(imin)=tv
            col=simplex(:,i);simplex(:,i)=simplex(:,imin);simplex(:,imin)=col
         end if
      end do
   end subroutine sort_simplex

end module mts_optimize
