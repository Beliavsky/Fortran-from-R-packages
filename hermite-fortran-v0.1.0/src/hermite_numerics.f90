module hermite_numerics
   use hermite_kinds, only : dp, i64
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: normal_pdf, normal_cdf, normal_quantile
   public :: logaddexp, rand_poisson, set_hermite_seed
   public :: nelder_mead, invert_matrix, numerical_hessian
   public :: chi_square1_upper

   abstract interface
      function objective_nd(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_nd
   end interface

contains

   pure elemental real(dp) function normal_pdf(x) result(p)
      real(dp), intent(in) :: x
      p = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a1=-3.969683028665376e+01_dp
      real(dp), parameter :: a2= 2.209460984245205e+02_dp
      real(dp), parameter :: a3=-2.759285104469687e+02_dp
      real(dp), parameter :: a4= 1.383577518672690e+02_dp
      real(dp), parameter :: a5=-3.066479806614716e+01_dp
      real(dp), parameter :: a6= 2.506628277459239e+00_dp
      real(dp), parameter :: b1=-5.447609879822406e+01_dp
      real(dp), parameter :: b2= 1.615858368580409e+02_dp
      real(dp), parameter :: b3=-1.556989798598866e+02_dp
      real(dp), parameter :: b4= 6.680131188771972e+01_dp
      real(dp), parameter :: b5=-1.328068155288572e+01_dp
      real(dp), parameter :: c1=-7.784894002430293e-03_dp
      real(dp), parameter :: c2=-3.223964580411365e-01_dp
      real(dp), parameter :: c3=-2.400758277161838e+00_dp
      real(dp), parameter :: c4=-2.549732539343734e+00_dp
      real(dp), parameter :: c5= 4.374664141464968e+00_dp
      real(dp), parameter :: c6= 2.938163982698783e+00_dp
      real(dp), parameter :: d1= 7.784695709041462e-03_dp
      real(dp), parameter :: d2= 3.224671290700398e-01_dp
      real(dp), parameter :: d3= 2.445134137142996e+00_dp
      real(dp), parameter :: d4= 3.754408661907416e+00_dp
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q,r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p > phigh) then
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
              ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else
         q = p-0.5_dp
         r = q*q
         x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
             (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      end if
   end function normal_quantile

   pure elemental real(dp) function logaddexp(a,b) result(y)
      real(dp), intent(in) :: a,b
      real(dp) :: m
      if (a <= -0.5_dp*huge(1.0_dp)) then
         y = b
      else if (b <= -0.5_dp*huge(1.0_dp)) then
         y = a
      else
         m = max(a,b)
         y = m+log(exp(a-m)+exp(b-m))
      end if
   end function logaddexp

   real(dp) function rand_uniform() result(u)
      call random_number(u)
      if (u <= 0.0_dp) u = tiny(1.0_dp)
      if (u >= 1.0_dp) u = 1.0_dp-spacing(1.0_dp)
   end function rand_uniform

   integer(i64) function rand_poisson(lambda) result(kout)
      real(dp), intent(in) :: lambda
      real(dp) :: l,p,u,a,b,inv_alpha,vr,us,v,x,logv
      integer(i64) :: k
      integer :: n

      if (lambda < 0.0_dp) error stop "rand_poisson: lambda must be nonnegative"
      if (lambda <= 0.0_dp) then
         kout = 0_i64
         return
      end if

      if (lambda < 30.0_dp) then
         l = exp(-lambda)
         p = 1.0_dp
         k = 0_i64
         do
            k = k+1_i64
            p = p*rand_uniform()
            if (p <= l) exit
         end do
         kout = k-1_i64
         return
      end if

      b = 0.931_dp+2.53_dp*sqrt(lambda)
      a = -0.059_dp+0.02483_dp*b
      inv_alpha = 1.1239_dp+1.1328_dp/(b-3.4_dp)
      vr = 0.9277_dp-3.6224_dp/(b-2.0_dp)
      do
         u = rand_uniform()-0.5_dp
         v = rand_uniform()
         us = 0.5_dp-abs(u)
         if (us <= 0.0_dp) cycle
         x = (2.0_dp*a/us+b)*u+lambda+0.43_dp
         n = floor(x)
         if (n < 0) cycle
         if (us >= 0.07_dp .and. v <= vr) then
            kout = int(n,i64)
            return
         end if
         if (us < 0.013_dp .and. v > us) cycle
         logv = log(v*inv_alpha/(a/(us*us)+b))
         if (logv <= -lambda+real(n,dp)*log(lambda)-log_gamma(real(n+1,dp))) then
            kout = int(n,i64)
            return
         end if
      end do
   end function rand_poisson

   subroutine set_hermite_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: state(:)
      call random_seed(size=n)
      allocate(state(n))
      do i = 1, n
         state(i) = mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=state)
   end subroutine set_hermite_seed

   subroutine order_simplex(f,idx)
      real(dp), intent(in) :: f(:)
      integer, intent(out) :: idx(size(f))
      integer :: i,j,t
      idx = [(i,i=1,size(f))]
      do i = 1, size(f)-1
         do j = i+1, size(f)
            if (f(idx(j)) < f(idx(i))) then
               t=idx(i)
               idx(i)=idx(j)
               idx(j)=t
            end if
         end do
      end do
   end subroutine order_simplex

   subroutine nelder_mead(fn,start,xbest,fbest,iterations,status,tol,max_iter,step)
      procedure(objective_nd) :: fn
      real(dp), intent(in) :: start(:)
      real(dp), allocatable, intent(out) :: xbest(:)
      real(dp), intent(out) :: fbest
      integer, intent(out) :: iterations,status
      real(dp), intent(in), optional :: tol,step
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: x(:,:),f(:),cent(:),xr(:),xe(:),xc(:)
      integer, allocatable :: idx(:)
      real(dp) :: fr,fe,fc,eps,st,spread
      integer :: n,j,best,worst,second,imax

      n = size(start)
      eps = 1.0e-8_dp
      st = 0.15_dp
      imax = max(1000,300*n)
      if (present(tol)) eps=tol
      if (present(step)) st=step
      if (present(max_iter)) imax=max_iter

      allocate(x(n,n+1),f(n+1),cent(n),xr(n),xe(n),xc(n),idx(n+1))
      x(:,1)=start
      do j=2,n+1
         x(:,j)=start
         x(j-1,j)=x(j-1,j)+st*max(1.0_dp,abs(start(j-1)))
      end do
      do j=1,n+1
         f(j)=fn(x(:,j))
      end do
      status=1
      do iterations=1,imax
         call order_simplex(f,idx)
         best=idx(1)
         second=idx(n)
         worst=idx(n+1)
         cent=0.0_dp
         do j=1,n
            cent=cent+x(:,idx(j))
         end do
         cent=cent/real(n,dp)

         xr=cent+(cent-x(:,worst))
         fr=fn(xr)
         if (fr < f(best)) then
            xe=cent+2.0_dp*(xr-cent)
            fe=fn(xe)
            if (fe < fr) then
               x(:,worst)=xe
               f(worst)=fe
            else
               x(:,worst)=xr
               f(worst)=fr
            end if
         else if (fr < f(second)) then
            x(:,worst)=xr
            f(worst)=fr
         else
            if (fr < f(worst)) then
               xc=cent+0.5_dp*(xr-cent)
            else
               xc=cent+0.5_dp*(x(:,worst)-cent)
            end if
            fc=fn(xc)
            if (fc < min(fr,f(worst))) then
               x(:,worst)=xc
               f(worst)=fc
            else
               do j=2,n+1
                  x(:,idx(j))=x(:,best)+0.5_dp*(x(:,idx(j))-x(:,best))
                  f(idx(j))=fn(x(:,idx(j)))
               end do
            end if
         end if

         spread=0.0_dp
         do j=2,n+1
            spread=max(spread,maxval(abs(x(:,idx(j))-x(:,best))))
         end do
         if (spread <= eps*(1.0_dp+maxval(abs(x(:,best))))) then
            status=0
            exit
         end if
      end do
      call order_simplex(f,idx)
      allocate(xbest(n))
      xbest=x(:,idx(1))
      fbest=f(idx(1))
   end subroutine nelder_mead

   subroutine invert_matrix(a,ainv,status)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: aug(:,:)
      real(dp) :: piv,tmp,scale
      integer :: n,i,j,k,p

      status=0
      if (size(a,1) /= size(a,2)) then
         allocate(ainv(0,0))
         status=1
         return
      end if
      n=size(a,1)
      allocate(aug(n,2*n),ainv(n,n))
      aug(:,1:n)=a
      aug(:,n+1:2*n)=0.0_dp
      do i=1,n
         aug(i,n+i)=1.0_dp
      end do
      scale=max(1.0_dp,maxval(abs(a)))

      do k=1,n
         p=k
         do i=k+1,n
            if (abs(aug(i,k)) > abs(aug(p,k))) p=i
         end do
         if (abs(aug(p,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            ainv=huge(1.0_dp)
            status=2
            return
         end if
         if (p /= k) then
            do j=1,2*n
               tmp=aug(k,j)
               aug(k,j)=aug(p,j)
               aug(p,j)=tmp
            end do
         end if
         piv=aug(k,k)
         aug(k,:)=aug(k,:)/piv
         do i=1,n
            if (i==k) cycle
            tmp=aug(i,k)
            aug(i,:)=aug(i,:)-tmp*aug(k,:)
         end do
      end do
      ainv=aug(:,n+1:2*n)
   end subroutine invert_matrix

   subroutine numerical_hessian(fn,x,hess,status,rel_step)
      procedure(objective_nd) :: fn
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: hess(:,:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: rel_step
      real(dp), allocatable :: xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:)
      real(dp) :: f0,hi,hj,rs
      integer :: n,i,j

      n=size(x)
      allocate(hess(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n))
      rs=1.0e-4_dp
      if (present(rel_step)) rs=rel_step
      f0=fn(x)
      if (f0 >= 0.1_dp*huge(1.0_dp)) then
         hess=huge(1.0_dp)
         status=1
         return
      end if

      do i=1,n
         hi=rs*max(1.0_dp,abs(x(i)))
         xp=x
         xm=x
         xp(i)=xp(i)+hi
         xm(i)=xm(i)-hi
         hess(i,i)=(fn(xp)-2.0_dp*f0+fn(xm))/(hi*hi)
         do j=i+1,n
            hj=rs*max(1.0_dp,abs(x(j)))
            xpp=x
            xpm=x
            xmp=x
            xmm=x
            xpp(i)=xpp(i)+hi
            xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi
            xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi
            xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi
            xmm(j)=xmm(j)-hj
            hess(i,j)=(fn(xpp)-fn(xpm)-fn(xmp)+fn(xmm))/(4.0_dp*hi*hj)
            hess(j,i)=hess(i,j)
         end do
      end do
      status=0
   end subroutine numerical_hessian

   pure elemental real(dp) function chi_square1_upper(x) result(p)
      real(dp), intent(in) :: x
      if (x <= 0.0_dp) then
         p=1.0_dp
      else
         p=erfc(sqrt(0.5_dp*x))
      end if
   end function chi_square1_upper

end module hermite_numerics
