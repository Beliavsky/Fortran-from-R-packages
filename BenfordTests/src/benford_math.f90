module benford_math
   use benford_kinds, only: dp
   implicit none
   private
   public :: normal_cdf, normal_quantile, chi_square_sf, symmetric_eigen

contains

   pure elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: q, r
      real(dp), parameter :: a1=-3.969683028665376d1, a2=2.209460984245205d2
      real(dp), parameter :: a3=-2.759285104469687d2, a4=1.383577518672690d2
      real(dp), parameter :: a5=-3.066479806614716d1, a6=2.506628277459239d0
      real(dp), parameter :: b1=-5.447609879822406d1, b2=1.615858368580409d2
      real(dp), parameter :: b3=-1.556989798598866d2, b4=6.680131188771972d1
      real(dp), parameter :: b5=-1.328068155288572d1
      real(dp), parameter :: c1=-7.784894002430293d-3, c2=-3.223964580411365d-1
      real(dp), parameter :: c3=-2.400758277161838d0, c4=-2.549732539343734d0
      real(dp), parameter :: c5=4.374664141464968d0, c6=2.938163982698783d0
      real(dp), parameter :: d1=7.784695709041462d-3, d2=3.224671290700398d-1
      real(dp), parameter :: d3=2.445134137142996d0, d4=3.754408661907416d0
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
             (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
              ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      end if
   end function normal_quantile

   pure real(dp) function gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax=1000
      real(dp), parameter :: eps=3.0d-14, fpmin=1.0d-300
      integer :: i
      real(dp) :: ap, del, sumv, b, c, d, h, an
      if (x < 0.0_dp .or. a <= 0.0_dp) then
         q = 0.0_dp
         return
      end if
      if (x == 0.0_dp) then
         q = 1.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do i=1,itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         q = 1.0_dp - sumv*exp(-x + a*log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/max(abs(b),fpmin)
         if (b < 0.0_dp) d = -d
         h = d
         do i=1,itmax
            an = -real(i,dp)*(real(i,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         q = exp(-x + a*log(x) - log_gamma(a))*h
      end if
      q = max(0.0_dp,min(1.0_dp,q))
   end function gamma_q

   pure real(dp) function chi_square_sf(x, df) result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: df
      if (x <= 0.0_dp) then
         p = 1.0_dp
      else
         p = gamma_q(0.5_dp*real(df,dp),0.5_dp*x)
      end if
   end function chi_square_sf

   subroutine symmetric_eigen(a, values, vectors, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: b(:,:)
      integer :: n, i, j, p, q, iter, maxiter, k
      real(dp) :: maxoff, app, aqq, apq, phi, c, s, bip, biq, vip, viq, tmp
      n = size(a,1)
      allocate(b(n,n),vectors(n,n),values(n))
      b = a
      vectors = 0.0_dp
      do i=1,n
         vectors(i,i)=1.0_dp
      end do
      maxiter = max(100,50*n*n)
      do iter=1,maxiter
         maxoff=0.0_dp; p=1; q=min(2,n)
         do j=2,n
            do i=1,j-1
               if (abs(b(i,j)) > maxoff) then
                  maxoff=abs(b(i,j)); p=i; q=j
               end if
            end do
         end do
         if (maxoff <= 1.0d-13*max(1.0_dp,maxval(abs(b)))) exit
         app=b(p,p); aqq=b(q,q); apq=b(p,q)
         phi=0.5_dp*atan2(2.0_dp*apq,aqq-app)
         c=cos(phi); s=sin(phi)
         do k=1,n
            if (k /= p .and. k /= q) then
               bip=b(k,p); biq=b(k,q)
               b(k,p)=c*bip-s*biq; b(p,k)=b(k,p)
               b(k,q)=s*bip+c*biq; b(q,k)=b(k,q)
            end if
         end do
         b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
         b(p,q)=0.0_dp; b(q,p)=0.0_dp
         do k=1,n
            vip=vectors(k,p); viq=vectors(k,q)
            vectors(k,p)=c*vip-s*viq
            vectors(k,q)=s*vip+c*viq
         end do
      end do
      do i=1,n
         values(i)=b(i,i)
      end do
      ! Match R's eigen(..., symmetric=TRUE): descending eigenvalues.
      do i=1,n-1
         p=i
         do j=i+1,n
            if (values(j) > values(p)) p=j
         end do
         if (p /= i) then
            tmp=values(i); values(i)=values(p); values(p)=tmp
            do k=1,n
               tmp=vectors(k,i); vectors(k,i)=vectors(k,p); vectors(k,p)=tmp
            end do
         end if
      end do
      if (present(info)) then
         if (iter > maxiter) then
            info=1
         else
            info=0
         end if
      end if
   end subroutine symmetric_eigen
end module benford_math
