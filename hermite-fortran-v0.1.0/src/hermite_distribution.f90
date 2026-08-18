module hermite_distribution
   use hermite_kinds, only : dp, i64
   use hermite_numerics, only : normal_pdf, normal_cdf, normal_quantile, &
      logaddexp, rand_poisson
   implicit none
   private

   public :: int_hermite, hermite_logpmf_exact, dhermite_exact
   public :: cofi, edg, dhermite, phermite, qhermite, rhermite
   public :: hermite_mean, hermite_variance

contains

   pure logical function valid_parameters(a,b,m) result(ok)
      real(dp), intent(in) :: a,b
      integer, intent(in) :: m
      ok = a >= 0.0_dp .and. b >= 0.0_dp .and. m >= 2
   end function valid_parameters

   pure elemental real(dp) function hermite_mean(a,b,m) result(mu)
      real(dp), intent(in) :: a,b
      integer, intent(in) :: m
      mu=a+real(m,dp)*b
   end function hermite_mean

   pure elemental real(dp) function hermite_variance(a,b,m) result(v)
      real(dp), intent(in) :: a,b
      integer, intent(in) :: m
      v=a+real(m*m,dp)*b
   end function hermite_variance

   real(dp) function hermite_logpmf_exact(x,a,b,m) result(lp)
      integer(i64), intent(in) :: x
      real(dp), intent(in) :: a,b
      integer, intent(in) :: m
      integer(i64) :: j,k,jmax
      real(dp) :: term,loga,logb

      if (.not. valid_parameters(a,b,m) .or. x < 0_i64) then
         lp=-huge(1.0_dp)
         return
      end if
      if (a <= 0.0_dp .and. b <= 0.0_dp) then
         if (x == 0_i64) then
            lp=0.0_dp
         else
            lp=-huge(1.0_dp)
         end if
         return
      end if

      loga=0.0_dp
      logb=0.0_dp
      if (a > 0.0_dp) loga=log(a)
      if (b > 0.0_dp) logb=log(b)
      lp=-huge(1.0_dp)
      jmax=x/int(m,i64)
      do j=0_i64,jmax
         k=x-int(m,i64)*j
         if (a <= 0.0_dp .and. k > 0_i64) cycle
         if (b <= 0.0_dp .and. j > 0_i64) cycle
         term=-(a+b)-log_gamma(real(k+1_i64,dp))-log_gamma(real(j+1_i64,dp))
         if (k > 0_i64) term=term+real(k,dp)*loga
         if (j > 0_i64) term=term+real(j,dp)*logb
         lp=logaddexp(lp,term)
      end do
   end function hermite_logpmf_exact

   real(dp) function dhermite_exact(x,a,b,m,log_p) result(p)
      real(dp), intent(in) :: x,a,b
      integer, intent(in), optional :: m
      logical, intent(in), optional :: log_p
      integer :: mm
      logical :: lpflag
      real(dp) :: lp

      mm=2
      lpflag=.false.
      if (present(m)) mm=m
      if (present(log_p)) lpflag=log_p
      if (x < 0.0_dp .or. abs(x-anint(x)) > 10.0_dp*epsilon(1.0_dp)) then
         if (lpflag) then
            p=-huge(1.0_dp)
         else
            p=0.0_dp
         end if
         return
      end if
      lp=hermite_logpmf_exact(int(nint(x),i64),a,b,mm)
      if (lpflag) then
         p=lp
      else if (lp <= -0.5_dp*huge(1.0_dp)) then
         p=0.0_dp
      else
         p=exp(lp)
      end if
   end function dhermite_exact

   real(dp) function int_hermite(x,a,b,m) result(p)
      integer(i64), intent(in) :: x
      real(dp), intent(in) :: a,b
      integer, intent(in) :: m
      real(dp), allocatable :: probs(:)
      real(dp) :: mu,d,p0
      integer :: k,nmax

      if (.not. valid_parameters(a,b,m) .or. x < 0_i64) then
         p=0.0_dp
         return
      end if
      mu=a+real(m,dp)*b
      if (mu <= 0.0_dp) then
         if (x==0_i64) then
            p=1.0_dp
         else
            p=0.0_dp
         end if
         return
      end if
      d=(a+real(m*m,dp)*b)/mu
      p0=exp(mu*(-1.0_dp+(d-1.0_dp)/real(m,dp)))
      if (x==0_i64) then
         p=p0
         return
      end if
      nmax=max(int(x),m)
      allocate(probs(nmax))
      probs=0.0_dp
      do k=1,min(m-1,int(x))
         probs(k)=p0*mu**k/exp(log_gamma(real(k+1,dp))) * &
                  ((real(m,dp)-d)/real(m-1,dp))**k
      end do
      if (x >= int(m,i64)) then
         do k=m,int(x)
            if (k==m) then
               probs(k)=mu*(p0*(d-1.0_dp)+probs(k-1)*(real(m,dp)-d)) / &
                        real(k*(m-1),dp)
            else
               probs(k)=mu*(probs(k-m)*(d-1.0_dp)+probs(k-1)*(real(m,dp)-d)) / &
                        real(k*(m-1),dp)
            end if
         end do
      end if
      p=max(0.0_dp,probs(int(x)))
   end function int_hermite

   pure real(dp) function edg(y,a,b,m) result(p)
      real(dp), intent(in) :: y,a,b
      integer, intent(in), optional :: m
      integer :: mm
      real(dp) :: v,r3,r4,x,he2,he3,he5

      mm=2
      if (present(m)) mm=m
      v=a+b*real(mm*mm,dp)
      if (v <= 0.0_dp) then
         if (y>=0.0_dp) then
            p=1.0_dp
         else
            p=0.0_dp
         end if
         return
      end if
      r3=(a+b*real(mm**3,dp))/v**1.5_dp
      r4=(a+b*real(mm**4,dp))/(v*v)
      x=(1.0_dp+1.0_dp/(24.0_dp*v))*(y+0.5_dp-a-b*real(mm,dp))/sqrt(v)
      he2=x*x-1.0_dp
      he3=x**3-3.0_dp*x
      he5=x**5-10.0_dp*x**3+15.0_dp*x
      p=normal_cdf(x)-normal_pdf(x)*(r3*he2/6.0_dp+r4*he3/24.0_dp+ &
        r3*r3*he5/72.0_dp)
   end function edg

   pure real(dp) function cofi(p,a,b,m) result(q)
      real(dp), intent(in) :: p,a,b
      integer, intent(in), optional :: m
      integer :: mm
      real(dp) :: v,r3,r4,u,y

      mm=2
      if (present(m)) mm=m
      if (p <= 0.0_dp) then
         q=-huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         q=huge(1.0_dp)
         return
      end if
      v=a+b*real(mm*mm,dp)
      if (v <= 0.0_dp) then
         q=0.0_dp
         return
      end if
      r3=(a+b*real(mm**3,dp))/v**1.5_dp
      r4=(a+b*real(mm**4,dp))/(v*v)
      u=normal_quantile(p)
      y=u+(u*u-1.0_dp)*r3/6.0_dp+(u**3-3.0_dp*u)*r4/24.0_dp - &
        (2.0_dp*u**3-5.0_dp*u)*r3*r3/36.0_dp
      q=y*sqrt(v)+a+b*real(mm,dp)
   end function cofi

   real(dp) function dhermite(x,a,b,m,log_p,exact) result(p)
      real(dp), intent(in) :: x,a,b
      integer, intent(in), optional :: m
      logical, intent(in), optional :: log_p,exact
      integer :: mm
      logical :: lpflag,ex
      real(dp) :: val

      mm=2
      lpflag=.false.
      ex=.false.
      if (present(m)) mm=m
      if (present(log_p)) lpflag=log_p
      if (present(exact)) ex=exact
      if (.not. valid_parameters(a,b,mm)) then
         if (lpflag) then
            p=-huge(1.0_dp)
         else
            p=0.0_dp
         end if
         return
      end if
      if (x < 0.0_dp .or. abs(x-anint(x)) > 10.0_dp*epsilon(1.0_dp)) then
         if (lpflag) then
            p=-huge(1.0_dp)
         else
            p=0.0_dp
         end if
         return
      end if

      if (ex .or. (a <= 20.0_dp .and. b <= 20.0_dp)) then
         val=dhermite_exact(x,a,b,mm,.false.)
      else
         val=edg(x,a,b,mm)-edg(x-1.0_dp,a,b,mm)
         if (val <= 0.0_dp) val=dhermite_exact(x,a,b,mm,.false.)
      end if
      if (lpflag) then
         if (val > 0.0_dp) then
            p=log(val)
         else
            p=-huge(1.0_dp)
         end if
      else
         p=val
      end if
   end function dhermite

   real(dp) function phermite(q,a,b,m,lower_tail,exact) result(p)
      real(dp), intent(in) :: q,a,b
      integer, intent(in), optional :: m
      logical, intent(in), optional :: lower_tail,exact
      integer :: mm,k
      logical :: lt,ex
      real(dp) :: s

      mm=2
      lt=.true.
      ex=.false.
      if (present(m)) mm=m
      if (present(lower_tail)) lt=lower_tail
      if (present(exact)) ex=exact
      if (.not. valid_parameters(a,b,mm)) then
         p=0.0_dp
         return
      end if
      if (q < 0.0_dp) then
         if (lt) then
            p=0.0_dp
         else
            p=1.0_dp
         end if
         return
      end if

      if (.not. ex .and. (a > 20.0_dp .or. b > 20.0_dp)) then
         s=min(1.0_dp,max(0.0_dp,edg(real(floor(q),dp),a,b,mm)))
      else
         s=0.0_dp
         do k=0,floor(q)
            s=s+dhermite_exact(real(k,dp),a,b,mm,.false.)
         end do
         s=min(1.0_dp,max(0.0_dp,s))
      end if
      if (lt) then
         p=s
      else
         p=1.0_dp-s
      end if
   end function phermite

   integer(i64) function qhermite(p,a,b,m,lower_tail,exact) result(q)
      real(dp), intent(in) :: p,a,b
      integer, intent(in), optional :: m
      logical, intent(in), optional :: lower_tail,exact
      integer :: mm,lo,hi,mid
      logical :: lt,ex
      real(dp) :: target,approx,mu,var

      mm=2
      lt=.true.
      ex=.false.
      if (present(m)) mm=m
      if (present(lower_tail)) lt=lower_tail
      if (present(exact)) ex=exact
      if (.not. valid_parameters(a,b,mm) .or. p < 0.0_dp .or. p > 1.0_dp) then
         q=-huge(1_i64)
         return
      end if
      if (lt) then
         target=p
      else
         target=1.0_dp-p
      end if
      if (target <= 0.0_dp) then
         q=0_i64
         return
      else if (target >= 1.0_dp) then
         q=huge(1_i64)
         return
      end if

      if (.not. ex .and. (a > 20.0_dp .or. b > 20.0_dp)) then
         approx=cofi(target,a,b,mm)
         q=max(0_i64,int(floor(approx),i64)+1_i64)
         return
      end if

      mu=hermite_mean(a,b,mm)
      var=hermite_variance(a,b,mm)
      lo=0
      hi=max(1,int(ceiling(mu+8.0_dp*sqrt(max(var,0.0_dp))+10.0_dp)))
      do while (phermite(real(hi,dp),a,b,mm,.true.,.true.) < target)
         if (hi > ishft(huge(1),-1)) then
            q=huge(1_i64)
         return
         end if
         hi=2*hi
      end do
      do while (lo < hi)
         mid=lo+(hi-lo)/2
         if (phermite(real(mid,dp),a,b,mm,.true.,.true.) >= target) then
            hi=mid
         else
            lo=mid+1
         end if
      end do
      q=int(lo,i64)
   end function qhermite

   subroutine rhermite(x,a,b,m)
      integer(i64), intent(out) :: x(:)
      real(dp), intent(in) :: a,b
      integer, intent(in), optional :: m
      integer :: mm,i
      integer(i64) :: p1,p2

      mm=2
      if (present(m)) mm=m
      if (.not. valid_parameters(a,b,mm)) error stop "rhermite: invalid parameters"
      do i=1,size(x)
         p1=rand_poisson(a)
         p2=rand_poisson(b)
         x(i)=p1+int(mm,i64)*p2
      end do
   end subroutine rhermite

end module hermite_distribution
