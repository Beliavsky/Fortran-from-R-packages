! SPDX-License-Identifier: GPL-3.0-only
module smoots_smoothing
   use smoots_kinds, only : dp
   use smoots_status, only : sm_ok, sm_invalid_input
   use smoots_linalg, only : solve_linear_system
   use smoots_stats, only : factorial_real, mean_value
   implicit none
   private
   public :: gsmooth, knsmooth, local_polynomial_smooth, kernel_smooth
   public :: lag_window_variance, rescale_derivative
   public :: trend_kernel_constants, derivative_kernel_constants
   public :: inflation_bandwidth, derivative_inflation_bandwidth
   public :: variance_bandwidth, trend_exponents, derivative_exponents
contains
   subroutine gsmooth(y, v, p, mu, bandwidth, boundary_adaptive, estimate, weights, status)
      real(dp), intent(in) :: y(:), bandwidth
      integer, intent(in) :: v, p, mu, boundary_adaptive
      real(dp), allocatable, intent(out) :: estimate(:), weights(:,:)
      integer, intent(out) :: status
      call local_polynomial_smooth(y,v,p,mu,bandwidth,boundary_adaptive,estimate,weights,status)
   end subroutine gsmooth

   subroutine knsmooth(y, mu, bandwidth, boundary_adaptive, estimate, status)
      real(dp), intent(in) :: y(:), bandwidth
      integer, intent(in) :: mu, boundary_adaptive
      real(dp), allocatable, intent(out) :: estimate(:)
      integer, intent(out) :: status
      call kernel_smooth(y,mu,bandwidth,boundary_adaptive,estimate,status)
   end subroutine knsmooth

   subroutine local_polynomial_smooth(y, v, p, mu, bandwidth, boundary_adaptive, estimate, weights, status)
      real(dp), intent(in) :: y(:), bandwidth
      integer, intent(in) :: v, p, mu, boundary_adaptive
      real(dp), allocatable, intent(out) :: estimate(:), weights(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: normal(:,:), rhs(:), solution(:), w(:), xpos(:)
      real(dp) :: u, kernel, scale, power
      integer :: n, hh, htm, i0, hr, nwin, j, r, s, row, source, center

      status = sm_ok
      n = size(y)
      allocate(estimate(n))
      estimate = 0.0_dp
      if (n < 3 .or. v < 0 .or. p < 1 .or. v >= p .or. mod(p-v,2) == 0 .or. &
          mu < 0 .or. bandwidth <= 0.0_dp .or. bandwidth >= 0.5_dp .or. &
          (boundary_adaptive /= 0 .and. boundary_adaptive /= 1)) then
         status = sm_invalid_input
         allocate(weights(0,0))
         return
      end if
      hh = int(real(n,dp)*bandwidth + 0.5_dp)
      hh = min(hh,(n-1)/2)
      if (hh < max(1,p/2+1)) then
         hh = max(1,p/2+1)
         if (2*hh+1 > n) then
            status = sm_invalid_input
            allocate(weights(0,0))
            return
         end if
      end if
      htm = 2*hh+1
      allocate(weights(htm,htm)); weights = 0.0_dp
      allocate(normal(p+1,p+1),rhs(p+1),solution(p+1),w(htm),xpos(htm))

      do i0 = 0, hh
         hr = hh + boundary_adaptive*(hh-i0)
         nwin = i0 + hr + 1
         w = 0.0_dp; xpos = 0.0_dp
         do j = 1, nwin
            xpos(j) = real(j-1-i0,dp)/real(hh,dp)
            u = real(j-1-i0,dp)/real(hr+1,dp)
            kernel = max(0.0_dp,1.0_dp-u*u)
            if (mu == 0) then
               w(j) = 1.0_dp
            else
               w(j) = kernel**mu
            end if
         end do
         normal = 0.0_dp
         do r = 0, p
            do s = 0, p
               normal(r+1,s+1) = sum(w(1:nwin)*xpos(1:nwin)**(r+s))
            end do
         end do
         rhs = 0.0_dp; rhs(v+1) = 1.0_dp
         call solve_linear_system(normal,rhs,solution,status)
         if (status /= sm_ok) return
         do j = 1, nwin
            power = 0.0_dp
            do r = 0, p
               power = power + solution(r+1)*xpos(j)**r
            end do
            weights(i0+1,j) = w(j)*power
         end do
      end do

      do i0 = 1, hh
         row = hh+1+i0
         source = hh+1-i0
         weights(row,:) = merge(-1.0_dp,1.0_dp,mod(v,2)==1)*weights(source,htm:1:-1)
      end do
      scale = factorial_real(v)*(real(n,dp)/real(hh,dp))**v
      weights = scale*weights

      if (hh > 0) then
         do row = 1, hh
            estimate(row) = dot_product(weights(row,:),y(1:htm))
         end do
      end if
      center = hh+1
      do row = center, n-hh
         estimate(row) = dot_product(weights(center,:),y(row-hh:row+hh))
      end do
      if (hh > 0) then
         do i0 = 1, hh
            row = n-hh+i0
            estimate(row) = dot_product(weights(hh+1+i0,:),y(n-htm+1:n))
         end do
      end if
   end subroutine local_polynomial_smooth

   subroutine kernel_smooth(y, mu, bandwidth, boundary_adaptive, estimate, status)
      real(dp), intent(in) :: y(:), bandwidth
      integer, intent(in) :: mu, boundary_adaptive
      real(dp), allocatable, intent(out) :: estimate(:)
      integer, intent(out) :: status
      real(dp), allocatable :: w(:)
      real(dp) :: u, sw
      integer :: n, hh, i0, hr, nwin, j, row
      n = size(y)
      allocate(estimate(n)); estimate = 0.0_dp
      status = sm_ok
      if (n < 3 .or. mu < 0 .or. bandwidth <= 0.0_dp .or. bandwidth >= 0.5_dp .or. &
          (boundary_adaptive /= 0 .and. boundary_adaptive /= 1)) then
         status = sm_invalid_input
         return
      end if
      hh = int(real(n,dp)*bandwidth+0.5_dp)
      hh = min(hh,(n-1)/2)
      if (hh < 1) hh = 1
      allocate(w(2*hh+1))
      do i0 = 0, hh
         hr = hh + boundary_adaptive*(hh-i0)
         nwin = i0+hr+1
         w = 0.0_dp
         do j = 1, nwin
            u = real(j-1-i0,dp)/(real(hr,dp)+0.5_dp)
            w(j) = 0.5_dp*max(0.0_dp,1.0_dp-u*u)**mu
         end do
         sw = sum(w(1:nwin))
         if (sw > 0.0_dp) w(1:nwin)=w(1:nwin)/sw
         row = i0+1
         if (row <= hh) then
            estimate(row)=dot_product(w(1:nwin),y(1:nwin))
            estimate(n-row+1)=dot_product(w(1:nwin),y(n:n-nwin+1:-1))
         end if
         if (row == hh+1) then
            do j = row, n-hh
               estimate(j)=dot_product(w(1:nwin),y(j-hh:j+hh))
            end do
         end if
      end do
   end subroutine kernel_smooth

   subroutine lag_window_variance(x, cf0, l0_opt, lg_opt, status)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: cf0
      integer, intent(out) :: l0_opt, lg_opt, status
      real(dp), allocatable :: gamma(:), lseq(:)
      real(dp) :: mu, c1, c2, c10, c20, xj, w, sum0, sum1
      integer :: n, lag, i, l1, next_l, nit, max_l
      n = size(x)
      cf0 = 0.0_dp; l0_opt=0; lg_opt=0; status=sm_ok
      if (n < 4) then
         status=sm_invalid_input; return
      end if
      allocate(gamma(n)); mu=mean_value(x)
      do lag=0,n-1
         gamma(lag+1)=dot_product(x(1:n-lag)-mu,x(1+lag:n)-mu)/real(n,dp)
      end do
      c1=(gamma(1)**2+2.0_dp*sum(gamma(2:n)**2))/(4.0_dp*acos(-1.0_dp))
      if (c1 <= tiny(1.0_dp)) then
         cf0=gamma(1); return
      end if
      nit=20; allocate(lseq(nit+1)); lseq=0.0_dp
      lseq(1)=int(real(n,dp)/2.0_dp+0.5_dp)
      lg_opt=int(lseq(1))
      do i=1,nit
         l1=int(lseq(i)/real(n,dp)**(2.0_dp/21.0_dp))+1
         l1=min(max(1,l1),n)
         c2=0.0_dp
         do lag=0,l1-1
            xj=real(lag,dp)/real(l1,dp)
            c2=c2+(real(lag,dp)*gamma(lag+1)*(1.0_dp-xj))**2
         end do
         c2=3.0_dp*(2.0_dp*c2)/(2.0_dp*acos(-1.0_dp))
         next_l=int(real(n,dp)**(1.0_dp/3.0_dp)*(max(c2/c1,0.0_dp))**(1.0_dp/3.0_dp))+1
         next_l=min(max(1,next_l),n-1)
         lseq(i+1)=real(next_l,dp)
         lg_opt=next_l
         if (next_l==int(lseq(i))) exit
      end do
      l1=int(real(lg_opt,dp)/real(n,dp)**(2.0_dp/21.0_dp))+1
      l1=min(max(1,l1),n)
      sum1=0.0_dp; sum0=0.0_dp
      do lag=0,l1-1
         xj=real(lag,dp)/real(l1,dp)
         sum1=sum1+real(lag,dp)*gamma(lag+1)*(1.0_dp-xj)
         sum0=sum0+gamma(lag+1)*(1.0_dp+cos(acos(-1.0_dp)*xj))/2.0_dp
      end do
      c20=3.0_dp*(2.0_dp*sum1)**2/(2.0_dp*acos(-1.0_dp))
      c10=(2.0_dp*sum0-gamma(1))**2/(2.0_dp*acos(-1.0_dp))
      if (c10 <= tiny(1.0_dp)) then
         l0_opt=1
      else
         l0_opt=int(real(n,dp)**(1.0_dp/3.0_dp)*((c20/c10)/2.0_dp)**(1.0_dp/3.0_dp))+1
      end if
      l0_opt=min(max(1,l0_opt),n-1)
      max_l=l0_opt
      cf0=-gamma(1)
      do lag=0,max_l
         w=real(max_l+1-lag,dp)/real(max_l+1,dp)
         cf0=cf0+2.0_dp*gamma(lag+1)*w
      end do
   end subroutine lag_window_variance

   pure function rescale_derivative(value, first_x, second_x, last_x, order) result(output)
      real(dp), intent(in) :: value(:), first_x, second_x, last_x
      integer, intent(in) :: order
      real(dp) :: output(size(value)), spacing, scale
      spacing=second_x-first_x
      scale=last_x-first_x+spacing
      output=value/scale**order
   end function rescale_derivative

   subroutine trend_kernel_constants(p,mu,boundary_cut,rp,muk,status)
      integer,intent(in)::p,mu
      real(dp),intent(in)::boundary_cut
      real(dp),intent(out)::rp,muk
      integer,intent(out)::status
      integer,parameter::m=20000
      integer::i,k
      real(dp)::u,h,kv,s1,s2,coef
      status=sm_ok; rp=0.0_dp; muk=0.0_dp
      if (.not.(p==1 .or. p==3) .or. mu<0 .or. mu>3 .or. boundary_cut<0.0_dp .or. boundary_cut>=0.5_dp) then
         status=sm_invalid_input; return
      end if
      k=p+1; h=2.0_dp/real(m,dp); s1=0.0_dp; s2=0.0_dp
      do i=0,m
         u=-1.0_dp+h*real(i,dp)
         kv=trend_kernel_value(p,mu,u)
         if (i==0 .or. i==m) then; coef=1.0_dp
         else if (mod(i,2)==0) then; coef=2.0_dp
         else; coef=4.0_dp; end if
         s1=s1+coef*kv*kv
         s2=s2+coef*u**k*kv
      end do
      rp=(h/3.0_dp)*s1
      muk=(h/3.0_dp)*s2
   end subroutine trend_kernel_constants

   pure function trend_kernel_value(p,mu,u) result(kv)
      integer,intent(in)::p,mu
      real(dp),intent(in)::u
      real(dp)::kv,z
      z=u*u
      if (p==1) then
         kv=max(0.0_dp,1.0_dp-z)**mu
      else
         select case(mu)
         case(0); kv=3.0_dp/8.0_dp*(3.0_dp-5.0_dp*z)
         case(1); kv=15.0_dp/32.0_dp*(3.0_dp-10.0_dp*z+7.0_dp*z*z)
         case(2); kv=105.0_dp/64.0_dp*(1.0_dp-5.0_dp*z+7.0_dp*z*z-3.0_dp*z**3)
         case default; kv=315.0_dp/512.0_dp*(3.0_dp-20.0_dp*z+42.0_dp*z*z-36.0_dp*z**3+11.0_dp*z**4)
         end select
      end if
   end function trend_kernel_value

   subroutine derivative_kernel_constants(d,mu,rp,muk,status)
      integer,intent(in)::d,mu
      real(dp),intent(out)::rp,muk
      integer,intent(out)::status
      integer,parameter::m=20000
      integer::i,k,p
      real(dp)::u,h,kv,s1,s2,coef
      status=sm_ok; rp=0.0_dp; muk=0.0_dp
      if ((d/=1 .and. d/=2) .or. mu<0 .or. mu>3) then
         status=sm_invalid_input; return
      end if
      p=d+1; k=p+1; h=2.0_dp/real(m,dp); s1=0.0_dp; s2=0.0_dp
      do i=0,m
         u=-1.0_dp+h*real(i,dp)
         kv=derivative_kernel_value(d,mu,u)
         if (i==0 .or. i==m) then; coef=1.0_dp
         else if (mod(i,2)==0) then; coef=2.0_dp
         else; coef=4.0_dp; end if
         s1=s1+coef*kv*kv
         s2=s2+coef*u**k*kv
      end do
      rp=(h/3.0_dp)*s1; muk=(h/3.0_dp)*s2
   end subroutine derivative_kernel_constants

   pure function derivative_kernel_value(d,mu,u) result(kv)
      integer,intent(in)::d,mu
      real(dp),intent(in)::u
      real(dp)::kv,z
      z=u*u
      if (d==1) then
         select case(mu)
         case(0); kv=-1.5_dp*u
         case(1); kv=15.0_dp/4.0_dp*(u-u**3)
         case(2); kv=105.0_dp/16.0_dp*(-u+2.0_dp*u**3-u**5)
         case default; kv=315.0_dp/32.0_dp*(-u+3.0_dp*u**3-3.0_dp*u**5+u**7)
         end select
      else
         select case(mu)
         case(0); kv=15.0_dp/4.0_dp*(1.0_dp-3.0_dp*z)
         case(1); kv=105.0_dp/16.0_dp*(-1.0_dp+6.0_dp*z-5.0_dp*z*z)
         case(2); kv=315.0_dp/32.0_dp*(-1.0_dp+9.0_dp*z-15.0_dp*z*z+7.0_dp*z**3)
         case default; kv=3465.0_dp/256.0_dp*(-1.0_dp+12.0_dp*z-30.0_dp*z*z+28.0_dp*z**3-9.0_dp*z**4)
         end select
      end if
   end function derivative_kernel_value

   pure function inflation_bandwidth(p,inflation,b) result(value)
      integer,intent(in)::p,inflation
      real(dp),intent(in)::b
      real(dp)::value,exponent
      select case(inflation)
      case(1); exponent=merge(5.0_dp/7.0_dp,9.0_dp/11.0_dp,p==1)
      case(2); exponent=merge(5.0_dp/9.0_dp,9.0_dp/13.0_dp,p==1)
      case default; exponent=0.5_dp
      end select
      value=b**exponent
   end function inflation_bandwidth

   pure function derivative_inflation_bandwidth(d,inflation,b) result(value)
      integer,intent(in)::d,inflation
      real(dp),intent(in)::b
      real(dp)::value,exponent
      select case(inflation)
      case(1); exponent=merge(7.0_dp/9.0_dp,9.0_dp/11.0_dp,d==1)
      case(2); exponent=merge(7.0_dp/11.0_dp,9.0_dp/13.0_dp,d==1)
      case default; exponent=0.5_dp
      end select
      value=b**exponent
   end function derivative_inflation_bandwidth

   pure function variance_bandwidth(p,mu,b) result(value)
      integer,intent(in)::p,mu
      real(dp),intent(in)::b
      real(dp)::value,mult
      if (p==1) then
         select case(mu); case(0); mult=1.3195_dp; case(1); mult=1.4310_dp; case(2); mult=1.4541_dp; case default; mult=1.4640_dp; end select
      else
         select case(mu); case(0); mult=1.2599_dp; case(1); mult=1.2913_dp; case(2); mult=1.3006_dp; case default; mult=1.3052_dp; end select
      end if
      value=mult*b
   end function variance_bandwidth

   pure subroutine trend_exponents(p,expo1,expo2)
      integer,intent(in)::p
      real(dp),intent(out)::expo1,expo2
      if (p==1) then; expo1=1.0_dp/5.0_dp; expo2=-5.0_dp/7.0_dp
      else; expo1=1.0_dp/9.0_dp; expo2=-9.0_dp/11.0_dp; end if
   end subroutine trend_exponents

   pure subroutine derivative_exponents(d,expo1,expo2)
      integer,intent(in)::d
      real(dp),intent(out)::expo1,expo2
      if (d==1) then; expo1=1.0_dp/7.0_dp; expo2=-7.0_dp/9.0_dp
      else; expo1=1.0_dp/9.0_dp; expo2=-9.0_dp/11.0_dp; end if
   end subroutine derivative_exponents
end module smoots_smoothing
