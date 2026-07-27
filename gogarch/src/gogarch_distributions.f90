! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gogarch_kinds, only : dp
   use gogarch_rng, only : random_uniform, random_normal, random_gamma, random_student_t
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: distribution_is_valid, innovation_logpdf, innovation_pdf
   public :: random_innovation, innovation_asym_power_moment
   public :: symmetric_absolute_moment, fs_location_scale

contains

   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code+32)
      end do
   end function lower_string

   pure function distribution_is_valid(distribution, shape, skew) result(ok)
      character(len=*), intent(in) :: distribution
      real(dp), intent(in), optional :: shape, skew
      logical :: ok
      character(len=:), allocatable :: dist
      real(dp) :: sh, sk
      sh = 8.0_dp
      sk = 1.0_dp
      if (present(shape)) sh = shape
      if (present(skew)) sk = skew
      dist = trim(adjustl(lower_string(distribution)))
      select case (dist)
      case ('norm')
         ok = .true.
      case ('snorm')
         ok = sk > 0.0_dp
      case ('std')
         ok = sh > 2.0_dp
      case ('sstd')
         ok = sh > 2.0_dp .and. sk > 0.0_dp
      case ('ged')
         ok = sh > 0.0_dp
      case ('sged')
         ok = sh > 0.0_dp .and. sk > 0.0_dp
      case default
         ok = .false.
      end select
   end function distribution_is_valid

   pure function innovation_pdf(z, distribution, shape, skew) result(value)
      real(dp), intent(in) :: z, shape, skew
      character(len=*), intent(in) :: distribution
      real(dp) :: value, lp
      lp = innovation_logpdf(z,distribution,shape,skew)
      if (lp <= log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(lp)
      end if
   end function innovation_pdf

   pure function innovation_logpdf(z, distribution, shape, skew) result(value)
      real(dp), intent(in) :: z, shape, skew
      character(len=*), intent(in) :: distribution
      real(dp) :: value, mu, sigma, raw
      character(len=:), allocatable :: dist
      dist = trim(adjustl(lower_string(distribution)))
      if (.not. distribution_is_valid(dist,shape,skew)) then
         value = -huge(1.0_dp)
         return
      end if
      select case (dist)
      case ('norm')
         value = normal_logpdf(z)
      case ('std')
         value = student_logpdf(z,shape)
      case ('ged')
         value = ged_logpdf(z,shape)
      case ('snorm')
         call fs_location_scale('norm',shape,skew,mu,sigma)
         raw = mu+sigma*z
         value = log(sigma)+fs_raw_logpdf(raw,'norm',shape,skew)
      case ('sstd')
         call fs_location_scale('std',shape,skew,mu,sigma)
         raw = mu+sigma*z
         value = log(sigma)+fs_raw_logpdf(raw,'std',shape,skew)
      case ('sged')
         call fs_location_scale('ged',shape,skew,mu,sigma)
         raw = mu+sigma*z
         value = log(sigma)+fs_raw_logpdf(raw,'ged',shape,skew)
      case default
         value = -huge(1.0_dp)
      end select
   end function innovation_logpdf

   function random_innovation(distribution, shape, skew) result(z)
      character(len=*), intent(in) :: distribution
      real(dp), intent(in) :: shape, skew
      real(dp) :: z, raw, mu, sigma, probability_positive, magnitude
      character(len=:), allocatable :: dist
      z = 0.0_dp
      raw = 0.0_dp
      mu = 0.0_dp
      sigma = 1.0_dp
      probability_positive = 0.5_dp
      magnitude = 0.0_dp
      dist = trim(adjustl(lower_string(distribution)))
      if (.not. distribution_is_valid(dist,shape,skew)) error stop 'random_innovation: invalid distribution parameters'
      select case (dist)
      case ('norm')
         z = random_normal()
      case ('std')
         z = sqrt((shape-2.0_dp)/shape)*random_student_t(shape)
      case ('ged')
         z = random_symmetric_ged(shape)
      case ('snorm','sstd','sged')
         select case (dist)
         case ('snorm')
            magnitude = abs(random_normal())
            call fs_location_scale('norm',shape,skew,mu,sigma)
         case ('sstd')
            magnitude = abs(sqrt((shape-2.0_dp)/shape)*random_student_t(shape))
            call fs_location_scale('std',shape,skew,mu,sigma)
         case default
            magnitude = abs(random_symmetric_ged(shape))
            call fs_location_scale('ged',shape,skew,mu,sigma)
         end select
         probability_positive = skew*skew/(1.0_dp+skew*skew)
         if (random_uniform() < probability_positive) then
            raw = skew*magnitude
         else
            raw = -magnitude/skew
         end if
         z = (raw-mu)/sigma
      end select
   end function random_innovation

   pure function symmetric_absolute_moment(distribution, order, shape) result(moment)
      character(len=*), intent(in) :: distribution
      real(dp), intent(in) :: order, shape
      real(dp) :: moment, lambda
      character(len=:), allocatable :: dist
      dist = trim(adjustl(lower_string(distribution)))
      if (order <= -1.0_dp) then
         moment = huge(1.0_dp)
         return
      end if
      select case (dist)
      case ('norm')
         moment = 2.0_dp**(0.5_dp*order)*exp(log_gamma(0.5_dp*(order+1.0_dp))-0.5_dp*log(pi))
      case ('std')
         if (shape <= max(2.0_dp,order)) then
            moment = huge(1.0_dp)
         else
            moment = (shape-2.0_dp)**(0.5_dp*order)*exp(log_gamma(0.5_dp*(order+1.0_dp))+ &
               log_gamma(0.5_dp*(shape-order))-0.5_dp*log(pi)-log_gamma(0.5_dp*shape))
         end if
      case ('ged')
         lambda = sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
         moment = lambda**order*exp(log_gamma((order+1.0_dp)/shape)-log_gamma(1.0_dp/shape))
      case default
         moment = huge(1.0_dp)
      end select
   end function symmetric_absolute_moment

   pure subroutine fs_location_scale(base_distribution, shape, skew, mean_raw, sd_raw)
      character(len=*), intent(in) :: base_distribution
      real(dp), intent(in) :: shape, skew
      real(dp), intent(out) :: mean_raw, sd_raw
      real(dp) :: first_abs, second_raw
      first_abs = symmetric_absolute_moment(base_distribution,1.0_dp,shape)
      mean_raw = first_abs*(skew-1.0_dp/skew)
      second_raw = (skew**3+skew**(-3))/(skew+1.0_dp/skew)
      sd_raw = sqrt(max(second_raw-mean_raw*mean_raw,1.0e-14_dp))
   end subroutine fs_location_scale

   function innovation_asym_power_moment(distribution, delta, gamma, shape, skew) result(moment)
      character(len=*), intent(in) :: distribution
      real(dp), intent(in) :: delta, gamma, shape, skew
      real(dp) :: moment, h, z, weight, bound
      integer, parameter :: ngrid = 2000
      integer :: i
      character(len=:), allocatable :: dist
      dist = trim(adjustl(lower_string(distribution)))
      if (.not. distribution_is_valid(dist,shape,skew) .or. delta <= 0.0_dp .or. abs(gamma) >= 1.0_dp) then
         moment = huge(1.0_dp)
         return
      end if
      if ((dist == 'norm' .or. dist == 'std' .or. dist == 'ged') .and. &
          (dist /= 'std' .or. shape > delta)) then
         moment = symmetric_absolute_moment(dist,delta,shape)*0.5_dp* &
            ((1.0_dp-gamma)**delta+(1.0_dp+gamma)**delta)
         return
      end if
      bound = 16.0_dp
      if (dist == 'std' .or. dist == 'sstd') bound = 40.0_dp
      h = 2.0_dp*bound/real(ngrid,dp)
      moment = 0.0_dp
      do i = 0, ngrid
         z = -bound+h*real(i,dp)
         if (i == 0 .or. i == ngrid) then
            weight = 1.0_dp
         else if (mod(i,2) == 0) then
            weight = 2.0_dp
         else
            weight = 4.0_dp
         end if
         moment = moment+weight*(abs(z)-gamma*z)**delta*innovation_pdf(z,dist,shape,skew)
      end do
      moment = moment*h/3.0_dp
      if (.not. ieee_is_finite(moment) .or. moment <= 0.0_dp) moment = huge(1.0_dp)
   end function innovation_asym_power_moment

   pure function normal_logpdf(z) result(value)
      real(dp), intent(in) :: z
      real(dp) :: value
      value = -0.5_dp*(log(2.0_dp*pi)+z*z)
   end function normal_logpdf

   pure function student_logpdf(z, shape) result(value)
      real(dp), intent(in) :: z, shape
      real(dp) :: value
      value = log_gamma(0.5_dp*(shape+1.0_dp))-log_gamma(0.5_dp*shape)- &
         0.5_dp*log(pi*(shape-2.0_dp))-0.5_dp*(shape+1.0_dp)*log(1.0_dp+z*z/(shape-2.0_dp))
   end function student_logpdf

   pure function ged_logpdf(z, shape) result(value)
      real(dp), intent(in) :: z, shape
      real(dp) :: value, lambda
      lambda = sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
      value = log(shape)-log(2.0_dp)-log(lambda)-log_gamma(1.0_dp/shape)-(abs(z)/lambda)**shape
   end function ged_logpdf

   pure function fs_raw_logpdf(x, base_distribution, shape, skew) result(value)
      real(dp), intent(in) :: x, shape, skew
      character(len=*), intent(in) :: base_distribution
      real(dp) :: value, transformed
      if (x >= 0.0_dp) then
         transformed = x/skew
      else
         transformed = x*skew
      end if
      select case (trim(adjustl(lower_string(base_distribution))))
      case ('norm')
         value = log(2.0_dp)-log(skew+1.0_dp/skew)+normal_logpdf(transformed)
      case ('std')
         value = log(2.0_dp)-log(skew+1.0_dp/skew)+student_logpdf(transformed,shape)
      case ('ged')
         value = log(2.0_dp)-log(skew+1.0_dp/skew)+ged_logpdf(transformed,shape)
      case default
         value = -huge(1.0_dp)
      end select
   end function fs_raw_logpdf

   function random_symmetric_ged(shape) result(z)
      real(dp), intent(in) :: shape
      real(dp) :: z, lambda, magnitude
      lambda = sqrt(exp(log_gamma(1.0_dp/shape)-log_gamma(3.0_dp/shape)))
      magnitude = lambda*random_gamma(1.0_dp/shape)**(1.0_dp/shape)
      if (random_uniform() < 0.5_dp) then
         z = magnitude
      else
         z = -magnitude
      end if
   end function random_symmetric_ged

end module gogarch_distributions
