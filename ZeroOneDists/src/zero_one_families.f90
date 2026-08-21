! SPDX-License-Identifier: MIT
module zero_one_families
   use zero_one_kinds, only : dp
   use zero_one_distributions
   use zero_one_special, only : quiet_nan
   implicit none
   private
   integer, parameter, public :: zod_ber=1, zod_ber2=2, zod_uhlg=3, zod_umb=4, zod_uphn=5
   integer, parameter, public :: link_identity=1, link_log=2, link_logit=3, link_inverse=4
   public :: family_npar, family_default_links, linkfun, linkinv, mu_eta
   public :: family_logpdf, family_cdf, family_quantile, family_score, family_working_hessian
   public :: family_deviance_increment, family_initial, family_mean
   public :: family_valid_parameters, family_valid_y
contains
   pure integer function family_npar(family) result(n)
      integer, intent(in) :: family
      select case (family)
      case (zod_ber,zod_ber2)
         n = 3
      case (zod_uphn)
         n = 2
      case (zod_uhlg,zod_umb)
         n = 1
      case default
         n = 0
      end select
   end function family_npar

   pure subroutine family_default_links(family, links)
      integer, intent(in) :: family
      integer, intent(out) :: links(3)
      links = link_identity
      select case (family)
      case (zod_ber,zod_ber2)
         links = [link_logit,link_log,link_logit]
      case (zod_uhlg,zod_umb)
         links(1) = link_log
      case (zod_uphn)
         links(1:2) = link_log
      end select
   end subroutine family_default_links

   elemental real(dp) function linkfun(mu, link) result(v)
      real(dp), intent(in) :: mu
      integer, intent(in) :: link
      select case (link)
      case (link_identity)
         v = mu
      case (link_log)
         v = log(mu)
      case (link_logit)
         v = log(mu/(1.0_dp-mu))
      case (link_inverse)
         v = 1.0_dp/mu
      case default
         v = quiet_nan()
      end select
   end function linkfun

   elemental real(dp) function linkinv(eta, link) result(v)
      real(dp), intent(in) :: eta
      integer, intent(in) :: link
      real(dp) :: e
      select case (link)
      case (link_identity)
         v = eta
      case (link_log)
         v = max(tiny(1.0_dp),exp(min(700.0_dp,eta)))
      case (link_logit)
         if (eta >= 0.0_dp) then
            v = 1.0_dp/(1.0_dp+exp(-min(eta,700.0_dp)))
         else
            e = exp(max(eta,-700.0_dp))
            v = e/(1.0_dp+e)
         end if
      case (link_inverse)
         v = 1.0_dp/eta
      case default
         v = quiet_nan()
      end select
   end function linkinv

   elemental real(dp) function mu_eta(eta, link) result(v)
      real(dp), intent(in) :: eta
      integer, intent(in) :: link
      real(dp) :: m
      select case (link)
      case (link_identity)
         v = 1.0_dp
      case (link_log)
         v = max(tiny(1.0_dp),exp(min(700.0_dp,eta)))
      case (link_logit)
         m = linkinv(eta,link_logit)
         v = max(tiny(1.0_dp),m*(1.0_dp-m))
      case (link_inverse)
         v = -1.0_dp/(eta*eta)
      case default
         v = quiet_nan()
      end select
   end function mu_eta

   real(dp) function family_logpdf(family, y, par) result(v)
      integer, intent(in) :: family
      real(dp), intent(in) :: y, par(3)
      select case (family)
      case (zod_ber)
         v = dber(y,par(1),par(2),par(3),.true.)
      case (zod_ber2)
         v = dber2(y,par(1),par(2),par(3),.true.)
      case (zod_uhlg)
         v = duhlg(y,par(1),.true.)
      case (zod_umb)
         v = dumb(y,par(1),.true.)
      case (zod_uphn)
         v = duphn(y,par(1),par(2),.true.)
      case default
         v = quiet_nan()
      end select
   end function family_logpdf

   real(dp) function family_cdf(family, y, par) result(v)
      integer, intent(in) :: family
      real(dp), intent(in) :: y, par(3)
      select case (family)
      case (zod_ber)
         v = pber(y,par(1),par(2),par(3))
      case (zod_ber2)
         v = pber2(y,par(1),par(2),par(3))
      case (zod_uhlg)
         v = puhlg(y,par(1))
      case (zod_umb)
         v = pumb(y,par(1))
      case (zod_uphn)
         v = puphn(y,par(1),par(2))
      case default
         v = quiet_nan()
      end select
   end function family_cdf

   real(dp) function family_quantile(family, p, par) result(v)
      integer, intent(in) :: family
      real(dp), intent(in) :: p, par(3)
      select case (family)
      case (zod_ber)
         v = qber(p,par(1),par(2),par(3))
      case (zod_ber2)
         v = qber2(p,par(1),par(2),par(3))
      case (zod_uhlg)
         v = quhlg(p,par(1))
      case (zod_umb)
         v = qumb(p,par(1))
      case (zod_uphn)
         v = quphn(p,par(1),par(2))
      case default
         v = quiet_nan()
      end select
   end function family_quantile

   subroutine family_score(family, y, par, score)
      integer, intent(in) :: family
      real(dp), intent(in) :: y, par(3)
      real(dp), intent(out) :: score(3)
      real(dp) :: h, pp(3), pm(3), lplus, lminus, l
      integer :: j, np
      score = 0.0_dp
      select case (family)
      case (zod_uhlg)
         score(1) = 1.0_dp/par(1)-2.0_dp*(1.0_dp-y)/(par(1)+(2.0_dp-par(1))*y)
         return
      case (zod_umb)
         l = log(1.0_dp/y)
         score(1) = -3.0_dp/par(1)+l*l/par(1)**3
         return
      case (zod_ber,zod_ber2)
         h = 1.0e-5_dp
      case (zod_uphn)
         h = 1.0e-2_dp
      case default
         score = quiet_nan()
         return
      end select
      np = family_npar(family)
      do j = 1, np
         pp = par
         pm = par
         pp(j) = par(j)+h
         pm(j) = par(j)-h
         if (.not.family_valid_parameters(family,pm)) then
            pm = par
            l = family_logpdf(family,y,par)
            lplus = family_logpdf(family,y,pp)
            score(j) = (lplus-l)/h
         else
            lplus = family_logpdf(family,y,pp)
            lminus = family_logpdf(family,y,pm)
            score(j) = (lplus-lminus)/(2.0_dp*h)
         end if
      end do
   end subroutine family_score

   subroutine family_working_hessian(family, y, par, hess)
      integer, intent(in) :: family
      real(dp), intent(in) :: y, par(3)
      real(dp), intent(out) :: hess(3,3)
      real(dp) :: score(3), l
      integer :: i, j, np
      hess = 0.0_dp
      if (family == zod_umb) then
         l = log(1.0_dp/y)
         hess(1,1) = 3.0_dp/par(1)**2-3.0_dp*l*l/par(1)**4
         return
      end if
      call family_score(family,y,par,score)
      np = family_npar(family)
      do i = 1, np
         do j = 1, np
            hess(i,j) = min(-score(i)*score(j),-1.0e-15_dp)
         end do
      end do
   end subroutine family_working_hessian

   real(dp) function family_deviance_increment(family, y, par) result(v)
      integer, intent(in) :: family
      real(dp), intent(in) :: y, par(3)
      v = -2.0_dp*family_logpdf(family,y,par)
   end function family_deviance_increment

   subroutine family_initial(family, y, par)
      integer, intent(in) :: family
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: par(3)
      par = 1.0_dp
      select case (family)
      case (zod_ber,zod_ber2)
         par = [0.5_dp,1.5_dp,0.5_dp]
      case (zod_uhlg)
         par(1) = 0.5_dp
      case (zod_umb)
         par(1) = sqrt(sum(log(y)**2)/(3.0_dp*real(size(y),dp)))
      case (zod_uphn)
         par(1:2) = 1.0_dp
      end select
   end subroutine family_initial

   real(dp) function family_mean(family, par) result(v)
      integer, intent(in) :: family
      real(dp), intent(in) :: par(3)
      select case (family)
      case (zod_ber)
         v = 0.5_dp*par(3)+(1.0_dp-par(3))*par(1)
      case (zod_ber2)
         v = par(1)
      case default
         v = quiet_nan()
      end select
   end function family_mean

   pure logical function family_valid_parameters(family, par) result(ok)
      integer, intent(in) :: family
      real(dp), intent(in) :: par(3)
      select case (family)
      case (zod_ber)
         ok = par(1)>0.0_dp .and. par(1)<1.0_dp .and. par(2)>0.0_dp .and. par(3)>0.0_dp .and. par(3)<1.0_dp
      case (zod_ber2)
         ok = par(1)>=0.0_dp .and. par(1)<=1.0_dp .and. par(2)>0.0_dp .and. par(3)>=0.0_dp .and. par(3)<=1.0_dp
      case (zod_uhlg,zod_umb)
         ok = par(1)>0.0_dp
      case (zod_uphn)
         ok = par(1)>0.0_dp .and. par(2)>0.0_dp
      case default
         ok = .false.
      end select
   end function family_valid_parameters

   pure logical function family_valid_y(family, y) result(ok)
      integer, intent(in) :: family
      real(dp), intent(in) :: y(:)
      select case (family)
      case (zod_ber,zod_ber2,zod_uhlg,zod_umb,zod_uphn)
         ok = all(y>0.0_dp .and. y<1.0_dp)
      case default
         ok = .false.
      end select
   end function family_valid_y
end module zero_one_families
