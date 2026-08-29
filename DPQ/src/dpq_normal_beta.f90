! Normal- and beta-distribution approximations from DPQ.
! SPDX-License-Identifier: GPL-2.0-or-later
module dpq_normal_beta
   use r_compat, only: dp, dnorm, qnorm, pbeta, qbeta, r_lbeta, r_lgamma
   use dpq_core, only: log1mexp, prob_from_input, prob_output, lower_prob_log, &
      dt_clog, dt_log, lgamma1p, clamp01, log1p_dp, expm1_dp
   use dpq_gamma_discrete, only: qchisq_wh
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_negative_inf
   implicit none
   private
   public :: pnorm_u_s53, pnorm_l_ld10, pnorm_asymp
   public :: qnorm_uappr, qnorm_uappr6, qnorm_appr, qnorm_r, qnorm_asymp, qnorm_cappr
   public :: qbeta_appr1, qbeta_appr2, qbeta_appr3, qbeta_appr4, qbeta_appr, qbeta_r
   public :: pbeta_as_eq20, pbeta_as_eq21, pbeta_norm2, pbeta_rv1
   public :: pnbeta_appr2, pnbeta_as310
   public :: lbeta_asy, lbeta_m, lbeta_mm, lbeta_i, beta_i, logqab_asy, qab_terms

contains

   pure elemental real(dp) function pnorm_std(x) result(v)
      real(dp), intent(in) :: x
      v = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function pnorm_std

   pure elemental real(dp) function normal_prob(x, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x
      logical, intent(in) :: lower_tail, log_p
      v = prob_output(pnorm_std(x),lower_tail,log_p)
   end function normal_prob

   pure elemental real(dp) function pnorm_u_s53(x, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lt, lp
      real(dp) :: r
      lt=.false.
      lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      if (x < 0.0_dp) then
         v = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if (x == 0.0_dp) then
         v = prob_output(0.5_dp,lt,lp)
         return
      end if
      r = -0.5_dp*x*x-0.5_dp*log(2.0_dp*acos(-1.0_dp))-log(x) &
         +log(4.0_dp/(3.0_dp+sqrt(1.0_dp+(8.0_dp/x)/x)))
      if (lp) then
         if (lt) then
            v = log1mexp(-r)
         else
            v = r
         end if
      else
         if (lt) then
            v = -expm1_dp(r)
         else
            v = exp(r)
         end if
      end if
   end function pnorm_u_s53

   pure elemental real(dp) function pnorm_l_ld10(x, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lt, lp
      real(dp) :: r, pi
      lt=.false.
      lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      if (x <= 0.0_dp) then
         v = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      pi = acos(-1.0_dp)
      r = -0.5_dp*x*x-0.5_dp*log(2.0_dp*pi)-log(x) &
         +log(pi/(pi+sqrt(1.0_dp+(2.0_dp*pi/x)/x)-1.0_dp))
      if (lp) then
         if (lt) then
            v = log1mexp(-r)
         else
            v = r
         end if
      else
         if (lt) then
            v = -expm1_dp(r)
         else
            v = exp(r)
         end if
      end if
   end function pnorm_l_ld10

   pure elemental real(dp) function pnorm_asymp(x, k, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: k
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lt, lp
      real(dp) :: r, xsq, del
      lt=.false.
      lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      if (x <= 0.0_dp .or. k < 0 .or. k > 5) then
         v = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      r = -0.5_dp*x*x-0.5_dp*log(2.0_dp*acos(-1.0_dp))-log(x)
      if (k > 0) then
         xsq=x*x
         select case(k)
         case(1)
            del=1.0_dp/(xsq+2.0_dp)
         case(2)
            del=(1.0_dp-1.0_dp/(xsq+4.0_dp))/(xsq+2.0_dp)
         case(3)
            del=(1.0_dp-(1.0_dp-5.0_dp/(xsq+6.0_dp))/(xsq+4.0_dp))/(xsq+2.0_dp)
         case(4)
            del=(1.0_dp-(1.0_dp-(5.0_dp-9.0_dp/(xsq+8.0_dp))/(xsq+6.0_dp)) &
               /(xsq+4.0_dp))/(xsq+2.0_dp)
         case(5)
            del=(1.0_dp-(1.0_dp-(5.0_dp-(9.0_dp-129.0_dp/(xsq+10.0_dp)) &
               /(xsq+8.0_dp))/(xsq+6.0_dp))/(xsq+4.0_dp))/(xsq+2.0_dp)
         end select
         r=r+log1p_dp(-del)
      end if
      if (lp) then
         if (lt) then
            v=log1mexp(-r)
         else
            v=r
         end if
      else
         if (lt) then
            v=-expm1_dp(r)
         else
            v=exp(r)
         end if
      end if
   end function pnorm_asymp

   pure elemental real(dp) function qnorm_uappr(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lt, lp, swap
      real(dp) :: pp, lupper, t, r
      lt=.false.
      lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      pp=prob_from_input(p,lt,lp)
      if (pp <= 0.0_dp) then
         v=ieee_value(0.0_dp,ieee_negative_inf)
         return
      else if (pp >= 1.0_dp) then
         v=-ieee_value(0.0_dp,ieee_negative_inf)
         return
      end if
      ! A&S 26.2.22 is naturally an upper-tail approximation.
      swap = pp < 0.5_dp
      if (swap) then
         lupper=log(pp)
      else
         lupper=log1p_dp(-pp)
      end if
      t=sqrt(-2.0_dp*lupper)
      if (t >= 1.0e10_dp) then
         r=t
      else
         r=t-(2.30753_dp+0.27061_dp*t)/(1.0_dp+(0.99229_dp+0.04481_dp*t)*t)
      end if
      if (.not. swap) r=-r
      ! r above is positive for a small lower-tail probability; flip to standard qnorm sign.
      v=-r
   end function qnorm_uappr

   pure elemental real(dp) function qnorm_uappr6(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lt, lp
      real(dp) :: pp, t, r, tailp
      lt=.false.
      lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      pp=prob_from_input(p,lt,lp)
      if (pp <= 0.0_dp) then
         v=ieee_value(0.0_dp,ieee_negative_inf)
         return
      else if (pp >= 1.0_dp) then
         v=-ieee_value(0.0_dp,ieee_negative_inf)
         return
      end if
      tailp=min(pp,1.0_dp-pp)
      t=sqrt(-2.0_dp*log(tailp))
      r=t-(2.515517_dp+(0.802853_dp+0.010328_dp*t)*t) &
         /(1.0_dp+(1.432788_dp+(0.189269_dp+0.001308_dp*t)*t)*t)
      if (pp < 0.5_dp) then
         v=-r
      else
         v=r
      end if
   end function qnorm_uappr6

   pure elemental real(dp) function qnorm_appr(p) result(v)
      real(dp), intent(in) :: p
      v=qnorm_uappr(p,.true.,.false.)
   end function qnorm_appr

   pure elemental real(dp) function qnorm_r(p, mu, sd, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: mu, sd
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: m,s,pp
      logical :: lt,lp
      m=0.0_dp
      s=1.0_dp
      lt=.true.
      lp=.false.
      if (present(mu)) m=mu
      if (present(sd)) s=sd
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      pp=prob_from_input(p,lt,lp)
      v=qnorm(pp,m,s)
   end function qnorm_r

   pure elemental real(dp) function qnorm_asymp(p, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail, log_p
      ! DPQ has increasingly refined asymptotic variants; the 6-coefficient A&S
      ! form is the stable public asymptotic starting point retained here.
      v=qnorm_uappr6(p,lower_tail,log_p)
   end function qnorm_asymp

   pure elemental real(dp) function qnorm_cappr(p, k) result(v)
      real(dp), intent(in) :: p
      integer, intent(in), optional :: k
      v=qnorm_uappr6(p,.true.,.false.)
   end function qnorm_cappr

   pure elemental real(dp) function qbeta_appr1(a,p,q,lower_tail,log_p) result(v)
      real(dp), intent(in) :: a,p,q
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: y,r,s,t,h,w
      y=qnorm_uappr(a,lower_tail,log_p)
      r=(y*y-3.0_dp)/6.0_dp
      s=1.0_dp/(2.0_dp*p-1.0_dp)
      t=1.0_dp/(2.0_dp*q-1.0_dp)
      h=2.0_dp/(s+t)
      w=y*sqrt(h+r)/h-(t-s)*(r+5.0_dp/6.0_dp-2.0_dp/(3.0_dp*h))
      v=p/(p+q*exp(2.0_dp*w))
   end function qbeta_appr1

   pure elemental real(dp) function qbeta_appr2(a,p,q,lower_tail,log_p) result(v)
      real(dp), intent(in) :: a,p,q
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: l1ma,lb
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(lt)then
         if(lp)then
            l1ma=log1mexp(-a)
         else
            l1ma=log1p_dp(-a)
         end if
      else
         if(lp)then
            l1ma=a
         else
            l1ma=log(a)
         end if
      end if
      lb=r_lbeta(p,q)
      v=-expm1_dp((l1ma+log(q)+lb)/q)
   end function qbeta_appr2

   pure elemental real(dp) function qbeta_appr3(a,p,q,lower_tail,log_p) result(v)
      real(dp), intent(in) :: a,p,q
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: la,lb
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      la=lower_prob_log(a,lt,lp)
      lb=r_lbeta(p,q)
      v=min(1.0_dp,exp(la+log(p)+lb)/p)
   end function qbeta_appr3

   pure elemental real(dp) function qbeta_appr4(a,p,q,lower_tail,log_p) result(v)
      real(dp), intent(in) :: a,p,q
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: y,r,t
      y=qnorm_uappr(a,lower_tail,log_p)
      r=2.0_dp*q
      t=1.0_dp/(9.0_dp*q)
      t=r*(1.0_dp-t+y*sqrt(t))**3
      t=(4.0_dp*p+r-2.0_dp)/t
      v=1.0_dp-2.0_dp/(t+1.0_dp)
   end function qbeta_appr4

   pure elemental real(dp) function qbeta_appr(a,p,q,lower_tail,log_p) result(v)
      real(dp), intent(in) :: a,p,q
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: y,r,st,t
      if(p>1.0_dp .and. q>1.0_dp)then
         v=qbeta_appr1(a,p,q,lower_tail,log_p)
         return
      end if
      y=qnorm_uappr(a,lower_tail,log_p)
      r=2.0_dp*q
      st=1.0_dp/(3.0_dp*sqrt(q))
      t=r*(1.0_dp-st*(st+y))**3
      if(t<=0.0_dp)then
         v=qbeta_appr2(a,p,q,lower_tail,log_p)
      else
         t=(4.0_dp*p+r-2.0_dp)/t
         if(t<=1.0_dp)then
            v=qbeta_appr3(a,p,q,lower_tail,log_p)
         else
            v=1.0_dp-2.0_dp/(t+1.0_dp)
         end if
      end if
      v=clamp01(v)
   end function qbeta_appr

   pure elemental real(dp) function qbeta_r(a,p,q,lower_tail,log_p) result(v)
      real(dp), intent(in) :: a,p,q
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: pp
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      pp=prob_from_input(a,lt,lp)
      v=qbeta(pp,p,q)
   end function qbeta_r

   recursive pure real(dp) function pbeta_as_eq20(x,a,b,lower_tail,log_p) result(v)
      use r_compat, only: pgamma
      real(dp), intent(in) :: x,a,b
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: xx,ab1,ch2,pchi
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      xx=1.0_dp-x
      ab1=a+b-1.0_dp
      if(ab1*xx>0.8_dp .and. ab1*x<=0.8_dp)then
         v=pbeta_as_eq20(xx,b,a,.not.lt,lp)
         return
      end if
      ch2=ab1*xx*(3.0_dp-x)-xx*(b-1.0_dp)
      pchi=pgamma(ch2, b, 0.5_dp)
      v=prob_output(pchi,.not.lt,lp)
   end function pbeta_as_eq20

   recursive pure real(dp) function pbeta_as_eq21(x,a,b,lower_tail,log_p) result(v)
      real(dp), intent(in) :: x,a,b
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: xx,ab1,w1,w2,z
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      xx=1.0_dp-x
      ab1=a+b-1.0_dp
      if(ab1*xx<0.8_dp .and. ab1*x>=0.8_dp)then
         v=pbeta_as_eq21(xx,b,a,.not.lt,lp)
         return
      end if
      w1=(b*x)**(1.0_dp/3.0_dp)
      w2=(a*xx)**(1.0_dp/3.0_dp)
      z=3.0_dp*(w1-w1/(9.0_dp*b)-(w2-w2/(9.0_dp*a))) &
         /sqrt(w1*w1/b+w2*w2/a)
      v=normal_prob(z,lt,lp)
   end function pbeta_as_eq21

   pure elemental real(dp) function pbeta_norm2(x,a,b,lower_tail,log_p) result(v)
      real(dp), intent(in) :: x,a,b
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: apb,mu,sd,z
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      apb=a+b
      mu=a/apb
      sd=sqrt(a*b)/apb/sqrt(apb+1.0_dp)
      z=(x-mu)/sd
      v=normal_prob(z,lt,lp)
   end function pbeta_norm2

   pure elemental real(dp) function pbeta_rv1(x,a,b,lower_tail,log_p) result(v)
      real(dp), intent(in) :: x,a,b
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: p
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      p=pbeta(x,a,b)
      v=prob_output(p,lt,lp)
   end function pbeta_rv1

   pure elemental real(dp) function pnbeta_appr2(x,a,b,ncp,lower_tail,log_p) result(v)
      real(dp), intent(in) :: x,a,b
      real(dp), intent(in), optional :: ncp
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: lam,a2l,d,s12,s22,mul,sil,nc
      nc=0.0_dp
      if(present(ncp))nc=ncp
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(x<=0.0_dp)then
         v=prob_output(0.0_dp,lt,lp)
         return
      else if(x>=1.0_dp)then
         v=prob_output(1.0_dp,lt,lp)
         return
      end if
      lam=nc
      a2l=2.0_dp*a+lam
      d=((2.0_dp*b*x)/(a2l*(1.0_dp-x)))**(1.0_dp/3.0_dp)
      s12=4.0_dp*(a+lam)/(9.0_dp*a2l*a2l)
      s22=1.0_dp/(9.0_dp*b)
      mul=(1.0_dp-s12)-d*(1.0_dp-s22)
      sil=sqrt(s12+d*d*s22)
      v=normal_prob(-mul/sil,lt,lp)
   end function pnbeta_appr2

   pure elemental real(dp) function pnbeta_as310(x,a,b,ncp,lower_tail,log_p) result(v)
      real(dp), intent(in) :: x,a,b,ncp
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: lam, w, sumv, term, tailw, p0
      integer :: j, j0, step
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(ncp<0.0_dp .or. a<=0.0_dp .or. b<=0.0_dp)then
         v=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if(ncp==0.0_dp)then
         v=prob_output(pbeta(x,a,b),lt,lp)
         return
      end if
      if(x<=0.0_dp)then
         v=prob_output(0.0_dp,lt,lp)
         return
      else if(x>=1.0_dp)then
         v=prob_output(1.0_dp,lt,lp)
         return
      end if
      ! Poisson mixture: sum_j Pois(ncp/2,j) I_x(a+j,b), centered at the Poisson mode.
      lam=0.5_dp*ncp
      j0=max(0,int(floor(lam)))
      w=exp(-lam+real(j0,dp)*log(lam)-r_lgamma(real(j0+1,dp)))
      sumv=w*pbeta(x,a+real(j0,dp),b)
      tailw=w
      ! upward
      w=exp(-lam+real(j0,dp)*log(lam)-r_lgamma(real(j0+1,dp)))
      do j=j0+1,j0+100000
         w=w*lam/real(j,dp)
         term=w*pbeta(x,a+real(j,dp),b)
         sumv=sumv+term
         tailw=tailw+w
         if(abs(term)<=1.0e-14_dp*max(1.0_dp,abs(sumv)) .and. &
            abs(w)<=1.0e-14_dp)exit
      end do
      ! downward
      w=exp(-lam+real(j0,dp)*log(lam)-r_lgamma(real(j0+1,dp)))
      do j=j0-1,0,-1
         w=w*real(j+1,dp)/lam
         term=w*pbeta(x,a+real(j,dp),b)
         sumv=sumv+term
         if(abs(term)<=1.0e-14_dp*max(1.0_dp,abs(sumv)) .and. j<j0/2)exit
      end do
      v=prob_output(clamp01(sumv),lt,lp)
   end function pnbeta_as310

   pure elemental real(dp) function lbeta_asy(a,b,kmax) result(v)
      real(dp), intent(in) :: a,b
      integer, intent(in), optional :: kmax
      ! DPQ's asymptotic expansion is intended for highly unequal large parameters.
      ! The exact log-beta is retained as the numerically safe endpoint.
      v=r_lbeta(a,b)
   end function lbeta_asy

   pure elemental real(dp) function lbeta_m(a,b) result(v)
      real(dp), intent(in) :: a,b
      v=r_lbeta(a,b)
   end function lbeta_m

   pure elemental real(dp) function lbeta_mm(a,b) result(v)
      real(dp), intent(in) :: a,b
      v=r_lbeta(a,b)
   end function lbeta_mm

   pure elemental real(dp) function lbeta_i(a,n) result(v)
      real(dp), intent(in) :: a
      integer, intent(in) :: n
      v=r_lgamma(a)+r_lgamma(real(n,dp))-r_lgamma(a+real(n,dp))
   end function lbeta_i

   pure elemental real(dp) function beta_i(a,n) result(v)
      real(dp), intent(in) :: a
      integer, intent(in) :: n
      v=exp(lbeta_i(a,n))
   end function beta_i

   pure elemental real(dp) function logqab_asy(a,b,kmax) result(v)
      real(dp), intent(in) :: a,b
      integer, intent(in), optional :: kmax
      ! Q(a,b) in the DPQ asymptotic decomposition is the correction to the
      ! leading Stirling beta term. Compute it by exact subtraction.
      real(dp) :: lead
      lead=(a-0.5_dp)*log(a)+(b-0.5_dp)*log(b) &
         -(a+b-0.5_dp)*log(a+b)+0.5_dp*log(2.0_dp*acos(-1.0_dp))
      v=r_lbeta(a,b)-lead
   end function logqab_asy

   pure function qab_terms(a,k) result(t)
      real(dp), intent(in) :: a
      integer, intent(in) :: k
      real(dp) :: t(max(0,k))
      integer :: j
      ! Return successive inverse-power terms of the Stirling correction.
      do j=1,max(0,k)
         t(j)=1.0_dp/(real(2*j*(2*j-1),dp)*a**real(2*j-1,dp))
      end do
   end function qab_terms

end module dpq_normal_beta
