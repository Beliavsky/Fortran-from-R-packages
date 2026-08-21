module circstats_distributions
    use circstats_kinds, only: dp, pi, twopi
    use circstats_special, only: i0e
    use circstats_utils, only: wrap_2pi, wrap_pi, randu, randn, randexp, randcauchy
    implicit none
    private
    public :: dcard, dtri, dvm, dmixedvm, dwrpcauchy, dwrpnorm, pvm
    public :: rcard, rtri, rvm, rmixedvm, rwrpcauchy, rwrpnorm
    public :: rstable, rwrpstab

contains

    pure elemental real(dp) function dcard(theta,mu,r) result(f)
        real(dp), intent(in) :: theta,mu,r
        if (abs(r) > 0.5_dp) then
            f = 0.0_dp
        else
            f = (1.0_dp+2.0_dp*r*cos(theta-mu))/twopi
        end if
    end function dcard

    pure elemental real(dp) function dtri(theta,r) result(f)
        real(dp), intent(in) :: theta,r
        real(dp) :: t
        if (r < 0.0_dp .or. r > 4.0_dp/(pi*pi)) then
            f = 0.0_dp
            return
        end if
        t = wrap_2pi(theta)
        f = (4.0_dp-pi*pi*r+2.0_dp*pi*r*abs(pi-t))/(8.0_dp*pi)
    end function dtri

    pure elemental real(dp) function dvm(theta,mu,kappa) result(f)
        real(dp), intent(in) :: theta,mu,kappa
        if (kappa < 0.0_dp) then
            f = 0.0_dp
        else
            f = exp(kappa*(cos(theta-mu)-1.0_dp))/(twopi*i0e(kappa))
        end if
    end function dvm

    pure elemental real(dp) function dmixedvm(theta,mu1,mu2,kappa1,kappa2,p) result(f)
        real(dp), intent(in) :: theta,mu1,mu2,kappa1,kappa2,p
        if (p < 0.0_dp .or. p > 1.0_dp) then
            f = 0.0_dp
        else
            f = p*dvm(theta,mu1,kappa1)+(1.0_dp-p)*dvm(theta,mu2,kappa2)
        end if
    end function dmixedvm

    pure elemental real(dp) function dwrpcauchy(theta,mu,rho) result(f)
        real(dp), intent(in) :: theta,mu,rho
        real(dp) :: den
        if (rho < 0.0_dp .or. rho >= 1.0_dp) then
            f = 0.0_dp
            return
        end if
        den = 1.0_dp+rho*rho-2.0_dp*rho*cos(theta-mu)
        f = (1.0_dp-rho*rho)/(twopi*den)
    end function dwrpcauchy

    real(dp) function dwrpnorm(theta,mu,rho,sd,tol) result(f)
        real(dp), intent(in) :: theta,mu
        real(dp), intent(in), optional :: rho,sd,tol
        real(dp) :: rr, ss, var, epsv, termp, termm, old
        integer :: k
        if (present(rho)) then
            rr = rho
        else
            ss = 1.0_dp
            if (present(sd)) ss = sd
            rr = exp(-0.5_dp*ss*ss)
        end if
        if (rr < 0.0_dp .or. rr > 1.0_dp) error stop "dwrpnorm: rho must be in [0,1]"
        if (rr <= tiny(1.0_dp)) then
            f = 1.0_dp/twopi
            return
        end if
        if (rr >= 1.0_dp-epsilon(1.0_dp)) then
            f = 0.0_dp
            return
        end if
        var = -2.0_dp*log(rr)
        epsv = 1.0e-10_dp
        if (present(tol)) epsv = tol
        f = normal_term(theta,mu,var,0)
        do k = 1, 100000
            old = f
            termp = normal_term(theta,mu,var,k)
            termm = normal_term(theta,mu,var,-k)
            f = f+termp+termm
            if (abs(f-old) <= epsv) exit
        end do
    contains
        pure real(dp) function normal_term(t,m,v,kk) result(z)
            real(dp), intent(in) :: t,m,v
            integer, intent(in) :: kk
            real(dp) :: d
            d = t-m+twopi*real(kk,dp)
            z = exp(-0.5_dp*d*d/v)/sqrt(twopi*v)
        end function normal_term
    end function dwrpnorm

    real(dp) function pvm(theta,mu,kappa,acc) result(p)
        real(dp), intent(in) :: theta,mu,kappa
        real(dp), intent(in), optional :: acc
        real(dp) :: t, tol
        if (kappa < 0.0_dp) error stop "pvm: kappa must be nonnegative"
        t = wrap_2pi(theta)
        if (abs(t) <= 10.0_dp*epsilon(1.0_dp) .and. theta > 0.0_dp) t = twopi
        tol = 1.0e-12_dp
        if (present(acc)) tol = max(acc,1.0e-14_dp)
        if (t <= 0.0_dp) then
            p = 0.0_dp
        else if (t >= twopi) then
            p = 1.0_dp
        else
            p = adaptive_simpson_vm(0.0_dp,t,mu,kappa,tol,20)
            p = max(0.0_dp,min(1.0_dp,p))
        end if
    end function pvm

    recursive real(dp) function adaptive_simpson_vm(a,b,mu,kappa,tol,depth) result(v)
        real(dp), intent(in) :: a,b,mu,kappa,tol
        integer, intent(in) :: depth
        real(dp) :: c, fa, fb, fc, whole
        c = 0.5_dp*(a+b)
        fa = dvm(a,mu,kappa)
        fb = dvm(b,mu,kappa)
        fc = dvm(c,mu,kappa)
        whole = (b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
        v = adaptive_step(a,b,fa,fb,fc,whole,mu,kappa,tol,depth)
    end function adaptive_simpson_vm

    recursive real(dp) function adaptive_step(a,b,fa,fb,fc,whole,mu,kappa,tol,depth) result(v)
        real(dp), intent(in) :: a,b,fa,fb,fc,whole,mu,kappa,tol
        integer, intent(in) :: depth
        real(dp) :: c,d,e,fd,fe,left,right,delta
        c = 0.5_dp*(a+b)
        d = 0.5_dp*(a+c)
        e = 0.5_dp*(c+b)
        fd = dvm(d,mu,kappa)
        fe = dvm(e,mu,kappa)
        left = (c-a)*(fa+4.0_dp*fd+fc)/6.0_dp
        right = (b-c)*(fc+4.0_dp*fe+fb)/6.0_dp
        delta = left+right-whole
        if (depth <= 0 .or. abs(delta) <= 15.0_dp*tol) then
            v = left+right+delta/15.0_dp
        else
            v = adaptive_step(a,c,fa,fc,fd,left,mu,kappa,0.5_dp*tol,depth-1) + &
                adaptive_step(c,b,fc,fb,fe,right,mu,kappa,0.5_dp*tol,depth-1)
        end if
    end function adaptive_step

    subroutine rcard(n,mu,r,x)
        integer, intent(in) :: n
        real(dp), intent(in) :: mu,r
        real(dp), intent(out) :: x(n)
        real(dp) :: cand,y,f,fmax
        integer :: i
        if (abs(r) > 0.5_dp) error stop "rcard: abs(r) must not exceed 0.5"
        fmax = (1.0_dp+2.0_dp*abs(r))/twopi
        i = 1
        do while (i <= n)
            cand = twopi*randu()
            y = fmax*randu()
            f = dcard(cand,mu,r)
            if (y <= f) then
                x(i) = cand
                i = i+1
            end if
        end do
    end subroutine rcard

    subroutine rtri(n,r,x)
        integer, intent(in) :: n
        real(dp), intent(in) :: r
        real(dp), intent(out) :: x(n)
        real(dp) :: u,a,b,c,disc,t1,t2,t
        integer :: i
        if (r < 0.0_dp .or. r > 4.0_dp/(pi*pi)) error stop "rtri: invalid r"
        if (r <= 100.0_dp*epsilon(1.0_dp)) then
            do i = 1,n
                x(i) = wrap_pi(twopi*randu())
            end do
            return
        end if
        a = pi*r
        do i = 1,n
            u = randu()
            if (u < 0.5_dp) then
                b = -(4.0_dp+pi*pi*r)
                c = 8.0_dp*pi*u
                disc = max(0.0_dp,b*b-4.0_dp*a*c)
                t1 = (-b+sqrt(disc))/(2.0_dp*a)
                t2 = (-b-sqrt(disc))/(2.0_dp*a)
                t = min(t1,t2)
            else
                b = 4.0_dp-3.0_dp*pi*pi*r
                c = 2.0_dp*pi**3*r-8.0_dp*pi*u
                disc = max(0.0_dp,b*b-4.0_dp*a*c)
                t1 = (-b+sqrt(disc))/(2.0_dp*a)
                t2 = (-b-sqrt(disc))/(2.0_dp*a)
                t = max(t1,t2)
            end if
            if (t > pi) t = t-twopi
            x(i) = t
        end do
    end subroutine rtri

    subroutine rvm(n,mean,kappa,x)
        integer, intent(in) :: n
        real(dp), intent(in) :: mean,kappa
        real(dp), intent(out) :: x(n)
        real(dp) :: a,b,r,z,f,c,u1,u2,u3,sgn
        integer :: obs
        if (kappa < 0.0_dp) error stop "rvm: kappa must be nonnegative"
        if (kappa <= 100.0_dp*epsilon(1.0_dp)) then
            do obs = 1,n
                x(obs) = twopi*randu()
            end do
            return
        end if
        a = 1.0_dp+sqrt(1.0_dp+4.0_dp*kappa*kappa)
        b = (a-sqrt(2.0_dp*a))/(2.0_dp*kappa)
        r = (1.0_dp+b*b)/(2.0_dp*b)
        obs = 1
        do while (obs <= n)
            u1 = randu()
            z = cos(pi*u1)
            f = (1.0_dp+r*z)/(r+z)
            c = kappa*(r-f)
            u2 = randu()
            if (c*(2.0_dp-c)-u2 > 0.0_dp .or. log(c/u2)+1.0_dp-c >= 0.0_dp) then
                u3 = randu()
                if (u3 >= 0.5_dp) then
                    sgn = 1.0_dp
                else
                    sgn = -1.0_dp
                end if
                x(obs) = wrap_2pi(sgn*acos(max(-1.0_dp,min(1.0_dp,f)))+mean)
                obs = obs+1
            end if
        end do
    end subroutine rvm

    subroutine rmixedvm(n,mu1,mu2,kappa1,kappa2,p,x)
        integer, intent(in) :: n
        real(dp), intent(in) :: mu1,mu2,kappa1,kappa2,p
        real(dp), intent(out) :: x(n)
        real(dp) :: tmp(1)
        integer :: i
        if (p < 0.0_dp .or. p > 1.0_dp) error stop "rmixedvm: p must be in [0,1]"
        do i=1,n
            if (randu() < p) then
                call rvm(1,mu1,kappa1,tmp)
            else
                call rvm(1,mu2,kappa2,tmp)
            end if
            x(i)=tmp(1)
        end do
    end subroutine rmixedvm

    subroutine rwrpcauchy(n,location,rho,x)
        integer, intent(in) :: n
        real(dp), intent(in) :: location,rho
        real(dp), intent(out) :: x(n)
        real(dp) :: scale
        integer :: i
        if (rho < 0.0_dp .or. rho > 1.0_dp) error stop "rwrpcauchy: rho must be in [0,1]"
        if (rho <= tiny(1.0_dp)) then
            do i=1,n
                x(i)=twopi*randu()
            end do
        else if (rho >= 1.0_dp-epsilon(1.0_dp)) then
            x=wrap_2pi(location)
        else
            scale=-log(rho)
            do i=1,n
                x(i)=wrap_2pi(randcauchy(location,scale))
            end do
        end if
    end subroutine rwrpcauchy

    subroutine rwrpnorm(n,mu,rho,x,sd)
        integer, intent(in) :: n
        real(dp), intent(in) :: mu
        real(dp), intent(in), optional :: rho,sd
        real(dp), intent(out) :: x(n)
        real(dp) :: rr,ss
        integer :: i
        if (present(rho)) then
            rr=rho
        else
            ss=1.0_dp
            if (present(sd)) ss=sd
            rr=exp(-0.5_dp*ss*ss)
        end if
        if (rr < 0.0_dp .or. rr > 1.0_dp) error stop "rwrpnorm: rho must be in [0,1]"
        if (rr <= tiny(1.0_dp)) then
            do i=1,n
                x(i)=twopi*randu()
            end do
        else if (rr >= 1.0_dp-epsilon(1.0_dp)) then
            x=wrap_2pi(mu)
        else
            ss=sqrt(-2.0_dp*log(rr))
            do i=1,n
                x(i)=wrap_2pi(mu+ss*randn())
            end do
        end if
    end subroutine rwrpnorm

    subroutine rstable(n,scale,index,skewness,x)
        integer, intent(in) :: n
        real(dp), intent(in) :: scale,index,skewness
        real(dp), intent(out) :: x(n)
        real(dp) :: alpha,beta,u,v,t,s,bang,ss,x0
        integer :: i
        alpha=index
        beta=skewness
        if (alpha <= 0.0_dp .or. alpha > 2.0_dp) error stop "rstable: index must be in (0,2]"
        if (abs(beta) > 1.0_dp) error stop "rstable: skewness must be in [-1,1]"
        if (scale <= 0.0_dp) error stop "rstable: scale must be positive"
        if (abs(alpha-2.0_dp) <= 20.0_dp*epsilon(1.0_dp)) then
            do i=1,n
                x(i)=sqrt(2.0_dp)*scale*randn()
            end do
            return
        end if
        if (abs(beta) <= 20.0_dp*epsilon(1.0_dp) .and. &
            abs(alpha-1.0_dp) <= 20.0_dp*epsilon(1.0_dp)) then
            do i=1,n
                x(i)=randcauchy(0.0_dp,scale)
            end do
            return
        end if
        do i=1,n
            u=pi*(randu()-0.5_dp)
            v=randexp()
            if (abs(beta) <= 20.0_dp*epsilon(1.0_dp)) then
                t=sin(alpha*u)/(cos(u)**(1.0_dp/alpha))
                s=(cos((1.0_dp-alpha)*u)/v)**((1.0_dp-alpha)/alpha)
                x(i)=scale*t*s
            else if (abs(alpha-1.0_dp) <= 20.0_dp*epsilon(1.0_dp)) then
                x0=(((pi/2.0_dp)+beta*u)*tan(u)-beta*log((pi/2.0_dp)*v*cos(u)/ &
                    ((pi/2.0_dp)+beta*u)))/(pi/2.0_dp)
                x(i)=scale*(x0+beta*log(scale)/(pi/2.0_dp))
            else
                t=beta*tan((pi/2.0_dp)*alpha)
                bang=atan(t)/alpha
                ss=(1.0_dp+t*t)**(1.0_dp/(2.0_dp*alpha))
                x0=ss*sin(alpha*(u+bang))/(cos(u)**(1.0_dp/alpha))* &
                    (cos(u-alpha*(u+bang))/v)**((1.0_dp-alpha)/alpha)
                x(i)=scale*x0
            end if
        end do
    end subroutine rstable

    subroutine rwrpstab(n,index,skewness,scale,x)
        integer, intent(in) :: n
        real(dp), intent(in) :: index,skewness,scale
        real(dp), intent(out) :: x(n)
        call rstable(n,scale,index,skewness,x)
        x=wrap_2pi(x)
    end subroutine rwrpstab
end module circstats_distributions
