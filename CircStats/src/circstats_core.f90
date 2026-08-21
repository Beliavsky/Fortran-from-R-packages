module circstats_core
    use circstats_kinds, only: dp, pi, twopi
    use circstats_types, only: circ_dispersion_result, circ_summary_result, trig_moment_result, circ_cor_result, change_point_result
    use circstats_utils, only: wrap_2pi, sort_real
    use circstats_special, only: i0e, i1e, log_i0, normal_cdf
    implicit none
    private
    public :: a1, a1inv, est_rho, circ_mean, circ_disp, circ_summary
    public :: circ_range, circ_cor, trig_moment, rad, deg, nck, est_kappa, change_pt

contains

    pure elemental real(dp) function a1(kappa) result(r)
        real(dp), intent(in) :: kappa
        if (abs(kappa) <= tiny(1.0_dp)) then
            r = 0.0_dp
        else
            r = i1e(kappa)/i0e(kappa)
        end if
    end function a1

    pure elemental real(dp) function a1inv(x) result(kappa)
        real(dp), intent(in) :: x
        real(dp) :: xx
        xx = max(0.0_dp,min(x,1.0_dp-epsilon(1.0_dp)))
        if (xx < 0.53_dp) then
            kappa = 2.0_dp*xx + xx**3 + 5.0_dp*xx**5/6.0_dp
        else if (xx < 0.85_dp) then
            kappa = -0.4_dp + 1.39_dp*xx + 0.43_dp/(1.0_dp-xx)
        else
            kappa = 1.0_dp/(xx**3 - 4.0_dp*xx**2 + 3.0_dp*xx)
        end if
    end function a1inv

    pure real(dp) function est_rho(x) result(rho)
        real(dp), intent(in) :: x(:)
        rho = hypot(sum(sin(x)),sum(cos(x)))/real(size(x),dp)
    end function est_rho

    pure real(dp) function circ_mean(x) result(mu)
        real(dp), intent(in) :: x(:)
        mu = atan2(sum(sin(x)),sum(cos(x)))
    end function circ_mean

    pure function circ_disp(x) result(res)
        real(dp), intent(in) :: x(:)
        type(circ_dispersion_result) :: res
        real(dp) :: c, s
        res%n = size(x)
        c = sum(cos(x))
        s = sum(sin(x))
        res%r = hypot(c,s)
        res%rbar = res%r/real(res%n,dp)
        res%variance = 1.0_dp-res%rbar
    end function circ_disp

    pure function circ_summary(x) result(res)
        real(dp), intent(in) :: x(:)
        type(circ_summary_result) :: res
        res%n = size(x)
        res%mean_dir = circ_mean(x)
        res%rho = est_rho(x)
    end function circ_summary

    function circ_range(x, p_value) result(rangev)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out), optional :: p_value
        real(dp) :: rangev, gap
        real(dp), allocatable :: w(:)
        real(dp) :: frac, term, pv
        integer :: n, i, k, stopk
        n = size(x)
        allocate(w(n))
        w = wrap_2pi(x)
        call sort_real(w)
        gap = w(1)-w(n)+twopi
        do i = 1, n-1
            gap = max(gap,w(i+1)-w(i))
        end do
        rangev = twopi-gap
        if (present(p_value)) then
            frac = 1.0_dp-rangev/twopi
            if (frac <= 0.0_dp) then
                p_value = 1.0_dp
                return
            end if
            stopk = floor(1.0_dp/frac)
            stopk = min(stopk,n)
            pv = 0.0_dp
            do k = 1, stopk
                if (1.0_dp-real(k,dp)*frac < 0.0_dp) cycle
                term = (-1.0_dp)**(k-1)*nck(n,k)*(1.0_dp-real(k,dp)*frac)**(n-1)
                pv = pv+term
            end do
            p_value = max(0.0_dp,min(1.0_dp,pv))
        end if
    end function circ_range

    function circ_cor(alpha,beta,do_test) result(res)
        real(dp), intent(in) :: alpha(:), beta(:)
        logical, intent(in), optional :: do_test
        type(circ_cor_result) :: res
        real(dp) :: am, bm, num, den, l20, l02, l22
        integer :: n
        logical :: test
        n = min(size(alpha),size(beta))
        am = circ_mean(alpha(1:n))
        bm = circ_mean(beta(1:n))
        num = sum(sin(alpha(1:n)-am)*sin(beta(1:n)-bm))
        den = sqrt(sum(sin(alpha(1:n)-am)**2)*sum(sin(beta(1:n)-bm)**2))
        res%r = num/den
        test = .false.
        if (present(do_test)) test = do_test
        if (test) then
            l20 = sum(sin(alpha(1:n)-am)**2)/real(n,dp)
            l02 = sum(sin(beta(1:n)-bm)**2)/real(n,dp)
            l22 = sum((sin(alpha(1:n)-am)**2)*(sin(beta(1:n)-bm)**2))/real(n,dp)
            res%statistic = sqrt(real(n,dp)*l20*l02/l22)*res%r
            res%p_value = 2.0_dp*(1.0_dp-normal_cdf(abs(res%statistic)))
        end if
    end function circ_cor

    pure function trig_moment(x,p,center) result(res)
        real(dp), intent(in) :: x(:)
        integer, intent(in), optional :: p
        logical, intent(in), optional :: center
        type(trig_moment_result) :: res
        integer :: pp
        real(dp) :: shift
        logical :: ctr
        pp = 1
        if (present(p)) pp = p
        ctr = .false.
        if (present(center)) ctr = center
        shift = 0.0_dp
        if (ctr) shift = circ_mean(x)
        res%sine = sum(sin(real(pp,dp)*(x-shift)))/real(size(x),dp)
        res%cosine = sum(cos(real(pp,dp)*(x-shift)))/real(size(x),dp)
        res%mu = atan2(res%sine,res%cosine)
        res%rho = hypot(res%sine,res%cosine)
    end function trig_moment

    pure elemental real(dp) function rad(degree) result(radian)
        real(dp), intent(in) :: degree
        radian = degree*pi/180.0_dp
    end function rad

    pure elemental real(dp) function deg(radian) result(degree)
        real(dp), intent(in) :: radian
        degree = radian*180.0_dp/pi
    end function deg

    pure real(dp) function nck(n,k) result(v)
        integer, intent(in) :: n,k
        if (k < 0 .or. k > n) then
            v = 0.0_dp
        else
            v = exp(log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp)))
        end if
    end function nck

    pure real(dp) function est_kappa(x,bias) result(kappa)
        real(dp), intent(in) :: x(:)
        logical, intent(in), optional :: bias
        real(dp) :: ml, mu
        integer :: n
        logical :: correct
        mu = circ_mean(x)
        ml = a1inv(sum(cos(x-mu))/real(size(x),dp))
        kappa = ml
        correct = .false.
        if (present(bias)) correct = bias
        if (correct) then
            n = size(x)
            if (ml < 2.0_dp .and. ml > 0.0_dp) then
                kappa = max(ml-2.0_dp/(real(n,dp)*ml),0.0_dp)
            else if (ml >= 2.0_dp) then
                kappa = real((n-1)**3,dp)*ml/real(n**3+n,dp)
            end if
        end if
    end function est_kappa

    function change_pt(x) result(res)
        real(dp), intent(in) :: x(:)
        type(change_point_result) :: res
        real(dp), allocatable :: r1(:), r2(:), v(:), rdiff(:)
        real(dp) :: rho
        integer :: n, k
        n = size(x)
        if (n < 4) error stop "change_pt: sample size must be at least 4"
        allocate(r1(n),r2(n),v(n),rdiff(n))
        rho = est_rho(x)
        v = -huge(1.0_dp)
        do k = 1, n-1
            r1(k) = est_rho(x(1:k))*real(k,dp)
            r2(k) = est_rho(x(k+1:n))*real(n-k,dp)
            if (k >= 2 .and. k <= n-2) then
                v(k) = real(k,dp)/real(n,dp)*phi_cp(r1(k)/real(k,dp)) + &
                    real(n-k,dp)/real(n,dp)*phi_cp(r2(k)/real(n-k,dp))
            end if
        end do
        r1(n) = rho*real(n,dp)
        r2(n) = 0.0_dp
        rdiff = r1+r2-rho*real(n,dp)
        res%n = n
        res%rho = rho
        res%rmax = maxval(rdiff)
        res%rave = sum(rdiff)/real(n,dp)
        res%k_r = maxloc(rdiff,dim=1)
        res%tmax = maxval(v(2:n-2))
        res%tave = sum(v(2:n-2))/real(n-3,dp)
        res%k_t = maxloc(v(2:n-2),dim=1)+1
    contains
        pure real(dp) function phi_cp(r) result(vv)
            real(dp), intent(in) :: r
            real(dp) :: arg
            arg = a1inv(r)
            vv = r*arg-log_i0(arg)
        end function phi_cp
    end function change_pt
end module circstats_core
