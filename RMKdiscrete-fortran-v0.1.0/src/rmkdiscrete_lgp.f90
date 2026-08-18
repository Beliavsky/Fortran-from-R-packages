! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
module rmkdiscrete_lgp
  use rmkdiscrete_kinds, only : dp
  use rmkdiscrete_math, only : qnan, pinf, ninf, real_equal, kahan_add, dpois, rpois, rnorm_std
  implicit none
  private

  type, public :: lgp_summary
    real(dp) :: mean = 0.0_dp
    real(dp) :: median = 0.0_dp
    real(dp) :: mode = 0.0_dp
    real(dp) :: variance = 0.0_dp
    real(dp) :: sd = 0.0_dp
    real(dp) :: third_central_moment = 0.0_dp
    real(dp) :: fourth_central_moment = 0.0_dp
    real(dp) :: pearson_skewness = 0.0_dp
    real(dp) :: skewness = 0.0_dp
    real(dp) :: kurtosis = 0.0_dp
  end type lgp_summary

  public :: lgp_from_mu_sigma2, lgp_from_theta_lambda
  public :: lgp_from_mu_theta, lgp_from_sigma2_lambda
  public :: lgp_from_sigma2_theta, lgp_from_mu_lambda
  public :: lgp_findmax, lgp_get_nc, dlgp, plgp, qlgp, rlgp, rlgp_sample, slgp

contains

  pure subroutine lgp_from_mu_sigma2(mu, sigma2, theta, lambda, status)
    real(dp), intent(in) :: mu, sigma2
    real(dp), intent(out) :: theta, lambda
    integer, intent(out), optional :: status
    if (mu < 0.0_dp .or. sigma2 <= 0.0_dp) then
      theta = qnan()
      lambda = qnan()
      if (present(status)) status = 1
      return
    end if
    theta = sqrt(mu**3/sigma2)
    if (real_equal(mu,0.0_dp)) then
      lambda = 0.0_dp
    else
      lambda = (mu-theta)/mu
    end if
    if (present(status)) status = merge(0, 1, theta >= 0.0_dp .and. abs(lambda) <= 1.0_dp)
  end subroutine lgp_from_mu_sigma2

  pure subroutine lgp_from_mu_theta(mu, theta, sigma2, lambda, status)
    real(dp), intent(in) :: mu, theta
    real(dp), intent(out) :: sigma2, lambda
    integer, intent(out), optional :: status
    if (mu <= 0.0_dp .or. theta < 0.0_dp) then
      sigma2 = qnan()
      lambda = qnan()
      if (present(status)) status = 1
      return
    end if
    lambda = (mu-theta)/mu
    sigma2 = theta/(1.0_dp-lambda)**3
    if (present(status)) status = merge(0,1,abs(lambda)<=1.0_dp)
  end subroutine lgp_from_mu_theta

  pure subroutine lgp_from_sigma2_lambda(sigma2, lambda, mu, theta, status)
    real(dp), intent(in) :: sigma2, lambda
    real(dp), intent(out) :: mu, theta
    integer, intent(out), optional :: status
    if (sigma2 < 0.0_dp .or. abs(lambda) > 1.0_dp .or. real_equal(lambda,1.0_dp)) then
      mu = qnan()
      theta = qnan()
      if (present(status)) status = 1
      return
    end if
    theta = sigma2*(1.0_dp-lambda)**3
    mu = theta/(1.0_dp-lambda)
    if (present(status)) status = 0
  end subroutine lgp_from_sigma2_lambda

  pure subroutine lgp_from_sigma2_theta(sigma2, theta, mu, lambda, status)
    real(dp), intent(in) :: sigma2, theta
    real(dp), intent(out) :: mu, lambda
    integer, intent(out), optional :: status
    if (sigma2 <= 0.0_dp .or. theta < 0.0_dp) then
      mu = qnan()
      lambda = qnan()
      if (present(status)) status = 1
      return
    end if
    mu = (sigma2*theta**2)**(1.0_dp/3.0_dp)
    if (real_equal(mu,0.0_dp)) then
      lambda = 0.0_dp
    else
      lambda = (mu-theta)/mu
    end if
    if (present(status)) status = merge(0,1,abs(lambda)<=1.0_dp)
  end subroutine lgp_from_sigma2_theta

  pure subroutine lgp_from_mu_lambda(mu, lambda, sigma2, theta, status)
    real(dp), intent(in) :: mu, lambda
    real(dp), intent(out) :: sigma2, theta
    integer, intent(out), optional :: status
    if (mu < 0.0_dp .or. abs(lambda) > 1.0_dp .or. real_equal(lambda,1.0_dp)) then
      sigma2 = qnan()
      theta = qnan()
      if (present(status)) status = 1
      return
    end if
    theta = mu*(1.0_dp-lambda)
    sigma2 = theta/(1.0_dp-lambda)**3
    if (present(status)) status = 0
  end subroutine lgp_from_mu_lambda

  pure subroutine lgp_from_theta_lambda(theta, lambda, mu, sigma2, status)
    real(dp), intent(in) :: theta, lambda
    real(dp), intent(out) :: mu, sigma2
    integer, intent(out), optional :: status
    if (theta < 0.0_dp .or. abs(lambda) > 1.0_dp .or. real_equal(lambda,1.0_dp)) then
      mu = qnan()
      sigma2 = qnan()
      if (present(status)) status = 1
      return
    end if
    mu = theta/(1.0_dp-lambda)
    sigma2 = theta/(1.0_dp-lambda)**3
    if (present(status)) status = 0
  end subroutine lgp_from_theta_lambda

  pure real(dp) function lgp_findmax(theta, lambda) result(mx)
    real(dp), intent(in) :: theta, lambda
    real(dp) :: ratio, fr
    if (theta < 0.0_dp .or. abs(lambda) > 1.0_dp) then
      mx = qnan()
      return
    end if
    if (real_equal(theta,0.0_dp)) then
      mx = 0.0_dp
      return
    end if
    if (lambda >= 0.0_dp) then
      mx = pinf()
      return
    end if
    ratio = -theta/lambda
    if (ratio <= 1.0_dp) then
      mx = 0.0_dp
    else
      fr = floor(ratio)
      mx = fr
      if (real_equal(fr,ratio)) mx = mx - 1.0_dp
    end if
  end function lgp_findmax

  pure real(dp) function dlgp_unnormalized_log(x, theta, lambda, mx) result(v)
    integer, intent(in) :: x
    real(dp), intent(in) :: theta, lambda, mx
    real(dp) :: a
    if (x < 0) then
      v = ninf()
      return
    end if
    if (real_equal(theta,0.0_dp)) then
      v = merge(0.0_dp, ninf(), x == 0)
      return
    end if
    if (real_equal(lambda,0.0_dp)) then
      v = dpois(x, theta, .true.)
      return
    end if
    if (lambda < 0.0_dp .and. real(x,dp) > mx) then
      v = ninf()
      return
    end if
    a = theta + lambda*real(x,dp)
    if (a <= 0.0_dp) then
      v = ninf()
      return
    end if
    v = log(theta) + real(x-1,dp)*log(a) - theta - lambda*real(x,dp) &
        - log_gamma(real(x+1,dp))
  end function dlgp_unnormalized_log

  real(dp) function lgp_get_nc(theta, lambda, tol) result(nc)
    real(dp), intent(in) :: theta, lambda
    real(dp), intent(in), optional :: tol
    real(dp) :: mx, s, c, term, eps
    integer :: x, imax
    if (theta < 0.0_dp .or. abs(lambda) > 1.0_dp) then
      nc = qnan()
      return
    end if
    if (lambda >= 0.0_dp) then
      nc = 1.0_dp
      return
    end if
    mx = lgp_findmax(theta, lambda)
    imax = int(mx)
    eps = 0.0_dp
    if (present(tol)) eps = max(0.0_dp, tol)
    s = 0.0_dp
    c = 0.0_dp
    do x = 0, imax
      term = exp(dlgp_unnormalized_log(x, theta, lambda, mx))
      call kahan_add(s, c, term)
      if (eps > 0.0_dp .and. imax > 200000 .and. abs(1.0_dp-s) <= eps) exit
    end do
    nc = s
  end function lgp_get_nc

  real(dp) function dlgp(x, theta, lambda, nc, give_log) result(v)
    integer, intent(in) :: x
    real(dp), intent(in) :: theta, lambda
    real(dp), intent(in), optional :: nc
    logical, intent(in), optional :: give_log
    logical :: gl
    real(dp) :: z, lv, mx
    gl = .false.
    if (present(give_log)) gl = give_log
    if (theta < 0.0_dp .or. abs(lambda) > 1.0_dp) then
      v = qnan()
      return
    end if
    if (x < 0) then
      v = merge(ninf(), 0.0_dp, gl)
      return
    end if
    if (lambda < 0.0_dp) then
      mx = lgp_findmax(theta, lambda)
    else
      mx = pinf()
    end if
    if (present(nc)) then
      z = nc
    else
      z = lgp_get_nc(theta, lambda)
    end if
    if (.not. (z > 0.0_dp)) then
      v = qnan()
      return
    end if
    lv = dlgp_unnormalized_log(x, theta, lambda, mx) - log(z)
    if (real_equal(lv,ninf())) then
      v = merge(ninf(), 0.0_dp, gl)
    else
      v = merge(lv, exp(lv), gl)
    end if
  end function dlgp

  real(dp) function plgp(q, theta, lambda, nc, lower_tail, log_p) result(v)
    real(dp), intent(in) :: q, theta, lambda
    real(dp), intent(in), optional :: nc
    logical, intent(in), optional :: lower_tail, log_p
    logical :: lower, lp
    real(dp) :: z, s, c, term, mx
    integer :: x, iq, imax
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    lp = .false.
    if (present(log_p)) lp = log_p
    if (theta < 0.0_dp .or. abs(lambda) > 1.0_dp) then
      v = qnan()
      return
    end if
    if (q < 0.0_dp) then
      v = merge(0.0_dp, 1.0_dp, lower)
      if (lp) v = log(v)
      return
    end if
    if (present(nc)) then
      z = nc
    else
      z = lgp_get_nc(theta, lambda)
    end if
    iq = int(floor(q))
    if (lambda < 0.0_dp) then
      mx = lgp_findmax(theta, lambda)
      imax = int(mx)
      if (iq >= imax) then
        v = merge(1.0_dp, 0.0_dp, lower)
        if (lp) v = merge(0.0_dp, ninf(), lower)
        return
      end if
    end if
    s = 0.0_dp
    c = 0.0_dp
    do x = 0, iq
      term = dlgp(x, theta, lambda, z)
      call kahan_add(s, c, term)
    end do
    s = min(1.0_dp, max(0.0_dp, s))
    if (.not. lower) s = max(0.0_dp, 1.0_dp-s)
    if (lp) then
      if (real_equal(s,0.0_dp)) then
        v = ninf()
      else
        v = log(s)
      end if
    else
      v = s
    end if
  end function plgp

  real(dp) function qlgp(p, theta, lambda, nc, lower_tail, log_p, status) result(q)
    real(dp), intent(in) :: p, theta, lambda
    real(dp), intent(in), optional :: nc
    logical, intent(in), optional :: lower_tail, log_p
    integer, intent(out), optional :: status
    logical :: lower, lp
    real(dp) :: pp, z, s, c, term, mx
    integer :: x, imax
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    lp = .false.
    if (present(log_p)) lp = log_p
    pp = p
    if (lp) pp = exp(pp)
    if (.not. lower) pp = 1.0_dp-pp
    if (theta < 0.0_dp .or. abs(lambda) > 1.0_dp .or. pp < 0.0_dp .or. pp > 1.0_dp) then
      q = qnan()
      if (present(status)) status = 1
      return
    end if
    if (real_equal(pp,0.0_dp)) then
      q = 0.0_dp
      if (present(status)) status = 0
      return
    end if
    if (lambda < 0.0_dp) then
      mx = lgp_findmax(theta, lambda)
      imax = int(mx)
      if (real_equal(pp,1.0_dp)) then
        q = mx
        if (present(status)) status = 0
        return
      end if
    else
      imax = 10000000
      if (real_equal(pp,1.0_dp)) then
        q = pinf()
        if (present(status)) status = 0
        return
      end if
    end if
    if (present(nc)) then
      z = nc
    else
      z = lgp_get_nc(theta, lambda)
    end if
    s = 0.0_dp
    c = 0.0_dp
    do x = 0, imax
      term = dlgp(x, theta, lambda, z)
      call kahan_add(s, c, term)
      if (s >= pp) then
        q = real(x,dp)
        if (present(status)) status = 0
        return
      end if
      if (real_equal(term,0.0_dp) .and. x > int(max(1.0_dp,theta/(max(1.0e-12_dp,1.0_dp-lambda))))) exit
    end do
    q = qnan()
    if (present(status)) status = 2
  end function qlgp

  integer function rlgp(theta, lambda) result(x)
    real(dp), intent(in) :: theta, lambda
    real(dp) :: u, y, z, k, m, sd, omega, s, pterm, cc, mx
    if (theta < 0.0_dp .or. abs(lambda) > 1.0_dp) then
      x = -huge(1)
      return
    end if
    if (real_equal(theta,0.0_dp)) then
      x = 0
      return
    end if
    if (real_equal(lambda,0.0_dp)) then
      x = rpois(theta)
      return
    end if
    if ((theta >= 10.0_dp .and. lambda < 0.0_dp) .or. &
        (theta >= 30.0_dp .and. lambda > 0.0_dp .and. lambda < 0.2_dp)) then
      m = theta/(1.0_dp-lambda)
      sd = sqrt(theta/(1.0_dp-lambda)**3)
      y = rnorm_std()
      x = max(0, int(floor(m+sd*y+0.5_dp)))
      return
    end if
    if (lambda < 0.0_dp) then
      mx = lgp_findmax(theta,lambda)
      if (mx <= 5.0_dp) then
        call random_number(u)
        x = int(qlgp(u,theta,lambda))
        return
      end if
      omega = exp(-lambda)
      x = 0
      s = exp(-theta)
      pterm = s
      call random_number(u)
      do while (u > s .and. real(x,dp) < mx)
        x = x + 1
        cc = theta - lambda + lambda*real(x,dp)
        if (cc <= 0.0_dp) then
          call random_number(u)
          x = int(qlgp(u,theta,lambda))
          return
        end if
        pterm = omega*cc*(1.0_dp+lambda/cc)**real(x-1,dp)*pterm/real(x,dp)
        if (.not. (pterm >= 0.0_dp)) then
          call random_number(u)
          x = int(qlgp(u,theta,lambda))
          return
        end if
        s = s + pterm
      end do
      return
    end if
    y = real(rpois(theta),dp)
    x = int(y)
    do while (y > 0.0_dp)
      k = lambda*y
      z = real(rpois(k),dp)
      x = x + int(z)
      y = z
    end do
  end function rlgp

  function rlgp_sample(n,theta,lambda) result(out)
    integer, intent(in) :: n
    real(dp), intent(in) :: theta,lambda
    integer, allocatable :: out(:)
    integer :: i
    allocate(out(max(0,n)))
    do i=1,n
      out(i)=rlgp(theta,lambda)
    end do
  end function rlgp_sample

  function slgp(theta, lambda, nc, do_numerically) result(out)
    real(dp), intent(in) :: theta, lambda
    real(dp), intent(in), optional :: nc
    logical, intent(in), optional :: do_numerically
    type(lgp_summary) :: out
    logical :: num
    real(dp) :: z, mx, p, c, raw1, raw2, raw3, raw4, lp, prev, term4
    integer :: x, imax
    logical :: median_found, mode_found
    num = .false.
    if (present(do_numerically)) num = do_numerically
    if (theta < 0.0_dp .or. abs(lambda) > 1.0_dp) then
      out%mean = qnan()
      out%median = qnan()
      out%mode = qnan()
      out%variance = qnan()
      out%sd = qnan()
      out%third_central_moment = qnan()
      out%fourth_central_moment = qnan()
      out%pearson_skewness = qnan()
      out%skewness = qnan()
      out%kurtosis = qnan()
      return
    end if
    if (present(nc)) then
      z = nc
    else
      z = lgp_get_nc(theta,lambda)
    end if
    if (lambda < 0.0_dp .and. num) then
      mx = lgp_findmax(theta,lambda)
      imax = int(mx)
      p=0.0_dp
      c=0.0_dp
      raw1=0.0_dp
      raw2=0.0_dp
      raw3=0.0_dp
      raw4=0.0_dp
      prev = -huge(1.0_dp)
      out%mode = mx
      median_found = .false.
      mode_found = .false.
      do x=0,imax
        lp = dlgp(x,theta,lambda,z,.true.)
        if (x > 0 .and. lp < prev .and. .not. mode_found) then
          out%mode = real(x-1,dp)
          mode_found = .true.
        end if
        call kahan_add(p,c,exp(lp))
        if (p >= 0.5_dp .and. .not. median_found) then
          out%median = real(x,dp)
          median_found = .true.
        end if
        raw1 = raw1 + exp(lp)*real(x,dp)
        raw2 = raw2 + exp(lp)*real(x,dp)**2
        raw3 = raw3 + exp(lp)*real(x,dp)**3
        term4 = exp(lp)*real(x,dp)**4
        raw4 = raw4 + term4
        prev = lp
      end do
      out%mean = raw1
      out%variance = raw2-raw1**2
      out%third_central_moment = raw3-3.0_dp*raw2*raw1+2.0_dp*raw1**3
      out%fourth_central_moment = raw4-4.0_dp*raw3*raw1+6.0_dp*raw2*raw1**2-3.0_dp*raw1**4
    else
      out%mean = theta/(1.0_dp-lambda)
      out%median = qlgp(0.5_dp,theta,lambda,z)
      out%mode = lgp_mode(theta,lambda,z)
      out%variance = theta/(1.0_dp-lambda)**3
      out%third_central_moment = theta*(1.0_dp+2.0_dp*lambda)/(1.0_dp-lambda)**5
      out%fourth_central_moment = 3.0_dp*theta**2/(1.0_dp-lambda)**6 + &
        theta*(1.0_dp+8.0_dp*lambda+6.0_dp*lambda**2)/(1.0_dp-lambda)**7
    end if
    out%sd = sqrt(max(0.0_dp,out%variance))
    if (out%sd > 0.0_dp) then
      out%pearson_skewness = (out%mean-out%mode)/out%sd
      out%skewness = out%third_central_moment/out%variance**1.5_dp
      out%kurtosis = out%fourth_central_moment/out%variance**2 - 3.0_dp
    else
      out%pearson_skewness = 0.0_dp
      out%skewness=0.0_dp
      out%kurtosis=0.0_dp
    end if
  end function slgp

  real(dp) function lgp_mode(theta,lambda,nc) result(mode)
    real(dp), intent(in) :: theta,lambda,nc
    real(dp) :: start_r, oldp, newp, mx
    integer :: x, imax
    if (theta*exp(-lambda) < 1.0_dp) then
      mode=0.0_dp
      return
    else if (real_equal(theta*exp(-lambda),1.0_dp)) then
      mode=0.5_dp
      return
    end if
    mx=lgp_findmax(theta,lambda)
    if (lambda < 0.0_dp) then
      imax=int(mx)
    else
      imax=10000000
    end if
    start_r=floor((theta-exp(lambda))/(exp(lambda)-2.0_dp*lambda))
    x=max(0,int(start_r))
    oldp=dlgp(x,theta,lambda,nc,.true.)
    do while(x < imax)
      x=x+1
      newp=dlgp(x,theta,lambda,nc,.true.)
      if (newp < oldp) then
        mode=real(x-1,dp)
        return
      end if
      oldp=newp
    end do
    mode=real(x,dp)
  end function lgp_mode

end module rmkdiscrete_lgp
