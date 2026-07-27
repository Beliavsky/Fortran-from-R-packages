! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from PeerPerformance 2.4.0, copyright 2012-2023 David Ardia and Kris Boudt.
module peerperformance_stats
  use peerperformance_kinds, only: dp
  use peerperformance_math, only: finite_value, missing_value, normal_quantile, &
                                  two_sided_normal_pvalue
  use peerperformance_linalg, only: long_run_covariance, sample_covariance, &
                                    ols_fit, outer_product
  use peerperformance_types, only: test_result
  implicit none
  private
  public :: sharpe, modified_sharpe, alpha_coefficients
  public :: alpha_testing, sharpe_testing_asymptotic, modified_sharpe_testing_asymptotic
  public :: sharpe_difference, modified_sharpe_difference
  public :: sharpe_standard_error, modified_sharpe_standard_error
  public :: covariance_for_blocks

contains

  subroutine sharpe(x, values, nobs, status)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: values(size(x,2))
    integer, intent(out), optional :: nobs(size(x,2))
    integer, intent(out), optional :: status
    real(dp), allocatable :: z(:)
    real(dp) :: mu, variance
    integer :: j, n
    logical :: mask(size(x,1))
    values = missing_value()
    if (present(status)) status = 0
    do j = 1, size(x,2)
      mask = finite_value(x(:,j))
      n = count(mask)
      if (present(nobs)) nobs(j) = n
      if (n <= 1) cycle
      z = pack(x(:,j),mask)
      mu = sum(z)/real(n,dp)
      variance = sum((z-mu)**2)/real(n-1,dp)
      if (variance <= tiny(1.0_dp)) cycle
      values(j) = mu/sqrt(variance)
    end do
  end subroutine sharpe

  subroutine modified_sharpe(x, level, values, nobs, na_negative, status)
    real(dp), intent(in) :: x(:,:), level
    real(dp), intent(out) :: values(size(x,2))
    integer, intent(out), optional :: nobs(size(x,2))
    logical, intent(in), optional :: na_negative
    integer, intent(out), optional :: status
    real(dp), allocatable :: z(:)
    real(dp) :: m1, m2, m3, m4, skew, kurt, mvar, za
    integer :: j, n
    logical :: mask(size(x,1)), reject_negative
    values = missing_value()
    if (present(status)) status = 0
    reject_negative = .true.
    if (present(na_negative)) reject_negative = na_negative
    if (level <= 0.0_dp .or. level >= 1.0_dp) then
      if (present(status)) status = 1
      return
    end if
    za = normal_quantile(1.0_dp-level)
    do j = 1, size(x,2)
      mask = finite_value(x(:,j))
      n = count(mask)
      if (present(nobs)) nobs(j) = n
      if (n <= 1) cycle
      z = pack(x(:,j),mask)
      m1 = sum(z)/real(n,dp)
      z = z-m1
      m2 = sum(z**2)/real(n,dp)
      if (m2 <= tiny(1.0_dp)) cycle
      m3 = sum(z**3)/real(n,dp)
      m4 = sum(z**4)/real(n,dp)
      skew = m3/m2**1.5_dp
      kurt = m4/(m2*m2)-3.0_dp
      mvar = -m1+sqrt(m2)*(-za-(za*za-1.0_dp)*skew/6.0_dp- &
             (za**3-3.0_dp*za)*kurt/24.0_dp+ &
             (2.0_dp*za**3-5.0_dp*za)*skew*skew/36.0_dp)
      if (reject_negative .and. mvar < 0.0_dp) cycle
      if (abs(mvar) <= tiny(1.0_dp)) cycle
      values(j) = m1/mvar
    end do
  end subroutine modified_sharpe

  subroutine alpha_coefficients(x, coefficients, nobs, factors, status)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: coefficients(:,:)
    integer, allocatable, intent(out), optional :: nobs(:)
    real(dp), intent(in), optional :: factors(:,:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: design(:,:), y(:), beta(:), resid(:), cov(:,:), se(:), ts(:), pv(:)
    logical, allocatable :: mask(:)
    logical :: ok
    integer :: ncoef, n, j, k
    if (present(status)) status = 0
    if (present(factors)) then
      if (size(factors,1) /= size(x,1)) then
        allocate(coefficients(0,0))
        if (present(status)) status = 1
        return
      end if
      ncoef = 1+size(factors,2)
    else
      ncoef = 1
    end if
    allocate(coefficients(ncoef,size(x,2)))
    coefficients = missing_value()
    if (present(nobs)) allocate(nobs(size(x,2)))
    allocate(mask(size(x,1)))
    do j = 1, size(x,2)
      mask = finite_value(x(:,j))
      if (present(factors)) then
        do k = 1, size(factors,2)
          mask = mask .and. finite_value(factors(:,k))
        end do
      end if
      n = count(mask)
      if (present(nobs)) nobs(j) = n
      if (n <= ncoef) cycle
      allocate(design(n,ncoef),y(n),beta(ncoef),resid(n),cov(ncoef,ncoef), &
               se(ncoef),ts(ncoef),pv(ncoef))
      design(:,1) = 1.0_dp
      if (present(factors)) then
        do k = 1, size(factors,2)
          design(:,k+1) = pack(factors(:,k),mask)
        end do
      end if
      y = pack(x(:,j),mask)
      call ols_fit(design,y,beta,resid,cov,se,ts,pv,.false.,ok)
      if (ok) coefficients(:,j) = beta
      deallocate(design,y,beta,resid,cov,se,ts,pv)
    end do
  end subroutine alpha_coefficients

  subroutine alpha_testing(x, y, result, factors, hac, screen_beta, min_obs)
    real(dp), intent(in) :: x(:), y(:)
    type(test_result), intent(out) :: result
    real(dp), intent(in), optional :: factors(:,:)
    logical, intent(in), optional :: hac, screen_beta
    integer, intent(in), optional :: min_obs
    logical, allocatable :: mask(:)
    real(dp), allocatable :: design(:,:), d(:), xx(:), yy(:)
    real(dp), allocatable :: beta(:), residuals(:), covariance(:,:), se(:), ts(:), pv(:)
    real(dp), allocatable :: bx(:), by(:), rx(:), ry(:), cx(:,:), cy(:,:), sx(:), sy(:), tx(:), ty(:), px(:), py(:)
    integer :: ncoef, n, k, minimum
    logical :: use_hac, all_coef, ok, okx, oky

    result%status = 0
    result%message = ''
    use_hac = .false.
    if (present(hac)) use_hac = hac
    all_coef = .false.
    if (present(screen_beta)) all_coef = screen_beta
    minimum = 10
    if (present(min_obs)) minimum = min_obs
    if (size(x) /= size(y)) then
      result%status = 1
      result%message = 'x and y must have equal length'
      return
    end if
    if (present(factors)) then
      if (size(factors,1) /= size(x)) then
        result%status = 1
        result%message = 'factors must have one row per return observation'
        return
      end if
      ncoef = 1+size(factors,2)
    else
      ncoef = 1
      all_coef = .false.
    end if
    allocate(mask(size(x)))
    mask = finite_value(x) .and. finite_value(y)
    if (present(factors)) then
      do k = 1, size(factors,2)
        mask = mask .and. finite_value(factors(:,k))
      end do
    end if
    n = count(mask)
    result%n = n
    if (n < minimum .or. n <= ncoef) then
      result%status = 2
      result%message = 'too few complete observations'
      return
    end if
    allocate(design(n,ncoef),d(n),xx(n),yy(n),beta(ncoef),residuals(n), &
             covariance(ncoef,ncoef),se(ncoef),ts(ncoef),pv(ncoef), &
             bx(ncoef),by(ncoef),rx(n),ry(n),cx(ncoef,ncoef),cy(ncoef,ncoef), &
             sx(ncoef),sy(ncoef),tx(ncoef),ty(ncoef),px(ncoef),py(ncoef))
    design(:,1) = 1.0_dp
    if (present(factors)) then
      do k = 1, size(factors,2)
        design(:,k+1) = pack(factors(:,k),mask)
      end do
    end if
    xx = pack(x,mask)
    yy = pack(y,mask)
    d = xx-yy
    call ols_fit(design,d,beta,residuals,covariance,se,ts,pv,use_hac,ok)
    call ols_fit(design,xx,bx,rx,cx,sx,tx,px,use_hac,okx)
    call ols_fit(design,yy,by,ry,cy,sy,ty,py,use_hac,oky)
    if (.not. ok .or. .not. okx .or. .not. oky) then
      result%status = 3
      result%message = 'singular or deterministic regression'
      return
    end if
    if (.not. all_coef) ncoef = 1
    allocate(result%estimate(ncoef,2),result%difference(ncoef), &
             result%standard_error(ncoef),result%tstat(ncoef),result%pvalue(ncoef))
    result%estimate(:,1) = bx(1:ncoef)
    result%estimate(:,2) = by(1:ncoef)
    result%difference = beta(1:ncoef)
    result%standard_error = se(1:ncoef)
    result%tstat = ts(1:ncoef)
    result%pvalue = pv(1:ncoef)
  end subroutine alpha_testing

  pure real(dp) function sharpe_difference(x, y, ttype) result(value)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: ttype
    real(dp) :: m1, m2, s1, s2
    integer :: n
    n = size(x)
    value = missing_value()
    if (size(y) /= n .or. n <= 1) return
    m1 = sum(x)/real(n,dp)
    m2 = sum(y)/real(n,dp)
    s1 = sqrt(sum((x-m1)**2)/real(n-1,dp))
    s2 = sqrt(sum((y-m2)**2)/real(n-1,dp))
    if (s1 <= tiny(1.0_dp) .or. s2 <= tiny(1.0_dp)) return
    if (ttype == 1) then
      value = m1/s1-m2/s2
    else
      value = m1*s2-m2*s1
    end if
  end function sharpe_difference

  subroutine sharpe_gradient_and_scores(x, y, ttype, gradient, scores, ok)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: ttype
    real(dp), intent(out) :: gradient(4), scores(size(x),4)
    logical, intent(out) :: ok
    real(dp) :: mu(2), raw2(2), central2(2), sx(2)
    real(dp) :: dm1x(4), dm1y(4), dsx(4), dsy(4)
    integer :: n
    n = size(x)
    gradient = 0.0_dp
    scores = 0.0_dp
    ok = size(y) == n .and. n > 1
    if (.not. ok) return
    mu = [sum(x),sum(y)]/real(n,dp)
    raw2 = [sum(x*x),sum(y*y)]/real(n,dp)
    central2 = raw2-mu*mu
    if (any(central2 <= tiny(1.0_dp))) then
      ok = .false.
      return
    end if
    sx = sqrt(central2)
    if (ttype == 1) then
      gradient(1) = raw2(1)/central2(1)**1.5_dp
      gradient(2) = -raw2(2)/central2(2)**1.5_dp
      gradient(3) = -0.5_dp*mu(1)/central2(1)**1.5_dp
      gradient(4) = 0.5_dp*mu(2)/central2(2)**1.5_dp
      scores(:,1) = x-mu(1)
      scores(:,2) = y-mu(2)
      scores(:,3) = x*x-raw2(1)
      scores(:,4) = y*y-raw2(2)
    else
      dm1x = [1.0_dp,0.0_dp,0.0_dp,0.0_dp]
      dm1y = [0.0_dp,0.0_dp,1.0_dp,0.0_dp]
      dsx = 0.5_dp/sx(1)*[-2.0_dp*mu(1),1.0_dp,0.0_dp,0.0_dp]
      dsy = 0.5_dp/sx(2)*[0.0_dp,0.0_dp,-2.0_dp*mu(2),1.0_dp]
      gradient = dm1x*sx(2)+dsy*mu(1)-dm1y*sx(1)-dsx*mu(2)
      scores(:,1) = x-mu(1)
      scores(:,2) = x*x-raw2(1)
      scores(:,3) = y-mu(2)
      scores(:,4) = y*y-raw2(2)
    end if
  end subroutine sharpe_gradient_and_scores

  real(dp) function sharpe_standard_error(x, y, hac, ttype, block_length) result(value)
    real(dp), intent(in) :: x(:), y(:)
    logical, intent(in) :: hac
    integer, intent(in) :: ttype
    integer, intent(in), optional :: block_length
    real(dp), allocatable :: scores(:,:), psi(:,:)
    real(dp) :: gradient(4)
    integer :: b, n
    logical :: ok
    n = size(x)
    value = missing_value()
    allocate(scores(n,4),psi(4,4))
    call sharpe_gradient_and_scores(x,y,ttype,gradient,scores,ok)
    if (.not. ok) return
    if (present(block_length)) then
      b = block_length
      psi = covariance_for_blocks(scores,b)
    else
      psi = long_run_covariance(scores,hac)
    end if
    value = sqrt(max(0.0_dp,dot_product(gradient,matmul(psi,gradient))/real(n,dp)))
  end function sharpe_standard_error

  subroutine sharpe_testing_asymptotic(x, y, result, hac, ttype, min_obs)
    real(dp), intent(in) :: x(:), y(:)
    type(test_result), intent(out) :: result
    logical, intent(in), optional :: hac
    integer, intent(in), optional :: ttype, min_obs
    logical, allocatable :: mask(:)
    real(dp), allocatable :: xx(:), yy(:), m(:,:)
    real(dp) :: d, se, vals(2)
    integer :: n, kind_test, minimum, counts(2)
    logical :: use_hac
    result%status = 0
    result%message = ''
    use_hac = .false.; if (present(hac)) use_hac = hac
    kind_test = 2; if (present(ttype)) kind_test = ttype
    minimum = 10; if (present(min_obs)) minimum = min_obs
    if (size(x) /= size(y)) then
      result%status = 1; result%message = 'x and y must have equal length'; return
    end if
    allocate(mask(size(x)))
    mask = finite_value(x) .and. finite_value(y)
    n = count(mask)
    result%n = n
    if (n < minimum) then
      result%status = 2; result%message = 'too few complete observations'; return
    end if
    xx = pack(x,mask); yy = pack(y,mask)
    d = sharpe_difference(xx,yy,kind_test)
    se = sharpe_standard_error(xx,yy,use_hac,kind_test)
    if (.not. finite_value(d) .or. .not. finite_value(se) .or. se <= tiny(1.0_dp)) then
      result%status = 3; result%message = 'degenerate Sharpe comparison'; return
    end if
    allocate(m(n,2)); m(:,1)=xx; m(:,2)=yy
    call sharpe(m,vals,counts)
    allocate(result%estimate(1,2),result%difference(1),result%standard_error(1), &
             result%tstat(1),result%pvalue(1))
    result%estimate(1,:) = vals
    result%difference(1) = vals(1)-vals(2)
    result%standard_error(1) = se
    result%tstat(1) = d/se
    result%pvalue(1) = two_sided_normal_pvalue(result%tstat(1))
  end subroutine sharpe_testing_asymptotic

  pure real(dp) function modified_sharpe_difference(x, y, level, na_negative, ttype) result(value)
    real(dp), intent(in) :: x(:), y(:), level
    logical, intent(in) :: na_negative
    integer, intent(in) :: ttype
    real(dp) :: m1(2), m2(2), m3(2), m4(2), sk(2), ku(2), mv(2), za
    integer :: n
    n = size(x)
    value = missing_value()
    if (size(y) /= n .or. n <= 1 .or. level <= 0.0_dp .or. level >= 1.0_dp) return
    m1 = [sum(x),sum(y)]/real(n,dp)
    m2 = [sum((x-m1(1))**2),sum((y-m1(2))**2)]/real(n,dp)
    if (any(m2 <= tiny(1.0_dp))) return
    m3 = [sum((x-m1(1))**3),sum((y-m1(2))**3)]/real(n,dp)
    m4 = [sum((x-m1(1))**4),sum((y-m1(2))**4)]/real(n,dp)
    sk = m3/m2**1.5_dp
    ku = m4/(m2*m2)-3.0_dp
    za = normal_quantile(1.0_dp-level)
    mv = -m1+sqrt(m2)*(-za-(za*za-1.0_dp)*sk/6.0_dp- &
         (za**3-3.0_dp*za)*ku/24.0_dp+ &
         (2.0_dp*za**3-5.0_dp*za)*sk*sk/36.0_dp)
    if (na_negative .and. any(mv < 0.0_dp)) return
    if (any(abs(mv) <= tiny(1.0_dp))) return
    if (ttype == 1) then
      value = m1(1)/mv(1)-m1(2)/mv(2)
    else
      value = m1(1)*mv(2)-m1(2)*mv(1)
    end if
  end function modified_sharpe_difference

  subroutine modified_gradient_and_scores(x, y, level, ttype, gradient, scores, mvar, ok)
    real(dp), intent(in) :: x(:), y(:), level
    integer, intent(in) :: ttype
    real(dp), intent(out) :: gradient(8), scores(size(x),8), mvar(2)
    logical, intent(out) :: ok
    real(dp) :: m1(2), m2(2), m3(2), m4(2), g2(2), g3(2), g4(2)
    real(dp) :: skew(2), kurt(2), za
    real(dp) :: dm1(8,2), dm2(8,2), dm3(8,2), dm4(8,2), ds(8,2), dk(8,2), dmv(8,2)
    real(dp) :: c1, c2
    integer :: n, j
    n = size(x)
    gradient = 0.0_dp; scores = 0.0_dp; mvar = 0.0_dp
    ok = size(y) == n .and. n > 1 .and. level > 0.0_dp .and. level < 1.0_dp
    if (.not. ok) return
    m1 = [sum(x),sum(y)]/real(n,dp)
    m2 = [sum((x-m1(1))**2),sum((y-m1(2))**2)]/real(n,dp)
    if (any(m2 <= tiny(1.0_dp))) then
      ok = .false.; return
    end if
    m3 = [sum((x-m1(1))**3),sum((y-m1(2))**3)]/real(n,dp)
    m4 = [sum((x-m1(1))**4),sum((y-m1(2))**4)]/real(n,dp)
    g2 = m2+m1*m1
    g3 = m3+3.0_dp*m1*g2-2.0_dp*m1**3
    g4 = m4+4.0_dp*m1*g3-6.0_dp*m1*m1*g2+3.0_dp*m1**4
    skew = m3/m2**1.5_dp
    kurt = m4/(m2*m2)-3.0_dp
    dm1 = 0.0_dp; dm2 = 0.0_dp; dm3 = 0.0_dp; dm4 = 0.0_dp
    dm1(1,1)=1.0_dp; dm1(5,2)=1.0_dp
    dm2(1,1)=-2.0_dp*m1(1); dm2(2,1)=1.0_dp
    dm2(5,2)=-2.0_dp*m1(2); dm2(6,2)=1.0_dp
    dm3(1,1)=-3.0_dp*g2(1)+6.0_dp*m1(1)**2; dm3(2,1)=-3.0_dp*m1(1); dm3(3,1)=1.0_dp
    dm3(5,2)=-3.0_dp*g2(2)+6.0_dp*m1(2)**2; dm3(6,2)=-3.0_dp*m1(2); dm3(7,2)=1.0_dp
    dm4(1,1)=-4.0_dp*g3(1)+12.0_dp*m1(1)*g2(1)-12.0_dp*m1(1)**3
    dm4(2,1)=6.0_dp*m1(1)**2; dm4(3,1)=-4.0_dp*m1(1); dm4(4,1)=1.0_dp
    dm4(5,2)=-4.0_dp*g3(2)+12.0_dp*m1(2)*g2(2)-12.0_dp*m1(2)**3
    dm4(6,2)=6.0_dp*m1(2)**2; dm4(7,2)=-4.0_dp*m1(2); dm4(8,2)=1.0_dp
    do j = 1, 2
      ds(:,j) = (m2(j)**1.5_dp*dm3(:,j)-1.5_dp*m3(j)*sqrt(m2(j))*dm2(:,j))/m2(j)**3
      dk(:,j) = (m2(j)**2*dm4(:,j)-2.0_dp*m4(j)*m2(j)*dm2(:,j))/m2(j)**4
    end do
    za = normal_quantile(1.0_dp-level)
    do j = 1, 2
      c1 = 0.5_dp/sqrt(m2(j)); c2 = sqrt(m2(j))
      dmv(:,j) = -dm1(:,j)-za*c1*dm2(:,j)
      dmv(:,j) = dmv(:,j)-(za*za-1.0_dp)/6.0_dp*(c1*skew(j)*dm2(:,j)+c2*ds(:,j))
      dmv(:,j) = dmv(:,j)+(2.0_dp*za**3-5.0_dp*za)/36.0_dp* &
                   (c1*skew(j)**2*dm2(:,j)+2.0_dp*c2*skew(j)*ds(:,j))
      dmv(:,j) = dmv(:,j)-(za**3-3.0_dp*za)/24.0_dp* &
                   (c1*kurt(j)*dm2(:,j)+c2*dk(:,j))
      mvar(j) = -m1(j)+c2*(-za-(za*za-1.0_dp)*skew(j)/6.0_dp- &
                  (za**3-3.0_dp*za)*kurt(j)/24.0_dp+ &
                  (2.0_dp*za**3-5.0_dp*za)*skew(j)**2/36.0_dp)
    end do
    if (any(abs(mvar) <= tiny(1.0_dp))) then
      ok = .false.; return
    end if
    if (ttype == 1) then
      gradient = (mvar(1)*dm1(:,1)-m1(1)*dmv(:,1))/mvar(1)**2- &
                 (mvar(2)*dm1(:,2)-m1(2)*dmv(:,2))/mvar(2)**2
    else
      gradient = mvar(2)*dm1(:,1)+m1(1)*dmv(:,2)- &
                 mvar(1)*dm1(:,2)-m1(2)*dmv(:,1)
    end if
    scores(:,1)=x-m1(1); scores(:,2)=x*x-g2(1)
    scores(:,3)=x**3-g3(1); scores(:,4)=x**4-g4(1)
    scores(:,5)=y-m1(2); scores(:,6)=y*y-g2(2)
    scores(:,7)=y**3-g3(2); scores(:,8)=y**4-g4(2)
  end subroutine modified_gradient_and_scores

  real(dp) function modified_sharpe_standard_error(x, y, level, hac, ttype, block_length) result(value)
    real(dp), intent(in) :: x(:), y(:), level
    logical, intent(in) :: hac
    integer, intent(in) :: ttype
    integer, intent(in), optional :: block_length
    real(dp), allocatable :: scores(:,:), psi(:,:)
    real(dp) :: gradient(8), mvar(2)
    integer :: b, n
    logical :: ok
    n = size(x)
    value = missing_value()
    allocate(scores(n,8),psi(8,8))
    call modified_gradient_and_scores(x,y,level,ttype,gradient,scores,mvar,ok)
    if (.not. ok) return
    if (present(block_length)) then
      b = block_length
      psi = covariance_for_blocks(scores,b)
    else
      psi = long_run_covariance(scores,hac)
    end if
    value = sqrt(max(0.0_dp,dot_product(gradient,matmul(psi,gradient))/real(n,dp)))
  end function modified_sharpe_standard_error

  subroutine modified_sharpe_testing_asymptotic(x, y, level, result, na_negative, hac, ttype, min_obs)
    real(dp), intent(in) :: x(:), y(:), level
    type(test_result), intent(out) :: result
    logical, intent(in), optional :: na_negative, hac
    integer, intent(in), optional :: ttype, min_obs
    logical, allocatable :: mask(:)
    real(dp), allocatable :: xx(:), yy(:), m(:,:)
    real(dp) :: d, se, vals(2)
    integer :: n, kind_test, minimum, counts(2)
    logical :: reject_negative, use_hac
    result%status = 0; result%message = ''
    reject_negative = .true.; if (present(na_negative)) reject_negative=na_negative
    use_hac = .false.; if (present(hac)) use_hac=hac
    kind_test = 2; if (present(ttype)) kind_test=ttype
    minimum = 10; if (present(min_obs)) minimum=min_obs
    if (size(x) /= size(y)) then
      result%status=1; result%message='x and y must have equal length'; return
    end if
    allocate(mask(size(x))); mask=finite_value(x) .and. finite_value(y)
    n=count(mask); result%n=n
    if (n < minimum) then
      result%status=2; result%message='too few complete observations'; return
    end if
    xx=pack(x,mask); yy=pack(y,mask)
    d=modified_sharpe_difference(xx,yy,level,reject_negative,kind_test)
    se=modified_sharpe_standard_error(xx,yy,level,use_hac,kind_test)
    if (.not. finite_value(d) .or. .not. finite_value(se) .or. se <= tiny(1.0_dp)) then
      result%status=3; result%message='degenerate modified-Sharpe comparison'; return
    end if
    allocate(m(n,2)); m(:,1)=xx; m(:,2)=yy
    call modified_sharpe(m,level,vals,counts,reject_negative)
    allocate(result%estimate(1,2),result%difference(1),result%standard_error(1), &
             result%tstat(1),result%pvalue(1))
    result%estimate(1,:)=vals
    result%difference(1)=vals(1)-vals(2)
    result%standard_error(1)=se
    result%tstat(1)=d/se
    result%pvalue(1)=two_sided_normal_pvalue(result%tstat(1))
  end subroutine modified_sharpe_testing_asymptotic

  function covariance_for_blocks(scores, block_length) result(psi)
    real(dp), intent(in) :: scores(:,:)
    integer, intent(in) :: block_length
    real(dp) :: psi(size(scores,2),size(scores,2))
    real(dp) :: zeta(size(scores,2))
    integer :: n, p, blocks, j, first, last
    n=size(scores,1); p=size(scores,2)
    psi=0.0_dp
    if (n <= 1) return
    if (block_length <= 1) then
      psi=sample_covariance(scores)
      return
    end if
    blocks=n/block_length
    if (blocks < 1) return
    do j=1,blocks
      first=(j-1)*block_length+1; last=j*block_length
      zeta=sqrt(real(block_length,dp))*sum(scores(first:last,:),dim=1)/real(block_length,dp)
      psi=psi+outer_product(zeta,zeta)
    end do
    psi=psi/real(blocks,dp)
  end function covariance_for_blocks

end module peerperformance_stats
