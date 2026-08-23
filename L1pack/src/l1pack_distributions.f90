module l1pack_distributions
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use l1pack_base, only: dp, pi, sqrt2, inv_sqrt2, gamma_rand_l1, normal_rand, &
    chol_lower, mahalanobis_one, logdet_chol
  implicit none
  private
  public :: dlaplace, plaplace, qlaplace, rlaplace
  public :: dlaplace_vec, plaplace_vec, qlaplace_vec, rlaplace_vec
  public :: dmlaplace, log_dmlaplace, rmlaplace

contains

  pure real(dp) function dlaplace(x,location,scale,log_density) result(y)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: location,scale
    logical, intent(in), optional :: log_density
    real(dp) :: mu,s,ly
    logical :: lg
    mu=0.0_dp; if(present(location))mu=location
    s=1.0_dp; if(present(scale))s=scale
    lg=.false.; if(present(log_density))lg=log_density
    if(s<=0.0_dp) then
      y=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    ly=-0.5_dp*log(2.0_dp)-log(s)-sqrt2*abs(x-mu)/s
    if(lg) then
      y=ly
    else
      y=exp(ly)
    end if
  end function dlaplace

  pure real(dp) function plaplace(q,location,scale,lower_tail,log_p) result(y)
    real(dp), intent(in) :: q
    real(dp), intent(in), optional :: location,scale
    logical, intent(in), optional :: lower_tail,log_p
    real(dp) :: mu,s,p
    logical :: lower,lp
    mu=0.0_dp; if(present(location))mu=location
    s=1.0_dp; if(present(scale))s=scale
    lower=.true.; if(present(lower_tail))lower=lower_tail
    lp=.false.; if(present(log_p))lp=log_p
    if(s<=0.0_dp) then
      y=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    if(q<mu) then
      p=0.5_dp*exp(sqrt2*(q-mu)/s)
    else
      p=1.0_dp-0.5_dp*exp(-sqrt2*(q-mu)/s)
    end if
    if(.not.lower)p=1.0_dp-p
    if(lp) then
      y=log(max(p,tiny(1.0_dp)))
    else
      y=p
    end if
  end function plaplace

  pure real(dp) function qlaplace(p,location,scale,lower_tail,log_p) result(x)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: location,scale
    logical, intent(in), optional :: lower_tail,log_p
    real(dp) :: mu,s,pp
    logical :: lower,lp
    mu=0.0_dp; if(present(location))mu=location
    s=1.0_dp; if(present(scale))s=scale
    lower=.true.; if(present(lower_tail))lower=lower_tail
    lp=.false.; if(present(log_p))lp=log_p
    if(s<=0.0_dp) then
      x=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    pp=p
    if(lp)pp=exp(pp)
    if(.not.lower)pp=1.0_dp-pp
    if(pp<0.0_dp .or. pp>1.0_dp) then
      x=ieee_value(0.0_dp,ieee_quiet_nan)
    else if(pp<=0.0_dp) then
      x=-huge(1.0_dp)
    else if(pp>=1.0_dp) then
      x=huge(1.0_dp)
    else if(pp<0.5_dp) then
      x=mu+s*inv_sqrt2*log(2.0_dp*pp)
    else
      x=mu-s*inv_sqrt2*log(2.0_dp*(1.0_dp-pp))
    end if
  end function qlaplace

  real(dp) function rlaplace(location,scale) result(x)
    real(dp), intent(in), optional :: location,scale
    real(dp) :: mu,s,u
    mu=0.0_dp; if(present(location))mu=location
    s=1.0_dp; if(present(scale))s=scale
    if(s<=0.0_dp) then
      x=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    call random_number(u)
    u=max(min(u,1.0_dp-epsilon(1.0_dp)),tiny(1.0_dp))
    x=qlaplace(u,mu,s)
  end function rlaplace

  pure function dlaplace_vec(x,location,scale,log_density) result(y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: location,scale
    logical, intent(in), optional :: log_density
    real(dp) :: y(size(x)),mu,s
    logical :: lg
    integer :: i
    mu=0.0_dp;if(present(location))mu=location
    s=1.0_dp;if(present(scale))s=scale
    lg=.false.;if(present(log_density))lg=log_density
    do i=1,size(x); y(i)=dlaplace(x(i),mu,s,lg); end do
  end function dlaplace_vec

  pure function plaplace_vec(x,location,scale,lower_tail,log_p) result(y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: location,scale
    logical, intent(in), optional :: lower_tail,log_p
    real(dp) :: y(size(x)),mu,s
    logical :: lower,lp
    integer :: i
    mu=0.0_dp;if(present(location))mu=location
    s=1.0_dp;if(present(scale))s=scale
    lower=.true.;if(present(lower_tail))lower=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    do i=1,size(x); y(i)=plaplace(x(i),mu,s,lower,lp); end do
  end function plaplace_vec

  pure function qlaplace_vec(p,location,scale,lower_tail,log_p) result(x)
    real(dp), intent(in) :: p(:)
    real(dp), intent(in), optional :: location,scale
    logical, intent(in), optional :: lower_tail,log_p
    real(dp) :: x(size(p)),mu,s
    logical :: lower,lp
    integer :: i
    mu=0.0_dp;if(present(location))mu=location
    s=1.0_dp;if(present(scale))s=scale
    lower=.true.;if(present(lower_tail))lower=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    do i=1,size(p); x(i)=qlaplace(p(i),mu,s,lower,lp); end do
  end function qlaplace_vec

  function rlaplace_vec(n,location,scale) result(x)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: location,scale
    real(dp) :: x(n),mu,s
    integer :: i
    mu=0.0_dp;if(present(location))mu=location
    s=1.0_dp;if(present(scale))s=scale
    do i=1,n;x(i)=rlaplace(mu,s);end do
  end function rlaplace_vec

  real(dp) function log_dmlaplace(x,center,scatter) result(y)
    real(dp), intent(in) :: x(:),center(:),scatter(:,:)
    real(dp) :: d2,ld
    integer :: p,info
    p=size(x)
    if(size(center)/=p .or. size(scatter,1)/=p .or. size(scatter,2)/=p) then
      y=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    ld=logdet_chol(scatter,info)
    if(info/=0) then
      y=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    d2=mahalanobis_one(x,center,scatter)
    y=log_gamma(0.5_dp*real(p,dp))-0.5_dp*real(p,dp)*log(pi)- &
      log_gamma(real(p,dp))-(real(p,dp)+1.0_dp)*log(2.0_dp)-ld-0.5_dp*sqrt(max(d2,0.0_dp))
  end function log_dmlaplace

  real(dp) function dmlaplace(x,center,scatter,log_density) result(y)
    real(dp), intent(in) :: x(:),center(:),scatter(:,:)
    logical, intent(in), optional :: log_density
    real(dp) :: ly
    logical :: lg
    lg=.false.;if(present(log_density))lg=log_density
    ly=log_dmlaplace(x,center,scatter)
    if(lg) then;y=ly;else;y=exp(ly);end if
  end function dmlaplace

  subroutine rmlaplace(n,center,scatter,y,info)
    integer, intent(in) :: n
    real(dp), intent(in) :: center(:),scatter(:,:)
    real(dp), intent(out) :: y(:,:)
    integer, intent(out), optional :: info
    integer :: p,i,j,ier
    real(dp) :: l(size(scatter,1),size(scatter,2)),z(size(center)),nr,radial
    p=size(center)
    if(present(info))info=0
    if(size(scatter,1)/=p.or.size(scatter,2)/=p.or.size(y,1)/=n.or.size(y,2)/=p) then
      if(present(info))info=-1
      return
    end if
    call chol_lower(scatter,l,ier)
    if(ier/=0) then
      if(present(info))info=ier
      return
    end if
    do i=1,n
      do j=1,p;z(j)=normal_rand();end do
      nr=sqrt(sum(z*z))
      if(nr<=tiny(1.0_dp)) then
        z=0.0_dp;z(1)=1.0_dp;nr=1.0_dp
      end if
      radial=gamma_rand_l1(real(p,dp),2.0_dp)
      y(i,:)=center+matmul(l,z*(radial/nr))
    end do
  end subroutine rmlaplace

end module l1pack_distributions
