! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_denoise
  use waveslim_kinds, only : dp
  use waveslim_types, only : wavelet_transform, wavelet_transform_2d
  use waveslim_math, only : median_value
  implicit none
  private
  public :: soft, bishrink, manual_thresh, universal_thresh
  public :: universal_thresh_modwt, sure_thresh, hybrid_thresh, da_thresh
  public :: denoise_dwt_2d, denoise_modwt_2d
contains
  pure elemental function soft(x,threshold) result(y)
    real(dp),intent(in)::x,threshold
    real(dp)::y,a
    a=max(abs(x)-threshold,0.0_dp)
    if(a+threshold>0.0_dp)then
    y=a/(a+threshold)*x
    else
    y=0.0_dp
    end if
  end function soft

  pure elemental function bishrink(y1,y2,threshold) result(w1)
    real(dp),intent(in)::y1,y2,threshold
    real(dp)::w1,r
    r=max(sqrt(y1*y1+y2*y2)-threshold,0.0_dp)
    if(r+threshold>0.0_dp)then
    w1=y1*r/(r+threshold)
    else
    w1=0.0_dp
    end if
  end function bishrink

  function manual_thresh(wt,max_level,value,hard) result(ans)
    type(wavelet_transform),intent(in)::wt
    integer,intent(in),optional::max_level
    real(dp),intent(in)::value
    logical,intent(in),optional::hard
    type(wavelet_transform)::ans
    integer::j,jmax
    real(dp)::sigma,t
    logical::h
    ans=wt
    jmax=min(4,wt%levels())
    if(present(max_level))jmax=min(max_level,wt%levels())
    h=.true.
    if(present(hard))h=hard
    sigma=median_value(abs(wt%detail(1)%values))/0.6745_dp
    t=sigma*value
    do j=1,jmax
    call apply_threshold(ans%detail(j)%values,t,h)
    end do
  end function manual_thresh

  function universal_thresh(wt,max_level,hard) result(ans)
    type(wavelet_transform),intent(in)::wt
    integer,intent(in),optional::max_level
    logical,intent(in),optional::hard
    type(wavelet_transform)::ans
    integer::n,j
    n=wt%original_length
    if(n<=0)n=sum([(size(wt%detail(j)%values),j=1,wt%levels())])+size(wt%smooth)
    ans=manual_thresh(wt,max_level,sqrt(2.0_dp*log(real(max(n,2),dp))),hard)
  end function universal_thresh

  function universal_thresh_modwt(wt,max_level,hard) result(ans)
    type(wavelet_transform),intent(in)::wt
    integer,intent(in),optional::max_level
    logical,intent(in),optional::hard
    type(wavelet_transform)::ans
    integer::j,jmax,n
    real(dp)::sigma,t
    logical::h
    ans=wt
    jmax=min(4,wt%levels())
    if(present(max_level))jmax=min(max_level,wt%levels())
    h=.true.
    if(present(hard))h=hard
    n=size(wt%detail(1)%values)
    sigma=sqrt(2.0_dp)*median_value(abs(wt%detail(1)%values))/0.6745_dp
    do j=1,jmax
    t=sigma*sqrt(2.0_dp*log(real(max(n,2),dp)))/sqrt(real(2**j,dp))
    call apply_threshold(ans%detail(j)%values,t,h)
    end do
  end function universal_thresh_modwt

  function sure_thresh(wt,max_level,hard) result(ans)
    type(wavelet_transform),intent(in)::wt
    integer,intent(in),optional::max_level
    logical,intent(in),optional::hard
    type(wavelet_transform)::ans
    integer::j,jmax,n,i
    real(dp)::sigma,t,best,risk,best_risk
    real(dp),allocatable::a(:),squares(:)
    logical::h
    ans=wt
    jmax=min(4,wt%levels())
    if(present(max_level))jmax=min(max_level,wt%levels())
    h=.true.
    if(present(hard))h=hard
    do j=1,jmax
      n=size(wt%detail(j)%values)
      sigma=median_value(abs(wt%detail(j)%values))/0.6745_dp
      if(sigma<=tiny(1.0_dp))cycle
      a=abs(wt%detail(j)%values)/sigma
      call sort_local(a)
      squares=a*a
      best_risk=huge(1.0_dp)
      best=a(n)
      do i=1,n
      t=a(i)
      risk=real(n-2*i,dp)+sum(min(squares,t*t))
      if(risk<best_risk)then
      best_risk=risk
      best=t
      end if
      end do
      call apply_threshold(ans%detail(j)%values,sigma*best,h)
    end do
  end function sure_thresh

  function hybrid_thresh(wt,max_level) result(ans)
    type(wavelet_transform),intent(in)::wt
    integer,intent(in),optional::max_level
    type(wavelet_transform)::ans
    integer::j,jmax,n
    real(dp)::energy,limit,sigma,t
    ans=wt
    jmax=min(4,wt%levels())
    if(present(max_level))jmax=min(max_level,wt%levels())
    do j=1,jmax
      n=size(wt%detail(j)%values)
      sigma=median_value(abs(wt%detail(j)%values))/0.6745_dp
      if(sigma<=tiny(1.0_dp))cycle
      energy=sum((wt%detail(j)%values/sigma)**2)
      limit=sqrt(real(j**3,dp)/real(2**j,dp))
      if((energy-real(n,dp))/real(n,dp)<=limit)then
      t=sigma*sqrt(2.0_dp*log(real(max(n,2),dp)))
      else
      t=sigma*sure_threshold_value(wt%detail(j)%values/sigma)
      end if
      call apply_threshold(ans%detail(j)%values,t,.false.)
    end do
  end function hybrid_thresh

  function da_thresh(wt,alpha,max_level,thresholds) result(ans)
    type(wavelet_transform),intent(in)::wt
    real(dp),intent(in),optional::alpha
    integer,intent(in),optional::max_level
    real(dp),allocatable,intent(out),optional::thresholds(:)
    type(wavelet_transform)::ans
    integer::j,jmax,n,keep
    real(dp)::a,t
    real(dp),allocatable::s(:)
    a=0.05_dp
    if(present(alpha))a=alpha
    ans=wt
    jmax=min(4,wt%levels())
    if(present(max_level))jmax=min(max_level,wt%levels())
    if(present(thresholds))allocate(thresholds(jmax))
    do j=1,jmax
      n=size(wt%detail(j)%values)
      s=wt%detail(j)%values**2
      call sort_local(s)
      keep=max(1,min(n,nint((1.0_dp-a)*real(n,dp))))
      t=sqrt(max(s(keep),0.0_dp))
      if(present(thresholds))thresholds(j)=t
      call apply_threshold(ans%detail(j)%values,t,.false.)
    end do
  end function da_thresh

  function denoise_dwt_2d(wt,method,hard) result(ans)
    type(wavelet_transform_2d),intent(in)::wt
    character(len=*),intent(in),optional::method
    logical,intent(in),optional::hard
    type(wavelet_transform_2d)::ans
    character(len=16)::m
    logical::h
    integer::j,n
    real(dp)::sigma,t
    ans=wt
    m='universal'
    if(present(method))m=trim(method)
    h=.true.
    if(present(hard))h=hard
    sigma=median_value(abs(reshape(wt%level(1)%hh,[size(wt%level(1)%hh)])))/0.6745_dp
    n=product(wt%original_shape)
    do j=1,size(wt%level)
    if(m=='sure')then
    t=sigma*sure_threshold_value(reshape(wt%level(j)%hh,[size(wt%level(j)%hh)]))
    else
    t=sigma*sqrt(2.0_dp*log(real(max(n,2),dp)))
    end if
    call apply_threshold_2d(ans%level(j)%lh,t,h)
    call apply_threshold_2d(ans%level(j)%hl,t,h)
    call apply_threshold_2d(ans%level(j)%hh,t,h)
    end do
  end function denoise_dwt_2d

  function denoise_modwt_2d(wt,method,hard) result(ans)
    type(wavelet_transform_2d),intent(in)::wt
    character(len=*),intent(in),optional::method
    logical,intent(in),optional::hard
    type(wavelet_transform_2d)::ans
    integer::j,n
    real(dp)::sigma,t
    logical::h
    character(len=16) :: m
    ans=wt
    m='universal'
    if(present(method))m=trim(method)
    h=.true.
    if(present(hard))h=hard
    sigma=sqrt(2.0_dp)*median_value(abs(reshape(wt%level(1)%hh,[size(wt%level(1)%hh)])))/0.6745_dp
    n=product(wt%original_shape)
    do j=1,size(wt%level)
    if(m=='sure')then
    t=sigma*sure_threshold_value(reshape(wt%level(j)%hh,[size(wt%level(j)%hh)]))/sqrt(real(4**j,dp))
    else
    t=sigma*sqrt(2.0_dp*log(real(max(n,2),dp)))/sqrt(real(4**j,dp))
    end if
    call apply_threshold_2d(ans%level(j)%lh,t,h)
    call apply_threshold_2d(ans%level(j)%hl,t,h)
    call apply_threshold_2d(ans%level(j)%hh,t,h)
    end do
  end function denoise_modwt_2d

  subroutine apply_threshold(x,t,hard)
    real(dp),intent(inout)::x(:)
    real(dp),intent(in)::t
    logical,intent(in)::hard
    if(hard)then
    where(abs(x)<=t)x=0.0_dp
    else
    x=sign(max(abs(x)-t,0.0_dp),x)
    end if
  end subroutine apply_threshold

  subroutine apply_threshold_2d(x,t,hard)
    real(dp),intent(inout)::x(:,:)
    real(dp),intent(in)::t
    logical,intent(in)::hard
    if(hard)then
    where(abs(x)<=t)x=0.0_dp
    else
    x=sign(max(abs(x)-t,0.0_dp),x)
    end if
  end subroutine apply_threshold_2d

  function sure_threshold_value(x) result(best)
    real(dp),intent(in)::x(:)
    real(dp)::best,best_risk,risk,t
    real(dp),allocatable::a(:),sq(:)
    integer::i,n
    n=size(x)
    a=abs(x)
    call sort_local(a)
    sq=a*a
    best=a(n)
    best_risk=huge(1.0_dp)
    do i=1,n
    t=a(i)
    risk=real(n-2*i,dp)+sum(min(sq,t*t))
    if(risk<best_risk)then
    best_risk=risk
    best=t
    end if
    end do
  end function sure_threshold_value

  subroutine sort_local(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::key
    do i=2,size(x)
    key=x(i)
    j=i-1
    do while(j>=1)
    if(x(j)<=key)exit
    x(j+1)=x(j)
    j=j-1
    end do
    x(j+1)=key
    end do
  end subroutine sort_local
end module waveslim_denoise
