module normalp_estimation
  use normalp_special, only: dp
  implicit none
  private
  public :: estimatep, paramp_fit, kurtosis_p, normalp_params
  type :: normalp_params
    real(dp) :: mean = 0.0_dp
    real(dp) :: mp = 0.0_dp
    real(dp) :: sd = 0.0_dp
    real(dp) :: sp = 0.0_dp
    real(dp) :: p = 2.0_dp
    integer :: no_conv = 0
  end type normalp_params
contains
  pure function meanv(x) result(v)
    real(dp), intent(in) :: x(:); real(dp) :: v
    v=sum(x)/real(size(x),dp)
  end function

  pure function estimatep(x, mu, p0, direct) result(pp)
    real(dp), intent(in) :: x(:), mu
    real(dp), intent(in), optional :: p0
    logical, intent(in), optional :: direct
    real(dp) :: pp, pstart, ssp, sp, sa, sb, vi, zi, tz, yy
    logical :: dir
    pstart=2.0_dp; if(present(p0)) pstart=p0
    dir=.false.; if(present(direct)) dir=direct
    ssp=sum(abs(x-mu)**pstart)/real(size(x),dp)
    sp=ssp**(1.0_dp/pstart)
    if(sp<=tiny(1.0_dp)) then; pp=2.0_dp; return; end if
    sa=sum(abs((x-mu)/sp)); sb=sum(((x-mu)/sp)**2)
    if(sa<=tiny(1.0_dp)) then; pp=2.0_dp; return; end if
    vi=sqrt(real(size(x),dp)*sb)/sa
    vi=vi+((vi-1.0_dp)/real(size(x),dp))*5.0_dp
    if(dir) then
      pp=solve_direct(vi)
      if(pp>10.0_dp) pp=11.5_dp
      return
    end if
    if(vi<1.1547005_dp) then
      pp=11.5_dp
    else
      zi=(1.4142135_dp-vi)/0.259513_dp
      if(zi<0.0_dp) then; pp=1.0_dp
      else if(zi<0.6200052_dp) then
        tz=zi/0.6200052_dp; yy=((0.4738581_dp*tz-0.4966873_dp)*tz+1.0532646_dp)*tz+1.2246159_dp; pp=1.0_dp+tz**yy
      else if(zi<0.7914632_dp) then
        tz=(zi-0.6200052_dp)/0.1714592_dp; yy=((0.5246979_dp*tz-0.8167733_dp)*tz+0.8805483_dp)*tz+1.0859246_dp; pp=2.0_dp+tz**yy
      else if(zi<0.8670333_dp) then
        tz=(zi-0.7914632_dp)/0.0755701_dp; yy=((0.0743092_dp*tz-0.1269859_dp)*tz+0.3588207_dp)*tz+1.1227837_dp; pp=3.0_dp+tz**yy
      else if(zi<0.9072536_dp) then
        tz=(zi-0.8670333_dp)/0.0402203_dp; yy=((0.1097723_dp*tz-0.2127039_dp)*tz+0.3529203_dp)*tz+1.0761256_dp; pp=4.0_dp+tz**yy
      else if(zi<0.9314555_dp) then
        tz=(zi-0.9072536_dp)/0.0242019_dp; yy=((0.0955441_dp*tz-0.1891569_dp)*tz+0.2961275_dp)*tz+1.0631784_dp; pp=5.0_dp+tz**yy
      else if(zi<0.9472072_dp) then
        tz=(zi-0.9314555_dp)/0.0157518_dp; yy=((0.0862627_dp*tz-0.1725326_dp)*tz+0.256885_dp)*tz+1.0540746_dp; pp=6.0_dp+tz**yy
      else if(zi<0.9580557_dp) then
        tz=(zi-0.9472072_dp)/0.0108484_dp; yy=((0.078785_dp*tz-0.1581388_dp)*tz+0.227011_dp)*tz+1.04735_dp; pp=7.0_dp+tz**yy
      else if(zi<0.9658545_dp) then
        tz=(zi-0.9580557_dp)/0.0077988_dp; yy=((0.0663921_dp*tz-0.1380841_dp)*tz+0.2010053_dp)*tz+1.0422984_dp; pp=8.0_dp+tz**yy
      else if(zi<0.9716534_dp) then
        tz=(zi-0.9658545_dp)/0.0057989_dp; yy=((0.0557199_dp*tz-0.1184033_dp)*tz+0.178176_dp)*tz+1.038582_dp; pp=9.0_dp+tz**yy
      else
        pp=10.0_dp+(zi-0.9716534_dp)/0.0283466_dp
      end if
    end if
  end function estimatep

  pure function solve_direct(vi) result(pbest)
    real(dp),intent(in)::vi; real(dp)::pbest,a,b,c,d,fc,fd
    integer::i
    a=1.0_dp; b=10.0_dp
    do i=1,100
      c=b-(b-a)/1.6180339887498948_dp; d=a+(b-a)/1.6180339887498948_dp
      fc=obj(c); fd=obj(d)
      if(fc<fd) then; b=d; else; a=c; end if
    end do
    pbest=0.5_dp*(a+b)
  contains
    pure function obj(p) result(v)
      real(dp),intent(in)::p; real(dp)::v,t
      t=sqrt(gamma(1.0_dp/p)*gamma(3.0_dp/p))/gamma(2.0_dp/p)
      v=(vi-t)**2
    end function
  end function solve_direct

  subroutine paramp_fit(x, fit, p_fixed)
    real(dp),intent(in)::x(:); type(normalp_params),intent(out)::fit
    real(dp),intent(in),optional::p_fixed
    real(dp)::p,pp,mp,oldmp,df
    integer::i,n
    n=size(x); fit%mean=meanv(x); fit%sd=sqrt(sum((x-fit%mean)**2)/real(n,dp)); mp=fit%mean
    if(present(p_fixed)) then
      p=p_fixed; df=1.0_dp; call lp_location(x,p,mp); fit%no_conv=0
    else
      df=2.0_dp; pp=2.0_dp; p=estimatep(x,mp,pp); fit%no_conv=0
      do i=1,100
        pp=p; oldmp=mp; call lp_location(x,pp,mp); p=estimatep(x,mp,pp)
        if(abs(p-pp)<=1.0e-4_dp .and. abs(mp-oldmp)<=1.0e-4_dp) exit
        if(i==100) fit%no_conv=1
      end do
    end if
    if(abs(p-1.0_dp)<1.0e-12_dp) mp=median_copy(x)
    if(p>=11.5_dp) then; mp=0.5_dp*(maxval(x)+minval(x)); fit%sp=0.5_dp*(maxval(x)-minval(x))
    else; fit%sp=(sum(abs(x-mp)**p)/real(n-int(df),dp))**(1.0_dp/p); end if
    fit%mp=mp; fit%p=p
  end subroutine paramp_fit

  subroutine lp_location(x,p,mp)
    real(dp),intent(in)::x(:),p; real(dp),intent(inout)::mp
    real(dp)::a,b,c,d,fc,fd
    integer::i
    if(abs(p-1.0_dp)<1.0e-12_dp) then; mp=median_copy(x); return; end if
    a=minval(x); b=maxval(x)
    do i=1,100
      c=b-(b-a)/1.6180339887498948_dp; d=a+(b-a)/1.6180339887498948_dp
      fc=sum(abs(x-c)**p); fd=sum(abs(x-d)**p)
      if(fc<fd) then; b=d; else; a=c; end if
    end do
    mp=0.5_dp*(a+b)
  end subroutine lp_location

  function median_copy(x) result(m)
    real(dp),intent(in)::x(:); real(dp)::m; real(dp),allocatable::z(:); real(dp)::tmp
    integer::i,j,n
    z=x; n=size(z)
    do i=2,n; tmp=z(i); j=i-1; do while(j>=1 .and. z(j)>tmp); z(j+1)=z(j); j=j-1; end do; z(j+1)=tmp; end do
    if(mod(n,2)==1) then; m=z((n+1)/2); else; m=0.5_dp*(z(n/2)+z(n/2+1)); end if
  end function median_copy

  subroutine kurtosis_p(stats, p, x, use_estimate)
    real(dp),intent(out)::stats(3); real(dp),intent(in),optional::p,x(:); logical,intent(in),optional::use_estimate
    real(dp)::pv,mp,vi,b2,bp,spv; integer::n; type(normalp_params)::fit; logical::uest
    if(.not.present(x)) then
      pv=2.0_dp; if(present(p)) pv=p
      vi=sqrt(gamma(1.0_dp/pv)*gamma(3.0_dp/pv))/gamma(2.0_dp/pv)
      b2=gamma(1.0_dp/pv)*gamma(5.0_dp/pv)/gamma(3.0_dp/pv)**2; bp=pv+1.0_dp
    else
      n=size(x); uest=.true.; if(present(use_estimate)) uest=use_estimate
      if(uest) then; call paramp_fit(x,fit); pv=fit%p; mp=fit%mp
      else; pv=2.0_dp; if(present(p)) pv=p; call paramp_fit(x,fit,pv); mp=fit%mp; end if
      vi=sqrt(real(n,dp)*sum((x-mp)**2))/sum(abs(x-mp)); b2=real(n,dp)*sum((x-mp)**4)/sum((x-mp)**2)**2
      spv=(sum(abs(x-mp)**pv)/real(n,dp))**(1.0_dp/pv); bp=sum(abs(x-mp)**(2.0_dp*pv))/(real(n,dp)*spv**(2.0_dp*pv))
    end if
    stats=[vi,b2,bp]
  end subroutine kurtosis_p
end module normalp_estimation
