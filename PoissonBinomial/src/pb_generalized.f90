! SPDX-License-Identifier: GPL-3.0-only
module pb_generalized
  use pb_kinds, only : dp, pi
  use pb_math, only : normal_cdf, normal_pdf
  use pb_numerics, only : normalize_pmf, convolve_fft, gcd_vector, expand_gpb, &
                          sample_from_pmf
  implicit none
  private

  type, public :: gpb_table
    integer :: lower = 0
    integer :: upper = -1
    real(dp), allocatable :: values(:)  ! values(i-lower+1) for outcome i
  end type gpb_table

  public :: dgpbinom, dgpbinom_at, dgpbinom_values
  public :: pgpbinom, pgpbinom_at, pgpbinom_values
  public :: qgpbinom, qgpbinom_values, rgpbinom
  public :: dgpb_convolve, dgpb_dividefft, dgpb_characteristic, pgpb_normal

contains

  subroutine prepare_gpb(probs, val_p, val_q, wts, complete_lo, complete_hi, &
                         inner_lo, inner_hi, phigh, steps)
    real(dp), intent(in) :: probs(:)
    integer, intent(in) :: val_p(:), val_q(:)
    integer, intent(in), optional :: wts(:)
    integer, intent(out) :: complete_lo, complete_hi, inner_lo, inner_hi
    real(dp), allocatable, intent(out) :: phigh(:)
    integer, allocatable, intent(out) :: steps(:)
    real(dp), allocatable :: p(:)
    integer, allocatable :: vp(:), vq(:)
    integer :: i, nr, k, lo, hi, sure, base
    real(dp) :: ph

    call expand_gpb(probs,val_p,val_q,wts,p,vp,vq)
    if (any(p < 0.0_dp) .or. any(p > 1.0_dp)) then
      error stop "probabilities must lie in [0,1]"
    end if
    complete_lo = 0
    complete_hi = 0
    do i = 1, size(p)
      complete_lo = complete_lo + min(vp(i),vq(i))
      complete_hi = complete_hi + max(vp(i),vq(i))
    end do

    nr = 0
    sure = 0
    base = 0
    do i = 1, size(p)
      lo = min(vp(i),vq(i)); hi = max(vp(i),vq(i))
      if (vp(i) >= vq(i)) then
        ph = p(i)
      else
        ph = 1.0_dp-p(i)
      end if
      if (hi == lo) then
        sure = sure + lo
      else if (ph <= 0.0_dp) then
        sure = sure + lo
      else if (ph >= 1.0_dp) then
        sure = sure + hi
      else
        nr = nr + 1
        base = base + lo
      end if
    end do
    allocate(phigh(nr),steps(nr))
    k = 0
    do i = 1, size(p)
      lo = min(vp(i),vq(i)); hi = max(vp(i),vq(i))
      if (vp(i) >= vq(i)) then
        ph = p(i)
      else
        ph = 1.0_dp-p(i)
      end if
      if (hi /= lo .and. ph > 0.0_dp .and. ph < 1.0_dp) then
        k = k + 1
        phigh(k) = ph
        steps(k) = hi-lo
      end if
    end do
    inner_lo = sure + base
    inner_hi = inner_lo + sum(steps)
  end subroutine prepare_gpb

  function conv_core(phigh, steps) result(pmf)
    real(dp), intent(in) :: phigh(:)
    integer, intent(in) :: steps(:)
    real(dp), allocatable :: pmf(:), old(:)
    integer :: i, j, last, d, total
    if (size(phigh) == 0) then
      allocate(pmf(0:0)); pmf(0)=1.0_dp; return
    end if
    total = sum(steps)
    allocate(pmf(0:total)); pmf=0.0_dp; pmf(0)=1.0_dp
    last=0
    do i=1,size(phigh)
      d=steps(i)
      old=pmf
      pmf=0.0_dp
      do j=0,last
        if(abs(old(j))<=tiny(1.0_dp)) cycle
        pmf(j)=pmf(j)+old(j)*(1.0_dp-phigh(i))
        pmf(j+d)=pmf(j+d)+old(j)*phigh(i)
      end do
      last=last+d
    end do
    call normalize_pmf(pmf)
  end function conv_core

  recursive function divide_core(phigh, steps) result(pmf)
    real(dp), intent(in) :: phigh(:)
    integer, intent(in) :: steps(:)
    real(dp), allocatable :: pmf(:), left(:), right(:)
    integer :: n,m,d
    n=size(phigh)
    if(n==0) then
      allocate(pmf(0:0)); pmf(0)=1.0_dp
    else if(n==1) then
      d=steps(1); allocate(pmf(0:d)); pmf=0.0_dp
      pmf(0)=1.0_dp-phigh(1); pmf(d)=phigh(1)
    else if(n<=20 .or. sum(steps)<=64) then
      pmf=conv_core(phigh,steps)
    else
      m=n/2
      left=divide_core(phigh(:m),steps(:m))
      right=divide_core(phigh(m+1:),steps(m+1:))
      pmf=convolve_fft(left,right)
      call normalize_pmf(pmf)
    end if
  end function divide_core

  function char_core(phigh, steps) result(pmf)
    real(dp), intent(in) :: phigh(:)
    integer, intent(in) :: steps(:)
    real(dp), allocatable :: pmf(:)
    complex(dp), allocatable :: phi(:)
    complex(dp) :: z, prod, accum
    integer :: s,nout,l,j,k
    real(dp) :: theta
    s=sum(steps); nout=s+1
    allocate(phi(0:s),pmf(0:s))
    do l=0,s
      prod=cmplx(1.0_dp,0.0_dp,kind=dp)
      do j=1,size(phigh)
        theta=2.0_dp*pi*real(l*steps(j),dp)/real(nout,dp)
        z=cmplx(cos(theta),sin(theta),kind=dp)
        prod=prod*((1.0_dp-phigh(j))+phigh(j)*z)
      end do
      phi(l)=prod
    end do
    do k=0,s
      accum=cmplx(0.0_dp,0.0_dp,kind=dp)
      do l=0,s
        theta=-2.0_dp*pi*real(k*l,dp)/real(nout,dp)
        accum=accum+phi(l)*cmplx(cos(theta),sin(theta),kind=dp)
      end do
      pmf(k)=real(accum,dp)/real(nout,dp)
    end do
    where(pmf<2.22e-16_dp) pmf=0.0_dp
    where(pmf>1.0_dp) pmf=1.0_dp
    call normalize_pmf(pmf)
  end function char_core

  function scaled_core(phigh, steps, which) result(full)
    real(dp), intent(in) :: phigh(:)
    integer, intent(in) :: steps(:)
    character(len=*), intent(in) :: which
    real(dp), allocatable :: full(:), small(:)
    integer, allocatable :: ds(:)
    integer :: g,i,total
    if(size(steps)==0) then
      allocate(full(0:0)); full(0)=1.0_dp; return
    end if
    g=gcd_vector(steps); if(g<=0) g=1
    allocate(ds(size(steps))); ds=steps/g
    allocate(small(0:0)); small = 0.0_dp
    select case(trim(which))
    case("Convolve"); small=conv_core(phigh,ds)
    case("DivideFFT"); small=divide_core(phigh,ds)
    case("Characteristic"); small=char_core(phigh,ds)
    case default; error stop "unknown exact generalized method"
    end select
    total=sum(steps)
    allocate(full(0:total)); full=0.0_dp
    do i=1,size(small)
      full((i-1)*g)=small(i)
    end do
    call normalize_pmf(full)
  end function scaled_core

  function dgpb_convolve(probs,val_p,val_q,wts) result(tab)
    real(dp), intent(in) :: probs(:)
    integer, intent(in) :: val_p(:),val_q(:)
    integer, intent(in), optional :: wts(:)
    type(gpb_table) :: tab
    integer :: clo,chi,ilo,ihi
    real(dp),allocatable::ph(:),inner(:)
    integer,allocatable::st(:)
    call prepare_gpb(probs,val_p,val_q,wts,clo,chi,ilo,ihi,ph,st)
    tab%lower=clo; tab%upper=chi
    allocate(tab%values(chi-clo+1)); tab%values=0.0_dp
    if(size(ph)==0) then
      tab%values(ilo-clo+1)=1.0_dp
    else
      inner=scaled_core(ph,st,"Convolve")
      tab%values(ilo-clo+1:ihi-clo+1)=inner
    end if
  end function dgpb_convolve

  function dgpb_dividefft(probs,val_p,val_q,wts) result(tab)
    real(dp), intent(in) :: probs(:)
    integer, intent(in) :: val_p(:),val_q(:)
    integer, intent(in), optional :: wts(:)
    type(gpb_table) :: tab
    integer :: clo,chi,ilo,ihi
    real(dp),allocatable::ph(:),inner(:)
    integer,allocatable::st(:)
    call prepare_gpb(probs,val_p,val_q,wts,clo,chi,ilo,ihi,ph,st)
    tab%lower=clo; tab%upper=chi
    allocate(tab%values(chi-clo+1)); tab%values=0.0_dp
    if(size(ph)==0) then
      tab%values(ilo-clo+1)=1.0_dp
    else
      inner=scaled_core(ph,st,"DivideFFT")
      tab%values(ilo-clo+1:ihi-clo+1)=inner
    end if
  end function dgpb_dividefft

  function dgpb_characteristic(probs,val_p,val_q,wts) result(tab)
    real(dp), intent(in) :: probs(:)
    integer, intent(in) :: val_p(:),val_q(:)
    integer, intent(in), optional :: wts(:)
    type(gpb_table) :: tab
    integer :: clo,chi,ilo,ihi
    real(dp),allocatable::ph(:),inner(:)
    integer,allocatable::st(:)
    call prepare_gpb(probs,val_p,val_q,wts,clo,chi,ilo,ihi,ph,st)
    tab%lower=clo; tab%upper=chi
    allocate(tab%values(chi-clo+1)); tab%values=0.0_dp
    if(size(ph)==0) then
      tab%values(ilo-clo+1)=1.0_dp
    else
      inner=scaled_core(ph,st,"Characteristic")
      tab%values(ilo-clo+1:ihi-clo+1)=inner
    end if
  end function dgpb_characteristic

  function pgpb_normal(probs,val_p,val_q,wts,refined,lower_tail) result(tab)
    real(dp), intent(in) :: probs(:)
    integer, intent(in) :: val_p(:),val_q(:)
    integer, intent(in), optional :: wts(:)
    logical, intent(in), optional :: refined,lower_tail
    type(gpb_table) :: tab
    integer :: clo,chi,ilo,ihi,g,i,z,rz
    real(dp),allocatable::ph(:)
    integer,allocatable::st(:),ds(:)
    real(dp)::mu,sigma,gamma3,x,val
    logical::r,lower
    r=.true.; if(present(refined)) r=refined
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    call prepare_gpb(probs,val_p,val_q,wts,clo,chi,ilo,ihi,ph,st)
    tab%lower=clo; tab%upper=chi
    allocate(tab%values(chi-clo+1))
    tab%values=merge(0.0_dp,1.0_dp,lower)
    if(size(ph)==0) then
      do i=clo,chi
        if(lower) then
          tab%values(i-clo+1)=merge(1.0_dp,0.0_dp,i>=ilo)
        else
          tab%values(i-clo+1)=merge(0.0_dp,1.0_dp,i>=ilo)
        end if
      end do
      return
    end if
    g=gcd_vector(st); if(g<=0) g=1
    allocate(ds(size(st))); ds=st/g
    mu=sum(ph*real(ds,dp))
    sigma=sqrt(sum(ph*(1.0_dp-ph)*real(ds,dp)**2))
    if(sigma>0.0_dp) then
      gamma3=sum(ph*(1.0_dp-ph)*(1.0_dp-2.0_dp*ph)*real(ds,dp)**3)/sigma**3
    else
      gamma3=0.0_dp
    end if
    do i=clo,chi
      if(i<ilo) then
        tab%values(i-clo+1)=merge(0.0_dp,1.0_dp,lower)
      else if(i>=ihi) then
        tab%values(i-clo+1)=merge(1.0_dp,0.0_dp,lower)
      else if(sigma<=tiny(1.0_dp)) then
        tab%values(i-clo+1)=merge(1.0_dp,0.0_dp,lower)
      else
        z=i-ilo
        rz=z/g
        x=(real(rz,dp)+0.5_dp-mu)/sigma
        val=normal_cdf(x)
        if(.not.lower) val=1.0_dp-val
        if(r) then
          if(lower) then
            val=val+gamma3*(1.0_dp-x*x)*normal_pdf(x)/6.0_dp
          else
            val=val-gamma3*(1.0_dp-x*x)*normal_pdf(x)/6.0_dp
          end if
        end if
        tab%values(i-clo+1)=max(0.0_dp,min(1.0_dp,val))
      end if
    end do
  end function pgpb_normal

  function normal_pmf(probs,val_p,val_q,wts,refined) result(tab)
    real(dp), intent(in) :: probs(:)
    integer, intent(in)::val_p(:),val_q(:)
    integer,intent(in),optional::wts(:)
    logical,intent(in)::refined
    type(gpb_table)::tab,cl,cu
    integer::i,mid
    real(dp),allocatable::pexp(:)
    integer,allocatable::vp(:),vq(:)
    call expand_gpb(probs,val_p,val_q,wts,pexp,vp,vq)
    cl=pgpb_normal(probs,val_p,val_q,wts,refined,.true.)
    cu=pgpb_normal(probs,val_p,val_q,wts,refined,.false.)
    tab%lower=cl%lower; tab%upper=cl%upper
    allocate(tab%values(tab%upper-tab%lower+1)); tab%values=0.0_dp
    mid=floor(sum(pexp*real(vp,dp)+(1.0_dp-pexp)*real(vq,dp))+0.5_dp)
    tab%values(1)=cl%values(1)
    do i=tab%lower+1,tab%upper
      if(i<=mid) then
        tab%values(i-tab%lower+1)=cl%values(i-tab%lower+1)-cl%values(i-tab%lower)
      else
        tab%values(i-tab%lower+1)=cu%values(i-tab%lower)-cu%values(i-tab%lower+1)
      end if
    end do
    where(tab%values<0.0_dp .and. tab%values>-1.0e-12_dp) tab%values=0.0_dp
    call normalize_pmf(tab%values)
  end function normal_pmf

  function dgpbinom(probs,val_p,val_q,method,wts) result(tab)
    real(dp), intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::wts(:)
    type(gpb_table)::tab
    character(len=32)::m
    m="DivideFFT"; if(present(method)) m=method
    select case(trim(adjustl(m)))
    case("DivideFFT","dividefft","DIVIDEFFT")
      tab=dgpb_dividefft(probs,val_p,val_q,wts)
    case("Convolve","convolve","CONVOLVE")
      tab=dgpb_convolve(probs,val_p,val_q,wts)
    case("Characteristic","characteristic","CHARACTERISTIC")
      tab=dgpb_characteristic(probs,val_p,val_q,wts)
    case("Normal","normal","NORMAL")
      tab=normal_pmf(probs,val_p,val_q,wts,.false.)
    case("RefinedNormal","refinednormal","REFINEDNORMAL")
      tab=normal_pmf(probs,val_p,val_q,wts,.true.)
    case default
      error stop "unknown generalized Poisson-binomial method"
    end select
  end function dgpbinom

  real(dp) function dgpbinom_at(x,probs,val_p,val_q,method,wts,log_p) result(p)
    integer,intent(in)::x
    real(dp),intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::wts(:)
    logical,intent(in),optional::log_p
    type(gpb_table)::tab
    logical::lp
    lp=.false.; if(present(log_p)) lp=log_p
    tab=dgpbinom(probs,val_p,val_q,method,wts)
    if(x<tab%lower .or. x>tab%upper) then
      p=0.0_dp
    else
      p=tab%values(x-tab%lower+1)
    end if
    if(lp) p=log(p)
  end function dgpbinom_at

  function pgpbinom(probs,val_p,val_q,method,wts,lower_tail) result(tab)
    real(dp),intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::wts(:)
    logical,intent(in),optional::lower_tail
    type(gpb_table)::tab,dtab
    logical::lower
    integer::i,n
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    if(present(method)) then
      if(trim(adjustl(method))=="Normal" .or. trim(adjustl(method))=="normal" .or. &
         trim(adjustl(method))=="NORMAL") then
        tab=pgpb_normal(probs,val_p,val_q,wts,.false.,lower); return
      else if(trim(adjustl(method))=="RefinedNormal" .or. &
              trim(adjustl(method))=="refinednormal" .or. &
              trim(adjustl(method))=="REFINEDNORMAL") then
        tab=pgpb_normal(probs,val_p,val_q,wts,.true.,lower); return
      end if
    end if
    dtab=dgpbinom(probs,val_p,val_q,method,wts)
    tab%lower=dtab%lower; tab%upper=dtab%upper
    n=tab%upper-tab%lower+1
    allocate(tab%values(n))
    if(lower) then
      tab%values(1)=dtab%values(1)
      do i=2,n; tab%values(i)=tab%values(i-1)+dtab%values(i); end do
      tab%values(n)=1.0_dp
    else
      tab%values(n)=0.0_dp
      do i=n-1,1,-1; tab%values(i)=tab%values(i+1)+dtab%values(i+1); end do
    end if
  end function pgpbinom

  real(dp) function pgpbinom_at(x,probs,val_p,val_q,method,wts,lower_tail,log_p) result(p)
    integer,intent(in)::x
    real(dp),intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::wts(:)
    logical,intent(in),optional::lower_tail,log_p
    type(gpb_table)::tab
    logical::lower,lp
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    tab=pgpbinom(probs,val_p,val_q,method,wts,lower)
    if(x<tab%lower) then
      p=merge(0.0_dp,1.0_dp,lower)
    else if(x>tab%upper) then
      p=merge(1.0_dp,0.0_dp,lower)
    else
      p=tab%values(x-tab%lower+1)
    end if
    if(lp) p=log(p)
  end function pgpbinom_at

  integer function qgpbinom(prob,probs,val_p,val_q,method,wts,lower_tail,log_p) result(q)
    real(dp),intent(in)::prob
    real(dp),intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::wts(:)
    logical,intent(in),optional::lower_tail,log_p
    type(gpb_table)::dtab,ctab
    logical::lower,lp
    real(dp)::target
    integer::i,lo,hi
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    target=prob; if(lp) target=exp(prob)
    if(target<0.0_dp .or. target>1.0_dp) error stop "quantile probability outside [0,1]"
    dtab=dgpbinom(probs,val_p,val_q,method,wts)
    lo=dtab%lower; hi=dtab%upper
    do while (lo < hi)
      if (abs(dtab%values(lo-dtab%lower+1)) > tiny(1.0_dp)) exit
      lo = lo + 1
    end do
    do while (hi > lo)
      if (abs(dtab%values(hi-dtab%lower+1)) > tiny(1.0_dp)) exit
      hi = hi - 1
    end do
    if(lower) then
      if(target<=0.0_dp) then; q=lo; return; end if
      if(target>=1.0_dp) then; q=hi; return; end if
      ctab=pgpbinom(probs,val_p,val_q,method,wts,.true.)
      q=hi
      do i=lo,hi
        if(ctab%values(i-ctab%lower+1)>=target) then; q=i; exit; end if
      end do
    else
      if(target>=1.0_dp) then; q=lo; return; end if
      if(target<=0.0_dp) then; q=hi; return; end if
      ctab=pgpbinom(probs,val_p,val_q,method,wts,.false.)
      q=hi
      do i=lo,hi
        if(ctab%values(i-ctab%lower+1)<=target) then; q=i; exit; end if
      end do
    end if
  end function qgpbinom

  function dgpbinom_values(x,probs,val_p,val_q,method,wts,log_p) result(v)
    integer,intent(in)::x(:)
    real(dp),intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::wts(:)
    logical,intent(in),optional::log_p
    real(dp),allocatable::v(:)
    type(gpb_table)::tab
    logical::lp
    integer::i
    lp=.false.; if(present(log_p)) lp=log_p
    tab=dgpbinom(probs,val_p,val_q,method,wts)
    allocate(v(size(x)))
    do i=1,size(x)
      if(x(i)<tab%lower .or. x(i)>tab%upper) then
        v(i)=0.0_dp
      else
        v(i)=tab%values(x(i)-tab%lower+1)
      end if
    end do
    if(lp) v=log(v)
  end function dgpbinom_values

  function pgpbinom_values(x,probs,val_p,val_q,method,wts,lower_tail,log_p) result(v)
    integer,intent(in)::x(:)
    real(dp),intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::wts(:)
    logical,intent(in),optional::lower_tail,log_p
    real(dp),allocatable::v(:)
    type(gpb_table)::tab
    logical::lower,lp
    integer::i
    lower=.true.; if(present(lower_tail)) lower=lower_tail
    lp=.false.; if(present(log_p)) lp=log_p
    tab=pgpbinom(probs,val_p,val_q,method,wts,lower)
    allocate(v(size(x)))
    do i=1,size(x)
      if(x(i)<tab%lower) then
        v(i)=merge(0.0_dp,1.0_dp,lower)
      else if(x(i)>tab%upper) then
        v(i)=merge(1.0_dp,0.0_dp,lower)
      else
        v(i)=tab%values(x(i)-tab%lower+1)
      end if
    end do
    if(lp) v=log(v)
  end function pgpbinom_values

  function qgpbinom_values(prob,probs,val_p,val_q,method,wts,lower_tail,log_p) result(q)
    real(dp),intent(in)::prob(:)
    real(dp),intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::wts(:)
    logical,intent(in),optional::lower_tail,log_p
    integer,allocatable::q(:)
    integer::i
    allocate(q(size(prob)))
    do i=1,size(prob)
      q(i)=qgpbinom(prob(i),probs,val_p,val_q,method,wts,lower_tail,log_p)
    end do
  end function qgpbinom_values

  subroutine rgpbinom(n,probs,val_p,val_q,draws,method,wts,generator)
    integer,intent(in)::n
    real(dp),intent(in)::probs(:)
    integer,intent(in)::val_p(:),val_q(:)
    integer,allocatable,intent(out)::draws(:)
    character(len=*),intent(in),optional::method,generator
    integer,intent(in),optional::wts(:)
    real(dp),allocatable::pexp(:)
    integer,allocatable::vp(:),vq(:)
    type(gpb_table)::tab
    character(len=32)::m,g
    integer::i,j,s
    real(dp)::u
    if(n<0) error stop "n must be nonnegative"
    m="DivideFFT"; if(present(method)) m=method
    g="Sample"; if(present(generator)) g=generator
    call expand_gpb(probs,val_p,val_q,wts,pexp,vp,vq)
    if(any(pexp<0.0_dp).or.any(pexp>1.0_dp)) error stop "probabilities must lie in [0,1]"
    allocate(draws(n))
    select case(trim(adjustl(g)))
    case("Bernoulli","bernoulli","BERNOULLI")
      do j=1,n
        s=0
        do i=1,size(pexp)
          call random_number(u)
          if(u<pexp(i)) then
            s=s+vp(i)
          else
            s=s+vq(i)
          end if
        end do
        draws(j)=s
      end do
    case("Sample","sample","SAMPLE")
      tab=dgpbinom(pexp,vp,vq,m)
      do j=1,n
        draws(j)=sample_from_pmf(tab%values,tab%lower)
      end do
    case default
      error stop "unknown generator"
    end select
  end subroutine rgpbinom

end module pb_generalized
