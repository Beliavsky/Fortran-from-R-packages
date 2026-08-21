! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_spectral
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tsa_kinds, only : dp
  use tsa_types, only : spectral_estimate, arimax_result
  use tsa_statistics, only : autocovariance
  use tsa_arimax, only : arima_fit
  use tseries_linalg, only : least_squares
  implicit none
  private

  public :: spec_pgram, spec_ar
  public :: modified_daniell_weights, spectral_kernel_df, spectral_kernel_bandwidth
  public :: tskernel_weights, daniell_kernel, modified_daniell_kernel
  public :: fejer_kernel, dirichlet_kernel

  interface spec_pgram
    module procedure spec_pgram_1d
    module procedure spec_pgram_2d
  end interface spec_pgram

  interface daniell_kernel
    module procedure daniell_kernel_scalar
    module procedure daniell_kernel_vector
  end interface daniell_kernel

  interface modified_daniell_kernel
    module procedure modified_daniell_kernel_scalar
    module procedure modified_daniell_kernel_vector
  end interface modified_daniell_kernel

  interface fejer_kernel
    module procedure fejer_kernel_real
    module procedure fejer_kernel_int
  end interface fejer_kernel

  interface dirichlet_kernel
    module procedure dirichlet_kernel_real
    module procedure dirichlet_kernel_int
  end interface dirichlet_kernel

contains

  function spec_pgram_1d(x, spans, kernel_weights, taper, pad, fast, demean, detrend, frequency, kernel_coef) result(res)
    ! Univariate numerical counterpart of stats::spec.pgram().
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: spans(:), pad
    real(dp), intent(in), optional :: kernel_weights(:), taper, frequency, kernel_coef(:)
    logical, intent(in), optional :: fast, demean, detrend
    type(spectral_estimate) :: res
    real(dp), allocatable :: y(:), work(:), pg(:), sm(:), kw(:)
    real(dp) :: tp, xfreq, u2, u4, tmid, sumt2, slope, ang, re, im
    logical :: do_fast, do_demean, do_detrend, have_kernel
    integer :: n0, n1, n, nspec, npad, k, j

    n0 = size(x)
    if (n0 < 2 .or. any(.not. ieee_is_finite(x))) then
      res%status = 1
      allocate(res%frequency(0),res%spectrum(0))
      return
    end if
    tp = 0.1_dp
    if (present(taper)) tp = taper
    if (tp < 0.0_dp .or. tp > 0.5_dp) then
      res%status = 2
      allocate(res%frequency(0),res%spectrum(0))
      return
    end if
    xfreq = 1.0_dp
    if (present(frequency)) xfreq = frequency
    if (xfreq <= 0.0_dp) then
      res%status = 3
      allocate(res%frequency(0),res%spectrum(0))
      return
    end if
    do_fast = .true.
    if (present(fast)) do_fast = fast
    do_demean = .false.
    if (present(demean)) do_demean = demean
    do_detrend = .true.
    if (present(detrend)) do_detrend = detrend

    allocate(y(n0))
    y = x
    if (do_detrend) then
      tmid = 0.5_dp*real(n0+1,dp)
      sumt2 = real(n0,dp)*(real(n0,dp)**2-1.0_dp)/12.0_dp
      slope = 0.0_dp
      do j = 1, n0
        slope = slope + y(j)*(real(j,dp)-tmid)
      end do
      slope = slope/sumt2
      y = y-sum(y)/real(n0,dp)
      do j = 1, n0
        y(j) = y(j)-slope*(real(j,dp)-tmid)
      end do
    else if (do_demean) then
      y = y-sum(y)/real(n0,dp)
    end if
    call apply_taper(y,tp)
    u2 = 1.0_dp-(5.0_dp/8.0_dp)*2.0_dp*tp
    u4 = 1.0_dp-(93.0_dp/128.0_dp)*2.0_dp*tp

    npad = 0
    if (present(pad)) npad = max(0,pad)
    n1 = n0*(1+npad)
    if (do_fast) then
      n = nextn235(n1)
    else
      n = n1
    end if
    allocate(work(n))
    work = 0.0_dp
    work(1:n0) = y

    allocate(pg(n))
    do k = 0, n-1
      re = 0.0_dp
      im = 0.0_dp
      do j = 0, n-1
        ang = 2.0_dp*acos(-1.0_dp)*real(k*j,dp)/real(n,dp)
        re = re+work(j+1)*cos(ang)
        im = im-work(j+1)*sin(ang)
      end do
      pg(k+1) = (re*re+im*im)/(real(n0,dp)*xfreq)
    end do
    pg(1) = 0.5_dp*(pg(2)+pg(n))

    have_kernel = .false.
    if (present(spans)) then
      if (size(spans) > 0) then
        call modified_daniell_weights(spans,kw,res%status)
        if (res%status /= 0) then
          allocate(res%frequency(0),res%spectrum(0))
          return
        end if
        have_kernel = .true.
      end if
    else if (present(kernel_coef)) then
      if (size(kernel_coef) > 0) then
        call tskernel_weights(kernel_coef,kw,res%status)
        if (res%status /= 0) then
          allocate(res%frequency(0),res%spectrum(0))
          return
        end if
        have_kernel = .true.
      end if
    else if (present(kernel_weights)) then
      if (size(kernel_weights) > 0) then
        if (mod(size(kernel_weights),2) == 0) then
          res%status = 4
          allocate(res%frequency(0),res%spectrum(0))
          return
        end if
        allocate(kw(size(kernel_weights)))
        kw = kernel_weights
        if (abs(sum(kw)-1.0_dp) > 1.0e-10_dp) then
          res%status = 5
          allocate(res%frequency(0),res%spectrum(0))
          return
        end if
        have_kernel = .true.
      end if
    end if

    if (have_kernel) then
      allocate(sm(n))
      call circular_smooth(pg,kw,sm)
      pg = sm
      res%degrees_freedom = spectral_kernel_df(kw)
      res%bandwidth = spectral_kernel_bandwidth(kw)
      res%method = 'Smoothed Periodogram'
    else
      res%degrees_freedom = 2.0_dp
      res%bandwidth = sqrt(1.0_dp/12.0_dp)
      res%method = 'Raw Periodogram'
    end if
    res%degrees_freedom = res%degrees_freedom/(u4/(u2*u2))
    res%degrees_freedom = res%degrees_freedom*real(n0,dp)/real(n,dp)
    res%bandwidth = res%bandwidth*xfreq/real(n,dp)

    nspec = n/2
    allocate(res%frequency(nspec),res%spectrum(nspec))
    do k = 1, nspec
      res%frequency(k) = real(k,dp)*xfreq/real(n,dp)
      res%spectrum(k) = pg(k+1)/u2
    end do
    res%n_used = n
    res%orig_n = n0
    res%taper = tp
    res%pad = npad
    res%detrend = do_detrend
    res%demean = do_demean
    res%status = 0
    res%n_series = 1
  end function spec_pgram_1d

  function spec_pgram_2d(x, spans, kernel_weights, taper, pad, fast, demean, &
      detrend, frequency, kernel_coef, taper_series) result(res)
    ! Multivariate stats::spec.pgram() numerical path.  This mirrors the
    ! mvfft/cross-periodogram calculation and R's pair ordering for squared
    ! coherency and phase: (1,2), (1,3), (2,3), (1,4), ... .
    real(dp), intent(in) :: x(:,:)
    integer, intent(in), optional :: spans(:), pad
    real(dp), intent(in), optional :: kernel_weights(:), taper, frequency, kernel_coef(:), taper_series(:)
    logical, intent(in), optional :: fast, demean, detrend
    type(spectral_estimate) :: res
    real(dp), allocatable :: y(:,:), work(:,:), kw(:)
    complex(dp), allocatable :: fft(:,:), pg(:,:,:), smc(:)
    real(dp) :: tp, xfreq, tmid, sumt2, slope, ang
    real(dp), allocatable :: tpv(:), u2v(:), u4v(:)
    real(dp) :: re, im, den
    logical :: do_fast, do_demean, do_detrend, have_kernel
    integer :: n0, n1, n, nspec, npad, nser, npair, i, j, k, t, ind

    n0 = size(x,1)
    nser = size(x,2)
    if (n0 < 2 .or. nser < 1 .or. any(.not. ieee_is_finite(x))) then
      res%status = 1
      allocate(res%frequency(0),res%spectrum(0),res%spectrum_matrix(0,0), &
        res%coherence(0,0),res%phase(0,0))
      return
    end if
    tp = 0.1_dp
    if (present(taper)) tp = taper
    allocate(tpv(nser),u2v(nser),u4v(nser))
    tpv = tp
    if (present(taper_series)) then
      if (size(taper_series) == 1) then
        tpv = taper_series(1)
      else if (size(taper_series) == nser) then
        tpv = taper_series
      else
        res%status = 2
        allocate(res%frequency(0),res%spectrum(0),res%spectrum_matrix(0,0), &
          res%coherence(0,0),res%phase(0,0))
        return
      end if
    end if
    if (any(tpv < 0.0_dp) .or. any(tpv > 0.5_dp)) then
      res%status = 2
      allocate(res%frequency(0),res%spectrum(0),res%spectrum_matrix(0,0), &
        res%coherence(0,0),res%phase(0,0))
      return
    end if
    xfreq = 1.0_dp
    if (present(frequency)) xfreq = frequency
    if (xfreq <= 0.0_dp) then
      res%status = 3
      allocate(res%frequency(0),res%spectrum(0),res%spectrum_matrix(0,0), &
        res%coherence(0,0),res%phase(0,0))
      return
    end if
    do_fast = .true.
    if (present(fast)) do_fast = fast
    do_demean = .false.
    if (present(demean)) do_demean = demean
    do_detrend = .true.
    if (present(detrend)) do_detrend = detrend

    allocate(y(n0,nser)); y=x
    if (do_detrend) then
      tmid = 0.5_dp*real(n0+1,dp)
      sumt2 = real(n0,dp)*(real(n0,dp)**2-1.0_dp)/12.0_dp
      do i=1,nser
        slope=0.0_dp
        do t=1,n0
          slope=slope+y(t,i)*(real(t,dp)-tmid)
        end do
        slope=slope/sumt2
        y(:,i)=y(:,i)-sum(y(:,i))/real(n0,dp)
        do t=1,n0
          y(t,i)=y(t,i)-slope*(real(t,dp)-tmid)
        end do
      end do
    else if (do_demean) then
      do i=1,nser
        y(:,i)=y(:,i)-sum(y(:,i))/real(n0,dp)
      end do
    end if
    do i=1,nser
      call apply_taper(y(:,i),tpv(i))
    end do
    u2v=1.0_dp-(5.0_dp/8.0_dp)*2.0_dp*tpv
    u4v=1.0_dp-(93.0_dp/128.0_dp)*2.0_dp*tpv

    npad=0
    if (present(pad)) npad=max(0,pad)
    n1=n0*(1+npad)
    if (do_fast) then
      n=nextn235(n1)
    else
      n=n1
    end if
    allocate(work(n,nser)); work=0.0_dp; work(1:n0,:)=y
    allocate(fft(n,nser)); fft=cmplx(0.0_dp,0.0_dp,kind=dp)
    do i=1,nser
      do k=0,n-1
        re=0.0_dp; im=0.0_dp
        do t=0,n-1
          ang=2.0_dp*acos(-1.0_dp)*real(k*t,dp)/real(n,dp)
          re=re+work(t+1,i)*cos(ang)
          im=im-work(t+1,i)*sin(ang)
        end do
        fft(k+1,i)=cmplx(re,im,kind=dp)
      end do
    end do
    allocate(pg(n,nser,nser))
    do i=1,nser
      do j=1,nser
        pg(:,i,j)=fft(:,i)*conjg(fft(:,j))/(real(n0,dp)*xfreq)
        pg(1,i,j)=0.5_dp*(pg(2,i,j)+pg(n,i,j))
      end do
    end do

    have_kernel=.false.
    if (present(spans)) then
      if (size(spans)>0) then
        call modified_daniell_weights(spans,kw,res%status)
        if (res%status/=0) return
        have_kernel=.true.
      end if
    else if (present(kernel_coef)) then
      if (size(kernel_coef)>0) then
        call tskernel_weights(kernel_coef,kw,res%status)
        if (res%status/=0) return
        have_kernel=.true.
      end if
    else if (present(kernel_weights)) then
      if (size(kernel_weights)>0) then
        if (mod(size(kernel_weights),2)==0) then
          res%status=4; return
        end if
        allocate(kw(size(kernel_weights))); kw=kernel_weights
        if (abs(sum(kw)-1.0_dp)>1.0e-10_dp) then
          res%status=5; return
        end if
        have_kernel=.true.
      end if
    end if
    if (have_kernel) then
      allocate(smc(n))
      do i=1,nser
        do j=1,nser
          call circular_smooth_complex(pg(:,i,j),kw,smc)
          pg(:,i,j)=smc
        end do
      end do
      res%degrees_freedom=spectral_kernel_df(kw)
      res%bandwidth=spectral_kernel_bandwidth(kw)
      res%method='Smoothed Periodogram'
    else
      res%degrees_freedom=2.0_dp
      res%bandwidth=sqrt(1.0_dp/12.0_dp)
      res%method='Raw Periodogram'
    end if
    allocate(res%degrees_freedom_series(nser))
    res%degrees_freedom_series=res%degrees_freedom/(u4v/(u2v*u2v))
    res%degrees_freedom_series=res%degrees_freedom_series*real(n0,dp)/real(n,dp)
    res%degrees_freedom=res%degrees_freedom_series(1)
    res%bandwidth=res%bandwidth*xfreq/real(n,dp)

    nspec=n/2
    npair=nser*(nser-1)/2
    allocate(res%frequency(nspec),res%spectrum_matrix(nspec,nser))
    allocate(res%coherence(nspec,npair),res%phase(nspec,npair))
    if (nser==1) allocate(res%spectrum(nspec))
    do k=1,nspec
      res%frequency(k)=real(k,dp)*xfreq/real(n,dp)
      do i=1,nser
        res%spectrum_matrix(k,i)=real(pg(k+1,i,i),dp)/u2v(i)
      end do
    end do
    if (nser==1) res%spectrum=res%spectrum_matrix(:,1)
    if (npair>0) then
      do j=2,nser
        do i=1,j-1
          ind=i+(j-1)*(j-2)/2
          do k=1,nspec
            den=res%spectrum_matrix(k,i)*res%spectrum_matrix(k,j)
            if (den>0.0_dp) then
              ! R computes coherency before taper correction.  Since den uses
              ! corrected autospectra, undo the two series-specific u2 factors.
              res%coherence(k,ind)=abs(pg(k+1,i,j))**2/(u2v(i)*u2v(j)*den)
            else
              res%coherence(k,ind)=0.0_dp
            end if
            res%phase(k,ind)=atan2(aimag(pg(k+1,i,j)),real(pg(k+1,i,j),dp))
          end do
        end do
      end do
    end if
    res%n_used=n
    res%orig_n=n0
    res%n_series=nser
    res%taper=tpv(1)
    allocate(res%taper_series(nser)); res%taper_series=tpv
    res%pad=npad
    res%detrend=do_detrend
    res%demean=do_demean
    res%status=0
  end function spec_pgram_2d

  function spec_ar(x, n_freq, order, order_max, demean, frequency, method, var_method, intercept) result(res)
    ! Numerical counterpart of stats::spec.ar(), including all univariate ar()
    ! fitting methods used by R: Yule-Walker, Burg/Burg2, OLS, and MLE.
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: n_freq, order, order_max, var_method
    logical, intent(in), optional :: demean, intercept
    real(dp), intent(in), optional :: frequency
    character(len=*), intent(in), optional :: method
    type(spectral_estimate) :: res
    real(dp), allocatable :: phi(:)
    real(dp) :: xfreq, varp, cs, sn, ang
    logical :: dm, inc
    integer :: n, nf, maxord, ord, j, k, vm, status
    character(len=24) :: meth, label

    n = size(x)
    if (n < 2 .or. any(.not. ieee_is_finite(x))) then
      res%status = 1
      allocate(res%frequency(0),res%spectrum(0))
      return
    end if
    xfreq = 1.0_dp
    if (present(frequency)) xfreq = frequency
    if (xfreq <= 0.0_dp) then
      res%status=4; allocate(res%frequency(0),res%spectrum(0)); return
    end if
    dm = .true.; if (present(demean)) dm = demean
    inc = dm; if (present(intercept)) inc = intercept
    vm = 1; if (present(var_method)) vm = var_method
    meth = 'yule-walker'; if (present(method)) meth = lowercase(trim(method))

    if (present(order)) then
      maxord=max(0,order)
    else if (present(order_max)) then
      maxord=max(0,min(order_max,n-1))
    else if (index(meth,'mle')==1) then
      maxord=min(n-1,min(12,int(floor(10.0_dp*log10(real(n,dp))))))
    else
      maxord=min(n-1,int(floor(10.0_dp*log10(real(n,dp)))))
    end if
    if (maxord >= n) then
      res%status=2; allocate(res%frequency(0),res%spectrum(0)); return
    end if

    if (index(meth,'burg')==1) then
      call fit_ar_burg(x,maxord,.not.present(order),dm,vm,phi,varp,ord,status)
      label=merge('Burg2 AR spectrum       ','Burg AR spectrum        ',vm==2)
    else if (index(meth,'ols')==1) then
      call fit_ar_ols_spectrum(x,maxord,.not.present(order),dm,inc,phi,varp,ord,status)
      label='OLS AR spectrum'
    else if (index(meth,'mle')==1) then
      call fit_ar_mle_spectrum(x,maxord,.not.present(order),dm,phi,varp,ord,status)
      label='MLE AR spectrum'
    else if (index(meth,'yw')==1 .or. index(meth,'yule')==1) then
      call fit_ar_yw_spectrum(x,maxord,.not.present(order),dm,phi,varp,ord,status)
      label='Yule-Walker AR spectrum'
    else
      status=5; allocate(phi(0)); varp=0.0_dp; ord=0
      label='AR spectrum'
    end if
    if (status /= 0) then
      res%status=status; allocate(res%frequency(0),res%spectrum(0)); return
    end if

    nf=500; if (present(n_freq)) nf=max(2,n_freq)
    allocate(res%frequency(nf),res%spectrum(nf))
    do k=1,nf
      res%frequency(k)=0.5_dp*xfreq*real(k-1,dp)/real(nf-1,dp)
      cs=0.0_dp; sn=0.0_dp
      do j=1,ord
        ang=2.0_dp*acos(-1.0_dp)*(res%frequency(k)/xfreq)*real(j,dp)
        cs=cs+phi(j)*cos(ang); sn=sn+phi(j)*sin(ang)
      end do
      res%spectrum(k)=varp/(xfreq*((1.0_dp-cs)**2+sn*sn))
    end do
    res%degrees_freedom=0.0_dp; res%bandwidth=0.0_dp
    res%n_used=n; res%orig_n=n; res%order=ord; res%method=label
    res%demean=dm; res%detrend=.false.; res%taper=0.0_dp; res%pad=0
    res%status=0
  end function spec_ar

  subroutine fit_ar_yw_spectrum(x,maxord,select_aic,demean,phi,varp,ord,status)
    real(dp),intent(in)::x(:); integer,intent(in)::maxord
    logical,intent(in)::select_aic,demean
    real(dp),allocatable,intent(out)::phi(:); real(dp),intent(out)::varp
    integer,intent(out)::ord,status
    real(dp),allocatable::y(:),acv(:),cur(:),nxt(:),vars(:),aic(:),coefmat(:,:)
    real(dp)::v,refl,best; integer::n,m,j,bestord
    n=size(x); allocate(y(n)); y=x
    if(demean)y=y-sum(y)/real(n,dp)
    allocate(acv(0:maxord)); call autocovariance(y,maxord,acv,demean=.false.)
    if(acv(0)<=tiny(1.0_dp))then; status=3; allocate(phi(0)); varp=0; ord=0; return; end if
    allocate(vars(0:maxord),aic(0:maxord),coefmat(maxord,maxord)); coefmat=0
    vars(0)=acv(0); aic(0)=real(n,dp)*log(vars(0))+merge(2.0_dp,0.0_dp,demean)
    allocate(cur(0)); v=acv(0)
    do m=1,maxord
      refl=acv(m); do j=1,m-1; refl=refl-cur(j)*acv(m-j); end do; refl=refl/v
      allocate(nxt(m)); nxt(m)=refl
      do j=1,m-1; nxt(j)=cur(j)-refl*cur(m-j); end do
      v=v*(1-refl*refl); vars(m)=v; coefmat(m,1:m)=nxt
      aic(m)=real(n,dp)*log(max(v,tiny(1.0_dp)))+2.0_dp*real(m,dp)+merge(2.0_dp,0.0_dp,demean)
      call move_alloc(nxt,cur)
    end do
    if(select_aic)then; best=minval(aic); bestord=minloc(aic,dim=1)-1; ord=bestord; else; ord=maxord; end if
    allocate(phi(ord)); if(ord>0)phi=coefmat(ord,1:ord)
    varp=vars(ord)*real(n,dp)/real(n-(ord+1),dp); status=0
  end subroutine fit_ar_yw_spectrum

  subroutine fit_ar_burg(x,maxord,select_aic,demean,var_method,phi,varp,ord,status)
    ! Direct computational translation of R stats/src/burg.c for univariate AR.
    real(dp),intent(in)::x(:)
    integer,intent(in)::maxord,var_method
    logical,intent(in)::select_aic,demean
    real(dp),allocatable,intent(out)::phi(:)
    real(dp),intent(out)::varp
    integer,intent(out)::ord,status
    real(dp),allocatable::y(:),u(:),v(:),u0(:),v1(:),v2(:),aic(:),cm(:,:)
    real(dp)::sm,d,phii
    integer::n,p,j,t

    n=size(x)
    allocate(y(n)); y=x
    if(demean)y=y-sum(y)/real(n,dp)
    if(sum(y*y)<=tiny(1.0_dp))then
      status=3; allocate(phi(0)); varp=0.0_dp; ord=0; return
    end if

    allocate(u(n),v(n),u0(n),v1(0:maxord),v2(0:maxord),aic(0:maxord),cm(maxord,maxord))
    cm=0.0_dp
    sm=0.0_dp
    do t=1,n
      u(t)=y(n+1-t)
      v(t)=u(t)
      sm=sm+y(t)*y(t)
    end do
    v1(0)=sm/real(n,dp)
    v2(0)=v1(0)
    aic(0)=real(n,dp)*log(v1(0))+merge(2.0_dp,0.0_dp,demean)

    do p=1,maxord
      sm=0.0_dp; d=0.0_dp
      do t=p+1,n
        sm=sm+v(t)*u(t-1)
        d=d+v(t)*v(t)+u(t-1)*u(t-1)
      end do
      if(d<=tiny(1.0_dp))then
        status=3; allocate(phi(0)); varp=0.0_dp; ord=0; return
      end if
      phii=2.0_dp*sm/d
      cm(p,p)=phii
      if(p>1)then
        do j=1,p-1
          cm(p,j)=cm(p-1,j)-phii*cm(p-1,p-j)
        end do
      end if

      u0=u
      do t=p+1,n
        u(t)=u0(t-1)-phii*v(t)
        v(t)=v(t)-phii*u0(t-1)
      end do
      v1(p)=v1(p-1)*(1.0_dp-phii*phii)
      d=0.0_dp
      do t=p+1,n
        d=d+v(t)*v(t)+u(t)*u(t)
      end do
      v2(p)=d/(2.0_dp*real(n-p,dp))
      aic(p)=real(n,dp)*log(max(merge(v1(p),v2(p),var_method==1),tiny(1.0_dp)))+ &
        2.0_dp*real(p,dp)+merge(2.0_dp,0.0_dp,demean)
    end do

    if(select_aic)then
      ord=minloc(aic,dim=1)-1
    else
      ord=maxord
    end if
    allocate(phi(ord))
    if(ord>0)phi=cm(ord,1:ord)
    varp=merge(v1(ord),v2(ord),var_method==1)
    status=0
  end subroutine fit_ar_burg

  subroutine fit_ar_ols_spectrum(x,maxord,select_aic,demean,intercept,phi,varp,ord,status)
    real(dp),intent(in)::x(:); integer,intent(in)::maxord
    logical,intent(in)::select_aic,demean,intercept
    real(dp),allocatable,intent(out)::phi(:); real(dp),intent(out)::varp
    integer,intent(out)::ord,status
    real(dp),allocatable::y(:),xx(:,:),yy(:),beta(:),rr(:),vars(:),aic(:),cm(:,:)
    integer::n,m,j,pcol,st; real(dp)::best
    n=size(x); allocate(y(n)); y=x; if(demean)y=y-sum(y)/real(n,dp)
    allocate(vars(0:maxord),aic(0:maxord),cm(maxord,maxord)); cm=0; aic=huge(1.0_dp)
    do m=0,maxord
      pcol=m+merge(1,0,intercept); allocate(xx(n-m,pcol),yy(n-m),beta(pcol),rr(n-m))
      yy=y(m+1:n); if(pcol>0)xx=0
      if(intercept)xx(:,1)=1.0_dp
      do j=1,m
        xx(:,merge(1,0,intercept)+j)=y(m+1-j:n-j)
      end do
      if(pcol>0)then
        call least_squares(xx,yy,beta,residuals=rr,status=st)
      else
        rr=yy; st=0
      end if
      if(st/=0)then; deallocate(xx,yy,beta,rr); exit; end if
      vars(m)=sum(rr*rr)/real(n-m,dp)
      aic(m)=real(n,dp)*log(max(vars(m),tiny(1.0_dp)))+2.0_dp*real(m+merge(1,0,intercept),dp)
      if(m>0)cm(m,1:m)=beta(merge(1,0,intercept)+1:pcol)
      deallocate(xx,yy,beta,rr)
    end do
    if(select_aic)then; best=minval(aic); ord=minloc(aic,dim=1)-1; else; ord=maxord; end if
    allocate(phi(ord)); if(ord>0)phi=cm(ord,1:ord); varp=vars(ord); status=0
  end subroutine fit_ar_ols_spectrum

  subroutine fit_ar_mle_spectrum(x,maxord,select_aic,demean,phi,varp,ord,status)
    real(dp),intent(in)::x(:); integer,intent(in)::maxord; logical,intent(in)::select_aic,demean
    real(dp),allocatable,intent(out)::phi(:); real(dp),intent(out)::varp
    integer,intent(out)::ord,status
    type(arimax_result)::fit
    real(dp),allocatable::aic(:),vars(:),cm(:,:)
    integer::m; real(dp)::xm,var0
    allocate(aic(0:maxord),vars(0:maxord),cm(maxord,maxord)); cm=0
    xm=merge(sum(x)/real(size(x),dp),0.0_dp,demean); var0=sum((x-xm)**2)/real(size(x),dp)
    aic(0)=real(size(x),dp)*log(max(var0,tiny(1.0_dp)))+2.0_dp*real(1+merge(1,0,demean),dp)+ &
      real(size(x),dp)*(1.0_dp+log(2.0_dp*acos(-1.0_dp)))
    vars(0)=var0
    do m=1,maxord
      fit=arima_fit(x,m,0,0,include_mean=demean,method='ML')
      if(fit%status/=0)then; aic(m)=huge(1.0_dp); vars(m)=huge(1.0_dp); cycle; end if
      aic(m)=fit%aic; vars(m)=fit%sigma2; cm(m,1:m)=fit%ar
    end do
    if(select_aic)then; ord=minloc(aic,dim=1)-1; else; ord=maxord; end if
    if(.not.ieee_is_finite(aic(ord)))then; status=6; allocate(phi(0)); varp=0; return; end if
    allocate(phi(ord)); if(ord>0)phi=cm(ord,1:ord); varp=vars(ord); status=0
  end subroutine fit_ar_mle_spectrum

  pure function lowercase(text) result(out)
    character(len=*),intent(in)::text; character(len=len(text))::out; integer::i,c
    out=text
    do i=1,len(text); c=iachar(out(i:i)); if(c>=iachar('A').and.c<=iachar('Z'))out(i:i)=achar(c+32); end do
  end function lowercase

  subroutine tskernel_weights(coef,weights,status)
    ! Expand R tskernel compact coefficients coef(0:m) into symmetric weights.
    real(dp),intent(in)::coef(:); real(dp),allocatable,intent(out)::weights(:); integer,intent(out)::status
    integer::m,j
    status=0
    if(size(coef)<1)then; allocate(weights(0)); status=1; return; end if
    if(abs(coef(1)+2.0_dp*sum(coef(2:))-1.0_dp)>1.0e-10_dp)then
      allocate(weights(0)); status=2; return
    end if
    m=size(coef)-1; allocate(weights(2*m+1)); weights(m+1)=coef(1)
    do j=1,m; weights(m+1-j)=coef(j+1); weights(m+1+j)=coef(j+1); end do
  end subroutine tskernel_weights

  function daniell_kernel_scalar(m) result(coef)
    integer,intent(in)::m; real(dp),allocatable::coef(:)
    if(m<0)then; allocate(coef(0)); return; end if
    allocate(coef(m+1)); coef=1.0_dp/real(2*m+1,dp)
  end function daniell_kernel_scalar

  function modified_daniell_kernel_scalar(m) result(coef)
    integer,intent(in)::m; real(dp),allocatable::coef(:)
    if(m<1)then; allocate(coef(0)); return; end if
    allocate(coef(m+1)); coef=1.0_dp/real(2*m,dp); coef(m+1)=1.0_dp/real(4*m,dp)
  end function modified_daniell_kernel_scalar

  function daniell_kernel_vector(m) result(coef)
    integer,intent(in)::m(:)
    real(dp),allocatable::coef(:),full(:),one_coef(:),one(:),tmp(:)
    integer::i,st,half
    if(size(m)<1.or.any(m<0))then; allocate(coef(0)); return; end if
    one_coef=daniell_kernel_scalar(m(1))
    call tskernel_weights(one_coef,full,st)
    if(st/=0)then; allocate(coef(0)); return; end if
    do i=2,size(m)
      one_coef=daniell_kernel_scalar(m(i))
      call tskernel_weights(one_coef,one,st)
      if(st/=0)then; allocate(coef(0)); return; end if
      call linear_convolution(full,one,tmp)
      call move_alloc(tmp,full)
    end do
    half=(size(full)-1)/2
    allocate(coef(half+1)); coef=full(half+1:size(full))
  end function daniell_kernel_vector

  function modified_daniell_kernel_vector(m) result(coef)
    integer,intent(in)::m(:)
    real(dp),allocatable::coef(:),full(:),one_coef(:),one(:),tmp(:)
    integer::i,st,half
    if(size(m)<1.or.any(m<1))then; allocate(coef(0)); return; end if
    one_coef=modified_daniell_kernel_scalar(m(1))
    call tskernel_weights(one_coef,full,st)
    if(st/=0)then; allocate(coef(0)); return; end if
    do i=2,size(m)
      one_coef=modified_daniell_kernel_scalar(m(i))
      call tskernel_weights(one_coef,one,st)
      if(st/=0)then; allocate(coef(0)); return; end if
      call linear_convolution(full,one,tmp)
      call move_alloc(tmp,full)
    end do
    half=(size(full)-1)/2
    allocate(coef(half+1)); coef=full(half+1:size(full))
  end function modified_daniell_kernel_vector

  function fejer_kernel_real(m,r) result(coef)
    integer,intent(in)::m
    real(dp),intent(in)::r
    real(dp),allocatable::coef(:)
    real(dp)::wj,den
    integer::j
    if(m<1.or.r<1.0_dp)then; allocate(coef(0)); return; end if
    allocate(coef(m+1)); coef(1)=r
    do j=1,m
      wj=2.0_dp*acos(-1.0_dp)*real(j,dp)/real(2*m+1,dp)
      coef(j+1)=sin(r*wj/2.0_dp)**2/(sin(wj/2.0_dp)**2*r)
    end do
    den=coef(1)+2.0_dp*sum(coef(2:)); coef=coef/den
  end function fejer_kernel_real

  function fejer_kernel_int(m,r) result(coef)
    integer,intent(in)::m,r
    real(dp),allocatable::coef(:)
    coef=fejer_kernel_real(m,real(r,dp))
  end function fejer_kernel_int

  function dirichlet_kernel_real(m,r) result(coef)
    integer,intent(in)::m
    real(dp),intent(in)::r
    real(dp),allocatable::coef(:)
    real(dp)::wj,den
    integer::j
    if(m<1.or.r<0.0_dp)then; allocate(coef(0)); return; end if
    allocate(coef(m+1)); coef(1)=2.0_dp*r+1.0_dp
    do j=1,m
      wj=2.0_dp*acos(-1.0_dp)*real(j,dp)/real(2*m+1,dp)
      coef(j+1)=sin((r+0.5_dp)*wj)/sin(wj/2.0_dp)
    end do
    den=coef(1)+2.0_dp*sum(coef(2:)); coef=coef/den
  end function dirichlet_kernel_real

  function dirichlet_kernel_int(m,r) result(coef)
    integer,intent(in)::m,r
    real(dp),allocatable::coef(:)
    coef=dirichlet_kernel_real(m,real(r,dp))
  end function dirichlet_kernel_int

  subroutine modified_daniell_weights(spans, weights, status)
    integer, intent(in) :: spans(:)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), allocatable :: one(:), tmp(:)
    integer :: i, m

    status = 0
    if (size(spans) == 0) then
      allocate(weights(1)); weights=1.0_dp
      return
    end if
    do i = 1, size(spans)
      m = spans(i)/2
      if (m < 1) then
        status = 1
        allocate(weights(0))
        return
      end if
      allocate(one(2*m+1))
      one = 1.0_dp/real(2*m,dp)
      one(1) = 1.0_dp/real(4*m,dp)
      one(2*m+1) = one(1)
      if (i == 1) then
        weights = one
      else
        call linear_convolution(weights,one,tmp)
        call move_alloc(tmp,weights)
      end if
      deallocate(one)
    end do
  end subroutine modified_daniell_weights

  pure real(dp) function spectral_kernel_df(weights) result(df)
    real(dp), intent(in) :: weights(:)
    df = 2.0_dp/sum(weights*weights)
  end function spectral_kernel_df

  pure real(dp) function spectral_kernel_bandwidth(weights) result(bw)
    real(dp), intent(in) :: weights(:)
    integer :: m, i, lag
    bw = 0.0_dp
    m = (size(weights)-1)/2
    do i = 1, size(weights)
      lag = i-m-1
      bw = bw+(1.0_dp/12.0_dp+real(lag*lag,dp))*weights(i)
    end do
    bw = sqrt(max(0.0_dp,bw))
  end function spectral_kernel_bandwidth

  subroutine apply_taper(x,p)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: p
    real(dp), allocatable :: w(:)
    integer :: m, i, n
    n = size(x)
    m = int(floor(real(n,dp)*p))
    if (m == 0) return
    allocate(w(m))
    do i = 1, m
      w(i)=0.5_dp*(1.0_dp-cos(acos(-1.0_dp)*real(2*i-1,dp)/real(2*m,dp)))
    end do
    do i = 1, m
      x(i)=w(i)*x(i)
      x(n-i+1)=w(i)*x(n-i+1)
    end do
  end subroutine apply_taper

  subroutine circular_smooth(x,w,y)
    real(dp), intent(in) :: x(:),w(:)
    real(dp), intent(out) :: y(:)
    integer :: n, m, i, j, lag, idx
    n = size(x)
    m = (size(w)-1)/2
    y = 0.0_dp
    do i = 1, n
      do j = 1, size(w)
        lag = j-m-1
        idx = modulo(i-1-lag,n)+1
        y(i)=y(i)+w(j)*x(idx)
      end do
    end do
  end subroutine circular_smooth

  subroutine circular_smooth_complex(x,w,y)
    complex(dp), intent(in) :: x(:)
    real(dp), intent(in) :: w(:)
    complex(dp), intent(out) :: y(:)
    integer :: n, m, i, j, lag, idx
    n=size(x)
    m=(size(w)-1)/2
    y=cmplx(0.0_dp,0.0_dp,kind=dp)
    do i=1,n
      do j=1,size(w)
        lag=j-m-1
        idx=modulo(i-1+lag,n)+1
        y(i)=y(i)+w(j)*x(idx)
      end do
    end do
  end subroutine circular_smooth_complex

  subroutine linear_convolution(a,b,c)
    real(dp), intent(in) :: a(:),b(:)
    real(dp), allocatable, intent(out) :: c(:)
    integer :: i,j
    allocate(c(size(a)+size(b)-1))
    c = 0.0_dp
    do i = 1, size(a)
      do j = 1, size(b)
        c(i+j-1)=c(i+j-1)+a(i)*b(j)
      end do
    end do
  end subroutine linear_convolution

  integer function nextn235(n0) result(n)
    integer, intent(in) :: n0
    integer :: q
    logical :: good
    n = max(1,n0)
    do
      q = n
      do while (mod(q,2)==0); q=q/2; end do
      do while (mod(q,3)==0); q=q/3; end do
      do while (mod(q,5)==0); q=q/5; end do
      good = q == 1
      if (good) exit
      n = n+1
    end do
  end function nextn235

end module tsa_spectral
