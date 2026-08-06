! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_transform_1d
  use waveslim_kinds, only : dp, sqrt2
  use waveslim_status, only : status_type, clear_status, set_status, &
    waveslim_invalid_input, waveslim_invalid_level
  use waveslim_types, only : wavelet_transform, mra_result, &
    complex_wavelet_transform, real_vector
  use waveslim_filters, only : wave_filter, hilbert_filter
  use waveslim_math, only : circular_shift, lower_string
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: dwt, idwt, dwt_nondyadic, modwt, imodwt, mra
  public :: dwt_step, idwt_step, modwt_step, imodwt_step
  public :: brick_wall, phase_shift, up_sample
  public :: dwt_hilbert, modwt_hilbert, idwt_hilbert, imodwt_hilbert
contains
  subroutine dwt_step(x, h, g, w, v)
    real(dp), intent(in) :: x(:), h(:), g(:)
    real(dp), allocatable, intent(out) :: w(:), v(:)
    integer :: m, l, t, n, u
    m = size(x)
    l = size(h)
    allocate(w(m/2), v(m/2))
    do t = 0, m/2 - 1
      u = 2*t + 2
      w(t+1) = h(1)*x(u)
      v(t+1) = g(1)*x(u)
      do n = 2, l
        u = u - 1
        if (u < 1) u = m
        w(t+1) = w(t+1) + h(n)*x(u)
        v(t+1) = v(t+1) + g(n)*x(u)
      end do
    end do
  end subroutine dwt_step

  subroutine idwt_step(w, v, h, g, x)
    real(dp), intent(in) :: w(:), v(:), h(:), g(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer :: m, l, t, u, k, even_i, odd_i
    m = size(w)
    l = size(h)
    allocate(x(2*m))
    x = 0.0_dp
    do t = 0, m-1
      u = t + 1
      do k = 0, l/2-1
        even_i = 2*t + 1
        odd_i = 2*t + 2
        even_i = modulo(even_i-1,2*m)+1
        odd_i = modulo(odd_i-1,2*m)+1
        x(even_i) = x(even_i) + h(2*k+2)*w(u) + g(2*k+2)*v(u)
        x(odd_i) = x(odd_i) + h(2*k+1)*w(u) + g(2*k+1)*v(u)
        u = modulo(u,m)+1
      end do
    end do
  end subroutine idwt_step

  subroutine modwt_step(x, level, h, g, w, v)
    real(dp), intent(in) :: x(:), h(:), g(:)
    integer, intent(in) :: level
    real(dp), allocatable, intent(out) :: w(:), v(:)
    integer :: n, l, t, k, j, stride
    n=size(x)
    l=size(h)
    stride=2**(level-1)
    allocate(w(n),v(n))
    do t=1,n
      k=t
      w(t)=h(1)*x(k)
      v(t)=g(1)*x(k)
      do j=2,l
        k=k-stride
        do while(k<1)
        k=k+n
        end do
        w(t)=w(t)+h(j)*x(k)
        v(t)=v(t)+g(j)*x(k)
      end do
    end do
  end subroutine modwt_step

  subroutine imodwt_step(w, v, level, h, g, x)
    real(dp),intent(in)::w(:),v(:),h(:),g(:)
    integer,intent(in)::level
    real(dp),allocatable,intent(out)::x(:)
    integer::n,l,t,k,j,stride
    n=size(w)
    l=size(h)
    stride=2**(level-1)
    allocate(x(n))
    do t=1,n
      k=t
      x(t)=h(1)*w(k)+g(1)*v(k)
      do j=2,l
        k=k+stride
        do while(k>n)
        k=k-n
        end do
        x(t)=x(t)+h(j)*w(k)+g(j)*v(k)
      end do
    end do
  end subroutine imodwt_step

  function dwt(x, wf, n_levels, boundary) result(result)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in), optional :: wf, boundary
    integer, intent(in), optional :: n_levels
    type(wavelet_transform) :: result
    character(len=16) :: wname, bname
    real(dp), allocatable :: work(:), w(:), v(:)
    integer :: j, levels, n
    type(status_type) :: st
    wname='la8'
    if(present(wf))wname=trim(wf)
    bname='periodic'
    if(present(boundary))bname=lower_string(trim(boundary))
    levels=4
    if(present(n_levels))levels=n_levels
    result%wavelet=wname
    result%boundary=bname
    result%method='dwt'
    result%original_length=size(x)
    if(bname=='reflection')then
      allocate(work(2*size(x)))
      work(1:size(x))=x
      work(size(x)+1:)=x(size(x):1:-1)
    else if(bname=='periodic')then
      work=x
    else
      call set_status(result%status,waveslim_invalid_input,'invalid boundary')
      return
    end if
    n=size(work)
    if(levels<1 .or. 2**levels>n .or. mod(n,2**levels)/=0)then
      call set_status(result%status,waveslim_invalid_level, &
        'DWT length must be divisible by 2**levels')
      return
    end if
    block
      use waveslim_types, only : wavelet_filter_type
      type(wavelet_filter_type) :: f
      f=wave_filter(wname,st)
      if(.not.st%ok())then
      result%status=st
      return
      end if
      allocate(result%detail(levels))
      do j=1,levels
        call dwt_step(work,f%hpf,f%lpf,w,v)
        result%detail(j)%values=w
        call move_alloc(v,work)
      end do
      result%smooth=work
    end block
    call clear_status(result%status)
  end function dwt

  function dwt_nondyadic(x, wf) result(result)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in),optional::wf
    type(wavelet_transform)::result
    real(dp),allocatable::xx(:)
    integer::n,j,m
    m=size(x)
    n=1
    do while(n<m)
    n=2*n
    end do
    allocate(xx(n))
    xx=0.0_dp
    xx(1:m)=x
    if(present(wf))then
      result=dwt(xx,wf,int(log(real(n,dp))/log(2.0_dp)))
    else
      result=dwt(xx,'la8',int(log(real(n,dp))/log(2.0_dp)))
    end if
    result%original_length=m
    do j=1,result%levels()
      if(m/2**j>0) result%detail(j)%values=result%detail(j)%values(1:m/2**j)
    end do
  end function dwt_nondyadic

  function idwt(y) result(x)
    type(wavelet_transform),intent(in)::y
    real(dp),allocatable::x(:)
    real(dp),allocatable::work(:),tmp(:)
    integer::j,nout
    block
      use waveslim_types, only : wavelet_filter_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(y%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0))
      return
      end if
      work=y%smooth
      do j=y%levels(),1,-1
        call idwt_step(y%detail(j)%values,work,f%hpf,f%lpf,tmp)
        call move_alloc(tmp,work)
      end do
    end block
    nout=size(work)
    if(trim(y%boundary)=='reflection')nout=y%original_length
    if(y%original_length>0 .and. y%original_length<nout)nout=y%original_length
    x=work(1:nout)
  end function idwt

  function modwt(x,wf,n_levels,boundary) result(result)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in),optional::wf,boundary
    integer,intent(in),optional::n_levels
    type(wavelet_transform)::result
    character(len=16)::wname,bname
    real(dp),allocatable::work(:),w(:),v(:),h(:),g(:)
    integer::j,levels
    type(status_type)::st
    wname='la8'
    if(present(wf))wname=trim(wf)
    bname='periodic'
    if(present(boundary))bname=lower_string(trim(boundary))
    levels=4
    if(present(n_levels))levels=n_levels
    result%wavelet=wname
    result%boundary=bname
    result%method='modwt'
    result%original_length=size(x)
    if(bname=='reflection')then
      allocate(work(2*size(x)))
      work(1:size(x))=x
      work(size(x)+1:)=x(size(x):1:-1)
    else if(bname=='periodic')then
    work=x
    else
    call set_status(result%status,waveslim_invalid_input,'invalid boundary')
    return
    end if
    if(levels<1 .or. 2**levels>size(work))then
      call set_status(result%status,waveslim_invalid_level,'MODWT exceeds sample size')
      return
    end if
    block
      use waveslim_types,only:wavelet_filter_type
      type(wavelet_filter_type)::f
      f=wave_filter(wname,st)
      if(.not.st%ok())then
      result%status=st
      return
      end if
      h=f%hpf/sqrt2
      g=f%lpf/sqrt2
      allocate(result%detail(levels))
      do j=1,levels
        call modwt_step(work,j,h,g,w,v)
        result%detail(j)%values=w
        call move_alloc(v,work)
      end do
    end block
    result%smooth=work
    call clear_status(result%status)
  end function modwt

  function imodwt(y) result(x)
    type(wavelet_transform),intent(in)::y
    real(dp),allocatable::x(:)
    real(dp),allocatable::work(:),tmp(:),h(:),g(:)
    integer::j,nout
    block
      use waveslim_types,only:wavelet_filter_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(y%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0))
      return
      end if
      h=f%hpf/sqrt2
      g=f%lpf/sqrt2
      work=y%smooth
      do j=y%levels(),1,-1
        call imodwt_step(y%detail(j)%values,work,j,h,g,tmp)
        call move_alloc(tmp,work)
      end do
    end block
    nout=size(work)
    if(trim(y%boundary)=='reflection')nout=y%original_length
    x=work(1:nout)
  end function imodwt

  function mra(x,wf,j_levels,method,boundary) result(ans)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in),optional::wf,method,boundary
    integer,intent(in),optional::j_levels
    type(mra_result)::ans
    type(wavelet_transform)::wt,z
    character(len=16)::wname,mname,bname
    integer::j,levels
    wname='la8'
    if(present(wf))wname=trim(wf)
    mname='modwt'
    if(present(method))mname=lower_string(trim(method))
    bname='periodic'
    if(present(boundary))bname=lower_string(trim(boundary))
    levels=4
    if(present(j_levels))levels=j_levels
    if(mname=='dwt')then
    wt=dwt(x,wname,levels,bname)
    else
    wt=modwt(x,wname,levels,bname)
    end if
    if(.not.wt%status%ok())then
    ans%status=wt%status
    return
    end if
    allocate(ans%detail(levels))
    ans%wavelet=wname
    ans%method=mname
    ans%boundary=bname
    z=zero_like(wt)
    z%smooth=wt%smooth
    if(mname=='dwt')then
    ans%smooth=idwt(z)
    else
    ans%smooth=imodwt(z)
    end if
    do j=1,levels
      z=zero_like(wt)
      z%detail(j)%values=wt%detail(j)%values
      if(mname=='dwt')then
      ans%detail(j)%values=idwt(z)
      else
      ans%detail(j)%values=imodwt(z)
      end if
    end do
    call clear_status(ans%status)
  end function mra

  function zero_like(wt) result(z)
    type(wavelet_transform),intent(in)::wt
    type(wavelet_transform)::z
    integer::j
    z%wavelet=wt%wavelet
    z%boundary=wt%boundary
    z%method=wt%method
    z%original_length=wt%original_length
    allocate(z%detail(wt%levels()))
    do j=1,wt%levels()
    allocate(z%detail(j)%values(size(wt%detail(j)%values)))
    z%detail(j)%values=0.0_dp
    end do
    allocate(z%smooth(size(wt%smooth)))
    z%smooth=0.0_dp
    call clear_status(z%status)
  end function zero_like

  subroutine brick_wall(wt,wf,method)
    type(wavelet_transform),intent(inout)::wt
    character(len=*),intent(in),optional::wf,method
    character(len=16)::wname,mname
    integer::j,n,m
    real(dp)::nan
    block
      use waveslim_types,only:wavelet_filter_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      wname=wt%wavelet
      if(present(wf))wname=trim(wf)
      mname=wt%method
      if(present(method))mname=lower_string(trim(method))
      f=wave_filter(wname,st)
      if(.not.st%ok())return
      m=f%length()
    end block
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    n=0
    do j=1,wt%levels()
      if(mname=='dwt')then
      n=ceiling(real(m-2,dp)*(1.0_dp-1.0_dp/real(2**j,dp)))
      else
      n=(2**j-1)*(m-1)
      end if
      n=min(n,size(wt%detail(j)%values))
      if(n>0)wt%detail(j)%values(1:n)=nan
    end do
    n=min(n,size(wt%smooth))
    if(n>0)wt%smooth(1:n)=nan
  end subroutine brick_wall

  subroutine phase_shift(wt,wf,inverse)
    type(wavelet_transform),intent(inout)::wt
    character(len=*),intent(in),optional::wf
    logical,intent(in),optional::inverse
    character(len=16)::wname
    logical::inv
    real(dp)::cg,ch
    integer::j,ph,n
    block
      use waveslim_types,only:wavelet_filter_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      wname=wt%wavelet
      if(present(wf))wname=trim(wf)
      inv=.false.
      if(present(inverse))inv=inverse
      f=wave_filter(wname,st)
      if(.not.st%ok())return
      cg=center_energy(f%lpf)
      ch=center_energy(f%hpf)
    end block
    do j=1,wt%levels()
      ph=nint(real(2**(j-1),dp)*(cg+ch)-cg)
      n=size(wt%detail(j)%values)
      if(inv)ph=-ph
      wt%detail(j)%values=circular_shift(wt%detail(j)%values,ph)
    end do
    ph=nint(real(2**wt%levels()-1,dp)*cg)
    if(inv)ph=-ph
    wt%smooth=circular_shift(wt%smooth,ph)
  end subroutine phase_shift

  pure function center_energy(g) result(c)
    real(dp),intent(in)::g(:)
    real(dp)::c
    integer::i
    c=0.0_dp
    do i=1,size(g)
    c=c+real(i-1,dp)*g(i)**2
    end do
    c=c/sum(g*g)
  end function center_energy

  function up_sample(x,f,fill_value) result(y)
    real(dp),intent(in)::x(:)
    integer,intent(in)::f
    real(dp),intent(in),optional::fill_value
    real(dp),allocatable::y(:)
    real(dp)::fill
    integer::i
    fill=0.0_dp
    if(present(fill_value))fill=fill_value
    if(f<1)then
    allocate(y(0))
    return
    end if
    allocate(y(size(x)*f))
    y=fill
    do i=1,size(x)
    y((i-1)*f+1)=x(i)
    end do
  end function up_sample

  function dwt_hilbert(x,wf,n_levels,boundary) result(ans)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in)::wf
    integer,intent(in),optional::n_levels
    character(len=*),intent(in),optional::boundary
    type(complex_wavelet_transform)::ans
    type(wavelet_transform)::a,b
    type(status_type)::st
    integer::levels,j
    character(len=16)::bname
    levels=4
    if(present(n_levels))levels=n_levels
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    block
      use waveslim_types,only:hilbert_filter_type
      type(hilbert_filter_type)::f
      f=hilbert_filter(wf,st)
      if(.not.st%ok())then
      ans%status=st
      return
      end if
      a=dwt_filters(x,f%h1,f%h0,levels,bname,wf)
      b=dwt_filters(x,f%g1,f%g0,levels,bname,wf)
    end block
    allocate(ans%detail(levels))
    do j=1,levels
    ans%detail(j)%values=cmplx(a%detail(j)%values,b%detail(j)%values,dp)
    end do
    ans%smooth=cmplx(a%smooth,b%smooth,dp)
    ans%wavelet=wf
    ans%boundary=bname
    ans%method='dwt_hilbert'
    ans%original_length=size(x)
    call clear_status(ans%status)
  end function dwt_hilbert

  function modwt_hilbert(x,wf,n_levels,boundary) result(ans)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in)::wf
    integer,intent(in),optional::n_levels
    character(len=*),intent(in),optional::boundary
    type(complex_wavelet_transform)::ans
    type(wavelet_transform)::a,b
    type(status_type)::st
    integer::levels,j
    character(len=16)::bname
    levels=4
    if(present(n_levels))levels=n_levels
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    block
      use waveslim_types,only:hilbert_filter_type
      type(hilbert_filter_type)::f
      f=hilbert_filter(wf,st)
      if(.not.st%ok())then
      ans%status=st
      return
      end if
      a=modwt_filters(x,f%h1/sqrt2,f%h0/sqrt2,levels,bname,wf)
      b=modwt_filters(x,f%g1/sqrt2,f%g0/sqrt2,levels,bname,wf)
    end block
    allocate(ans%detail(levels))
    do j=1,levels
    ans%detail(j)%values=cmplx(a%detail(j)%values,b%detail(j)%values,dp)
    end do
    ans%smooth=cmplx(a%smooth,b%smooth,dp)
    ans%wavelet=wf
    ans%boundary=bname
    ans%method='modwt_hilbert'
    ans%original_length=size(x)
    call clear_status(ans%status)
  end function modwt_hilbert

  function idwt_hilbert(y,tree) result(x)
    type(complex_wavelet_transform),intent(in)::y
    integer,intent(in),optional::tree
    real(dp),allocatable::x(:)
    type(wavelet_transform)::r
    integer::j,t
    t=1
    if(present(tree))t=tree
    r%wavelet=y%wavelet
    r%boundary=y%boundary
    r%method='dwt'
    r%original_length=y%original_length
    allocate(r%detail(y%levels()))
    do j=1,y%levels()
    if(t==1)then
    r%detail(j)%values=real(y%detail(j)%values,dp)
    else
    r%detail(j)%values=aimag(y%detail(j)%values)
    end if
    end do
    if(t==1)then
    r%smooth=real(y%smooth,dp)
    else
    r%smooth=aimag(y%smooth)
    end if
    x=idwt_hilbert_tree(r,t)
  end function idwt_hilbert

  function imodwt_hilbert(y,tree) result(x)
    type(complex_wavelet_transform),intent(in)::y
    integer,intent(in),optional::tree
    real(dp),allocatable::x(:)
    type(wavelet_transform)::r
    integer::j,t
    t=1
    if(present(tree))t=tree
    r%wavelet=y%wavelet
    r%boundary=y%boundary
    r%method='modwt'
    r%original_length=y%original_length
    allocate(r%detail(y%levels()))
    do j=1,y%levels()
    if(t==1)then
    r%detail(j)%values=real(y%detail(j)%values,dp)
    else
    r%detail(j)%values=aimag(y%detail(j)%values)
    end if
    end do
    if(t==1)then
    r%smooth=real(y%smooth,dp)
    else
    r%smooth=aimag(y%smooth)
    end if
    x=imodwt_hilbert_tree(r,t)
  end function imodwt_hilbert

  function dwt_filters(x,h,g,levels,boundary,wname) result(r)
    real(dp),intent(in)::x(:),h(:),g(:)
    integer,intent(in)::levels
    character(len=*),intent(in)::boundary,wname
    type(wavelet_transform)::r
    real(dp),allocatable::work(:),w(:),v(:)
    integer::j
    if(trim(boundary)=='reflection')then
    allocate(work(2*size(x)))
    work(1:size(x))=x
    work(size(x)+1:)=x(size(x):1:-1)
    else
    work=x
    end if
    allocate(r%detail(levels))
    do j=1,levels
    call dwt_step(work,h,g,w,v)
    r%detail(j)%values=w
    call move_alloc(v,work)
    end do
    r%smooth=work
    r%wavelet=wname
    r%boundary=boundary
    r%method='dwt'
    r%original_length=size(x)
    call clear_status(r%status)
  end function dwt_filters

  function modwt_filters(x,h,g,levels,boundary,wname) result(r)
    real(dp),intent(in)::x(:),h(:),g(:)
    integer,intent(in)::levels
    character(len=*),intent(in)::boundary,wname
    type(wavelet_transform)::r
    real(dp),allocatable::work(:),w(:),v(:)
    integer::j
    if(trim(boundary)=='reflection')then
    allocate(work(2*size(x)))
    work(1:size(x))=x
    work(size(x)+1:)=x(size(x):1:-1)
    else
    work=x
    end if
    allocate(r%detail(levels))
    do j=1,levels
    call modwt_step(work,j,h,g,w,v)
    r%detail(j)%values=w
    call move_alloc(v,work)
    end do
    r%smooth=work
    r%wavelet=wname
    r%boundary=boundary
    r%method='modwt'
    r%original_length=size(x)
    call clear_status(r%status)
  end function modwt_filters

  function idwt_hilbert_tree(r,tree) result(x)
    type(wavelet_transform),intent(in)::r
    integer,intent(in)::tree
    real(dp),allocatable::x(:),work(:),tmp(:)
    integer::j,nout
    type(status_type)::st
    block
      use waveslim_types,only:hilbert_filter_type
      type(hilbert_filter_type)::f
      f=hilbert_filter(r%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0))
      return
      end if
      work=r%smooth
      do j=r%levels(),1,-1
        if(tree==1)then
        call idwt_step(r%detail(j)%values,work,f%h1,f%h0,tmp)
        else
        call idwt_step(r%detail(j)%values,work,f%g1,f%g0,tmp)
        end if
        call move_alloc(tmp,work)
      end do
    end block
    nout=size(work)
    if(trim(r%boundary)=='reflection')nout=r%original_length
    x=work(1:nout)
  end function idwt_hilbert_tree

  function imodwt_hilbert_tree(r,tree) result(x)
    type(wavelet_transform),intent(in)::r
    integer,intent(in)::tree
    real(dp),allocatable::x(:),work(:),tmp(:)
    integer::j,nout
    type(status_type)::st
    block
      use waveslim_types,only:hilbert_filter_type
      type(hilbert_filter_type)::f
      f=hilbert_filter(r%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0))
      return
      end if
      work=r%smooth
      do j=r%levels(),1,-1
        if(tree==1)then
        call imodwt_step(r%detail(j)%values,work,j,f%h1/sqrt2,f%h0/sqrt2,tmp)
        else
        call imodwt_step(r%detail(j)%values,work,j,f%g1/sqrt2,f%g0/sqrt2,tmp)
        end if
        call move_alloc(tmp,work)
      end do
    end block
    nout=size(work)
    if(trim(r%boundary)=='reflection')nout=r%original_length
    x=work(1:nout)
  end function imodwt_hilbert_tree

end module waveslim_transform_1d
