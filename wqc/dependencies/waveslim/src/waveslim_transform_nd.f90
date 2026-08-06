! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_transform_nd
  use waveslim_kinds, only : dp, sqrt2
  use waveslim_status, only : clear_status, set_status, waveslim_invalid_input, waveslim_invalid_level
  use waveslim_types, only : wavelet_transform_2d, wavelet_transform_3d, mra_result, real_vector
  use waveslim_filters, only : wave_filter
  use waveslim_transform_1d, only : dwt_step, idwt_step, modwt_step, imodwt_step
  implicit none
  private
  public :: dwt_2d, idwt_2d, modwt_2d, imodwt_2d, mra_2d
  public :: dwt_3d, idwt_3d, modwt_3d, imodwt_3d, mra_3d
  public :: convolve_2d, shift_2d
  public :: dwt2_step, idwt2_step
contains
  function dwt_2d(x,wf,j_levels,boundary) result(ans)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in),optional::wf,boundary
    integer,intent(in),optional::j_levels
    type(wavelet_transform_2d)::ans
    character(len=16)::wname,bname
    integer::levels,j,nr,nc
    real(dp),allocatable::work(:,:),ll(:,:),lh(:,:),hl(:,:),hh(:,:)
    wname='la8'
    if(present(wf))wname=trim(wf)
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    levels=4
    if(present(j_levels))levels=j_levels
    ans%wavelet=wname
    ans%boundary=bname
    ans%method='dwt'
    ans%original_shape=shape(x)
    if(bname=='reflection')then
      allocate(work(2*size(x,1),2*size(x,2)))
      work(1:size(x,1),1:size(x,2))=x
      work(size(x,1)+1:,1:size(x,2))=x(size(x,1):1:-1,:)
      work(:,size(x,2)+1:)=work(:,size(x,2):1:-1)
    else if(bname=='periodic')then
      allocate(work(size(x,1),size(x,2)))
      work=x
    else
      call set_status(ans%status,waveslim_invalid_input,'invalid boundary')
      return
    end if
    nr=size(work,1)
    nc=size(work,2)
    if(levels<1.or.mod(nr,2**levels)/=0.or.mod(nc,2**levels)/=0)then
      call set_status(ans%status,waveslim_invalid_level,'2D DWT dimensions must be divisible by 2**levels')
      return
    end if
    allocate(ans%level(levels))
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wname,st)
      if(.not.st%ok())then
      ans%status=st
      return
      end if
      do j=1,levels
        call dwt2_step(work,f%hpf,f%lpf,ll,lh,hl,hh)
        ans%level(j)%lh=lh
        ans%level(j)%hl=hl
        ans%level(j)%hh=hh
        call move_alloc(ll,work)
      end do
    end block
    ans%smooth=work
    call clear_status(ans%status)
  end function dwt_2d

  function idwt_2d(wt) result(x)
    type(wavelet_transform_2d),intent(in)::wt
    real(dp),allocatable::x(:,:)
    real(dp),allocatable::work(:,:),tmp(:,:)
    integer::j,nr,nc
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wt%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0,0))
      return
      end if
      work=wt%smooth
      do j=size(wt%level),1,-1
        call idwt2_step(work,wt%level(j)%lh,wt%level(j)%hl,wt%level(j)%hh,f%hpf,f%lpf,tmp)
        call move_alloc(tmp,work)
      end do
    end block
    nr=size(work,1)
    nc=size(work,2)
    if(trim(wt%boundary)=='reflection')then
    nr=wt%original_shape(1)
    nc=wt%original_shape(2)
    end if
    x=work(1:nr,1:nc)
  end function idwt_2d

  function modwt_2d(x,wf,j_levels,boundary) result(ans)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in),optional::wf,boundary
    integer,intent(in),optional::j_levels
    type(wavelet_transform_2d)::ans
    character(len=16)::wname,bname
    integer::levels,j
    real(dp),allocatable::work(:,:),ll(:,:),lh(:,:),hl(:,:),hh(:,:),h(:),g(:)
    wname='la8'
    if(present(wf))wname=trim(wf)
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    levels=4
    if(present(j_levels))levels=j_levels
    ans%wavelet=wname
    ans%boundary=bname
    ans%method='modwt'
    ans%original_shape=shape(x)
    if(bname=='reflection')then
      allocate(work(2*size(x,1),2*size(x,2)))
      work(1:size(x,1),1:size(x,2))=x
      work(size(x,1)+1:,1:size(x,2))=x(size(x,1):1:-1,:)
      work(:,size(x,2)+1:)=work(:,size(x,2):1:-1)
    else if(bname=='periodic')then
    allocate(work(size(x,1),size(x,2)))
    work=x
    else
    call set_status(ans%status,waveslim_invalid_input,'invalid boundary')
    return
    end if
    if(levels<1.or.2**levels>min(size(work,1),size(work,2)))then
      call set_status(ans%status,waveslim_invalid_level,'2D MODWT exceeds dimensions')
      return
    end if
    allocate(ans%level(levels))
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wname,st)
      if(.not.st%ok())then
      ans%status=st
      return
      end if
      h=f%hpf/sqrt2
      g=f%lpf/sqrt2
      do j=1,levels
        call modwt2_step(work,j,h,g,ll,lh,hl,hh)
        ans%level(j)%lh=lh
        ans%level(j)%hl=hl
        ans%level(j)%hh=hh
        call move_alloc(ll,work)
      end do
    end block
    ans%smooth=work
    call clear_status(ans%status)
  end function modwt_2d

  function imodwt_2d(wt) result(x)
    type(wavelet_transform_2d),intent(in)::wt
    real(dp),allocatable::x(:,:)
    real(dp),allocatable::work(:,:),tmp(:,:),h(:),g(:)
    integer::j,nr,nc
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wt%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0,0))
      return
      end if
      h=f%hpf/sqrt2
      g=f%lpf/sqrt2
      work=wt%smooth
      do j=size(wt%level),1,-1
        call imodwt2_step(work,wt%level(j)%lh,wt%level(j)%hl,wt%level(j)%hh,j,h,g,tmp)
        call move_alloc(tmp,work)
      end do
    end block
    nr=size(work,1)
    nc=size(work,2)
    if(trim(wt%boundary)=='reflection')then
    nr=wt%original_shape(1)
    nc=wt%original_shape(2)
    end if
    x=work(1:nr,1:nc)
  end function imodwt_2d

  function mra_2d(x,wf,j_levels,method,boundary) result(parts)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in),optional::wf,method,boundary
    integer,intent(in),optional::j_levels
    type(wavelet_transform_2d),allocatable::parts(:)
    type(wavelet_transform_2d)::wt,blank
    character(len=16)::wname,mname,bname
    integer::levels,j
    wname='la8'
    if(present(wf))wname=trim(wf)
    mname='modwt'
    if(present(method))mname=trim(method)
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    levels=4
    if(present(j_levels))levels=j_levels
    if(mname=='dwt')then
    wt=dwt_2d(x,wname,levels,bname)
    else
    wt=modwt_2d(x,wname,levels,bname)
    end if
    blank=zero_2d_like(wt)
    allocate(parts(levels+1),source=blank)
    parts(levels+1)%smooth=wt%smooth
    do j=1,levels
      parts(j)%level(j)%lh=wt%level(j)%lh
      parts(j)%level(j)%hl=wt%level(j)%hl
      parts(j)%level(j)%hh=wt%level(j)%hh
    end do
  end function mra_2d

  function zero_2d_like(wt) result(z)
    type(wavelet_transform_2d),intent(in)::wt
    type(wavelet_transform_2d)::z
    integer::j
    z%wavelet=wt%wavelet
    z%boundary=wt%boundary
    z%method=wt%method
    z%original_shape=wt%original_shape
    allocate(z%level(size(wt%level)))
    do j=1,size(wt%level)
      allocate(z%level(j)%lh(size(wt%level(j)%lh,1),size(wt%level(j)%lh,2)))
      z%level(j)%lh=0.0_dp
      allocate(z%level(j)%hl(size(wt%level(j)%hl,1),size(wt%level(j)%hl,2)))
      z%level(j)%hl=0.0_dp
      allocate(z%level(j)%hh(size(wt%level(j)%hh,1),size(wt%level(j)%hh,2)))
      z%level(j)%hh=0.0_dp
    end do
    allocate(z%smooth(size(wt%smooth,1),size(wt%smooth,2)))
    z%smooth=0.0_dp
    call clear_status(z%status)
  end function zero_2d_like

  subroutine dwt2_step(x,h,g,ll,lh,hl,hh)
    real(dp),intent(in)::x(:,:),h(:),g(:)
    real(dp),allocatable,intent(out)::ll(:,:),lh(:,:),hl(:,:),hh(:,:)
    real(dp),allocatable::low_h(:,:),high_h(:,:),w(:),v(:)
    integer::i,j,nr,nc
    nr=size(x,1)
    nc=size(x,2)
    allocate(low_h(nr,nc/2),high_h(nr,nc/2))
    do i=1,nr
    call dwt_step(x(i,:),h,g,w,v)
    high_h(i,:)=w
    low_h(i,:)=v
    end do
    allocate(ll(nr/2,nc/2),lh(nr/2,nc/2),hl(nr/2,nc/2),hh(nr/2,nc/2))
    do j=1,nc/2
      call dwt_step(low_h(:,j),h,g,w,v)
      lh(:,j)=w
      ll(:,j)=v
      call dwt_step(high_h(:,j),h,g,w,v)
      hh(:,j)=w
      hl(:,j)=v
    end do
  end subroutine dwt2_step

  subroutine idwt2_step(ll,lh,hl,hh,h,g,x)
    real(dp),intent(in)::ll(:,:),lh(:,:),hl(:,:),hh(:,:),h(:),g(:)
    real(dp),allocatable,intent(out)::x(:,:)
    real(dp),allocatable::low_h(:,:),high_h(:,:),tmp(:)
    integer::i,j,nr,nc
    nr=size(ll,1)
    nc=size(ll,2)
    allocate(low_h(2*nr,nc),high_h(2*nr,nc))
    do j=1,nc
    call idwt_step(lh(:,j),ll(:,j),h,g,tmp)
    low_h(:,j)=tmp
    call idwt_step(hh(:,j),hl(:,j),h,g,tmp)
    high_h(:,j)=tmp
    end do
    allocate(x(2*nr,2*nc))
    do i=1,2*nr
    call idwt_step(high_h(i,:),low_h(i,:),h,g,tmp)
    x(i,:)=tmp
    end do
  end subroutine idwt2_step

  subroutine modwt2_step(x,level,h,g,ll,lh,hl,hh)
    real(dp),intent(in)::x(:,:),h(:),g(:)
    integer,intent(in)::level
    real(dp),allocatable,intent(out)::ll(:,:),lh(:,:),hl(:,:),hh(:,:)
    real(dp),allocatable::low_h(:,:),high_h(:,:),w(:),v(:)
    integer::i,j,nr,nc
    nr=size(x,1)
    nc=size(x,2)
    allocate(low_h(nr,nc),high_h(nr,nc))
    do i=1,nr
    call modwt_step(x(i,:),level,h,g,w,v)
    high_h(i,:)=w
    low_h(i,:)=v
    end do
    allocate(ll(nr,nc),lh(nr,nc),hl(nr,nc),hh(nr,nc))
    do j=1,nc
    call modwt_step(low_h(:,j),level,h,g,w,v)
    lh(:,j)=w
    ll(:,j)=v
    call modwt_step(high_h(:,j),level,h,g,w,v)
    hh(:,j)=w
    hl(:,j)=v
    end do
  end subroutine modwt2_step

  subroutine imodwt2_step(ll,lh,hl,hh,level,h,g,x)
    real(dp),intent(in)::ll(:,:),lh(:,:),hl(:,:),hh(:,:),h(:),g(:)
    integer,intent(in)::level
    real(dp),allocatable,intent(out)::x(:,:)
    real(dp),allocatable::low_h(:,:),high_h(:,:),tmp(:)
    integer::i,j,nr,nc
    nr=size(ll,1)
    nc=size(ll,2)
    allocate(low_h(nr,nc),high_h(nr,nc))
    do j=1,nc
    call imodwt_step(lh(:,j),ll(:,j),level,h,g,tmp)
    low_h(:,j)=tmp
    call imodwt_step(hh(:,j),hl(:,j),level,h,g,tmp)
    high_h(:,j)=tmp
    end do
    allocate(x(nr,nc))
    do i=1,nr
    call imodwt_step(high_h(i,:),low_h(i,:),level,h,g,tmp)
    x(i,:)=tmp
    end do
  end subroutine imodwt2_step

  function dwt_3d(x,wf,j_levels,boundary) result(ans)
    real(dp),intent(in)::x(:,:,:)
    character(len=*),intent(in),optional::wf,boundary
    integer,intent(in),optional::j_levels
    type(wavelet_transform_3d)::ans
    character(len=16)::wname,bname
    integer::levels,j
    real(dp),allocatable::work(:,:,:),bands(:,:,:,:)
    wname='la8'
    if(present(wf))wname=trim(wf)
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    levels=4
    if(present(j_levels))levels=j_levels
    ans%wavelet=wname
    ans%boundary=bname
    ans%method='dwt'
    ans%original_shape=shape(x)
    if(bname/='periodic')then
    call set_status(ans%status,waveslim_invalid_input,'3D port currently supports periodic boundaries')
    return
    end if
    if(any(mod(shape(x),2**levels)/=0))then
    call set_status(ans%status,waveslim_invalid_level,'3D dimensions must be divisible by 2**levels')
    return
    end if
    allocate(work(size(x,1),size(x,2),size(x,3)))
    work=x
    allocate(ans%level(levels))
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wname,st)
      if(.not.st%ok())then
      ans%status=st
      return
      end if
      do j=1,levels
      call dwt3_step(work,f%hpf,f%lpf,bands)
      ans%level(j)%band=bands(:,:,:,2:8)
      work=bands(:,:,:,1)
      end do
    end block
    ans%smooth=work
    call clear_status(ans%status)
  end function dwt_3d

  function idwt_3d(wt) result(x)
    type(wavelet_transform_3d),intent(in)::wt
    real(dp),allocatable::x(:,:,:),work(:,:,:),bands(:,:,:,:),tmp(:,:,:)
    integer::j
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wt%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0,0,0))
      return
      end if
      work=wt%smooth
      do j=size(wt%level),1,-1
      allocate(bands(size(work,1),size(work,2),size(work,3),8))
      bands(:,:,:,1)=work
      bands(:,:,:,2:8)=wt%level(j)%band
      call idwt3_step(bands,f%hpf,f%lpf,tmp)
      call move_alloc(tmp,work)
      deallocate(bands)
      end do
    end block
    x=work
  end function idwt_3d

  function modwt_3d(x,wf,j_levels,boundary) result(ans)
    real(dp),intent(in)::x(:,:,:)
    character(len=*),intent(in),optional::wf,boundary
    integer,intent(in),optional::j_levels
    type(wavelet_transform_3d)::ans
    character(len=16)::wname,bname
    integer::levels,j
    real(dp),allocatable::work(:,:,:),bands(:,:,:,:),h(:),g(:)
    wname='la8'
    if(present(wf))wname=trim(wf)
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    levels=4
    if(present(j_levels))levels=j_levels
    ans%wavelet=wname
    ans%boundary=bname
    ans%method='modwt'
    ans%original_shape=shape(x)
    if(bname/='periodic')then
    call set_status(ans%status,waveslim_invalid_input,'3D port currently supports periodic boundaries')
    return
    end if
    allocate(work(size(x,1),size(x,2),size(x,3)))
    work=x
    allocate(ans%level(levels))
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wname,st)
      if(.not.st%ok())then
      ans%status=st
      return
      end if
      h=f%hpf/sqrt2
      g=f%lpf/sqrt2
      do j=1,levels
      call modwt3_step(work,j,h,g,bands)
      ans%level(j)%band=bands(:,:,:,2:8)
      work=bands(:,:,:,1)
      end do
    end block
    ans%smooth=work
    call clear_status(ans%status)
  end function modwt_3d

  function imodwt_3d(wt) result(x)
    type(wavelet_transform_3d),intent(in)::wt
    real(dp),allocatable::x(:,:,:),work(:,:,:),bands(:,:,:,:),tmp(:,:,:),h(:),g(:)
    integer::j
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wt%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0,0,0))
      return
      end if
      h=f%hpf/sqrt2
      g=f%lpf/sqrt2
      work=wt%smooth
      do j=size(wt%level),1,-1
      allocate(bands(size(work,1),size(work,2),size(work,3),8))
      bands(:,:,:,1)=work
      bands(:,:,:,2:8)=wt%level(j)%band
      call imodwt3_step(bands,j,h,g,tmp)
      call move_alloc(tmp,work)
      deallocate(bands)
      end do
    end block
    x=work
  end function imodwt_3d

  function mra_3d(x,wf,j_levels,method,boundary) result(parts)
    real(dp),intent(in)::x(:,:,:)
    character(len=*),intent(in),optional::wf,method,boundary
    integer,intent(in),optional::j_levels
    real(dp),allocatable::parts(:,:,:,:)
    type(wavelet_transform_3d)::wt,z
    character(len=16)::wname,mname,bname
    integer::levels,j
    wname='la8'
    if(present(wf))wname=trim(wf)
    mname='modwt'
    if(present(method))mname=trim(method)
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    levels=4
    if(present(j_levels))levels=j_levels
    if(mname=='dwt')then
    wt=dwt_3d(x,wname,levels,bname)
    else
    wt=modwt_3d(x,wname,levels,bname)
    end if
    allocate(parts(size(x,1),size(x,2),size(x,3),levels+1))
    parts=0.0_dp
    do j=1,levels
      z=zero_3d_like(wt)
      z%level(j)%band=wt%level(j)%band
      if(mname=='dwt')then
      parts(:,:,:,j)=idwt_3d(z)
      else
      parts(:,:,:,j)=imodwt_3d(z)
      end if
    end do
    z=zero_3d_like(wt)
    z%smooth=wt%smooth
    if(mname=='dwt')then
    parts(:,:,:,levels+1)=idwt_3d(z)
    else
    parts(:,:,:,levels+1)=imodwt_3d(z)
    end if
  end function mra_3d

  function zero_3d_like(wt) result(z)
    type(wavelet_transform_3d),intent(in)::wt
    type(wavelet_transform_3d)::z
    integer::j
    z%wavelet=wt%wavelet
    z%boundary=wt%boundary
    z%method=wt%method
    z%original_shape=wt%original_shape
    allocate(z%level(size(wt%level)))
    do j=1,size(wt%level)
    allocate(z%level(j)%band(size(wt%level(j)%band,1),size(wt%level(j)%band,2),size(wt%level(j)%band,3),7))
    z%level(j)%band=0.0_dp
    end do
    allocate(z%smooth(size(wt%smooth,1),size(wt%smooth,2),size(wt%smooth,3)))
    z%smooth=0.0_dp
    call clear_status(z%status)
  end function zero_3d_like

  subroutine dwt3_step(x,h,g,bands)
    real(dp),intent(in)::x(:,:,:),h(:),g(:)
    real(dp),allocatable,intent(out)::bands(:,:,:,:)
    real(dp),allocatable::a(:,:,:,:),b(:,:,:,:),w(:),v(:)
    integer::i,j,k,n1,n2,n3,p
    n1=size(x,1)
    n2=size(x,2)
    n3=size(x,3)
    allocate(a(n1,n2,n3/2,2))
    do i=1,n1
    do j=1,n2
    call dwt_step(x(i,j,:),h,g,w,v)
    a(i,j,:,1)=v
    a(i,j,:,2)=w
    end do
    end do
    allocate(b(n1,n2/2,n3/2,4))
    do i=1,n1
    do k=1,n3/2
    do p=1,2
    call dwt_step(a(i,:,k,p),h,g,w,v)
    b(i,:,k,2*p-1)=v
    b(i,:,k,2*p)=w
    end do
    end do
    end do
    allocate(bands(n1/2,n2/2,n3/2,8))
    do j=1,n2/2
    do k=1,n3/2
    do p=1,4
    call dwt_step(b(:,j,k,p),h,g,w,v)
    bands(:,j,k,2*p-1)=v
    bands(:,j,k,2*p)=w
    end do
    end do
    end do
  end subroutine dwt3_step

  subroutine idwt3_step(bands,h,g,x)
    real(dp),intent(in)::bands(:,:,:,:),h(:),g(:)
    real(dp),allocatable,intent(out)::x(:,:,:)
    real(dp),allocatable::b(:,:,:,:),a(:,:,:,:),tmp(:)
    integer::i,j,k,p,n1,n2,n3
    n1=size(bands,1)
    n2=size(bands,2)
    n3=size(bands,3)
    allocate(b(2*n1,n2,n3,4))
    do j=1,n2
    do k=1,n3
    do p=1,4
    call idwt_step(bands(:,j,k,2*p),bands(:,j,k,2*p-1),h,g,tmp)
    b(:,j,k,p)=tmp
    end do
    end do
    end do
    allocate(a(2*n1,2*n2,n3,2))
    do i=1,2*n1
    do k=1,n3
    do p=1,2
    call idwt_step(b(i,:,k,2*p),b(i,:,k,2*p-1),h,g,tmp)
    a(i,:,k,p)=tmp
    end do
    end do
    end do
    allocate(x(2*n1,2*n2,2*n3))
    do i=1,2*n1
    do j=1,2*n2
    call idwt_step(a(i,j,:,2),a(i,j,:,1),h,g,tmp)
    x(i,j,:)=tmp
    end do
    end do
  end subroutine idwt3_step

  subroutine modwt3_step(x,level,h,g,bands)
    real(dp),intent(in)::x(:,:,:),h(:),g(:)
    integer,intent(in)::level
    real(dp),allocatable,intent(out)::bands(:,:,:,:)
    real(dp),allocatable::a(:,:,:,:),b(:,:,:,:),w(:),v(:)
    integer::i,j,k,n1,n2,n3,p
    n1=size(x,1)
    n2=size(x,2)
    n3=size(x,3)
    allocate(a(n1,n2,n3,2))
    do i=1,n1
    do j=1,n2
    call modwt_step(x(i,j,:),level,h,g,w,v)
    a(i,j,:,1)=v
    a(i,j,:,2)=w
    end do
    end do
    allocate(b(n1,n2,n3,4))
    do i=1,n1
    do k=1,n3
    do p=1,2
    call modwt_step(a(i,:,k,p),level,h,g,w,v)
    b(i,:,k,2*p-1)=v
    b(i,:,k,2*p)=w
    end do
    end do
    end do
    allocate(bands(n1,n2,n3,8))
    do j=1,n2
    do k=1,n3
    do p=1,4
    call modwt_step(b(:,j,k,p),level,h,g,w,v)
    bands(:,j,k,2*p-1)=v
    bands(:,j,k,2*p)=w
    end do
    end do
    end do
  end subroutine modwt3_step

  subroutine imodwt3_step(bands,level,h,g,x)
    real(dp),intent(in)::bands(:,:,:,:),h(:),g(:)
    integer,intent(in)::level
    real(dp),allocatable,intent(out)::x(:,:,:)
    real(dp),allocatable::b(:,:,:,:),a(:,:,:,:),tmp(:)
    integer::i,j,k,p,n1,n2,n3
    n1=size(bands,1)
    n2=size(bands,2)
    n3=size(bands,3)
    allocate(b(n1,n2,n3,4))
    do j=1,n2
    do k=1,n3
    do p=1,4
    call imodwt_step(bands(:,j,k,2*p),bands(:,j,k,2*p-1),level,h,g,tmp)
    b(:,j,k,p)=tmp
    end do
    end do
    end do
    allocate(a(n1,n2,n3,2))
    do i=1,n1
    do k=1,n3
    do p=1,2
    call imodwt_step(b(i,:,k,2*p),b(i,:,k,2*p-1),level,h,g,tmp)
    a(i,:,k,p)=tmp
    end do
    end do
    end do
    allocate(x(n1,n2,n3))
    do i=1,n1
    do j=1,n2
    call imodwt_step(a(i,j,:,2),a(i,j,:,1),level,h,g,tmp)
    x(i,j,:)=tmp
    end do
    end do
  end subroutine imodwt3_step

  function convolve_2d(x,y,conjugate_kernel,circular) result(z)
    real(dp),intent(in)::x(:,:),y(:,:)
    logical,intent(in),optional::conjugate_kernel,circular
    real(dp),allocatable::z(:,:)
    logical::conj,circ
    integer::i,j,a,b,ii,jj,nr,nc,kr,kc
    conj=.true.
    if(present(conjugate_kernel))conj=conjugate_kernel
    circ=.true.
    if(present(circular))circ=circular
    nr=size(x,1)
    nc=size(x,2)
    kr=size(y,1)
    kc=size(y,2)
    allocate(z(nr,nc))
    z=0.0_dp
    do i=1,nr
    do j=1,nc
    do a=1,kr
    do b=1,kc
      if(conj)then
      ii=i-(a-1)
      jj=j-(b-1)
      else
      ii=i+(a-1)
      jj=j+(b-1)
      end if
      if(circ)then
      ii=modulo(ii-1,nr)+1
      jj=modulo(jj-1,nc)+1
      z(i,j)=z(i,j)+x(ii,jj)*y(a,b)
      else if(ii>=1.and.ii<=nr.and.jj>=1.and.jj<=nc)then
      z(i,j)=z(i,j)+x(ii,jj)*y(a,b)
      end if
    end do
    end do
    end do
    end do
  end function convolve_2d

  function shift_2d(x,row_shift,column_shift) result(y)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::row_shift,column_shift
    real(dp),allocatable::y(:,:)
    integer::i,j,nr,nc
    nr=size(x,1)
    nc=size(x,2)
    allocate(y(nr,nc))
    do i=1,nr
    do j=1,nc
    y(i,j)=x(modulo(i-1+row_shift,nr)+1,modulo(j-1+column_shift,nc)+1)
    end do
    end do
  end function shift_2d
end module waveslim_transform_nd
