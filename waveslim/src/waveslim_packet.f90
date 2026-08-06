! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_packet
  use waveslim_kinds, only : dp, sqrt2
  use waveslim_status, only : clear_status, set_status, waveslim_invalid_input, waveslim_invalid_level
  use waveslim_types, only : packet_transform, real_vector, test_result
  use waveslim_filters, only : wave_filter
  use waveslim_transform_1d, only : dwt_step, idwt_step, modwt_step
  use waveslim_math, only : autocorrelation, chi_square_cdf, fft_complex, mean_value
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: dwpt, modwpt, idwpt, packet_basis, ortho_basis
  public :: dwpt_brick_wall, css_test, entropy_test, cpgram_test, portmanteau_test
contains
  function dwpt(x,wf,n_levels,boundary) result(ans)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in),optional::wf,boundary
    integer,intent(in),optional::n_levels
    type(packet_transform)::ans
    character(len=16)::wname,bname
    integer::levels,j,n,parent
    real(dp),allocatable::work(:),w(:),v(:)
    wname='la8'
    if(present(wf))wname=trim(wf)
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    levels=4
    if(present(n_levels))levels=n_levels
    ans%wavelet=wname
    ans%boundary=bname
    ans%method='dwpt'
    ans%original_length=size(x)
    if(bname=='reflection')then
      allocate(work(2*size(x)))
      work(1:size(x))=x
      work(size(x)+1:)=x(size(x):1:-1)
    else if(bname=='periodic')then
    allocate(work(size(x)))
    work=x
    else
    call set_status(ans%status,waveslim_invalid_input,'invalid boundary')
    return
    end if
    if(levels<1 .or. 2**levels>size(work) .or. mod(size(work),2**levels)/=0)then
      call set_status(ans%status,waveslim_invalid_level,'invalid packet depth')
      return
    end if
    allocate(ans%level(0:levels))
    allocate(ans%level(0)%node(1))
    ans%level(0)%node(1)%values=work
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
        allocate(ans%level(j)%node(2**j))
        do parent=1,2**(j-1)
          call dwt_step(ans%level(j-1)%node(parent)%values,f%hpf,f%lpf,w,v)
          n=2*parent-1
          if(mod(parent-1,2)==0)then
            ans%level(j)%node(n)%values=v
            ans%level(j)%node(n+1)%values=w
          else
            ans%level(j)%node(n)%values=w
            ans%level(j)%node(n+1)%values=v
          end if
        end do
      end do
    end block
    call clear_status(ans%status)
  end function dwpt

  function modwpt(x,wf,n_levels,boundary) result(ans)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in),optional::wf,boundary
    integer,intent(in),optional::n_levels
    type(packet_transform)::ans
    character(len=16)::wname,bname
    integer::levels,j,n,parent
    real(dp),allocatable::work(:),w(:),v(:),h(:),g(:)
    wname='la8'
    if(present(wf))wname=trim(wf)
    bname='periodic'
    if(present(boundary))bname=trim(boundary)
    levels=4
    if(present(n_levels))levels=n_levels
    ans%wavelet=wname
    ans%boundary=bname
    ans%method='modwpt'
    ans%original_length=size(x)
    if(bname=='reflection')then
      allocate(work(2*size(x)))
      work(1:size(x))=x
      work(size(x)+1:)=x(size(x):1:-1)
    else if(bname=='periodic')then
    allocate(work(size(x)))
    work=x
    else
    call set_status(ans%status,waveslim_invalid_input,'invalid boundary')
    return
    end if
    if(levels<1 .or. 2**levels>size(work))then
      call set_status(ans%status,waveslim_invalid_level,'invalid packet depth')
      return
    end if
    allocate(ans%level(0:levels))
    allocate(ans%level(0)%node(1))
    ans%level(0)%node(1)%values=work
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
        allocate(ans%level(j)%node(2**j))
        do parent=1,2**(j-1)
          call modwt_step(ans%level(j-1)%node(parent)%values,j,h,g,w,v)
          n=2*parent-1
          if(mod(parent-1,2)==0)then
            ans%level(j)%node(n)%values=v
            ans%level(j)%node(n+1)%values=w
          else
            ans%level(j)%node(n)%values=w
            ans%level(j)%node(n+1)%values=v
          end if
        end do
      end do
    end block
    call clear_status(ans%status)
  end function modwpt

  function idwpt(tree,basis_mask) result(x)
    type(packet_transform),intent(in)::tree
    logical,intent(in),optional::basis_mask(:)
    real(dp),allocatable::x(:)
    type(real_vector),allocatable::nodes(:),parents(:)
    logical,allocatable::mask(:),pmask(:)
    integer::j,n,k,offset,total,m,nout
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(tree%wavelet,st)
      if(.not.st%ok())then
      allocate(x(0))
      return
      end if
      j=tree%levels()
      allocate(nodes(2**j),mask(2**j))
      mask=.true.
      do n=1,2**j
      nodes(n)%values=tree%level(j)%node(n)%values
      end do
      if(present(basis_mask))then
        total=2**(j+1)-2
        if(size(basis_mask)==total)then
          offset=2**j-2
          mask=basis_mask(offset+1:offset+2**j)
        end if
      end if
      do j=tree%levels(),1,-1
        allocate(parents(2**(j-1)),pmask(2**(j-1)))
        pmask=.false.
        do n=1,2**(j-1)
          k=2*n-1
          if(mask(k).or.mask(k+1))then
            m=size(nodes(k)%values)
            if(mod(n-1,2)==0)then
              call idwt_step(nodes(k+1)%values,nodes(k)%values,f%hpf,f%lpf,parents(n)%values)
            else
              call idwt_step(nodes(k)%values,nodes(k+1)%values,f%hpf,f%lpf,parents(n)%values)
            end if
            pmask(n)=.true.
          else if(j>1)then
            parents(n)%values=tree%level(j-1)%node(n)%values
          end if
        end do
        call move_alloc(parents,nodes)
        call move_alloc(pmask,mask)
      end do
      x=nodes(1)%values
    end block
    nout=size(x)
    if(trim(tree%boundary)=='reflection')nout=tree%original_length
    if(tree%original_length>0 .and. tree%original_length<nout)nout=tree%original_length
    x=x(1:nout)
  end function idwpt

  function packet_basis(tree,levels,nodes) result(mask)
    type(packet_transform),intent(in)::tree
    integer,intent(in)::levels(:),nodes(:)
    logical,allocatable::mask(:)
    integer::total,i,idx
    total=2**(tree%levels()+1)-2
    allocate(mask(total))
    mask=.false.
    do i=1,min(size(levels),size(nodes))
      if(levels(i)>=1 .and. levels(i)<=tree%levels())then
        idx=2**levels(i)-2+nodes(i)+1
        if(nodes(i)>=0 .and. nodes(i)<2**levels(i))mask(idx)=.true.
      end if
    end do
  end function packet_basis

  function ortho_basis(tree_test) result(basis_mask)
    logical,intent(in)::tree_test(:)
    logical,allocatable::basis_mask(:)
    integer::total,jmax,j,n,parent,idx,pidx
    total=size(tree_test)
    jmax=nint(log(real(total+2,dp))/log(2.0_dp))-1
    allocate(basis_mask(total))
    basis_mask=.false.
    do j=jmax,1,-1
      do n=0,2**j-1
        idx=2**j-1+n
        if(tree_test(idx))then
          basis_mask(idx)=.true.
          parent=n/2
          do while(j>1)
            pidx=2**(j-1)-1+parent
            basis_mask(pidx)=.false.
            parent=parent/2
            exit
          end do
        end if
      end do
    end do
    do n=0,2**jmax-1
      idx=2**jmax-1+n
      if(.not.any_selected_ancestor(basis_mask,jmax,n))basis_mask(idx)=.true.
    end do
  end function ortho_basis

  logical function any_selected_ancestor(mask,level,node)
    logical,intent(in)::mask(:)
    integer,intent(in)::level,node
    integer::j,n,idx
    any_selected_ancestor=.false.
    n=node
    do j=level,1,-1
      idx=2**j-1+n
      if(mask(idx))then
      any_selected_ancestor=.true.
      return
      end if
      n=n/2
    end do
  end function any_selected_ancestor

  subroutine dwpt_brick_wall(tree,wf,method)
    type(packet_transform),intent(inout)::tree
    character(len=*),intent(in),optional::wf,method
    character(len=16)::wname,mname
    integer::j,n,l,m
    real(dp)::nan
    block
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      wname=tree%wavelet
      if(present(wf))wname=trim(wf)
      mname=tree%method
      if(present(method))mname=trim(method)
      f=wave_filter(wname,st)
      if(.not.st%ok())return
      m=f%length()
    end block
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    do j=1,tree%levels()
      if(mname=='dwpt')then
        select case(j)
        case(1)
        l=(m-2)/2
        case(2)
        l=(m-2)/2+m/4
        case(3)
        l=(m-2)/2+(m/2+m/4)/2
        case default
        l=m-2
        end select
      else
      l=min((2**j-1)*(m-1),size(tree%level(j)%node(1)%values))
      end if
      do n=1,size(tree%level(j)%node)
      if(l>0)tree%level(j)%node(n)%values(1:min(l,size(tree%level(j)%node(n)%values)))=nan
      end do
    end do
  end subroutine dwpt_brick_wall

  function css_test(tree) result(pass)
    type(packet_transform),intent(in)::tree
    logical,allocatable::pass(:)
    integer::total,j,n,k,nn,i
    real(dp),allocatable::x(:)
    real(dp)::ss,d,plus,minus,crit,cum
    total=2**(tree%levels()+1)-2
    allocate(pass(total))
    pass=.false.
    k=0
    do j=1,tree%levels()
    do n=1,2**j
    k=k+1
    x=pack(tree%level(j)%node(n)%values,ieee_is_finite(tree%level(j)%node(n)%values))
    nn=size(x)
    if(nn<2)cycle
    ss=sum(x*x)
    if(ss<=0)cycle
    d=0
    cum=0
    do i=1,nn
    cum=cum+x(i)**2
    plus=real(i,dp)/real(nn-1,dp)-cum/ss
    minus=cum/ss-real(i-1,dp)/real(nn-1,dp)
    d=max(d,abs(plus),abs(minus))
    end do
    crit=1.224_dp/(sqrt(real(nn,dp))+0.12_dp+0.11_dp/sqrt(real(nn,dp)))
    pass(k)=d<crit
    end do
    end do
  end function css_test

  function entropy_test(tree) result(entropy)
    type(packet_transform),intent(in)::tree
    real(dp),allocatable::entropy(:)
    integer::total,j,n,k,i
    real(dp)::v
    total=2**(tree%levels()+1)-2
    allocate(entropy(total))
    k=0
    do j=1,tree%levels()
    do n=1,2**j
    k=k+1
    entropy(k)=0.0_dp
    do i=1,size(tree%level(j)%node(n)%values)
    v=tree%level(j)%node(n)%values(i)**2
    if(v>0 .and. ieee_is_finite(v))entropy(k)=entropy(k)+v*log(v)
    end do
    end do
    end do
  end function entropy_test

  function cpgram_test(tree,p,taper) result(pass)
    type(packet_transform),intent(in)::tree
    real(dp),intent(in),optional::p,taper
    logical,allocatable::pass(:)
    integer::total,j,n,k,nn,m,i
    real(dp),allocatable::x(:),per(:)
    complex(dp),allocatable::z(:)
    real(dp)::alpha,tap,crit,cum,ss,d,w
    alpha=0.05_dp
    if(present(p))alpha=p
    tap=0.1_dp
    if(present(taper))tap=taper
    total=2**(tree%levels()+1)-2
    allocate(pass(total))
    pass=.false.
    k=0
    do j=1,tree%levels()
    do n=1,2**j
    k=k+1
    x=pack(tree%level(j)%node(n)%values,ieee_is_finite(tree%level(j)%node(n)%values))
    nn=size(x)
    if(nn<8)cycle
    x=x-mean_value(x)
    m=nint(tap*real(nn,dp))
    if(m>0)then
    do i=1,m
    w=0.5_dp*(1.0_dp-cos(acos(-1.0_dp)*real(i,dp)/real(m+1,dp)))
    x(i)=x(i)*w
    x(nn-i+1)=x(nn-i+1)*w
    end do
    end if
    allocate(z(nn))
    z=cmplx(x,0.0_dp,dp)
    if(iand(nn,nn-1)/=0)then
    deallocate(z)
    cycle
    end if
    call fft_complex(z)
    allocate(per(nn/2))
    per=abs(z(2:nn/2+1))**2/real(nn,dp)
    ss=sum(per)
    if(ss<=0)then
    deallocate(z,per)
    cycle
    end if
    d=0
    cum=0
    do i=1,size(per)
    cum=cum+per(i)
    d=max(d,abs(cum/ss-real(i,dp)/real(size(per),dp)))
    end do
    if(alpha<=0.01_dp)then
    crit=1.628_dp/(sqrt(real(size(per),dp))+0.12_dp+0.11_dp/sqrt(real(size(per),dp)))
    else
    crit=1.358_dp/(sqrt(real(size(per),dp))+0.12_dp+0.11_dp/sqrt(real(size(per),dp)))
    end if
    pass(k)=d<crit
    deallocate(z,per)
    end do
    end do
  end function cpgram_test

  function portmanteau_test(tree,p,ljung_box) result(pass)
    type(packet_transform),intent(in)::tree
    real(dp),intent(in),optional::p
    logical,intent(in),optional::ljung_box
    logical,allocatable::pass(:)
    integer::total,j,n,k,nn,h,l
    real(dp),allocatable::x(:)
    real(dp)::alpha,q,crit
    logical::lb
    alpha=0.05_dp
    if(present(p))alpha=p
    lb=.false.
    if(present(ljung_box))lb=ljung_box
    total=2**(tree%levels()+1)-2
    allocate(pass(total))
    pass=.false.
    k=0
    do j=1,tree%levels()
    do n=1,2**j
    k=k+1
    x=pack(tree%level(j)%node(n)%values,ieee_is_finite(tree%level(j)%node(n)%values))
    nn=size(x)
    h=nn/2
    if(h<1)cycle
    q=0.0_dp
    do l=1,h
    if(lb)then
    q=q+autocorrelation(x,l)**2/real(nn-l,dp)
    else
    q=q+autocorrelation(x,l)**2
    end if
    end do
    if(lb)then
    q=real(nn*(nn+2),dp)*q
    else
    q=real(nn,dp)*q
    end if
    crit=chi_square_quantile(1.0_dp-alpha,real(h,dp))
    pass(k)=q<=crit
    end do
    end do
  end function portmanteau_test

  function chi_square_quantile(p,df) result(x)
    real(dp),intent(in)::p,df
    real(dp)::x,lo,hi,mid
    integer::i
    lo=0.0_dp
    hi=max(df+10.0_dp*sqrt(2.0_dp*df),10.0_dp)
    do while(chi_square_cdf(hi,df)<p)
    hi=2.0_dp*hi
    end do
    do i=1,100
    mid=0.5_dp*(lo+hi)
    if(chi_square_cdf(mid,df)<p)then
    lo=mid
    else
    hi=mid
    end if
    end do
    x=0.5_dp*(lo+hi)
  end function chi_square_quantile
end module waveslim_packet
