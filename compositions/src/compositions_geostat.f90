! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_geostat
  use compositions_kinds, only: dp, pi
  use compositions_geometry, only: closure
  use compositions_linalg, only: invert_matrix
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: variogram_result, logratio_variogram, vg_to_lrvg
  public :: vgram_spherical, vgram_exponential, vgram_gaussian, vgram_linear, vgram_power, vgram_nugget, vgram_cardsin
  public :: compositional_ordinary_kriging, compositional_general_kriging

  type :: variogram_result
    real(dp), allocatable :: vg(:,:,:)
    real(dp), allocatable :: h(:,:,:)
    integer, allocatable :: count(:,:,:)
  end type
contains
  function logratio_variogram(comp,loc,bins,azimuth,azimuth_tol) result(res)
    real(dp), intent(in) :: comp(:,:),loc(:,:),bins(:,:)
    real(dp), intent(in), optional :: azimuth(:),azimuth_tol
    type(variogram_result) :: res
    integer :: n,d,s,nb,i,j,k,l,b
    real(dp) :: distv,delta(size(loc,2)),proj,ctol,a1,a2,b1,b2
    logical :: directional
    n=size(comp,1); d=size(comp,2); s=size(loc,2); nb=size(bins,1)
    if(size(loc,1)/=n.or.size(bins,2)/=2) error stop 'logratio_variogram: dimension mismatch'
    if(any(comp<=0.0_dp)) error stop 'logratio_variogram: positive complete compositions required'
    allocate(res%vg(nb,d,d),res%h(nb,d,d),res%count(nb,d,d)); res%vg=0.0_dp; res%h=0.0_dp; res%count=0
    directional=present(azimuth)
    ctol=-1.0_dp
    if(present(azimuth_tol)) ctol=cos(azimuth_tol*pi/180.0_dp)
    do i=1,n
      do j=1,n
        delta=loc(i,:)-loc(j,:); distv=sqrt(sum(delta*delta))
        if(directional.and.distv>0.0_dp) then
          if(size(azimuth)/=s) error stop 'logratio_variogram: azimuth dimension mismatch'
          proj=dot_product(delta,azimuth)/distv
          if(proj<ctol) cycle
        end if
        do b=1,nb
          if(distv>bins(b,1).and.distv<=bins(b,2)) then
            do k=1,d
              a1=log(comp(i,k)); a2=log(comp(j,k))
              do l=1,d
                b1=log(comp(i,l)); b2=log(comp(j,l))
                res%count(b,k,l)=res%count(b,k,l)+1
                res%vg(b,k,l)=res%vg(b,k,l)+((a1-b1)-(a2-b2))**2
                res%h(b,k,l)=res%h(b,k,l)+distv
              end do
            end do
          end if
        end do
      end do
    end do
    do b=1,nb; do k=1,d; do l=1,d
      if(res%count(b,k,l)>0) then
        res%vg(b,k,l)=res%vg(b,k,l)/real(res%count(b,k,l),dp)
        res%h(b,k,l)=res%h(b,k,l)/real(res%count(b,k,l),dp)
      else
        res%vg(b,k,l)=0.0_dp
        res%h(b,k,l)=0.0_dp
      end if
    end do; end do; end do
  end function logratio_variogram

  function vg_to_lrvg(vg) result(lr)
    real(dp), intent(in) :: vg(:,:,:)
    real(dp) :: lr(size(vg,1),size(vg,2),size(vg,3))
    integer :: b,i,j
    if(size(vg,2)/=size(vg,3)) error stop 'vg_to_lrvg: last dimensions must match'
    do b=1,size(vg,1); do i=1,size(vg,2); do j=1,size(vg,3)
      lr(b,i,j)=vg(b,i,i)+vg(b,j,j)-vg(b,i,j)-vg(b,j,i)
    end do; end do; end do
  end function vg_to_lrvg

  elemental real(dp) function vgram_spherical(h,nugget,sill,range) result(v)
    real(dp), intent(in) :: h,nugget,sill,range
    real(dp) :: s,q
    if(h<=range*1.0e-8_dp) then; v=0.0_dp; return; end if
    s=sill-nugget; q=h/range
    if(h>range) then; v=nugget+s; else; v=nugget+s*(1.5_dp*q-0.5_dp*q**3); end if
  end function vgram_spherical

  elemental real(dp) function vgram_exponential(h,nugget,sill,range) result(v)
    real(dp), intent(in) :: h,nugget,sill,range
    real(dp) :: s,r
    if(h<=range*1.0e-8_dp) then; v=0.0_dp; return; end if
    s=sill-nugget; r=-range/log(0.05_dp); v=nugget+s*(1.0_dp-exp(-h/r))
  end function vgram_exponential

  elemental real(dp) function vgram_gaussian(h,nugget,sill,range) result(v)
    real(dp), intent(in) :: h,nugget,sill,range
    real(dp) :: s,r
    if(h<=range*1.0e-8_dp) then; v=0.0_dp; return; end if
    s=sill-nugget; r=range/sqrt(-log(0.05_dp)); v=nugget+s*(1.0_dp-exp(-(h/r)**2))
  end function vgram_gaussian

  elemental real(dp) function vgram_linear(h,nugget,sill,range) result(v)
    real(dp), intent(in) :: h,nugget,sill,range
    if(h<=range*1.0e-8_dp) then; v=0.0_dp; else; v=nugget+h*(sill-nugget)/range; end if
  end function vgram_linear

  elemental real(dp) function vgram_power(h,nugget,sill,range) result(v)
    real(dp), intent(in) :: h,nugget,sill,range
    if(h<=range*1.0e-8_dp) then; v=0.0_dp; else; v=nugget+sill*h**range; end if
  end function vgram_power

  elemental real(dp) function vgram_nugget(h,nugget,tol) result(v)
    real(dp), intent(in) :: h,nugget,tol
    if(h>tol) then; v=nugget; else; v=0.0_dp; end if
  end function vgram_nugget

  elemental real(dp) function vgram_cardsin(h,nugget,sill,range) result(v)
    real(dp), intent(in) :: h,nugget,sill,range
    real(dp) :: s,r
    if(h<=range*1.0e-8_dp) then; v=0.0_dp; return; end if
    s=sill-nugget; r=-range/log(0.05_dp); v=nugget+s*(1.0_dp-r*sin(h/r)/h)
  end function vgram_cardsin

  function compositional_ordinary_kriging(z,gamma_obs,gamma_new) result(pred)
    real(dp), intent(in) :: z(:,:),gamma_obs(:,:,:,:),gamma_new(:,:,:,:)
    real(dp), allocatable :: pred(:,:)
    real(dp), allocatable :: nmat(:,:),ninv(:,:),w(:),coef(:),rhs(:,:),lp(:)
    integer :: n,d,np,m,i,j,p,q,row,col,tr
    n=size(z,1); d=size(z,2); np=size(gamma_new,2)
    if(any(z<=0.0_dp)) error stop 'compositional_ordinary_kriging: complete positive compositions required'
    if(any(shape(gamma_obs)/=[n,n,d,d])) error stop 'compositional_ordinary_kriging: gamma_obs shape mismatch'
    if(size(gamma_new,1)/=n.or.size(gamma_new,3)/=d.or.size(gamma_new,4)/=d) &
      error stop 'compositional_ordinary_kriging: gamma_new shape mismatch'
    m=(d-1)*(n+1); allocate(nmat(m,m)); nmat=0.0_dp
    do i=1,n; do p=1,d-1
      row=(i-1)*(d-1)+p
      do j=1,n; do q=1,d-1
        col=(j-1)*(d-1)+q
        nmat(row,col)=gamma_obs(i,j,p,q)+gamma_obs(i,j,d,d)-gamma_obs(i,j,d,q)-gamma_obs(i,j,p,d)
      end do; end do
      do q=1,d-1
        col=n*(d-1)+q; nmat(row,col)=merge(2.0_dp,1.0_dp,p==q); nmat(col,row)=nmat(row,col)
      end do
    end do; end do
    call invert_matrix(nmat,ninv)
    allocate(w(m)); w=0.0_dp
    do i=1,n; do p=1,d-1; row=(i-1)*(d-1)+p; w(row)=log(z(i,p)/z(i,d)); end do; end do
    coef=matmul(ninv,w); allocate(pred(np,d),rhs(m,d),lp(d))
    do j=1,np
      rhs=0.0_dp
      do i=1,n; do p=1,d-1
        row=(i-1)*(d-1)+p
        do q=1,d-1
          rhs(row,q)=gamma_new(i,j,p,q)+gamma_new(i,j,d,d)-gamma_new(i,j,d,q)-gamma_new(i,j,p,d)
        end do
      end do; end do
      do p=1,d-1; row=n*(d-1)+p; do q=1,d-1; rhs(row,q)=merge(2.0_dp,1.0_dp,p==q); end do; end do
      lp=matmul(transpose(rhs),coef); pred(j,:)=closure(exp(lp-maxval(lp)))
    end do
  end function compositional_ordinary_kriging

  function compositional_general_kriging(z,f,gamma_obs,fnew,gamma_new,krig_var,source_error) result(pred)
    !! Source-style gsiCGSkriging with arbitrary trends and partially observed compositions.
    !! A component is usable exactly when it is finite and strictly positive; BDL/zero,
    !! structural-zero and R-style missing encodings are omitted from that row's ALR.
    real(dp), intent(in) :: z(:,:),f(:,:),gamma_obs(:,:,:,:),fnew(:,:),gamma_new(:,:,:,:)
    real(dp), allocatable, intent(out), optional :: krig_var(:,:,:)
    logical, intent(in), optional :: source_error
    real(dp), allocatable :: pred(:,:),nmat(:,:),ninv(:,:),w(:),rhs(:,:),coef(:),ealr(:,:),tmat(:,:)
    integer, allocatable :: nmv(:),ref(:),idx(:,:)
    integer :: n,d,fd,np,i,j,p,q,row,col,k,m,iref,jref,iidx,jidx,info
    real(dp) :: lr,meanerr
    logical :: srcerr

    n=size(z,1); d=size(z,2); fd=size(f,2); np=size(fnew,1)
    if(size(f,1)/=n) error stop 'compositional_general_kriging: F row mismatch'
    if(size(fnew,2)/=fd) error stop 'compositional_general_kriging: trend column mismatch'
    if(any(shape(gamma_obs)/=[n,n,d,d])) error stop 'compositional_general_kriging: gamma_obs shape mismatch'
    if(any(shape(gamma_new)/=[n,np,d,d])) error stop 'compositional_general_kriging: gamma_new shape mismatch'
    allocate(nmv(n),ref(n),idx(n,d)); nmv=0; ref=1; idx=0
    do i=1,n
      do p=1,d
        if(ieee_is_finite(z(i,p)).and.z(i,p)>0.0_dp) then
          nmv(i)=nmv(i)+1; idx(i,nmv(i))=p; ref(i)=p
        end if
      end do
    end do
    m=(d-1)*fd+sum(max(0,nmv-1))
    if(m<=0) error stop 'compositional_general_kriging: no estimable log-ratios'
    allocate(nmat(m,m)); nmat=0.0_dp
    row=0
    do i=1,n
      iref=ref(i)
      do p=1,max(0,nmv(i)-1)
        iidx=idx(i,p); row=row+1; col=0
        do j=1,n
          jref=ref(j)
          do q=1,max(0,nmv(j)-1)
            jidx=idx(j,q); col=col+1
            nmat(row,col)=gamma_obs(i,j,iidx,jidx)+gamma_obs(i,j,iref,jref) &
              -gamma_obs(i,j,iref,jidx)-gamma_obs(i,j,iidx,jref)
          end do
        end do
        do j=1,fd
          jref=d
          do q=1,d-1
            jidx=q; col=col+1
            nmat(row,col)=f(i,j)*(indicator(iidx,jidx)+indicator(iref,jref) &
              -indicator(iidx,jref)-indicator(iref,jidx))
          end do
        end do
      end do
    end do
    do i=1,fd
      iref=d
      do p=1,d-1
        iidx=p; row=row+1; col=0
        do j=1,n
          jref=ref(j)
          do q=1,max(0,nmv(j)-1)
            jidx=idx(j,q); col=col+1
            nmat(row,col)=f(j,i)*(indicator(iidx,jidx)+indicator(iref,jref) &
              -indicator(iidx,jref)-indicator(iref,jidx))
          end do
        end do
        do j=1,fd
          do q=1,d-1
            col=col+1; nmat(row,col)=0.0_dp
          end do
        end do
      end do
    end do
    call invert_matrix(nmat,ninv,info)
    if(info/=0) error stop 'compositional_general_kriging: kriging matrix is singular'
    allocate(w(m)); w=0.0_dp; k=0
    do i=1,n
      if(nmv(i)>0) lr=log(z(i,ref(i)))
      do p=1,max(0,nmv(i)-1)
        k=k+1; w(k)=log(z(i,idx(i,p)))-lr
      end do
    end do
    coef=matmul(ninv,w)
    allocate(pred(np,d),rhs(m,d)); pred=0.0_dp
    if(present(krig_var)) then
      allocate(krig_var(np,d,d)); krig_var=0.0_dp
      allocate(tmat(d,d-1)); tmat=-1.0_dp/real(d,dp)
      do p=1,d-1; tmat(p,p)=tmat(p,p)+1.0_dp; end do
    end if
    srcerr=.false.; if(present(source_error)) srcerr=source_error
    do j=1,np
      rhs=0.0_dp; k=0; jref=d
      do i=1,n
        iref=ref(i)
        do p=1,max(0,nmv(i)-1)
          iidx=idx(i,p); k=k+1
          do q=1,d-1
            jidx=q
            rhs(k,q)=gamma_new(i,j,iidx,jidx)+gamma_new(i,j,iref,jref) &
              -gamma_new(i,j,iref,jidx)-gamma_new(i,j,iidx,jref)
          end do
          rhs(k,d)=0.0_dp
        end do
      end do
      do i=1,fd
        iref=d
        do p=1,d-1
          iidx=p; k=k+1
          do q=1,d-1
            jidx=q
            rhs(k,q)=fnew(j,i)*(indicator(iidx,jidx)+indicator(iref,jref) &
              -indicator(iidx,jref)-indicator(iref,jidx))
          end do
          rhs(k,d)=0.0_dp
        end do
      end do
      pred(j,:)=matmul(transpose(rhs),coef)
      pred(j,:)=closure(exp(pred(j,:)-maxval(pred(j,:))))
      if(present(krig_var)) then
        ealr=matmul(transpose(rhs),matmul(ninv,rhs))
        if(srcerr) then
          meanerr=sum(ealr)/real(d*d,dp)
          krig_var(j,:,:)=ealr+meanerr
        else
          krig_var(j,:,:)=matmul(tmat,matmul(ealr(1:d-1,1:d-1),transpose(tmat)))
        end if
      end if
    end do
  contains
    pure real(dp) function indicator(a,b) result(v)
      integer, intent(in) :: a,b
      if(a==b) then; v=1.0_dp; else; v=0.0_dp; end if
    end function indicator
  end function compositional_general_kriging

end module compositions_geostat
