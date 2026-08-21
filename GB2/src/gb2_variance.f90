! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_variance
  use gb2_kinds, only : dp
  use gb2_likelihood, only : dlogf_gb2, d2logf_gb2
  use gb2_optimizer, only : invert_matrix
  use gb2_indicators, only : main_gb2
  use gb2_compound, only : vofp_cgb2, pofv_cgb2
  use gb2_compound_indicators, only : main_cgb2
  use survey_types, only : survey_design_t, svystat_t
  use survey_estimators, only : svy_total
  implicit none
  private
  public :: varscore_gb2, vepar_gb2, derivind_gb2, veind_gb2
  public :: varscore_mixture, survey_score_variance, vepar_mixture
  public :: derivind_cgb2, veind_cgb2, veind_cavgb2_groups
contains
  subroutine varscore_gb2(x,par,vsc,w,hs)
    real(dp), intent(in) :: x(:),par(4)
    real(dp), intent(out) :: vsc(4,4)
    real(dp), intent(in), optional :: w(:),hs(:)
    real(dp) :: s(4),wi,hi
    integer :: i,j,k
    vsc=0.0_dp
    do k=1,size(x)
      wi=1.0_dp
      hi=1.0_dp
      if(present(w)) wi=w(k)
      if(present(hs)) hi=hs(k)
      call dlogf_gb2(x(k),par(1),par(2),par(3),par(4),s)
      do j=1,4
      do i=1,4
      vsc(i,j)=vsc(i,j)+(wi*hi)**2*s(i)*s(j)
      end do
      end do
    end do
  end subroutine varscore_gb2

  subroutine vepar_gb2(x,vsc,par,vcov,w,hs,wdm,ok)
    real(dp), intent(in) :: x(:),vsc(4,4),par(4)
    real(dp), intent(out) :: vcov(4,4)
    real(dp), intent(in), optional :: w(:),hs(:)
    real(dp), intent(out), optional :: wdm(4,4)
    logical, intent(out), optional :: ok
    real(dp) :: a(4,4),ai(4,4),d(4,4),wi,hi
    integer :: k
    logical :: good
    a=0.0_dp
    do k=1,size(x)
      wi=1.0_dp
      hi=1.0_dp
      if(present(w)) wi=w(k)
      if(present(hs)) hi=hs(k)
      call d2logf_gb2(x(k),par(1),par(2),par(3),par(4),d)
      a=a-wi*hi*d
    end do
    call invert_matrix(a,ai,good)
    if(good) then
    vcov=matmul(ai,matmul(vsc,transpose(ai)))
    else
    vcov=0.0_dp
    end if
    if(present(wdm)) wdm=a
    if(present(ok)) ok=good
  end subroutine vepar_gb2

  subroutine derivind_gb2(par,jac,prop)
    real(dp), intent(in) :: par(4)
    real(dp), intent(out) :: jac(6,4)
    real(dp), intent(in), optional :: prop
    real(dp) :: pp(4),pm(4),yp(6),ym(6),h,pr
    integer :: j
    pr=0.6_dp
    if(present(prop)) pr=prop
    do j=1,4
      h=1.0e-5_dp*max(1.0_dp,abs(par(j)))
      h=min(h,0.25_dp*par(j))
      pp=par
      pm=par
      pp(j)=par(j)+h
      pm(j)=par(j)-h
      call main_gb2(pr,pp(1),pp(2),pp(3),pp(4),yp)
      call main_gb2(pr,pm(1),pm(2),pm(3),pm(4),ym)
      jac(:,j)=(yp-ym)/(2.0_dp*h)
    end do
  end subroutine derivind_gb2

  subroutine veind_gb2(vpar,par,vcov,prop,jac)
    real(dp), intent(in) :: vpar(4,4),par(4)
    real(dp), intent(out) :: vcov(6,6)
    real(dp), intent(in), optional :: prop
    real(dp), intent(out), optional :: jac(6,4)
    real(dp) :: j(6,4)
    call derivind_gb2(par,j,prop)
    vcov=matmul(j,matmul(vpar,transpose(j)))
    if(present(jac)) jac=j
  end subroutine veind_gb2

  subroutine varscore_mixture(u,vsc,w)
    real(dp), intent(in) :: u(:,:)
    real(dp), intent(out) :: vsc(:,:)
    real(dp), intent(in), optional :: w(:)
    real(dp) :: wi
    integer :: i,j,k,p
    p=size(u,2)
    if(any(shape(vsc)/=[p,p])) error stop 'varscore_mixture: shape mismatch'
    vsc=0.0_dp
    do k=1,size(u,1)
      wi=1.0_dp
      if(present(w)) wi=w(k)
      do j=1,p
        do i=1,p
          vsc(i,j)=vsc(i,j)+(wi*u(k,i))*(wi*u(k,j))
        end do
      end do
    end do
  end subroutine varscore_mixture

  subroutine survey_score_variance(u,design,vsc)
    real(dp), intent(in) :: u(:,:)
    type(survey_design_t), intent(in) :: design
    real(dp), intent(out) :: vsc(:,:)
    type(svystat_t) :: stat
    stat=svy_total(u,design)
    if(any(shape(vsc)/=shape(stat%variance))) &
      error stop 'survey_score_variance: shape mismatch'
    vsc=stat%variance
  end subroutine survey_score_variance

  subroutine vepar_mixture(vsc,hess,vcov,ok)
    real(dp), intent(in) :: vsc(:,:),hess(:,:)
    real(dp), intent(out) :: vcov(:,:)
    logical, intent(out), optional :: ok
    real(dp), allocatable :: hi(:,:)
    logical :: good
    allocate(hi(size(hess,1),size(hess,2)))
    call invert_matrix(hess,hi,good)
    if(good) then
    vcov=matmul(hi,matmul(vsc,transpose(hi)))
    else
    vcov=0.0_dp
    end if
    if(present(ok)) ok=good
  end subroutine vepar_mixture

  subroutine derivind_cgb2(shape1,scale,shape2,shape3,pl0,pl,jac,prop,decomp)
    real(dp), intent(in) :: shape1,scale,shape2,shape3,pl0(:),pl(:)
    real(dp), intent(out) :: jac(:,:)
    real(dp), intent(in), optional :: prop
    character(len=*), intent(in), optional :: decomp
    real(dp), allocatable :: v(:),vp(:),vm(:),pp(:),pm(:)
    real(dp) :: yp(5),ym(5),h,pr
    integer :: j,p
    p=size(pl)-1
    if(any(shape(jac)/=[5,p])) error stop 'derivind_cgb2: shape mismatch'
    allocate(v(p),vp(p),vm(p),pp(size(pl)),pm(size(pl)))
    call vofp_cgb2(pl,v)
    pr=0.6_dp
    if(present(prop)) pr=prop
    do j=1,p
      h=1.0e-5_dp*max(1.0_dp,abs(v(j)))
      vp=v
      vm=v
      vp(j)=v(j)+h
      vm(j)=v(j)-h
      call pofv_cgb2(vp,pp)
      call pofv_cgb2(vm,pm)
      call main_cgb2(pr,shape1,scale,shape2,shape3,pl0,pp,yp,decomp)
      call main_cgb2(pr,shape1,scale,shape2,shape3,pl0,pm,ym,decomp)
      jac(:,j)=(yp-ym)/(2.0_dp*h)
    end do
  end subroutine derivind_cgb2

  subroutine veind_cgb2(vpar,shape1,scale,shape2,shape3,pl0,pl,vcov,prop,decomp,jac)
    real(dp), intent(in) :: vpar(:,:),shape1,scale,shape2,shape3,pl0(:),pl(:)
    real(dp), intent(out) :: vcov(5,5)
    real(dp), intent(in), optional :: prop
    character(len=*), intent(in), optional :: decomp
    real(dp), intent(out), optional :: jac(:,:)
    real(dp), allocatable :: j(:,:)
    allocate(j(5,size(pl)-1))
    call derivind_cgb2(shape1,scale,shape2,shape3,pl0,pl,j,prop,decomp)
    vcov=matmul(j,matmul(vpar,transpose(j)))
    if(present(jac)) jac=j
  end subroutine veind_cgb2
  subroutine veind_cavgb2_groups(vpar,group_pl,shape1,scale,shape2,shape3,pl0,estimates,stderr,vcovs,decomp)
    real(dp), intent(in) :: vpar(:,:),group_pl(:,:),shape1,scale,shape2,shape3,pl0(:)
    real(dp), intent(out) :: estimates(:,:),stderr(:,:),vcovs(:,:,:)
    character(len=*), intent(in), optional :: decomp
    real(dp), allocatable :: j(:,:),vg(:,:)
    real(dp) :: vals(5)
    integer :: k,l1,jj,ii,ng
    ng=size(group_pl,1)
    l1=size(group_pl,2)-1
    if(any(shape(estimates)/=[5,ng]) .or. any(shape(stderr)/=[5,ng])) error stop 'veind_cavgb2_groups: output shape mismatch'
    if(any(shape(vcovs)/=[5,5,ng])) error stop 'veind_cavgb2_groups: vcovs shape mismatch'
    if(any(shape(vpar)/=[ng*l1,ng*l1])) error stop 'veind_cavgb2_groups: parameter covariance assumes group-dummy auxiliary design'
    allocate(j(5,l1),vg(l1,l1))
    do k=1,ng
      call main_cgb2(0.6_dp,shape1,scale,shape2,shape3,pl0,group_pl(k,:),vals,decomp)
      estimates(:,k)=vals
      call derivind_cgb2(shape1,scale,shape2,shape3,pl0,group_pl(k,:),j,decomp=decomp)
      do jj=1,l1
      do ii=1,l1
      vg(ii,jj)=vpar(k+(ii-1)*ng,k+(jj-1)*ng)
      end do
      end do
      vcovs(:,:,k)=matmul(j,matmul(vg,transpose(j)))
      stderr(:,k)=sqrt(max(0.0_dp,[(vcovs(ii,ii,k),ii=1,5)]))
    end do
  end subroutine veind_cavgb2_groups

end module gb2_variance
