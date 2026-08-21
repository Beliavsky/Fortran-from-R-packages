! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_taylor
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, LONELY_FAIL, LONELY_REMOVE, LONELY_CERTAINTY, LONELY_ADJUST, LONELY_AVERAGE
  use survey_linalg, only : outer_product
  implicit none
  private
  public :: svyrecvar, onestage_variance
contains

  function svyrecvar(x,design) result(v)
    real(dp), intent(in) :: x(:,:)
    type(survey_design_t), intent(in) :: design
    real(dp) :: v(size(x,2),size(x,2))
    if(size(x,1)/=design%n) error stop 'svyrecvar: row count mismatch'
    v=multistage_core(x,design%cluster,design%strata,design%samp_size,design%pop_size, &
                      design%lonely_psu,design%ultimate_cluster,design%adjust_domain_lonely,1)
  end function svyrecvar

  recursive function multistage_core(x,clusters,stratas,npsus,popsizes,lonely,ultimate,adjust_domain,stage) result(v)
    real(dp), intent(in) :: x(:,:), npsus(:,:), popsizes(:,:)
    integer, intent(in) :: clusters(:,:), stratas(:,:)
    integer, intent(in) :: lonely,stage
    logical, intent(in) :: ultimate,adjust_domain
    real(dp) :: v(size(x,2),size(x,2))
    real(dp), allocatable :: subx(:,:), subn(:,:), subp(:,:), vsub(:,:)
    integer, allocatable :: subc(:,:), subs(:,:), cu(:), idx(:)
    integer :: h,i,m,n,k,p
    real(dp) :: fac
    n=size(x,1); p=size(x,2); k=size(clusters,2)
    v=onestage_variance(x,stratas(:,1),clusters(:,1),npsus(:,1),popsizes(:,1),lonely,adjust_domain,stage)
    if(ultimate .or. k<=1) return
    ! R survey only recurses when an FPC is supplied. Treat huge pop sizes as no FPC.
    if(all(popsizes(:,1)>=huge(1.0_dp)/4.0_dp)) return
    call unique_int(clusters(:,1),cu)
    do h=1,size(cu)
      m=count(clusters(:,1)==cu(h))
      allocate(idx(m)); idx=pack([(i,i=1,n)],clusters(:,1)==cu(h))
      allocate(subx(m,p),subc(m,k-1),subs(m,k-1),subn(m,k-1),subp(m,k-1),vsub(p,p))
      subx=x(idx,:); subc=clusters(idx,2:k); subs=stratas(idx,2:k)
      subn=npsus(idx,2:k); subp=popsizes(idx,2:k)
      vsub=multistage_core(subx,subc,subs,subn,subp,lonely,.false.,adjust_domain,stage+1)
      if(popsizes(idx(1),1)>0 .and. popsizes(idx(1),1)<huge(1.0_dp)/4.0_dp) then
        fac=npsus(idx(1),1)/popsizes(idx(1),1)
      else
        fac=0.0_dp
      end if
      v=v+fac*vsub
      deallocate(idx,subx,subc,subs,subn,subp,vsub)
    end do
  end function multistage_core

  function onestage_variance(x,strata,cluster,npsu,popsize,lonely,adjust_domain,stage) result(v)
    real(dp), intent(in) :: x(:,:), npsu(:), popsize(:)
    integer, intent(in) :: strata(:),cluster(:),lonely,stage
    logical, intent(in) :: adjust_domain
    real(dp) :: v(size(x,2),size(x,2))
    integer, allocatable :: su(:), cu(:), idxs(:), idxc(:)
    real(dp), allocatable :: psutot(:,:), center(:)
    real(dp) :: f,scale,nh,poph
    integer :: p,n,h,j,i,m,nobs,nstrata,nok
    logical :: singleton,ok
    p=size(x,2); n=size(x,1); v=0
    if(n==0) return
    call unique_int(strata,su); nstrata=size(su); nok=0
    allocate(center(p)); center=0
    if(lonely==LONELY_ADJUST) then
      nh=0
      do h=1,nstrata
        i=find_first(strata,su(h)); nh=nh+npsu(i)
      end do
      if(nh>0) center=sum(x,dim=1)/nh
    end if
    do h=1,nstrata
      m=count(strata==su(h)); allocate(idxs(m)); idxs=pack([(i,i=1,n)],strata==su(h))
      i=idxs(1); nh=npsu(i); poph=popsize(i)
      call unique_int(cluster(idxs),cu); nobs=size(cu)
      singleton=(nobs==1 .and. nh<=1.0_dp+epsilon(1.0_dp))
      if(adjust_domain .and. nobs==1 .and. nh>1.0_dp) singleton=.true.
      if(singleton .and. lonely==LONELY_FAIL) then
        write(*,'(a,i0)') 'lonely PSU at stage ',stage
        error stop 'survey_taylor: lonely PSU'
      end if
      if(singleton .and. lonely==LONELY_AVERAGE) then
        deallocate(idxs,cu); cycle
      end if
      allocate(psutot(max(nobs,nint(nh)),p)); psutot=0
      do j=1,nobs
        idxc=pack(idxs,cluster(idxs)==cu(j))
        psutot(j,:)=sum(x(idxc,:),dim=1)
        deallocate(idxc)
      end do
      if(singleton .and. lonely==LONELY_ADJUST) then
        do j=1,size(psutot,1); psutot(j,:)=psutot(j,:)-center; end do
      else
        center=sum(psutot,dim=1)/max(nh,1.0_dp)
        do j=1,size(psutot,1); psutot(j,:)=psutot(j,:)-center; end do
      end if
      if(poph>=huge(1.0_dp)/4.0_dp) then
        f=1.0_dp
      else if(poph>0) then
        f=max(0.0_dp,(poph-nh)/poph)
      else
        f=1.0_dp
      end if
      if(nh>1.0_dp) then
        scale=f*nh/(nh-1.0_dp)
      else
        scale=f
      end if
      ok=.true.
      if(singleton .and. (lonely==LONELY_REMOVE .or. lonely==LONELY_CERTAINTY)) then
        if(lonely==LONELY_CERTAINTY .and. poph>1.0_dp) ok=.false.
        if(lonely==LONELY_REMOVE) ok=.false.
      end if
      if(ok) then
        do j=1,size(psutot,1)
          v=v+scale*outer_product(psutot(j,:),psutot(j,:))
        end do
        nok=nok+1
      end if
      deallocate(idxs,cu,psutot)
    end do
    if(lonely==LONELY_AVERAGE .and. nok>0 .and. nok<nstrata) v=v*real(nstrata,dp)/real(nok,dp)
  end function onestage_variance

  integer function find_first(x,value) result(pos)
    integer, intent(in) :: x(:),value
    integer :: i
    pos=1
    do i=1,size(x)
      if(x(i)==value) then; pos=i; return; end if
    end do
  end function find_first

  subroutine unique_int(x,u)
    integer, intent(in) :: x(:)
    integer, allocatable, intent(out) :: u(:)
    integer, allocatable :: t(:)
    integer :: i,m
    allocate(t(size(x))); m=0
    do i=1,size(x)
      if(m==0 .or. .not.any(t(1:m)==x(i))) then; m=m+1; t(m)=x(i); end if
    end do
    allocate(u(m)); if(m>0) u=t(1:m)
  end subroutine unique_int
end module survey_taylor
