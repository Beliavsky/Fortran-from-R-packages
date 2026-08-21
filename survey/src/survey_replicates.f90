! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_replicates
  use survey_kinds, only : dp
  use survey_types, only : rep_design_t, svystat_t, ratio_result_t
  use survey_linalg, only : outer_product
  implicit none
  private
  public :: svr_var, make_jk1, make_jkn, make_brr, make_bootstrap, rep_mean, rep_total, rep_ratio
contains

  function svr_var(theta,scale,rscales,mse,coef) result(v)
    real(dp), intent(in) :: theta(:,:),scale,rscales(:)
    logical, intent(in) :: mse
    real(dp), intent(in) :: coef(:)
    real(dp) :: v(size(theta,2),size(theta,2))
    real(dp) :: center(size(theta,2)), d(size(theta,2)), sw
    integer :: r
    if(size(theta,1)/=size(rscales) .or. size(theta,2)/=size(coef)) error stop 'svr_var: shape mismatch'
    if(mse) then
      center=coef
    else
      sw=real(count(rscales>0),dp)
      if(sw<=0) error stop 'svr_var: no positive rscales'
      center=0
      do r=1,size(theta,1); if(rscales(r)>0) center=center+theta(r,:)/sw; end do
    end if
    v=0
    do r=1,size(theta,1)
      d=theta(r,:)-center
      v=v+scale*rscales(r)*outer_product(d,d)
    end do
  end function svr_var

  subroutine make_jk1(psu,base_weight,rep,fpc_correction,mse)
    integer, intent(in) :: psu(:)
    real(dp), intent(in) :: base_weight(:)
    type(rep_design_t), intent(out) :: rep
    real(dp), intent(in), optional :: fpc_correction
    logical, intent(in), optional :: mse
    integer, allocatable :: u(:)
    real(dp) :: fpc
    integer :: i,j,npsu
    if(size(psu)/=size(base_weight)) error stop 'make_jk1: size mismatch'
    call unique_int(psu,u); npsu=size(u)
    if(npsu<2) error stop 'make_jk1: need at least two PSUs'
    rep%n=size(psu); rep%r=npsu
    allocate(rep%weight(rep%n),rep%repweights(rep%n,rep%r),rep%rscales(rep%r))
    rep%weight=base_weight; rep%rscales=1
    fpc=1; if(present(fpc_correction)) fpc=fpc_correction
    rep%scale=fpc*real(npsu-1,dp)/real(npsu,dp); rep%mse=.false.; if(present(mse)) rep%mse=mse
    do j=1,npsu
      do i=1,rep%n
        if(psu(i)==u(j)) then
          rep%repweights(i,j)=0
        else
          rep%repweights(i,j)=base_weight(i)*real(npsu,dp)/real(npsu-1,dp)
        end if
      end do
    end do
  end subroutine make_jk1

  subroutine make_jkn(strata,psu,base_weight,rep,mse)
    integer, intent(in) :: strata(:),psu(:)
    real(dp), intent(in) :: base_weight(:)
    type(rep_design_t), intent(out) :: rep
    logical, intent(in), optional :: mse
    integer, allocatable :: repsu(:), repstr(:)
    integer :: i,k,nreps,nh
    if(size(strata)/=size(psu) .or. size(psu)/=size(base_weight)) error stop 'make_jkn: size mismatch'
    call unique_pairs(strata,psu,repsu,repstr); nreps=size(repsu)
    rep%n=size(psu); rep%r=nreps; allocate(rep%weight(rep%n),rep%repweights(rep%n,nreps),rep%rscales(nreps))
    rep%weight=base_weight; rep%repweights=spread(base_weight,2,nreps); rep%scale=1; rep%mse=.false.; if(present(mse)) rep%mse=mse
    do k=1,nreps
      nh=count_unique(pack(psu,strata==repstr(k)))
      if(nh<2) then
        rep%rscales(k)=0; cycle
      end if
      rep%rscales(k)=real(nh-1,dp)/real(nh,dp)
      do i=1,rep%n
        if(strata(i)==repstr(k)) then
          if(psu(i)==repsu(k)) then
            rep%repweights(i,k)=0
          else
            rep%repweights(i,k)=base_weight(i)*real(nh,dp)/real(nh-1,dp)
          end if
        end if
      end do
    end do
  end subroutine make_jkn

  subroutine make_brr(strata,psu,base_weight,rep,fay_rho,mse)
    integer, intent(in) :: strata(:),psu(:)
    real(dp), intent(in) :: base_weight(:)
    type(rep_design_t), intent(out) :: rep
    real(dp), intent(in), optional :: fay_rho
    logical, intent(in), optional :: mse
    integer, allocatable :: su(:), cu(:), psu1(:),psu2(:)
    integer, allocatable :: hmat(:,:)
    integer :: nstrata,h,i,rr,nreps
    real(dp) :: rho,sel
    call unique_int(strata,su); nstrata=size(su); allocate(psu1(nstrata),psu2(nstrata))
    do h=1,nstrata
      call unique_int(pack(psu,strata==su(h)),cu)
      if(size(cu)/=2) error stop 'make_brr: each stratum must have exactly two PSUs'
      psu1(h)=cu(1); psu2(h)=cu(2)
    end do
    call sylvester_hadamard(nstrata+1,hmat); nreps=size(hmat,1)
    rep%n=size(psu); rep%r=nreps; allocate(rep%weight(rep%n),rep%repweights(rep%n,nreps),rep%rscales(nreps))
    rep%weight=base_weight; rep%rscales=1
    rep%mse=.false.; if(present(mse)) rep%mse=mse
    rho=0; if(present(fay_rho)) rho=fay_rho; rep%scale=1.0_dp/(real(nreps,dp)*(1.0_dp-rho)**2)
    do rr=1,nreps
      do i=1,rep%n
        h=find_index(su,strata(i))
        if(hmat(rr,h+1)>0) then
          if(psu(i)==psu1(h)) then; sel=2-rho; else; sel=rho; end if
        else
          if(psu(i)==psu2(h)) then; sel=2-rho; else; sel=rho; end if
        end if
        rep%repweights(i,rr)=base_weight(i)*sel
      end do
    end do
  end subroutine make_brr

  subroutine make_bootstrap(strata,psu,base_weight,rep,nrep,mse)
    integer, intent(in) :: strata(:),psu(:),nrep
    real(dp), intent(in) :: base_weight(:)
    type(rep_design_t), intent(out) :: rep
    logical, intent(in), optional :: mse
    integer, allocatable :: su(:),cu(:), draw(:)
    integer :: nstrata,h,i,j,rr,nh,pick
    real(dp) :: u,hm
    call unique_int(strata,su); nstrata=size(su)
    rep%n=size(psu); rep%r=nrep; allocate(rep%weight(rep%n),rep%repweights(rep%n,nrep),rep%rscales(nrep))
    rep%weight=base_weight; rep%repweights=0; rep%rscales=1; rep%mse=.false.; if(present(mse)) rep%mse=mse
    hm=0
    do h=1,nstrata; nh=count_unique(pack(psu,strata==su(h))); hm=hm+1.0_dp/real(nh,dp); end do
    hm=real(nstrata,dp)/hm; rep%scale=hm/((hm-1.0_dp)*real(nrep-1,dp))
    do rr=1,nrep
      do h=1,nstrata
        call unique_int(pack(psu,strata==su(h)),cu); nh=size(cu); allocate(draw(nh)); draw=0
        do j=1,nh
          call random_number(u); pick=min(nh,1+int(u*nh)); draw(pick)=draw(pick)+1
        end do
        do i=1,rep%n
          if(strata(i)==su(h)) then
            j=find_index(cu,psu(i)); rep%repweights(i,rr)=base_weight(i)*real(draw(j),dp)
          end if
        end do
        deallocate(draw,cu)
      end do
    end do
  end subroutine make_bootstrap

  function rep_total(x,design) result(ans)
    real(dp), intent(in) :: x(:,:)
    type(rep_design_t), intent(in) :: design
    type(svystat_t) :: ans
    real(dp), allocatable :: theta(:,:)
    integer :: r,j
    if(size(x,1)/=design%n) error stop 'rep_total: row mismatch'
    allocate(ans%estimate(size(x,2)),theta(design%r,size(x,2)),ans%variance(size(x,2),size(x,2)))
    do j=1,size(x,2); ans%estimate(j)=dot_product(design%weight,x(:,j)); end do
    do r=1,design%r; do j=1,size(x,2); theta(r,j)=dot_product(design%repweights(:,r),x(:,j)); end do; end do
    ans%variance=svr_var(theta,design%scale,design%rscales,design%mse,ans%estimate)
  end function rep_total

  function rep_mean(x,design) result(ans)
    real(dp), intent(in) :: x(:,:)
    type(rep_design_t), intent(in) :: design
    type(svystat_t) :: ans
    real(dp), allocatable :: theta(:,:)
    real(dp) :: sw
    integer :: r,j
    if(size(x,1)/=design%n) error stop 'rep_mean: row mismatch'
    allocate(ans%estimate(size(x,2)),theta(design%r,size(x,2)),ans%variance(size(x,2),size(x,2)))
    sw=sum(design%weight); do j=1,size(x,2); ans%estimate(j)=dot_product(design%weight,x(:,j))/sw; end do
    do r=1,design%r
      sw=sum(design%repweights(:,r))
      do j=1,size(x,2); theta(r,j)=dot_product(design%repweights(:,r),x(:,j))/sw; end do
    end do
    ans%variance=svr_var(theta,design%scale,design%rscales,design%mse,ans%estimate)
  end function rep_mean

  function rep_ratio(numer,denom,design) result(ans)
    real(dp), intent(in) :: numer(:,:),denom(:,:)
    type(rep_design_t), intent(in) :: design
    type(ratio_result_t) :: ans
    real(dp), allocatable :: theta(:,:),coef(:),vv(:,:)
    integer :: i,j,k,r,nn,nd
    nn=size(numer,2); nd=size(denom,2); allocate(ans%ratio(nn,nd),ans%variance(nn,nd),theta(design%r,nn*nd),coef(nn*nd),vv(nn*nd,nn*nd))
    k=0
    do j=1,nd; do i=1,nn; k=k+1; ans%ratio(i,j)=dot_product(design%weight,numer(:,i))/dot_product(design%weight,denom(:,j)); coef(k)=ans%ratio(i,j); end do; end do
    do r=1,design%r
      k=0; do j=1,nd; do i=1,nn; k=k+1; theta(r,k)=dot_product(design%repweights(:,r),numer(:,i))/dot_product(design%repweights(:,r),denom(:,j)); end do; end do
    end do
    vv=svr_var(theta,design%scale,design%rscales,design%mse,coef); k=0
    do j=1,nd; do i=1,nn; k=k+1; ans%variance(i,j)=vv(k,k); end do; end do
  end function rep_ratio

  subroutine sylvester_hadamard(mincols,H)
    integer, intent(in) :: mincols
    integer, allocatable, intent(out) :: H(:,:)
    integer, allocatable :: A(:,:),B(:,:)
    integer :: n
    n=1; allocate(A(1,1)); A=1
    do while(n<mincols)
      allocate(B(2*n,2*n)); B(1:n,1:n)=A; B(1:n,n+1:)=A; B(n+1:,1:n)=A; B(n+1:,n+1:)=-A
      call move_alloc(B,A); n=2*n
    end do
    call move_alloc(A,H)
  end subroutine sylvester_hadamard

  integer function count_unique(x) result(n)
    integer,intent(in)::x(:); integer,allocatable::u(:); call unique_int(x,u); n=size(u)
  end function count_unique
  integer function find_index(x,v) result(k)
    integer,intent(in)::x(:),v; integer::i; k=0; do i=1,size(x); if(x(i)==v) then;k=i;return;end if;end do
  end function find_index
  subroutine unique_int(x,u)
    integer,intent(in)::x(:); integer,allocatable,intent(out)::u(:); integer,allocatable::t(:); integer::i,n
    allocate(t(size(x)));n=0;do i=1,size(x);if(n==0.or..not.any(t(1:n)==x(i)))then;n=n+1;t(n)=x(i);end if;end do;allocate(u(n));if(n>0)u=t(1:n)
  end subroutine unique_int
  subroutine unique_pairs(strata,psu,u_psu,u_strata)
    integer,intent(in)::strata(:),psu(:);integer,allocatable,intent(out)::u_psu(:),u_strata(:)
    integer,allocatable::tp(:),ts(:);integer::i,n
    allocate(tp(size(psu)),ts(size(psu)));n=0
    do i=1,size(psu)
      if(n==0.or..not.any((tp(1:n)==psu(i)).and.(ts(1:n)==strata(i))))then;n=n+1;tp(n)=psu(i);ts(n)=strata(i);end if
    end do
    allocate(u_psu(n),u_strata(n));u_psu=tp(1:n);u_strata=ts(1:n)
  end subroutine unique_pairs
end module survey_replicates
