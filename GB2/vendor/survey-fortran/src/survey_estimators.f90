! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_estimators
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, svystat_t, ratio_result_t
  use survey_taylor, only : svyrecvar
  implicit none
  private
  public :: svy_total, svy_mean, svy_ratio, svy_covariance, svy_cdf, svy_table1, svy_table2, svy_kappa, svy_cprod
contains

  function svy_total(x,design,influence) result(r)
    real(dp), intent(in) :: x(:,:)
    type(survey_design_t), intent(in) :: design
    logical, intent(in), optional :: influence
    type(svystat_t) :: r
    real(dp), allocatable :: z(:,:)
    logical :: infl
    integer :: j
    if(size(x,1)/=design%n) error stop 'svy_total: row count mismatch'
    allocate(z(design%n,size(x,2)))
    do j=1,size(x,2); z(:,j)=x(:,j)*design%weight; end do
    allocate(r%estimate(size(x,2)),r%variance(size(x,2),size(x,2)))
    r%estimate=sum(z,dim=1); r%variance=svyrecvar(z,design)
    infl=.false.; if(present(influence)) infl=influence
    if(infl) then; allocate(r%influence,source=z); end if
  end function svy_total

  function svy_mean(x,design,influence) result(r)
    real(dp), intent(in) :: x(:,:)
    type(survey_design_t), intent(in) :: design
    logical, intent(in), optional :: influence
    type(svystat_t) :: r
    real(dp), allocatable :: z(:,:)
    real(dp) :: sw
    logical :: infl
    integer :: j
    if(size(x,1)/=design%n) error stop 'svy_mean: row count mismatch'
    sw=sum(design%weight); if(sw<=0) error stop 'svy_mean: nonpositive total weight'
    allocate(r%estimate(size(x,2)),z(design%n,size(x,2)))
    do j=1,size(x,2); r%estimate(j)=dot_product(design%weight,x(:,j))/sw; end do
    do j=1,size(x,2); z(:,j)=(x(:,j)-r%estimate(j))*design%weight/sw; end do
    allocate(r%variance(size(x,2),size(x,2))); r%variance=svyrecvar(z,design)
    infl=.false.; if(present(influence)) infl=influence
    if(infl) then; allocate(r%influence,source=z); end if
  end function svy_mean

  function svy_ratio(numer,denom,design,influence) result(r)
    real(dp), intent(in) :: numer(:,:),denom(:,:)
    type(survey_design_t), intent(in) :: design
    logical, intent(in), optional :: influence
    type(ratio_result_t) :: r
    real(dp), allocatable :: nt(:),dt(:),z(:,:), cov(:,:)
    logical :: infl
    integer :: i,j,k,nn,nd
    if(size(numer,1)/=design%n .or. size(denom,1)/=design%n) error stop 'svy_ratio: row count mismatch'
    nn=size(numer,2); nd=size(denom,2)
    allocate(nt(nn),dt(nd),r%ratio(nn,nd),r%variance(nn,nd),z(design%n,nn*nd),cov(nn*nd,nn*nd))
    do i=1,nn; nt(i)=dot_product(design%weight,numer(:,i)); end do
    do j=1,nd; dt(j)=dot_product(design%weight,denom(:,j)); end do
    k=0
    do j=1,nd
      if(abs(dt(j))<=tiny(1.0_dp)) error stop 'svy_ratio: zero denominator total'
      do i=1,nn
        k=k+1; r%ratio(i,j)=nt(i)/dt(j)
        z(:,k)=(numer(:,i)-r%ratio(i,j)*denom(:,j))*design%weight/dt(j)
      end do
    end do
    cov=svyrecvar(z,design); k=0
    do j=1,nd; do i=1,nn; k=k+1; r%variance(i,j)=cov(k,k); end do; end do
    infl=.false.; if(present(influence)) infl=influence
    if(infl) allocate(r%influence,source=z)
  end function svy_ratio

  subroutine svy_covariance(x,design,covest,varcov)
    real(dp), intent(in) :: x(:,:)
    type(survey_design_t), intent(in) :: design
    real(dp), intent(out) :: covest(:,:)
    real(dp), intent(out), optional :: varcov(:,:)
    type(svystat_t) :: m, q
    real(dp), allocatable :: zz(:,:)
    integer :: p,n,i,j,k
    real(dp) :: kish
    n=size(x,1); p=size(x,2)
    if(size(covest,1)/=p .or. size(covest,2)/=p) error stop 'svy_covariance: shape mismatch'
    m=svy_mean(x,design)
    if(count(design%weight>0)>1) then
      kish=real(count(design%weight>0),dp)/real(count(design%weight>0)-1,dp)
    else
      covest=0; if(present(varcov)) varcov=0; return
    end if
    allocate(zz(n,p*p)); k=0
    do j=1,p
      do i=1,p
        k=k+1; zz(:,k)=(x(:,i)-m%estimate(i))*(x(:,j)-m%estimate(j))*kish
      end do
    end do
    q=svy_mean(zz,design); k=0
    do j=1,p; do i=1,p; k=k+1; covest(i,j)=q%estimate(k); end do; end do
    if(present(varcov)) then
      if(any(shape(varcov)/=[p*p,p*p])) error stop 'svy_covariance: varcov shape mismatch'
      varcov=q%variance
    end if
  end subroutine svy_covariance

  function svy_cprod(x,y,design) result(r)
    real(dp), intent(in) :: x(:,:),y(:,:)
    type(survey_design_t), intent(in) :: design
    type(svystat_t) :: r
    real(dp), allocatable :: z(:,:)
    integer :: i,j,k
    if(size(x,1)/=design%n .or. size(y,1)/=design%n) error stop 'svy_cprod: row count mismatch'
    allocate(z(design%n,size(x,2)*size(y,2))); k=0
    do j=1,size(y,2); do i=1,size(x,2); k=k+1; z(:,k)=x(:,i)*y(:,j); end do; end do
    r=svy_total(z,design)
  end function svy_cprod

  function svy_cdf(x,points,design) result(r)
    real(dp), intent(in) :: x(:),points(:)
    type(survey_design_t), intent(in) :: design
    type(svystat_t) :: r
    real(dp), allocatable :: ind(:,:)
    integer :: j
    if(size(x)/=design%n) error stop 'svy_cdf: row count mismatch'
    allocate(ind(design%n,size(points)))
    do j=1,size(points); ind(:,j)=merge(1.0_dp,0.0_dp,x<=points(j)); end do
    r=svy_mean(ind,design)
  end function svy_cdf

  function svy_table1(category,nlevels,design) result(r)
    integer, intent(in) :: category(:),nlevels
    type(survey_design_t), intent(in) :: design
    type(svystat_t) :: r
    real(dp), allocatable :: ind(:,:)
    integer :: j
    if(size(category)/=design%n) error stop 'svy_table1: row count mismatch'
    allocate(ind(design%n,nlevels)); ind=0
    do j=1,nlevels; where(category==j) ind(:,j)=1.0_dp; end do
    r=svy_total(ind,design)
  end function svy_table1

  function svy_table2(rowcat,nrowlev,colcat,ncollev,design) result(r)
    integer, intent(in) :: rowcat(:),colcat(:),nrowlev,ncollev
    type(survey_design_t), intent(in) :: design
    type(svystat_t) :: r
    real(dp), allocatable :: ind(:,:)
    integer :: i,j,k
    if(size(rowcat)/=design%n .or. size(colcat)/=design%n) error stop 'svy_table2: row count mismatch'
    allocate(ind(design%n,nrowlev*ncollev)); ind=0; k=0
    do j=1,ncollev; do i=1,nrowlev; k=k+1; where(rowcat==i .and. colcat==j) ind(:,k)=1; end do; end do
    r=svy_total(ind,design)
  end function svy_table2

  subroutine svy_kappa(rater1,rater2,nlevels,design,kappa,se)
    integer, intent(in) :: rater1(:),rater2(:),nlevels
    type(survey_design_t), intent(in) :: design
    real(dp), intent(out) :: kappa,se
    type(svystat_t) :: tab
    real(dp), allocatable :: p(:,:), rowp(:),colp(:), grad(:), gmat(:,:), covp(:,:)
    real(dp) :: total,po,pe
    integer :: i,j,k,L
    tab=svy_table2(rater1,nlevels,rater2,nlevels,design)
    L=nlevels; total=sum(tab%estimate); if(total<=0) error stop 'svy_kappa: empty table'
    allocate(p(L,L),rowp(L),colp(L),grad(L*L),gmat(1,L*L),covp(L*L,L*L))
    p=reshape(tab%estimate,[L,L])/total
    rowp=sum(p,dim=2); colp=sum(p,dim=1)
    po=sum([(p(i,i),i=1,L)]); pe=dot_product(rowp,colp)
    kappa=(po-pe)/(1.0_dp-pe)
    ! numerical delta gradient wrt cell proportions on simplex, then use total-table influence normalization
    grad=0; k=0
    do j=1,L
      do i=1,L
        k=k+1
        grad(k)=((merge(1.0_dp,0.0_dp,i==j)-(colp(i)+rowp(j)))*(1-pe) + (po-pe)*(colp(i)+rowp(j))) /(1-pe)**2
      end do
    end do
    ! covariance of estimated proportions via delta from totals
    covp=tab%variance/total**2
    do i=1,L*L
      do j=1,L*L
        covp(i,j)=covp(i,j) - tab%estimate(i)*sum(tab%variance(:,j))/total**3 &
                            - tab%estimate(j)*sum(tab%variance(i,:))/total**3 &
                            + tab%estimate(i)*tab%estimate(j)*sum(tab%variance)/total**4
      end do
    end do
    se=sqrt(max(0.0_dp,dot_product(grad,matmul(covp,grad))))
  end subroutine svy_kappa

end module survey_estimators
