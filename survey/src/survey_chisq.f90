! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_chisq
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, svystat_t, chisq_result_t
  use survey_estimators, only : svy_total, svy_mean
  use survey_design, only : design_degf
  use survey_linalg, only : sym_pinv
  use survey_special, only : f_survival, chisq_survival
  implicit none
  private
  public :: svy_chisq_wald, svy_chisq_rao_scott
contains

  subroutine svy_chisq_wald(rowcat,nr,colcat,nc,design,result,adjusted)
    integer, intent(in) :: rowcat(:),colcat(:),nr,nc
    type(survey_design_t), intent(in) :: design
    type(chisq_result_t), intent(out) :: result
    logical, intent(in), optional :: adjusted
    type(svystat_t) :: alltot
    real(dp), allocatable :: z(:,:), jac(:,:), ycell(:), vc(:,:), vu(:,:), yu(:), vinv(:,:)
    real(dp) :: grand, raw
    integer :: n,ncell,q,i,j,k,u,df,nu,rank,info
    logical :: adj
    n=design%n; ncell=nr*nc; q=ncell+nr+nc+1
    if(size(rowcat)/=n .or. size(colcat)/=n) error stop 'svy_chisq_wald: size mismatch'
    if(nr<2 .or. nc<2) error stop 'svy_chisq_wald: need at least 2x2 table'
    if(any(rowcat<1).or.any(rowcat>nr).or.any(colcat<1).or.any(colcat>nc)) error stop 'svy_chisq_wald: category out of range'
    allocate(z(n,q)); z=0.0_dp
    do k=1,n
      i=rowcat(k); j=colcat(k)
      z(k,i+(j-1)*nr)=1.0_dp
      z(k,ncell+i)=1.0_dp
      z(k,ncell+nr+j)=1.0_dp
      z(k,q)=1.0_dp
    end do
    alltot=svy_total(z,design)
    grand=alltot%estimate(q)
    if(grand<=0.0_dp) error stop 'svy_chisq_wald: empty table'
    allocate(jac(ncell,q),ycell(ncell),vc(ncell,ncell)); jac=0.0_dp
    do j=1,nc
      do i=1,nr
        k=i+(j-1)*nr
        ycell(k)=alltot%estimate(k)-alltot%estimate(ncell+i)*alltot%estimate(ncell+nr+j)/grand
        jac(k,k)=1.0_dp
        jac(k,ncell+i)=-alltot%estimate(ncell+nr+j)/grand
        jac(k,ncell+nr+j)=-alltot%estimate(ncell+i)/grand
        jac(k,q)=alltot%estimate(ncell+i)*alltot%estimate(ncell+nr+j)/(grand*grand)
      end do
    end do
    vc=matmul(jac,matmul(alltot%variance,transpose(jac)))
    df=(nr-1)*(nc-1); allocate(vu(df,df),yu(df),vinv(df,df)); u=0
    do j=2,nc
      do i=2,nr
        u=u+1; k=i+(j-1)*nr; yu(u)=ycell(k)
      end do
    end do
    u=0
    do j=2,nc
      do i=2,nr
        u=u+1; k=i+(j-1)*nr
        call fill_use_cov_row(vc,nr,nc,k,vu(u,:))
      end do
    end do
    call sym_pinv(vu,vinv,rank,info=info)
    raw=dot_product(yu,matmul(vinv,yu)); nu=max(1,design_degf(design)); adj=.false.; if(present(adjusted)) adj=adjusted
    result%rank=rank; result%numerator_df=real(df,dp)
    if(adj) then
      result%denominator_df=real(max(1,nu-df+1),dp)
      result%statistic=raw*result%denominator_df/(real(df,dp)*real(nu,dp))
    else
      result%denominator_df=real(nu,dp)
      result%statistic=raw/real(df,dp)
    end if
    result%p_value=f_survival(result%statistic,result%numerator_df,result%denominator_df)
  end subroutine svy_chisq_wald

  subroutine svy_chisq_rao_scott(rowcat,nr,colcat,nc,design,result,chisq_version)
    integer, intent(in) :: rowcat(:),colcat(:),nr,nc
    type(survey_design_t), intent(in) :: design
    type(chisq_result_t), intent(out) :: result
    logical, intent(in), optional :: chisq_version
    type(svystat_t) :: cellmean
    real(dp), allocatable :: mm(:,:), p(:), rowp(:),colp(:), x1(:,:), xi(:,:), xtx(:,:), xtxi(:,:), cmat(:,:), &
      id(:,:), denom(:,:), deninv(:,:), numr(:,:), delta(:,:)
    real(dp) :: pearson,tr,tr2,d0,nobs,expected
    integer :: n,ncell,df,p1,i,j,k,rank,info,nu
    logical :: ch
    n=design%n; ncell=nr*nc; df=(nr-1)*(nc-1); p1=nr+nc-1
    if(size(rowcat)/=n .or. size(colcat)/=n) error stop 'svy_chisq_rao_scott: size mismatch'
    allocate(mm(n,ncell)); mm=0.0_dp
    do k=1,n; mm(k,rowcat(k)+(colcat(k)-1)*nr)=1.0_dp; end do
    cellmean=svy_mean(mm,design); allocate(p(ncell),rowp(nr),colp(nc)); p=cellmean%estimate
    rowp=0.0_dp; colp=0.0_dp
    do j=1,nc; do i=1,nr; k=i+(j-1)*nr; rowp(i)=rowp(i)+p(k); colp(j)=colp(j)+p(k); end do; end do
    nobs=real(n,dp); pearson=0.0_dp
    do j=1,nc; do i=1,nr; k=i+(j-1)*nr; expected=rowp(i)*colp(j); if(expected>0) pearson=pearson+nobs*(p(k)-expected)**2/expected; end do; end do

    allocate(x1(ncell,p1),xi(ncell,df)); x1=0.0_dp; xi=0.0_dp
    do j=1,nc
      do i=1,nr
        k=i+(j-1)*nr; x1(k,1)=1.0_dp
        if(i>=2) x1(k,i)=1.0_dp
        if(j>=2) x1(k,nr+j-1)=1.0_dp
        if(i>=2 .and. j>=2) xi(k,(j-2)*(nr-1)+(i-1))=1.0_dp
      end do
    end do
    allocate(xtx(p1,p1),xtxi(p1,p1),cmat(ncell,df)); xtx=matmul(transpose(x1),x1)
    call sym_pinv(xtx,xtxi,rank,info=info); cmat=xi-matmul(x1,matmul(xtxi,matmul(transpose(x1),xi)))
    allocate(id(ncell,ncell)); id=0.0_dp
    do k=1,ncell; if(p(k)>0) id(k,k)=1.0_dp/p(k); end do
    allocate(denom(df,df),deninv(df,df),numr(df,df),delta(df,df))
    denom=matmul(transpose(cmat),matmul(id/nobs,cmat))
    numr=matmul(transpose(cmat),matmul(id,matmul(cellmean%variance,matmul(id,cmat))))
    call sym_pinv(denom,deninv,rank,info=info); delta=matmul(deninv,numr)
    tr=0.0_dp; tr2=0.0_dp
    do i=1,df; tr=tr+delta(i,i); end do
    do i=1,df; do j=1,df; tr2=tr2+delta(i,j)*delta(j,i); end do; end do
    if(tr<=tiny(1.0_dp) .or. tr2<=tiny(1.0_dp)) error stop 'svy_chisq_rao_scott: singular adjustment'
    d0=tr*tr/tr2; nu=max(1,design_degf(design)); ch=.false.; if(present(chisq_version)) ch=chisq_version
    result%rank=rank
    if(ch) then
      result%statistic=pearson/(tr/real(df,dp)); result%numerator_df=real(df,dp); result%denominator_df=0.0_dp
      result%p_value=chisq_survival(result%statistic,result%numerator_df)
    else
      result%statistic=pearson/tr; result%numerator_df=d0; result%denominator_df=d0*real(nu,dp)
      result%p_value=f_survival(result%statistic,result%numerator_df,result%denominator_df)
    end if
  end subroutine svy_chisq_rao_scott

  subroutine fill_use_cov_row(vc,nr,nc,source,row)
    real(dp), intent(in) :: vc(:,:)
    integer, intent(in) :: nr,nc,source
    real(dp), intent(out) :: row(:)
    integer :: i,j,u,k
    u=0
    do j=2,nc
      do i=2,nr
        u=u+1; k=i+(j-1)*nr; row(u)=vc(source,k)
      end do
    end do
  end subroutine fill_use_cov_row
end module survey_chisq
