! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_ivreg
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, rep_design_t, ivreg_result_t
  use survey_linalg, only : sym_pinv
  use survey_taylor, only : svyrecvar
  use survey_replicates, only : svr_var
  use survey_design, only : design_degf
  implicit none
  private
  public :: svy_ivreg, rep_ivreg
contains
  subroutine svy_ivreg(x,z,y,design,result,rescale_weights)
    real(dp), intent(in) :: x(:,:),z(:,:),y(:)
    type(survey_design_t), intent(in) :: design
    type(ivreg_result_t), intent(out) :: result
    logical, intent(in), optional :: rescale_weights
    real(dp), allocatable :: w(:), beta(:),bread(:,:),xhat(:,:),score(:,:),infl(:,:),fit(:),resid(:)
    real(dp) :: sigma2,sw
    integer :: n,p,rank
    logical :: resc
    n=size(x,1);p=size(x,2)
    if(size(z,1)/=n .or. size(y)/=n .or. design%n/=n) error stop 'svy_ivreg: shape mismatch'
    if(size(z,2)<p) error stop 'svy_ivreg: need at least as many instruments as regressors'
    allocate(w(n));w=design%weight;resc=.true.;if(present(rescale_weights))resc=rescale_weights
    if(resc) w=w/(sum(w)/real(n,dp))
    call weighted_2sls(x,z,y,w,beta,bread,xhat,rank)
    allocate(fit(n),resid(n),score(n,p),infl(n,p));fit=matmul(x,beta);resid=y-fit
    score=0.0_dp
    block
      integer :: i,j
      do i=1,n;do j=1,p;score(i,j)=w(i)*xhat(i,j)*resid(i);end do;end do
    end block
    infl=matmul(score,bread)
    allocate(result%coef(p),result%vcov(p,p),result%naive_vcov(p,p),result%fitted(n),result%residual(n))
    result%coef=beta;result%vcov=svyrecvar(infl,design);result%fitted=fit;result%residual=resid;result%rank=rank
    sw=sum(w);sigma2=sum(w*resid*resid)/max(1.0_dp,sw-real(rank,dp));result%naive_vcov=bread*sigma2
    result%df_residual=max(0,design_degf(design)+1-rank)
  end subroutine svy_ivreg

  subroutine rep_ivreg(x,z,y,design,result)
    real(dp), intent(in) :: x(:,:),z(:,:),y(:)
    type(rep_design_t), intent(in) :: design
    type(ivreg_result_t), intent(out) :: result
    real(dp), allocatable :: beta(:),betatmp(:),bread(:,:),xhat(:,:),betarep(:,:),vv(:,:),w(:),fit(:),resid(:)
    real(dp) :: sigma2,sw
    integer :: p,r,rank
    p=size(x,2);allocate(w(design%n));w=design%weight
    call weighted_2sls(x,z,y,w,beta,bread,xhat,rank)
    allocate(betarep(design%r,p),vv(p,p))
    do r=1,design%r
      w=design%repweights(:,r); call weighted_2sls(x,z,y,w,betatmp,bread,xhat,rank); betarep(r,:)=betatmp
    end do
    vv=svr_var(betarep,design%scale,design%rscales,design%mse,beta)
    allocate(fit(design%n),resid(design%n));fit=matmul(x,beta);resid=y-fit;w=design%weight;sw=sum(w);sigma2=sum(w*resid*resid)/max(1.0_dp,sw-real(rank,dp))
    allocate(result%coef(p),result%vcov(p,p),result%naive_vcov(p,p),result%fitted(design%n),result%residual(design%n))
    result%coef=beta;result%vcov=vv;result%naive_vcov=bread*sigma2;result%fitted=fit;result%residual=resid;result%rank=rank;result%df_residual=max(0,design%r-rank)
  end subroutine rep_ivreg

  subroutine weighted_2sls(x,z,y,w,beta,bread,xhat,rank)
    real(dp),intent(in)::x(:,:),z(:,:),y(:),w(:)
    real(dp),allocatable,intent(out)::beta(:),bread(:,:),xhat(:,:)
    integer,intent(out)::rank
    real(dp),allocatable::ztwz(:,:),zinv(:,:),ztwx(:,:),gamma(:,:),a(:,:),ainv(:,:),rhs(:)
    integer::n,p,q,i,j,k,rz,info
    n=size(x,1);p=size(x,2);q=size(z,2)
    allocate(ztwz(q,q),zinv(q,q),ztwx(q,p),gamma(q,p),xhat(n,p),a(p,p),ainv(p,p),rhs(p),beta(p),bread(p,p))
    ztwz=0.0_dp;ztwx=0.0_dp
    do i=1,n
      do j=1,q
        do k=1,q;ztwz(j,k)=ztwz(j,k)+w(i)*z(i,j)*z(i,k);end do
        do k=1,p;ztwx(j,k)=ztwx(j,k)+w(i)*z(i,j)*x(i,k);end do
      end do
    end do
    call sym_pinv(ztwz,zinv,rz,info=info);gamma=matmul(zinv,ztwx);xhat=matmul(z,gamma)
    a=0.0_dp;rhs=0.0_dp
    do i=1,n
      do j=1,p
        rhs(j)=rhs(j)+w(i)*xhat(i,j)*y(i)
        do k=1,p;a(j,k)=a(j,k)+w(i)*xhat(i,j)*x(i,k);end do
      end do
    end do
    ! Xhat' W X is symmetric up to roundoff for a weighted projection; symmetrize before pseudoinversion.
    a=0.5_dp*(a+transpose(a));call sym_pinv(a,ainv,rank,info=info);beta=matmul(ainv,rhs);bread=ainv
  end subroutine weighted_2sls
end module survey_ivreg
