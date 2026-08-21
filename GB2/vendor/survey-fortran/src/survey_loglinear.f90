! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_loglinear
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, loglin_result_t, loglin_test_t, glm_result_t, &
    FAMILY_POISSON, LINK_LOG
  use survey_design, only : make_design, design_degf
  use survey_estimators, only : svy_mean
  use survey_glm, only : svy_glm
  use survey_linalg, only : sym_pinv, symmetric_eigen, outer_product
  use survey_special, only : weighted_chisq_survival
  implicit none
  private
  public :: svy_loglin_cells, svy_loglin_prob, svy_loglin_compare
contains

  subroutine svy_loglin_cells(cell, ncell, x_cells, design, result)
    integer, intent(in) :: cell(:), ncell
    real(dp), intent(in) :: x_cells(:,:)
    type(survey_design_t), intent(in) :: design
    type(loglin_result_t), intent(out) :: result
    real(dp), allocatable :: ind(:,:)
    integer :: i
    ! The local block avoids exposing svystat_t in this public interface.
    block
      use survey_types, only : svystat_t
      type(svystat_t) :: phat
      if (size(cell) /= design%n .or. size(x_cells,1) /= ncell) error stop 'svy_loglin_cells: shape mismatch'
      if (minval(cell) < 1 .or. maxval(cell) > ncell) error stop 'svy_loglin_cells: cell index out of range'
      allocate(ind(design%n,ncell)); ind = 0.0_dp
      do i = 1, design%n
        ind(i,cell(i)) = 1.0_dp
      end do
      phat = svy_mean(ind,design)
      call svy_loglin_prob(phat%estimate,phat%variance,x_cells,real(design%n,dp),design_degf(design),result)
    end block
  end subroutine svy_loglin_cells

  subroutine svy_loglin_prob(prob, prob_vcov, x_cells, sample_n, df_null, result)
    real(dp), intent(in) :: prob(:), prob_vcov(:,:), x_cells(:,:), sample_n
    integer, intent(in) :: df_null
    type(loglin_result_t), intent(out) :: result
    type(survey_design_t) :: dcell
    type(glm_result_t) :: fit
    integer, allocatable :: cluster(:,:)
    real(dp), allocatable :: weight(:), counts(:), pfit(:), pmat(:,:), xc(:,:), a(:,:), ainv(:,:), b(:,:)
    integer :: m, p, i, j, rank, info
    m = size(prob); p = size(x_cells,2)
    if (m < 2 .or. p < 1 .or. size(x_cells,1) /= m .or. sample_n <= 0.0_dp) error stop 'svy_loglin_prob: invalid shape'
    if (any(shape(prob_vcov) /= [m,m])) error stop 'svy_loglin_prob: covariance shape mismatch'
    if (any(prob < 0.0_dp) .or. abs(sum(prob)-1.0_dp) > 1.0e-5_dp) error stop 'svy_loglin_prob: probabilities must sum to one'
    allocate(weight(m),cluster(m,1),counts(m)); weight = 1.0_dp; counts = sample_n*prob
    do i = 1, m; cluster(i,1) = i; end do
    call make_design(weight,cluster,dcell)
    call svy_glm(x_cells,counts,dcell,fit,FAMILY_POISSON,LINK_LOG,rescale_weights=.false.)
    allocate(pfit(m)); pfit = fit%fitted/sample_n
    if (p == 1) then
      allocate(result%coef(0),result%vcov(0,0))
    else
      allocate(xc(m,p-1),a(p-1,p-1),ainv(p-1,p-1),b(p-1,p-1),result%coef(p-1),result%vcov(p-1,p-1))
      xc = x_cells(:,2:p)
      do j = 1, p-1
        xc(:,j) = xc(:,j)-sum(xc(:,j))/real(m,dp)
      end do
      allocate(pmat(m,m)); pmat = -outer_product(pfit,pfit)/sample_n
      do i = 1, m; pmat(i,i) = pmat(i,i)+pfit(i)/sample_n; end do
      a = matmul(transpose(xc),matmul(pmat,xc))
      call sym_pinv(a,ainv,rank,info=info)
      b = matmul(transpose(xc),matmul(prob_vcov,xc))
      result%vcov = matmul(ainv,matmul(b,ainv))/(sample_n*sample_n)
      result%coef = fit%coef(2:p)
    end if
    result%intercept = fit%coef(1)
    allocate(result%fitted_prob(m),result%cell_prob(m),result%cell_vcov(m,m))
    result%fitted_prob = pfit
    result%cell_prob = prob
    result%cell_vcov = prob_vcov
    result%deviance = fit%deviance
    result%df_residual = max(0,m-p)
    result%df_null = df_null
    result%converged = fit%converged
  end subroutine svy_loglin_prob

  subroutine svy_loglin_compare(prob, prob_vcov, x_small, x_large, sample_n, df_null, test)
    real(dp), intent(in) :: prob(:), prob_vcov(:,:), x_small(:,:), x_large(:,:), sample_n
    integer, intent(in) :: df_null
    type(loglin_test_t), intent(out) :: test
    type(loglin_result_t) :: small, large
    real(dp), allocatable :: x1(:,:),x2(:,:),p1(:,:),psat(:,:),a11(:,:),a12(:,:),a11i(:,:),resid(:,:), &
      basis(:,:), ga(:,:), gb(:,:), eval(:)
    integer :: m,p1c,p2c,i,j,r,rank,info
    call svy_loglin_prob(prob,prob_vcov,x_small,sample_n,df_null,small)
    call svy_loglin_prob(prob,prob_vcov,x_large,sample_n,df_null,large)
    m = size(prob); p1c = size(x_small,2)-1; p2c = size(x_large,2)-1
    if (size(x_small,1) /= m .or. size(x_large,1) /= m .or. p2c <= p1c) error stop 'svy_loglin_compare: models must be nested in increasing dimension'
    allocate(x1(m,max(0,p1c)),x2(m,p2c))
    if (p1c > 0) then
      x1 = x_small(:,2:)
      do j=1,p1c; x1(:,j)=x1(:,j)-sum(x1(:,j))/real(m,dp); end do
    end if
    x2 = x_large(:,2:)
    do j=1,p2c; x2(:,j)=x2(:,j)-sum(x2(:,j))/real(m,dp); end do
    allocate(p1(m,m)); p1=-outer_product(large%fitted_prob,large%fitted_prob)/sample_n
    do i=1,m; p1(i,i)=p1(i,i)+large%fitted_prob(i)/sample_n; end do
    if (p1c > 0) then
      allocate(a11(p1c,p1c),a12(p1c,p2c),a11i(p1c,p1c))
      a11=matmul(transpose(x1),matmul(p1,x1)); a12=matmul(transpose(x1),matmul(p1,x2))
      call sym_pinv(a11,a11i,rank,info=info)
      allocate(resid(m,p2c)); resid=x2-matmul(x1,matmul(a11i,a12))
    else
      allocate(resid(m,p2c)); resid=x2
    end if
    call independent_basis(resid,basis,r)
    if (r < 1) error stop 'svy_loglin_compare: no additional independent terms'
    allocate(psat(m,m)); psat=-outer_product(prob,prob)/sample_n
    do i=1,m; psat(i,i)=psat(i,i)+prob(i)/sample_n; end do
    allocate(ga(r,r),gb(r,r)); ga=matmul(transpose(basis),matmul(psat,basis)); gb=matmul(transpose(basis),matmul(prob_vcov,basis))
    call generalized_symmetric_eigen(ga,gb,eval,info)
    if (info /= 0) error stop 'svy_loglin_compare: generalized eigenvalue failure'
    allocate(test%lambda(size(eval))); test%lambda=max(eval,0.0_dp); test%df=size(eval)
    test%deviance=max(0.0_dp,small%deviance-large%deviance)
    test%score=sample_n*sum((large%fitted_prob-small%fitted_prob)**2/max(small%fitted_prob,tiny(1.0_dp)))
    test%p_deviance_satterthwaite=weighted_chisq_survival(test%deviance,test%lambda,.false.)
    test%p_deviance_saddle=weighted_chisq_survival(test%deviance,test%lambda,.true.)
    test%p_score_satterthwaite=weighted_chisq_survival(test%score,test%lambda,.false.)
    test%p_score_saddle=weighted_chisq_survival(test%score,test%lambda,.true.)
  end subroutine svy_loglin_compare

  subroutine independent_basis(a,basis,r)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: basis(:,:)
    integer, intent(out) :: r
    real(dp), allocatable :: q(:,:),v(:)
    real(dp) :: nv,tol
    integer :: i,j,n,p
    n=size(a,1);p=size(a,2);allocate(q(n,p),v(n));q=0.0_dp;r=0
    tol=1.0e-10_dp*max(1.0_dp,maxval(abs(a)))
    do j=1,p
      v=a(:,j)
      do i=1,r;v=v-dot_product(q(:,i),v)*q(:,i);end do
      nv=sqrt(dot_product(v,v))
      if(nv>tol)then;r=r+1;q(:,r)=v/nv;end if
    end do
    allocate(basis(n,r));if(r>0)basis=q(:,1:r)
  end subroutine independent_basis

  subroutine generalized_symmetric_eigen(a,b,eval,info)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp),allocatable,intent(out)::eval(:)
    integer,intent(out)::info
    real(dp),allocatable::ae(:),av(:,:),invsqrt(:,:),c(:,:),cv(:,:)
    real(dp)::cutoff
    integer::n,i
    n=size(a,1);allocate(ae(n),av(n,n),invsqrt(n,n),c(n,n),cv(n,n),eval(n))
    call symmetric_eigen(a,ae,av,info);if(info/=0)return
    cutoff=1.0e-10_dp*max(1.0_dp,maxval(abs(ae)));invsqrt=0.0_dp
    do i=1,n
      if(ae(i)>cutoff)invsqrt=invsqrt+outer_product(av(:,i),av(:,i))/sqrt(ae(i))
    end do
    c=matmul(invsqrt,matmul(b,invsqrt));call symmetric_eigen(c,eval,cv,info)
  end subroutine generalized_symmetric_eigen

end module survey_loglinear
