! SPDX-License-Identifier: GPL-3.0-only
! Upstream-style compatibility entry points for HDShOP 0.1.7.
module hdshop_compat
  use hdshop_kinds, only: dp
  use hdshop_shrinkage, only: matrix_shrink_result, mean_shrink_result, &
    cov_shrink_bgp14, inv_cov_shrink_bgp16, nonlin_shrink_lw, &
    sigma_sample_estimator, mean_bs, mean_js, mean_bop19
  use hdshop_portfolio, only: portfolio_result, frontier_result, &
    mean_variance_portfolio, traditional_portfolio, shrinkage_mv_portfolio, &
    shrinkage_gmv_portfolio, mv_shrink_portfolio, bayesian_frontier
  use hdshop_inference, only: mvsp_test_result, test_mvsp
  use hdshop_random, only: random_covariance_matrix
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private
  public :: covshrinkbgp14, invcovshrinkbgp16, nonlin_shrinklw
  public :: covarestim, meanestim, randcovmtrx
  public :: new_meanvar_portfolio, new_mv_portfolio_traditional
  public :: new_mv_portfolio_traditional_pgn
  public :: new_mv_portfolio_weights_bdops21
  public :: new_mv_portfolio_weights_bdops21_pgn
  public :: new_gmv_portfolio_weights_bdps19
  public :: new_gmv_portfolio_weights_bdps19_pgn
  public :: mvshrinkportfolio, frontier_data, meanvar_portfolio
  public :: validate_meanvar_portfolio, plot_frontier

contains

  function covshrinkbgp14(n,tm,scm) result(res)
    integer,intent(in)::n
    real(dp),intent(in)::tm(:,:),scm(:,:)
    type(matrix_shrink_result)::res
    res=cov_shrink_bgp14(n,tm,scm)
  end function covshrinkbgp14

  function invcovshrinkbgp16(n,p,tm,iscm) result(res)
    integer,intent(in)::n,p
    real(dp),intent(in)::tm(:,:),iscm(:,:)
    type(matrix_shrink_result)::res
    if(p/=size(iscm,1))then
      allocate(res%matrix(size(iscm,1),size(iscm,2)));res%matrix=0.0_dp
      res%message='p does not match inverse covariance dimension'
    else
      res=inv_cov_shrink_bgp16(n,tm,iscm)
    end if
  end function invcovshrinkbgp16

  function nonlin_shrinklw(x) result(s)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable::s(:,:)
    s=nonlin_shrink_lw(x)
  end function nonlin_shrinklw

  function covarestim(x,kind,target) result(s)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in)::kind
    real(dp),intent(in),optional::target(:,:)
    real(dp),allocatable::s(:,:)
    type(matrix_shrink_result)::fit
    if(trim(kind)=='trad')then
      s=sigma_sample_estimator(x)
    else if(trim(kind)=='LW20')then
      s=nonlin_shrink_lw(x)
    else if(trim(kind)=='BGP14' .and. present(target))then
      fit=cov_shrink_bgp14(size(x,2),target,sigma_sample_estimator(x));s=fit%matrix
    else
      allocate(s(size(x,1),size(x,1)));s=0.0_dp
    end if
  end function covarestim

  function meanestim(x,kind,y0,mu0) result(means)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in)::kind
    real(dp),intent(in),optional::y0,mu0(:)
    real(dp),allocatable::means(:)
    type(mean_shrink_result)::fit
    if(trim(kind)=='trad')then
      means=sum(x,dim=2)/real(size(x,2),dp)
    else if(trim(kind)=='bs')then
      fit=mean_bs(x);means=fit%means
    else if(trim(kind)=='js')then
      fit=mean_js(x,y0);means=fit%means
    else if(trim(kind)=='BOP19' .and. present(mu0))then
      fit=mean_bop19(x,mu0);means=fit%means
    else
      allocate(means(size(x,1)));means=0.0_dp
    end if
  end function meanestim

  function randcovmtrx(p,eigenvalues,seed) result(s)
    integer,intent(in)::p
    real(dp),intent(in),optional::eigenvalues(:)
    integer(int64),intent(in),optional::seed
    real(dp),allocatable::s(:,:)
    s=random_covariance_matrix(p,eigenvalues,seed)
  end function randcovmtrx

  function new_meanvar_portfolio(mean_vec,cov_mtrx,gamma) result(res)
    real(dp),intent(in)::mean_vec(:),cov_mtrx(:,:),gamma
    type(portfolio_result)::res
    res=mean_variance_portfolio(mean_vec,cov_mtrx,gamma)
  end function new_meanvar_portfolio

  function new_mv_portfolio_traditional(x,gamma) result(res)
    real(dp),intent(in)::x(:,:),gamma
    type(portfolio_result)::res
    res=traditional_portfolio(x,gamma)
  end function new_mv_portfolio_traditional

  function new_mv_portfolio_traditional_pgn(x,gamma) result(res)
    real(dp),intent(in)::x(:,:),gamma
    type(portfolio_result)::res
    res=traditional_portfolio(x,gamma)
  end function new_mv_portfolio_traditional_pgn

  function new_mv_portfolio_weights_bdops21(x,gamma,b,beta) result(res)
    real(dp),intent(in)::x(:,:),gamma,b(:),beta
    type(portfolio_result)::res
    res=shrinkage_mv_portfolio(x,gamma,b,beta)
  end function new_mv_portfolio_weights_bdops21

  function new_mv_portfolio_weights_bdops21_pgn(x,gamma,b,beta) result(res)
    real(dp),intent(in)::x(:,:),gamma,b(:),beta
    type(portfolio_result)::res
    res=shrinkage_mv_portfolio(x,gamma,b,beta)
  end function new_mv_portfolio_weights_bdops21_pgn

  function new_gmv_portfolio_weights_bdps19(x,b,beta) result(res)
    real(dp),intent(in)::x(:,:),b(:),beta
    type(portfolio_result)::res
    res=shrinkage_gmv_portfolio(x,b,beta)
  end function new_gmv_portfolio_weights_bdps19

  function new_gmv_portfolio_weights_bdps19_pgn(x,b,beta) result(res)
    real(dp),intent(in)::x(:,:),b(:),beta
    type(portfolio_result)::res
    res=shrinkage_gmv_portfolio(x,b,beta)
  end function new_gmv_portfolio_weights_bdps19_pgn

  function mvshrinkportfolio(x,gamma,kind,b,beta) result(res)
    real(dp),intent(in)::x(:,:),gamma
    character(len=*),intent(in)::kind
    real(dp),intent(in),optional::b(:),beta
    type(portfolio_result)::res
    res=mv_shrink_portfolio(x,gamma,kind,b,beta)
  end function mvshrinkportfolio

  function frontier_data(x,weights,npoints) result(res)
    real(dp),intent(in)::x(:,:),weights(:,:)
    integer,intent(in),optional::npoints
    type(frontier_result)::res
    res=bayesian_frontier(x,weights,npoints)
  end function frontier_data

  function meanvar_portfolio(mean_vec,cov_mtrx,gamma) result(res)
    real(dp),intent(in)::mean_vec(:),cov_mtrx(:,:),gamma
    type(portfolio_result)::res
    res=mean_variance_portfolio(mean_vec,cov_mtrx,gamma)
  end function meanvar_portfolio

  logical function validate_meanvar_portfolio(res) result(valid)
    type(portfolio_result),intent(in)::res
    valid=res%ok .and. allocated(res%covariance) .and. allocated(res%inverse_covariance) .and. &
      allocated(res%means) .and. allocated(res%weights)
    if(valid)valid=size(res%covariance,1)==size(res%covariance,2) .and. &
      size(res%means)==size(res%weights) .and. size(res%covariance,1)==size(res%weights)
  end function validate_meanvar_portfolio

  function plot_frontier(x,weights,npoints) result(res)
    real(dp),intent(in)::x(:,:),weights(:,:)
    integer,intent(in),optional::npoints
    type(frontier_result)::res
    res=bayesian_frontier(x,weights,npoints)
  end function plot_frontier

end module hdshop_compat
