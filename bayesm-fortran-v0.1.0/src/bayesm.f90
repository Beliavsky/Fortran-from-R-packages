module bayesm
  use bayesm_kinds, only: dp
  use bayesm_types
  use bayesm_design, only: create_x
  use bayesm_densities, only: lnd_mvn, lnd_mvst, lnd_iwishart, lnd_ichisq
  use bayesm_utils, only: cond_mom, nmat, num_eff, log_marg_den_nr, cget_c
  use bayesm_rng, only: rng_seed, rand_dirichlet, rand_mvst, rwishart_draw
  use bayesm_regression, only: breg, rmultireg, runireg, runireg_gibbs, rsur_gibbs
  use bayesm_mixture, only: rmixture, mix_den, mix_den_bi, mom_mix, e_mix_marg_den, &
    rmix_gibbs, rnmix_gibbs, cluster_mix
  use bayesm_mnl, only: llmnl, llnhlogit, mnl_hess, rmnl_indep_metrop, simnhlogit
  use bayesm_probit, only: rtrun, rbprobit_gibbs, rmnp_gibbs, rmvp_gibbs, ghkvec, &
    llmnp, mnp_prob, rbinorm_gibbs, rordprobit_gibbs
  use bayesm_negbin, only: rnegbin_rw, rhier_negbin_rw
  use bayesm_hierarchical, only: rhier_linear_model, rhier_mnl_rw_mixture, rhier_bin_logit, &
    rhier_linear_mixture, rhier_mnl_dp
  use bayesm_iv, only: riv_gibbs, riv_dp
  use bayesm_dp, only: rdp_gibbs
  use bayesm_scale, only: rscale_usage
  use bayesm_blp, only: rbayes_blp, r2sigma, share2mu, log_jacob
  implicit none
  private
  public :: dp, rng_seed
  public :: create_x, lnd_mvn, lnd_mvst, lnd_iwishart, lnd_ichisq
  public :: cond_mom, nmat, num_eff, log_marg_den_nr, cget_c
  public :: breg, rmultireg, runireg, runireg_gibbs, rsur_gibbs
  public :: rmixture, mix_den, mix_den_bi, mom_mix, e_mix_marg_den, rmix_gibbs, rnmix_gibbs, cluster_mix
  public :: llmnl, llnhlogit, mnl_hess, rmnl_indep_metrop, simnhlogit
  public :: rtrun, rbprobit_gibbs, rmnp_gibbs, rmvp_gibbs, ghkvec, llmnp, mnp_prob, rbinorm_gibbs
  public :: rordprobit_gibbs, rnegbin_rw, rhier_negbin_rw
  public :: rhier_linear_model, rhier_mnl_rw_mixture, rhier_bin_logit, rhier_linear_mixture, rhier_mnl_dp
  public :: riv_gibbs, riv_dp, rdp_gibbs, rscale_usage, rbayes_blp, r2sigma, share2mu, log_jacob
  public :: rdirichlet, rmvst, rwishart
  public :: normal_component, normal_mixture, reg_data, mnl_data, unireg_result, multireg_draw
  public :: mixture_step_result, nmix_result, probit_result, matrix_mcmc_result, hier_linear_result
  public :: hier_mnl_result, hier_negbin_result, negbin_result, iv_result, dp_mixture_result
  public :: sur_result, ordprobit_result, mnl_metrop_result, scale_usage_result, blp_result, wishart_result
contains
  function rdirichlet(n,alpha) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha(:)
    real(dp) :: x(n,size(alpha)),tmp(size(alpha))
    integer :: i
    do i=1,n
      call rand_dirichlet(alpha,tmp)
      x(i,:)=tmp
    end do
  end function rdirichlet

  function rmvst(n,nu,mu,root) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: nu,mu(:),root(:,:)
    real(dp) :: x(n,size(mu)),tmp(size(mu))
    integer :: i
    do i=1,n
      call rand_mvst(nu,mu,root,tmp)
      x(i,:)=tmp
    end do
  end function rmvst

  function rwishart(nu,v) result(out)
    real(dp), intent(in) :: nu,v(:,:)
    type(wishart_result) :: out
    integer :: n,stat
    n=size(v,1)
    allocate(out%w(n,n),out%iw(n,n),out%c(n,n),out%ci(n,n))
    call rwishart_draw(nu,v,out%w,out%iw,out%c,out%ci,stat)
    if (stat/=0) error stop "rwishart: covariance is not positive definite"
  end function rwishart
end module bayesm
