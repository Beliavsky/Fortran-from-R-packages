program test_hmisc
  use hmisc
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none
  real(dp) :: x3(3), w3(3), q, c, dxy
  real(dp) :: mat(5,2), r(2,2), hd(2,2), aad(2,2), mad(2,2)
  real(dp) :: nrel, nconc, nuncert, ci, gamma, sd
  real(dp) :: bottom, top, pow, bn1, bn2, aa, bb, mu, sigma, nri_e, nri_ne, nri
  real(dp) :: fvd, fvr, fvlo, fmp
  real(dp) :: gammap, sdp, c12, c21, c1, c2, g1, g2
  real(dp) :: vals(4), sums(4), ranks(4), refs(3), ask(2), draws(2)
  integer :: np(2,2), nh(2,2), nu, groups(7), idx(2)
  integer :: y3(3)
  real(dp) :: po_power, po_eff, po_se, po_n, po_eff2
  real(dp) :: pmod(3), mv(2), sout(3), mres(2)
  real(dp) :: sm, ssd, slo, shi, med, jpow, cipow
  real(dp) :: aout(5)
  real(dp), allocatable :: ex(:), ef(:), jmean(:), jshr(:), jshrink(:)
  integer, allocatable :: jids(:), jn(:), ccat(:,:)
  real(dp), allocatable :: spline_basis(:,:), gx(:), gy(:), xuniq(:), yuniq(:)
  real(dp), allocatable :: ecx(:), ecy(:), nomv(:), nomm(:,:), filled(:)
  real(dp) :: lagout(4), qnan
  real(dp) :: drho, ddeff, rr2, rr2a, rr2e, rr2ae, best, blo, bhi
  real(dp) :: scorem(3), scores(3)
  real(dp), allocatable :: invout(:)
  integer :: gstatus, ftu_n1, ftu_n2, dn, dk
  real(dp) :: mh_rr, mh_lo, mh_hi
  integer :: mh_counts(2), rm_status, qr_status
  real(dp) :: lp(2), lpl(2), lpu(2), ln(2), lnl(2), lnu(2)
  real(dp) :: lpc(2), lpcl(2), lpcu(2), lnc(2), lncl(2), lncu(2)
  integer, allocatable :: nd_subs(:), nd_y(:)
  real(dp), allocatable :: nd_weights(:)
  integer :: rm_draws(2,3)
  real(dp) :: qx(4,2), qxt(4,2), qr(2,2), qri(2,2), qbar(2), qrec(4,2)
  logical :: ev(3)
  real(dp) :: sv_a(3,2), sv_b(3), sv_x(2), sv_inv(2,2), smn(3)
  integer :: sv_status, wk_status
  integer :: wk_idx(2)
  real(dp), allocatable :: pc_scores(:), pc_coef(:)
  real(dp) :: pc_frac
  integer :: pc_status
  real(dp), allocatable :: lm_coef(:), lm_res(:), lm_fit(:)
  real(dp) :: lm_r2, aep_diff(3,2), aep_rat(2,2)
  real(dp) :: rcs_c, rcs_d, rcs_ad, rcs_sd, rcs_z, rcs_p
  integer :: lm_status, aep_status, rcs_status
  real(dp) :: lr_chisq, lr_hr, sp_rho2, sp_f, sp_df1, sp_df2, sp_p, sp_adj
  integer :: lr_status, sp_n, sp_status, eb_status, jit_status
  integer, allocatable :: eb(:)
  real(dp) :: jit_out(5)
  real(dp) :: hdq(3), ct_stat, ct_df, ct_p, kw_stat, kw_df1, kw_df2, kw_p
  real(dp) :: tt_mean(2), tt_rho, tt_corr(2), tt_var, tt_vu, tt_de, tt_eff, tt_se
  real(dp) :: tt_ci(2), tt_z, tt_p
  integer :: hdq_status, ct_status, kw_status, tt_status
  integer :: tt_n(2), tt_nc(2)
  real(dp), allocatable :: km_t(:), km_s(:), km_q(:)
  integer, allocatable :: bs_labels(:), bs_n(:), bs_miss(:)
  integer, allocatable :: b2_vlab(:), b2_hlab(:), b2_n(:,:), b2_miss(:,:)
  real(dp), allocatable :: bs_mean(:,:), b2_mean(:,:,:)
  integer :: km_status, bs_status, b2_status
  integer, allocatable :: fm_match(:,:), fm_n(:), sf_assign(:), sf_order(:), sf_cond(:)
  real(dp), allocatable :: fm_dist(:,:)
  integer, allocatable :: pu_major(:), pu_minor(:)
  real(dp), allocatable :: pu_diff(:), pu_mid(:), pu_sd(:), pu_lo(:), pu_hi(:), pu_lm(:), pu_um(:)
  integer :: fm_status, sf_status, pu_status, boot_status, smear_status
  integer :: bps_status, sps_status, invf_status
  integer :: bz_status, cut_status, poc_status, rcf_status, rcr_status
  real(dp) :: bps_power, bps_lo, bps_hi, sps_power, sps_maxf, sps_maxc
  real(dp), allocatable :: invf_roots(:,:)
  real(dp), allocatable :: bz_x(:), bz_y(:), cut_eff(:), poc_or(:,:)
  real(dp), allocatable :: rcf_y(:), rcr_coef(:)
  real(dp) :: st_r2, st_f, st_df1, st_df2, st_p, cs_stat, cs_df, cs_excess, cs_p
  integer :: st_n, st_status, cs_status, cl_status, pm_status
  integer, allocatable :: cl_out(:), cl_map(:)
  real(dp), allocatable :: pm_scores(:,:), pm_sload(:,:), pm_oload(:,:), pm_var(:), pm_scale(:)
  real(dp), allocatable :: bkm_est(:), mc_score(:,:)
  integer, allocatable :: mc_matches(:,:), mc_nmatch(:)
  real(dp) :: le_center(2), le_rect(4), le_area
  integer :: bkm_status, mc_status, le_status
  integer, allocatable :: cut_group(:)
  integer, allocatable :: invf_nroots(:)
  real(dp) :: boot_mean, boot_lo, boot_hi, smear_out(2)
  real(dp) :: gb_mean_post, gb_var_post, gb_mean_pred, gb_var_pred
  real(dp) :: gb_pred_d, gb_pred_c, gb_post_d, gb_post_c, gb_post_m, gb_power
  real(dp) :: gb_mix_critical, gb_mix_power
  integer :: gb_status
  real(dp) :: gb2_prob
  integer :: gb2_status, ppo_status, ptr_status, yn_status
  real(dp), allocatable :: ppo_obs(:,:), ppo_pred(:,:), ptr_prop(:,:,:), yn_y(:,:), yn_prop(:)
  integer, allocatable :: ptr_count(:,:,:), yn_order(:)
  real(dp), allocatable :: boot_reps(:)
  integer, allocatable :: rf_order(:), rf_varmiss(:), rf_obsmiss(:)
  real(dp), allocatable :: ma_mean(:,:), so_prob(:,:)
  real(dp), allocatable :: og_grouped(:), clo_x(:), clo_y(:), cs_x(:), cs_y(:)
  integer, allocatable :: cs_id(:)
  integer, allocatable :: ma_levels(:), ma_count(:,:)
  integer :: rf_nimp, rf_status, ma_status, so_status, og_m, og_status, clo_status, cvs_status
  logical :: rf_missing(5,3), so_absorb(3)
  real(dp) :: ma_x(4,2), so_off(2,3,2)
  real(dp), allocatable :: wl_fit(:)
  real(dp) :: sm_off(2,2,3,2), sm_u(2,2)
  integer, allocatable :: sm_states(:,:)
  integer :: wl_status, smk_status
  logical :: sm_absorb(3)
  real(dp), allocatable :: imp_y(:), mov_x(:), mov_mean(:), mov_q1(:), mov_med(:), mov_q3(:)
  integer, allocatable :: imp_idx(:), mov_n(:), sum_levels(:,:), sum_counts(:,:)
  real(dp), allocatable :: sum_means(:,:)
  integer :: imp_status, mov_status, sum_status
  real(dp), allocatable :: at_out(:,:), al_coef(:), al_fit(:), al_res(:)
  real(dp), allocatable :: dfr_red(:,:), dr_unique(:,:), vc_sim(:,:), rd_r2(:)
  real(dp), allocatable :: gb_prob(:,:), gb_pm(:,:), gb_ps(:,:), mo_occ(:,:)
  integer, allocatable :: dfr_keep(:), dfr_reason(:), dr_counts(:), vc_np(:,:)
  integer, allocatable :: vc_left(:), vc_right(:), rd_keep(:), rd_remove(:)
  real(dp), allocatable :: vc_height(:)
  real(dp) :: al_r2, al_mae, al_medae
  integer :: at_status, al_status, dfr_status, dr_status, dr_n, vc_status, rd_status, gbs_status, mo_status
  integer :: failures

  failures = 0
  x3 = [1.0_dp, 2.0_dp, 3.0_dp]
  w3 = [1.0_dp, 1.0_dp, 1.0_dp]
  call check_close(gini_md(x3), 4.0_dp/3.0_dp, 1.0e-12_dp, 'GiniMd', failures)
  call check_close(trap_rule([0.0_dp,1.0_dp,2.0_dp],[0.0_dp,1.0_dp,4.0_dp]), 3.0_dp, 1.0e-12_dp, 'trap.rule', failures)
  call check_close(pseudomedian(x3), 2.0_dp, 1.0e-12_dp, 'pMedian', failures)
  call check_close(weighted_mean(x3,w3), 2.0_dp, 1.0e-12_dp, 'wtd.mean', failures)
  call check_close(weighted_variance(x3,w3,.false.,.false.), 1.0_dp, 1.0e-12_dp, 'wtd.var', failures)
  q = weighted_quantile(x3,w3,0.5_dp,.false.)
  call check_close(q, 2.0_dp, 1.0e-12_dp, 'wtd.quantile', failures)

  vals = 0.0_dp
  sums = 0.0_dp
  call weighted_table([1.0_dp,1.0_dp,2.0_dp,4.0_dp],[1.0_dp,2.0_dp,1.0_dp,1.0_dp],vals,sums,nu,.false.)
  call check_int(nu, 3, 'wtd.table n', failures)
  call check_close(sums(1), 3.0_dp, 1.0e-12_dp, 'wtd.table weight', failures)
  call weighted_rank([1.0_dp,1.0_dp,2.0_dp,4.0_dp],[1.0_dp,2.0_dp,1.0_dp,1.0_dp],ranks,.false.)
  call check_close(ranks(1), 2.0_dp, 1.0e-12_dp, 'wtd.rank tie', failures)

  y3 = [0,0,1]
  call somers2(x3,y3,w3,c,dxy)
  call check_close(c, 1.0_dp, 1.0e-12_dp, 'somers2 C', failures)
  call check_close(dxy, 1.0_dp, 1.0e-12_dp, 'somers2 Dxy', failures)

  mat(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  mat(:,2) = [2.0_dp,4.0_dp,6.0_dp,8.0_dp,10.0_dp]
  call rcorr(mat,.false.,r,np)
  call check_close(r(1,2),1.0_dp,1.0e-12_dp,'rcorr Pearson',failures)
  call rcorr(mat,.true.,r,np)
  call check_close(r(1,2),1.0_dp,1.0e-12_dp,'rcorr Spearman',failures)
  call hoeffd(mat,hd,aad,mad,nh)
  call check_int(nh(1,2),5,'hoeffd n',failures)
  if (ieee_is_nan(hd(1,2))) then
    print *, 'FAIL hoeffd finite'
    failures = failures + 1
  end if

  ev = [.true.,.true.,.true.]
  call rcorr_cens(x3,x3,ev,.false.,nrel,nconc,nuncert,ci,gamma,sd)
  call check_close(ci,1.0_dp,1.0e-12_dp,'rcorr.cens C',failures)
  call check_close(gamma,1.0_dp,1.0e-12_dp,'rcorr.cens Dxy',failures)
  call rcorrp_cens(x3,-x3,x3,ev,2,.false.,gammap,sdp,c12,c21,c1,c2,g1,g2,nrel,nuncert)
  if (ieee_is_nan(gammap)) then
    print *, 'FAIL rcorrp.cens finite'
    failures = failures + 1
  end if

  call cut_gn([1.0_dp,1.0_dp,2.0_dp,3.0_dp,3.0_dp,4.0_dp,5.0_dp],2,groups)
  call check_int(groups(1),1,'cutGn first',failures)
  call check_int(groups(7),3,'cutGn last',failures)

  refs = [0.0_dp, 10.0_dp, 20.0_dp]
  ask = [9.0_dp, 18.0_dp]
  call which_closest(refs,ask,idx)
  call check_int(idx(1),2,'whichClosest 1',failures)
  call check_int(idx(2),3,'whichClosest 2',failures)
  draws = [0.1_dp,0.8_dp]
  call which_close_pw(refs,ask,draws,1.0_dp,idx)
  if (any(idx < 1) .or. any(idx > 3)) then
    print *, 'FAIL whichClosePW bounds'
    failures = failures + 1
  end if

  if (.not. (samplesize_bin(0.025_dp,0.8_dp,0.7_dp,0.5_dp,0.5_dp) > 0.0_dp)) then
    print *, 'FAIL samplesize.bin'
    failures = failures + 1
  end if

  call dual_sd([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp], 3, 3.0_dp, bottom, top)
  call check_close(bottom, sqrt(2.5_dp), 1.0e-12_dp, 'dualSD bottom', failures)
  call check_close(top, sqrt(2.5_dp), 1.0e-12_dp, 'dualSD top', failures)

  pow = bpower(0.5_dp, 0.7_dp, 100.0_dp, 100.0_dp, 0.05_dp)
  call check_close(pow, 0.82810977_dp, 2.0e-7_dp, 'bpower', failures)
  call bsamsize(0.5_dp, 0.7_dp, 0.5_dp, 0.05_dp, 0.8_dp, bn1, bn2)
  call check_close(bn1, 92.998845_dp, 2.0e-5_dp, 'bsamsize n1', failures)
  call check_close(bn2, 92.998845_dp, 2.0e-5_dp, 'bsamsize n2', failures)

  call ballocation(0.5_dp, 0.7_dp, fvd, fvr, fvlo, 200.0_dp, fmp)
  call check_close(fvlo, 1.0_dp-fvd, 1.0e-12_dp, 'ballocation logodds', failures)
  if (fmp <= 0.0_dp .or. fmp >= 1.0_dp) then
    print *, 'FAIL ballocation maxpower'
    failures = failures + 1
  end if

  call weibull2_fit([1.0_dp,2.0_dp], [0.8_dp,0.5_dp], aa, bb)
  call check_close(weibull_survival(1.0_dp,aa,bb), 0.8_dp, 1.0e-12_dp, 'Weibull2 t1', failures)
  call check_close(weibull_survival(2.0_dp,aa,bb), 0.5_dp, 1.0e-12_dp, 'Weibull2 t2', failures)

  call gompertz2_fit([1.0_dp,2.0_dp], [0.8_dp,0.5_dp], aa, bb)
  call check_close(gompertz_survival(1.0_dp,aa,bb), 0.8_dp, 1.0e-12_dp, 'Gompertz2 t1', failures)
  call check_close(gompertz_survival(2.0_dp,aa,bb), 0.5_dp, 1.0e-12_dp, 'Gompertz2 t2', failures)

  call lognorm2_fit([1.0_dp,2.0_dp], [0.8_dp,0.5_dp], mu, sigma)
  call check_close(lognorm_survival(1.0_dp,mu,sigma), 0.8_dp, 2.0e-8_dp, 'Lognorm2 t1', failures)
  call check_close(lognorm_survival(2.0_dp,mu,sigma), 0.5_dp, 2.0e-8_dp, 'Lognorm2 t2', failures)

  call improve_prob([0.2_dp,0.6_dp,0.3_dp,0.7_dp], [0.3_dp,0.5_dp,0.2_dp,0.8_dp], &
                    [1,1,0,0], nri_e, nri_ne, nri)
  call check_close(nri_e, 0.0_dp, 1.0e-12_dp, 'improveProb events', failures)
  call check_close(nri_ne, 0.0_dp, 1.0e-12_dp, 'improveProb nonevents', failures)
  call check_close(nri, 0.0_dp, 1.0e-12_dp, 'improveProb total', failures)


  call popower([0.2_dp,0.3_dp,0.5_dp], 1.5_dp, 100.0_dp, 100.0_dp, 0.05_dp, &
               po_power, po_eff, po_se)
  call check_close(po_power, 0.3262109193_dp, 3.0e-8_dp, 'popower power', failures)
  call check_close(po_eff, 0.8400210005_dp, 2.0e-10_dp, 'popower efficiency', failures)
  call check_close(po_se, 0.2685975481_dp, 2.0e-10_dp, 'popower se', failures)
  call posamsize([0.2_dp,0.3_dp,0.5_dp], 1.5_dp, 0.5_dp, 0.05_dp, 0.8_dp, po_n, po_eff2)
  call check_close(po_n, 682.028715_dp, 3.0e-5_dp, 'posamsize n', failures)
  call check_close(po_eff2, 0.8400018058_dp, 3.0e-10_dp, 'posamsize efficiency', failures)
  call pomodm([0.2_dp,0.3_dp,0.5_dp], 1.5_dp, pmod)
  call check_close(pmod(1), 1.0_dp/7.0_dp, 1.0e-12_dp, 'pomodm p1', failures)
  call check_close(pmod(2), 0.257142857142857_dp, 1.0e-12_dp, 'pomodm p2', failures)
  call check_close(pmod(3), 0.6_dp, 1.0e-12_dp, 'pomodm p3', failures)

  call rcspline_eval([0.0_dp,1.0_dp,2.0_dp,3.0_dp], [0.0_dp,1.0_dp,2.0_dp], &
                     .true., .false., 0, spline_basis)
  call check_int(size(spline_basis,2), 2, 'rcspline.eval columns', failures)
  call check_close(spline_basis(4,1), 3.0_dp, 1.0e-12_dp, 'rcspline.eval x', failures)
  call check_close(spline_basis(2,2), 1.0_dp, 1.0e-12_dp, 'rcspline.eval basis1', failures)
  call check_close(spline_basis(3,2), 6.0_dp, 1.0e-12_dp, 'rcspline.eval basis2', failures)
  call check_close(spline_basis(4,2), 12.0_dp, 1.0e-12_dp, 'rcspline.eval basis3', failures)

  call groupn([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp], &
              [2.0_dp,4.0_dp,6.0_dp,8.0_dp,10.0_dp], 2, gx, gy, gstatus)
  call check_int(gstatus, 0, 'groupn status', failures)
  call check_int(size(gx), 3, 'groupn groups', failures)
  call check_close(gx(1), 1.5_dp, 1.0e-12_dp, 'groupn first x', failures)
  call check_close(gx(3), 4.5_dp, 1.0e-12_dp, 'groupn last x', failures)
  call check_close(gy(3), 9.0_dp, 1.0e-12_dp, 'groupn last y', failures)

  call check_close(spearman_rho([1.0_dp,2.0_dp,3.0_dp], [3.0_dp,1.0_dp,2.0_dp]), &
                   -0.5_dp, 1.0e-12_dp, 'spearman', failures)
  call stepfun_eval([1.0_dp,2.0_dp,3.0_dp], [10.0_dp,20.0_dp,30.0_dp], &
                    [1.5_dp,2.0_dp,2.5_dp], .false., sout)
  call check_close(sout(1), 10.0_dp, 1.0e-12_dp, 'stepfun.eval left 1', failures)
  call check_close(sout(2), 20.0_dp, 1.0e-12_dp, 'stepfun.eval exact', failures)
  call check_close(sout(3), 20.0_dp, 1.0e-12_dp, 'stepfun.eval left 2', failures)
  call stepfun_eval([1.0_dp,2.0_dp,3.0_dp], [10.0_dp,20.0_dp,30.0_dp], &
                    [1.5_dp,2.0_dp,2.5_dp], .true., sout)
  call check_close(sout(1), 20.0_dp, 1.0e-12_dp, 'stepfun.eval right 1', failures)
  call check_close(sout(3), 30.0_dp, 1.0e-12_dp, 'stepfun.eval right 2', failures)

  call xy_sort_no_dup_no_na([2.0_dp,1.0_dp,2.0_dp], [20.0_dp,10.0_dp,21.0_dp], xuniq, yuniq)
  call check_int(size(xuniq), 2, 'xySortNoDupNoNA length', failures)
  call check_close(xuniq(1), 1.0_dp, 1.0e-12_dp, 'xySortNoDupNoNA x1', failures)
  call check_close(yuniq(2), 20.0_dp, 1.0e-12_dp, 'xySortNoDupNoNA y2', failures)
  call check_int(n_coincident([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,4.0_dp,5.0_dp], &
                             [0.1_dp,0.2_dp,0.3_dp,0.4_dp,0.5_dp,0.4_dp,0.5_dp], 400), &
                 2, 'nCoincident', failures)

  mv = [2.0_dp, 3.0_dp]
  call matxv(reshape([1.0_dp,3.0_dp,2.0_dp,4.0_dp],[2,2]), mv, mres)
  call check_close(mres(1), 8.0_dp, 1.0e-12_dp, 'matxv row1', failures)
  call check_close(mres(2), 18.0_dp, 1.0e-12_dp, 'matxv row2', failures)
  call matxv(reshape([1.0_dp,3.0_dp,2.0_dp,4.0_dp],[2,2]), [5.0_dp,2.0_dp,3.0_dp], mres, 1)
  call check_close(mres(1), 13.0_dp, 1.0e-12_dp, 'matxv intercept row1', failures)


  jpow = cpower(5.0_dp, 1000.0_dp, 0.2_dp, 25.0_dp, 2.0_dp, 3.0_dp)
  call check_close(jpow, 0.4650106284_dp, 3.0e-8_dp, 'cpower', failures)
  cipow = ciapower(5.0_dp, 800.0_dp, 700.0_dp, 0.2_dp, 0.3_dp, 20.0_dp, 10.0_dp, 2.0_dp, 3.0_dp)
  call check_close(cipow, 0.0791755626_dp, 3.0e-8_dp, 'ciapower', failures)

  call wtd_ecdf([1.0_dp,1.0_dp,2.0_dp,4.0_dp], [1.0_dp,2.0_dp,1.0_dp,1.0_dp], 1, ex, ef)
  call check_int(size(ex), 4, 'wtd.Ecdf length', failures)
  call check_close(ex(1), 1.0_dp, 1.0e-12_dp, 'wtd.Ecdf leading x', failures)
  call check_close(ef(1), 0.0_dp, 1.0e-12_dp, 'wtd.Ecdf leading zero', failures)
  call check_close(ef(2), 0.6_dp, 1.0e-12_dp, 'wtd.Ecdf first cdf', failures)
  call check_close(ef(4), 1.0_dp, 1.0e-12_dp, 'wtd.Ecdf last cdf', failures)

  call smean_sd([1.0_dp,2.0_dp,3.0_dp,4.0_dp], sm, ssd)
  call check_close(sm, 2.5_dp, 1.0e-12_dp, 'smean.sd mean', failures)
  call check_close(ssd, sqrt(5.0_dp/3.0_dp), 1.0e-12_dp, 'smean.sd sd', failures)
  call smean_sdl([1.0_dp,2.0_dp,3.0_dp,4.0_dp], 1.0_dp, sm, slo, shi)
  call check_close(slo, 2.5_dp-sqrt(5.0_dp/3.0_dp), 1.0e-12_dp, 'smean.sdl lower', failures)
  call check_close(shi, 2.5_dp+sqrt(5.0_dp/3.0_dp), 1.0e-12_dp, 'smean.sdl upper', failures)
  call smedian_hilow([1.0_dp,2.0_dp,3.0_dp,4.0_dp], 0.5_dp, med, slo, shi)
  call check_close(med, 2.5_dp, 1.0e-12_dp, 'smedian.hilow median', failures)
  call check_close(slo, 1.75_dp, 1.0e-12_dp, 'smedian.hilow lower', failures)
  call check_close(shi, 3.25_dp, 1.0e-12_dp, 'smedian.hilow upper', failures)

  call approx_extrap([0.0_dp,1.0_dp,2.0_dp], [0.0_dp,1.0_dp,4.0_dp], &
                     [-1.0_dp,0.5_dp,1.5_dp,2.0_dp,3.0_dp], aout)
  call check_close(aout(1), -1.0_dp, 1.0e-12_dp, 'approxExtrap left', failures)
  call check_close(aout(2), 0.5_dp, 1.0e-12_dp, 'approxExtrap interpolate1', failures)
  call check_close(aout(3), 2.5_dp, 1.0e-12_dp, 'approxExtrap interpolate2', failures)
  call check_close(aout(5), 7.0_dp, 1.0e-12_dp, 'approxExtrap right', failures)

  call james_stein([1.0_dp,2.0_dp,3.0_dp,2.0_dp,4.0_dp,6.0_dp,5.0_dp,5.0_dp,5.0_dp,8.0_dp,9.0_dp,10.0_dp], &
                   [1,1,1,2,2,2,3,3,3,4,4,4], jids, jn, jmean, jshr, jshrink)
  call check_int(size(jids), 4, 'james.stein groups', failures)
  call check_close(jmean(4), 9.0_dp, 1.0e-12_dp, 'james.stein mean', failures)
  call check_close(jshrink(1), 77.0_dp/78.0_dp, 1.0e-12_dp, 'james.stein shrink1', failures)


  call check_close(ftupwr(0.4_dp, 0.6_dp, 300.0_dp, 1.0_dp, 0.05_dp), &
                   0.921505383459_dp, 2.0e-8_dp, 'ftupwr', failures)
  call ftuss(0.4_dp, 0.6_dp, 1.0_dp, 0.05_dp, 0.2_dp, ftu_n1, ftu_n2)
  call check_int(ftu_n1, 107, 'ftuss n1', failures)
  call check_int(ftu_n2, 107, 'ftuss n2', failures)

  qnan = ieee_value(0.0_dp, ieee_quiet_nan)
  call ecdf_steps([1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, qnan], ecx, ecy, do_extend=.false.)
  call check_int(size(ecx), 3, 'ecdfSteps unique count', failures)
  call check_close(ecx(2), 2.0_dp, 1.0e-12_dp, 'ecdfSteps x', failures)
  call check_close(ecy(1), 0.5_dp, 1.0e-12_dp, 'ecdfSteps first probability', failures)
  call check_close(ecy(3), 1.0_dp, 1.0e-12_dp, 'ecdfSteps last probability', failures)
  call ecdf_steps([0.0_dp, 10.0_dp], ecx, ecy)
  call check_close(ecx(1), -0.5_dp, 1.0e-12_dp, 'ecdfSteps default lower extension', failures)
  call check_close(ecx(4), 10.5_dp, 1.0e-12_dp, 'ecdfSteps default upper extension', failures)

  call lag_numeric([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 1, lagout)
  if (.not. ieee_is_nan(lagout(1))) then
    print *, 'FAIL Lag leading NaN'
    failures = failures + 1
  end if
  call check_close(lagout(4), 3.0_dp, 1.0e-12_dp, 'Lag positive shift', failures)
  call lag_numeric([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], -1, lagout)
  call check_close(lagout(1), 2.0_dp, 1.0e-12_dp, 'Lag negative shift', failures)
  if (.not. ieee_is_nan(lagout(4))) then
    print *, 'FAIL Lag trailing NaN'
    failures = failures + 1
  end if

  call nomiss([1.0_dp, qnan, 3.0_dp], nomv)
  call check_int(size(nomv), 2, 'nomiss vector size', failures)
  call check_close(nomv(2), 3.0_dp, 1.0e-12_dp, 'nomiss vector order', failures)
  call nomiss(reshape([1.0_dp, 2.0_dp, qnan, 4.0_dp, 5.0_dp, 6.0_dp], [3,2]), nomm)
  call check_int(size(nomm,1), 2, 'nomiss matrix rows', failures)
  call check_close(nomm(2,2), 5.0_dp, 1.0e-12_dp, 'nomiss matrix value', failures)

  allocate(filled(3))
  call fillin([1.0_dp, qnan, 3.0_dp], 9.0_dp, filled)
  call check_close(filled(2), 9.0_dp, 1.0e-12_dp, 'fillin scalar', failures)
  call fillin([qnan, 2.0_dp, qnan], [7.0_dp, 8.0_dp], filled)
  call check_close(filled(1), 7.0_dp, 1.0e-12_dp, 'fillin vector first', failures)
  call check_close(filled(3), 7.0_dp, 1.0e-12_dp, 'fillin vector recycle', failures)

  call cumcategory([1, 2, 3, 2], 3, ccat)
  call check_int(size(ccat,2), 2, 'cumcategory columns', failures)
  call check_int(ccat(1,1), 0, 'cumcategory low', failures)
  call check_int(ccat(3,1), 1, 'cumcategory middle indicator', failures)
  call check_int(ccat(3,2), 1, 'cumcategory high indicator', failures)
  call check_close(jshr(2), 4.051282051282051_dp, 1.0e-12_dp, 'james.stein shrunk2', failures)

  call deff([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp], [1,1,2,2,3,3], &
            dn, dk, drho, ddeff)
  call check_int(dn, 6, 'deff n', failures)
  call check_int(dk, 3, 'deff clusters', failures)
  call check_close(drho, 29.0_dp/35.0_dp, 1.0e-12_dp, 'deff rho', failures)
  call check_close(ddeff, 64.0_dp/35.0_dp, 1.0e-12_dp, 'deff design effect', failures)

  call r2_measures(138.6267_dp, 1, 100.0_dp, 1, rr2, rr2a, frequencies=[50.0_dp,50.0_dp], &
                   r2_ess=rr2e, r2adj_ess=rr2ae)
  call check_close(rr2, 0.7499931596264476_dp, 1.0e-12_dp, 'R2Measures R2', failures)
  call check_close(rr2a, 0.7474805491085085_dp, 1.0e-12_dp, 'R2Measures adjusted', failures)
  call check_close(rr2e, 0.8425041231831880_dp, 1.0e-12_dp, 'R2Measures effective', failures)
  call check_close(rr2ae, 0.8403901160969194_dp, 1.0e-12_dp, 'R2Measures effective adjusted', failures)

  call binconf_wilson(20.0_dp, 100.0_dp, 0.05_dp, best, blo, bhi)
  call check_close(best, 0.2_dp, 1.0e-12_dp, 'binconf Wilson estimate', failures)
  if (blo <= 0.0_dp .or. bhi >= 1.0_dp .or. blo >= best .or. bhi <= best) then
    print *, 'FAIL binconf Wilson bounds'
    failures = failures + 1
  end if
  call binconf_asymptotic(20.0_dp, 100.0_dp, 0.05_dp, best, blo, bhi)
  call check_close(best, 0.2_dp, 1.0e-12_dp, 'binconf asymptotic estimate', failures)
  call check_close((blo+bhi)/2.0_dp, 0.2_dp, 1.0e-12_dp, 'binconf asymptotic symmetry', failures)

  call invert_tabulated_linear([1.0_dp,2.0_dp,4.0_dp,6.0_dp], &
                               [0.0_dp,1.0_dp,1.0_dp,3.0_dp], &
                               [-1.0_dp,0.5_dp,1.0_dp,2.0_dp,4.0_dp], invout)
  call check_close(invout(1), 1.0_dp, 1.0e-12_dp, 'invertTabulated left clamp', failures)
  call check_close(invout(2), 2.0_dp, 1.0e-12_dp, 'invertTabulated interpolate', failures)
  call check_close(invout(3), 3.0_dp, 1.0e-12_dp, 'invertTabulated tie mean', failures)
  call check_close(invout(4), 4.5_dp, 1.0e-12_dp, 'invertTabulated upper interpolate', failures)
  call check_close(invout(5), 6.0_dp, 1.0e-12_dp, 'invertTabulated right clamp', failures)

  call score_binary_max(reshape([1.0_dp,0.0_dp,1.0_dp, 0.0_dp,1.0_dp,qnan], [3,2]), &
                        [1.0_dp,2.0_dp], scorem)
  call check_close(scorem(1), 1.0_dp, 1.0e-12_dp, 'score.binary max row1', failures)
  call check_close(scorem(2), 2.0_dp, 1.0e-12_dp, 'score.binary max row2', failures)
  call check_close(scorem(3), 1.0_dp, 1.0e-12_dp, 'score.binary max NA removal', failures)
  call score_binary_sum(reshape([1.0_dp,0.0_dp,1.0_dp, 0.0_dp,1.0_dp,qnan], [3,2]), &
                        [1.0_dp,2.0_dp], .true., scores)
  call check_close(scores(1), 1.0_dp, 1.0e-12_dp, 'score.binary sum row1', failures)
  call check_close(scores(2), 2.0_dp, 1.0e-12_dp, 'score.binary sum row2', failures)
  call check_close(scores(3), 1.0_dp, 1.0e-12_dp, 'score.binary sum NA removal', failures)

  call mhgr([1.0_dp,0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp], &
            [1,1,2,2,1,1,2,2], [1,1,1,1,2,2,2,2], 0.95_dp, &
            mh_rr, mh_lo, mh_hi, mh_counts)
  call check_close(mh_rr, 1.0_dp/3.0_dp, 1.0e-12_dp, 'mhgr rr', failures)
  call check_int(mh_counts(1), 4, 'mhgr group1 count', failures)
  call check_int(mh_counts(2), 4, 'mhgr group2 count', failures)
  if (mh_lo >= mh_rr .or. mh_hi <= mh_rr) then
    print *, 'FAIL mhgr confidence interval'
    failures = failures + 1
  end if

  call lrcum([10.0_dp,20.0_dp], [5.0_dp,10.0_dp], [2.0_dp,4.0_dp], [20.0_dp,40.0_dp], 0.95_dp, &
             lp, lpl, lpu, ln, lnl, lnu, lpc, lpcl, lpcu, lnc, lncl, lncu)
  call check_close(lp(1), 25.0_dp/6.0_dp, 1.0e-12_dp, 'lrcum lrpos', failures)
  call check_close(ln(1), 5.0_dp/24.0_dp, 1.0e-12_dp, 'lrcum lrneg', failures)
  call check_close(lpc(2), (25.0_dp/6.0_dp)**2, 1.0e-12_dp, 'lrcum cumulative positive', failures)
  call check_close(lnc(2), (5.0_dp/24.0_dp)**2, 1.0e-12_dp, 'lrcum cumulative negative', failures)
  if (lpl(1) >= lp(1) .or. lpu(1) <= lp(1) .or. lnl(1) >= ln(1) .or. lnu(1) <= ln(1)) then
    print *, 'FAIL lrcum confidence intervals'
    failures = failures + 1
  end if

  call num_denom_setup([2.0_dp,0.0_dp,3.0_dp,qnan], [5.0_dp,4.0_dp,3.0_dp,2.0_dp], &
                       nd_subs, nd_weights, nd_y)
  call check_int(size(nd_subs), 4, 'num.denom.setup size', failures)
  call check_int(nd_subs(1), 1, 'num.denom.setup positive sub1', failures)
  call check_int(nd_subs(2), 3, 'num.denom.setup positive sub2', failures)
  call check_int(nd_subs(3), 1, 'num.denom.setup complement sub1', failures)
  call check_int(nd_subs(4), 2, 'num.denom.setup complement sub2', failures)
  call check_close(nd_weights(1), 2.0_dp, 1.0e-12_dp, 'num.denom.setup numerator weight', failures)
  call check_close(nd_weights(3), 3.0_dp, 1.0e-12_dp, 'num.denom.setup complement weight', failures)
  call check_int(nd_y(1), 1, 'num.denom.setup event', failures)
  call check_int(nd_y(4), 0, 'num.denom.setup nonevent', failures)

  call rmultinom(transpose(reshape([0.2_dp,0.5_dp,0.3_dp, 0.6_dp,0.1_dp,0.3_dp], [3,2])), &
                  transpose(reshape([0.1_dp,0.4_dp,0.85_dp, 0.7_dp,0.95_dp,0.2_dp], [3,2])), &
                  rm_draws, rm_status)
  call check_int(rm_status, 0, 'rMultinom status', failures)
  call check_int(rm_draws(1,1), 1, 'rMultinom draw11', failures)
  call check_int(rm_draws(1,2), 2, 'rMultinom draw12', failures)
  call check_int(rm_draws(1,3), 3, 'rMultinom draw13', failures)
  call check_int(rm_draws(2,1), 2, 'rMultinom draw21', failures)
  call check_int(rm_draws(2,2), 3, 'rMultinom draw22', failures)
  call check_int(rm_draws(2,3), 1, 'rMultinom draw23', failures)

  qx = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, 2.0_dp,1.0_dp,4.0_dp,3.0_dp], [4,2])
  call qrxcenter(qx, qxt, qr, qri, qbar, qr_status)
  call check_int(qr_status, 0, 'qrxcenter status', failures)
  call check_close(qbar(1), 2.5_dp, 1.0e-12_dp, 'qrxcenter mean1', failures)
  call check_close(qbar(2), 2.5_dp, 1.0e-12_dp, 'qrxcenter mean2', failures)
  call check_close(sum(qxt(:,1)), 0.0_dp, 1.0e-12_dp, 'qrxcenter centered1', failures)
  call check_close(sum(qxt(:,2)), 0.0_dp, 1.0e-12_dp, 'qrxcenter centered2', failures)
  call check_close(dot_product(qxt(:,1), qxt(:,1)), 3.0_dp, 1.0e-12_dp, 'qrxcenter scale1', failures)
  call check_close(dot_product(qxt(:,2), qxt(:,2)), 3.0_dp, 1.0e-12_dp, 'qrxcenter scale2', failures)
  call check_close(dot_product(qxt(:,1), qxt(:,2)), 0.0_dp, 1.0e-12_dp, 'qrxcenter orthogonal', failures)
  qrec = matmul(qxt, qri)
  call check_close(maxval(abs(qrec(:,1) - (qx(:,1) - qbar(1)))), 0.0_dp, 1.0e-12_dp, &
                   'qrxcenter inverse transform1', failures)
  call check_close(maxval(abs(qrec(:,2) - (qx(:,2) - qbar(2)))), 0.0_dp, 1.0e-12_dp, &
                   'qrxcenter inverse transform2', failures)
  qrec = matmul(qx - spread(qbar, 1, 4), qr)
  call check_close(maxval(abs(qrec - qxt)), 0.0_dp, 1.0e-12_dp, 'qrxcenter forward transform', failures)

  sv_a = reshape([1.0_dp,1.0_dp,1.0_dp, 1.0_dp,2.0_dp,3.0_dp], [3,2])
  sv_b = [1.0_dp,2.0_dp,2.0_dp]
  call solvet(sv_a, sv_b, sv_x, status=sv_status)
  call check_int(sv_status, 0, 'solvet vector status', failures)
  call check_close(sv_x(1), 2.0_dp/3.0_dp, 1.0e-12_dp, 'solvet intercept', failures)
  call check_close(sv_x(2), 0.5_dp, 1.0e-12_dp, 'solvet slope', failures)
  call solvet_inverse(reshape([2.0_dp,1.0_dp,1.0_dp,3.0_dp],[2,2]), sv_inv, status=sv_status)
  call check_int(sv_status, 0, 'solvet inverse status', failures)
  call check_close(sv_inv(1,1), 0.6_dp, 1.0e-12_dp, 'solvet inverse 11', failures)
  call check_close(sv_inv(1,2), -0.2_dp, 1.0e-12_dp, 'solvet inverse 12', failures)
  call check_close(sv_inv(2,2), 0.4_dp, 1.0e-12_dp, 'solvet inverse 22', failures)

  call smean_cl_normal([1.0_dp,2.0_dp,3.0_dp,4.0_dp], 2.0_dp, smn)
  call check_close(smn(1), 2.5_dp, 1.0e-12_dp, 'smean.cl.normal mean', failures)
  call check_close(smn(2), 2.5_dp - 2.0_dp*sqrt(5.0_dp/12.0_dp), 1.0e-12_dp, &
                   'smean.cl.normal lower', failures)
  call check_close(smn(3), 2.5_dp + 2.0_dp*sqrt(5.0_dp/12.0_dp), 1.0e-12_dp, &
                   'smean.cl.normal upper', failures)

  call which_close_k([0.0_dp,10.0_dp,20.0_dp], [9.0_dp,18.0_dp], 2, &
                     [0.5_dp,0.5_dp,0.5_dp], [0.0_dp,0.99_dp], wk_idx, wk_status)
  call check_int(wk_status, 0, 'whichClosek status', failures)
  call check_int(wk_idx(1), 2, 'whichClosek first', failures)
  call check_int(wk_idx(2), 2, 'whichClosek second', failures)

  call pc1(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, 2.0_dp,3.0_dp,5.0_dp,7.0_dp], [4,2]), &
           pc_scores, pc_coef, pc_frac, pc_status, 10.0_dp)
  call check_int(pc_status, 0, 'pc1 status', failures)
  call check_close(pc_frac, 0.9948891332786447_dp, 1.0e-10_dp, 'pc1 fraction variance', failures)
  call check_close(minval(pc_scores), 0.0_dp, 1.0e-10_dp, 'pc1 rescaled minimum', failures)
  call check_close(maxval(pc_scores), 10.0_dp, 1.0e-10_dp, 'pc1 rescaled maximum', failures)


  call lm_fit_qr_bare(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[4,1]), &
                      [3.0_dp,5.0_dp,7.0_dp,9.0_dp], .true., &
                      lm_coef, lm_res, lm_fit, lm_r2, lm_status)
  call check_int(lm_status, 0, 'lm.fit.qr.bare status', failures)
  call check_close(lm_coef(1), 1.0_dp, 1.0e-12_dp, 'lm.fit.qr.bare intercept', failures)
  call check_close(lm_coef(2), 2.0_dp, 1.0e-12_dp, 'lm.fit.qr.bare slope', failures)
  call check_close(maxval(abs(lm_res)), 0.0_dp, 1.0e-12_dp, 'lm.fit.qr.bare residuals', failures)
  call check_close(lm_r2, 1.0_dp, 1.0e-12_dp, 'lm.fit.qr.bare R2', failures)

  call abs_error_pred([1.0_dp,2.0_dp,4.0_dp], [1.0_dp,3.0_dp,3.0_dp], &
                      aep_diff, aep_rat, aep_status)
  call check_int(aep_status, 0, 'abs.error.pred status', failures)
  call check_close(aep_diff(1,1), 2.0_dp/3.0_dp, 1.0e-12_dp, &
                   'abs.error.pred prediction spread mean', failures)
  call check_close(aep_diff(2,1), 2.0_dp/3.0_dp, 1.0e-12_dp, &
                   'abs.error.pred prediction error mean', failures)
  call check_close(aep_diff(3,1), 1.0_dp, 1.0e-12_dp, &
                   'abs.error.pred outcome spread mean', failures)
  call check_close(aep_diff(1,2), 0.0_dp, 1.0e-12_dp, &
                   'abs.error.pred prediction spread median', failures)
  call check_close(aep_diff(2,2), 1.0_dp, 1.0e-12_dp, &
                   'abs.error.pred prediction error median', failures)
  call check_close(aep_rat(2,2), 1.0_dp, 1.0e-12_dp, &
                   'abs.error.pred median ratio', failures)

  call rcorrcens_summary([1.0_dp,3.0_dp,2.0_dp,4.0_dp], &
                         [1.0_dp,2.0_dp,3.0_dp,4.0_dp], &
                         [.true.,.true.,.true.,.true.], .false., &
                         rcs_c, rcs_d, rcs_ad, rcs_sd, rcs_z, rcs_p, rcs_status)
  call check_int(rcs_status, 0, 'rcorrcens status', failures)
  call check_close(rcs_c, 5.0_dp/6.0_dp, 1.0e-12_dp, 'rcorrcens C', failures)
  call check_close(rcs_d, 2.0_dp/3.0_dp, 1.0e-12_dp, 'rcorrcens Dxy', failures)
  call check_close(rcs_sd, 1.0_dp/3.0_dp, 1.0e-12_dp, 'rcorrcens SD', failures)
  call check_close(rcs_z, 2.0_dp, 1.0e-12_dp, 'rcorrcens Z', failures)
  call check_close(rcs_p, 0.0455002638963584_dp, 1.0e-12_dp, 'rcorrcens P', failures)

  call logrank([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp], &
               [.true.,.true.,.true.,.true.,.true.,.true.], [1,1,1,2,2,2], &
               lr_chisq, lr_hr, lr_status)
  call check_int(lr_status, 0, 'logrank status', failures)
  call check_close(lr_chisq, 5.051660516605167_dp, 1.0e-12_dp, 'logrank chisq', failures)
  call check_close(lr_hr, 0.2371134020618557_dp, 1.0e-12_dp, 'logrank hazard ratio', failures)

  call spearman2_numeric([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp], &
                         [1.0_dp,3.0_dp,2.0_dp,5.0_dp,4.0_dp], 1, &
                         sp_rho2, sp_f, sp_df1, sp_df2, sp_p, sp_adj, sp_n, sp_status)
  call check_int(sp_status, 0, 'spearman2 status', failures)
  call check_int(sp_n, 5, 'spearman2 n', failures)
  call check_close(sp_rho2, 0.64_dp, 1.0e-12_dp, 'spearman2 rho2', failures)
  call check_close(sp_f, 16.0_dp/3.0_dp, 1.0e-12_dp, 'spearman2 F', failures)
  call check_close(sp_df1, 1.0_dp, 1.0e-12_dp, 'spearman2 df1', failures)
  call check_close(sp_df2, 3.0_dp, 1.0e-12_dp, 'spearman2 df2', failures)
  call check_close(sp_p, 0.10408803866182785_dp, 2.0e-12_dp, 'spearman2 P', failures)
  call check_close(sp_adj, 0.52_dp, 1.0e-12_dp, 'spearman2 adjusted', failures)

  call equal_bins([10,7], reshape([2,4,3,2,0,1],[2,3]), [2,3], eb, eb_status)
  call check_int(eb_status, 0, 'equalBins status', failures)
  call check_int(size(eb), 5, 'equalBins size', failures)
  call check_int(eb(1), 5, 'equalBins group1 first', failures)
  call check_int(eb(2), 4, 'equalBins group1 second', failures)
  call check_int(eb(3), 4, 'equalBins group2 first', failures)
  call check_int(eb(4), 2, 'equalBins group2 second', failures)
  call check_int(eb(5), 1, 'equalBins group2 third', failures)

  call jitter2_numeric([1.0_dp,1.0_dp,2.0_dp,2.0_dp,4.0_dp], 1.0_dp/3.0_dp, &
                       0.0_dp, 0.0_dp, [0.9_dp,0.1_dp,0.2_dp,0.8_dp], jit_out, jit_status)
  call check_int(jit_status, 0, 'jitter2 status', failures)
  call check_close(jit_out(1), 7.0_dp/6.0_dp, 1.0e-12_dp, 'jitter2 tie1a', failures)
  call check_close(jit_out(2), 5.0_dp/6.0_dp, 1.0e-12_dp, 'jitter2 tie1b', failures)
  call check_close(jit_out(3), 11.0_dp/6.0_dp, 1.0e-12_dp, 'jitter2 tie2a', failures)
  call check_close(jit_out(4), 13.0_dp/6.0_dp, 1.0e-12_dp, 'jitter2 tie2b', failures)
  call check_close(jit_out(5), 4.0_dp, 1.0e-12_dp, 'jitter2 singleton', failures)

  call hdquantile([1.0_dp,2.0_dp,3.0_dp], [0.0_dp,0.5_dp,1.0_dp], hdq, hdq_status)
  call check_int(hdq_status, 0, 'hdquantile status', failures)
  call check_close(hdq(1), 1.0_dp, 1.0e-12_dp, 'hdquantile minimum', failures)
  call check_close(hdq(2), 2.0_dp, 1.0e-12_dp, 'hdquantile median', failures)
  call check_close(hdq(3), 3.0_dp, 1.0e-12_dp, 'hdquantile maximum', failures)

  call cat_test_chisq(reshape([10.0_dp,30.0_dp,20.0_dp,40.0_dp],[2,2]), &
                      ct_stat, ct_df, ct_p, ct_status)
  call check_int(ct_status, 0, 'catTestchisq status', failures)
  call check_close(ct_stat, 0.7936507936507936_dp, 1.0e-12_dp, 'catTestchisq statistic', failures)
  call check_close(ct_df, 1.0_dp, 1.0e-12_dp, 'catTestchisq df', failures)
  call check_close(ct_p, 0.3729984836134869_dp, 2.0e-12_dp, 'catTestchisq P', failures)

  call con_test_kw([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp], &
                   [1.0_dp,3.0_dp,2.0_dp,5.0_dp,4.0_dp], kw_stat, kw_df1, kw_df2, kw_p, kw_status)
  call check_int(kw_status, 0, 'conTestkw status', failures)
  call check_close(kw_stat, 16.0_dp/3.0_dp, 1.0e-12_dp, 'conTestkw F', failures)
  call check_close(kw_df1, 1.0_dp, 1.0e-12_dp, 'conTestkw df1', failures)
  call check_close(kw_df2, 3.0_dp, 1.0e-12_dp, 'conTestkw df2', failures)
  call check_close(kw_p, 0.10408803866182785_dp, 2.0e-12_dp, 'conTestkw P', failures)

  call t_test_cluster([1.0_dp,2.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,5.0_dp,6.0_dp], &
                      [1,1,2,2,3,3,4,4], [1,1,1,1,2,2,2,2], 0.95_dp, &
                      tt_n, tt_nc, tt_mean, tt_rho, tt_corr, tt_var, tt_vu, tt_de, &
                      tt_eff, tt_se, tt_ci, tt_z, tt_p, tt_status)
  call check_int(tt_status, 0, 't.test.cluster status', failures)
  call check_int(tt_n(1), 4, 't.test.cluster n1', failures)
  call check_int(tt_nc(2), 2, 't.test.cluster clusters2', failures)
  call check_close(tt_mean(1), 2.0_dp, 1.0e-12_dp, 't.test.cluster mean1', failures)
  call check_close(tt_mean(2), 5.0_dp, 1.0e-12_dp, 't.test.cluster mean2', failures)
  call check_close(tt_rho, 3.0_dp/7.0_dp, 1.0e-12_dp, 't.test.cluster rho', failures)
  call check_close(tt_corr(1), 10.0_dp/7.0_dp, 1.0e-12_dp, 't.test.cluster correction', failures)
  call check_close(tt_var, 5.0_dp/14.0_dp, 1.0e-12_dp, 't.test.cluster variance', failures)
  call check_close(tt_de, 10.0_dp/7.0_dp, 1.0e-12_dp, 't.test.cluster design effect', failures)
  call check_close(tt_eff, 3.0_dp, 1.0e-12_dp, 't.test.cluster effect', failures)
  call check_close(tt_se, sqrt(5.0_dp/14.0_dp), 1.0e-12_dp, 't.test.cluster SE', failures)
  call check_close(tt_z, 5.019960159204453_dp, 2.0e-12_dp, 't.test.cluster Z', failures)
  call check_close(tt_p, 5.168219753811347e-7_dp, 2.0e-18_dp, 't.test.cluster P', failures)


  call km_quick([1.0_dp,2.0_dp,2.0_dp,3.0_dp], [.true.,.true.,.false.,.true.], &
                km_t, km_s, km_status, [0.5_dp,1.0_dp,2.0_dp,2.5_dp,3.0_dp], km_q)
  call check_int(km_status, 0, 'km.quick status', failures)
  call check_int(size(km_t), 3, 'km.quick event times', failures)
  call check_close(km_s(1), 0.75_dp, 1.0e-12_dp, 'km.quick S1', failures)
  call check_close(km_s(2), 0.5_dp, 1.0e-12_dp, 'km.quick S2', failures)
  call check_close(km_s(3), 0.0_dp, 1.0e-12_dp, 'km.quick S3', failures)
  call check_close(km_q(1), 1.0_dp, 1.0e-12_dp, 'km.quick before first', failures)
  call check_close(km_q(2), 0.75_dp, 1.0e-12_dp, 'km.quick exact event', failures)
  call check_close(km_q(3), 0.5_dp, 1.0e-12_dp, 'km.quick tied time', failures)

  call bystats_mean(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
                             10.0_dp,20.0_dp,30.0_dp,40.0_dp],[4,2]), &
                    [2,1,2,1], bs_labels, bs_n, bs_miss, bs_mean, bs_status)
  call check_int(bs_status, 0, 'bystats status', failures)
  call check_int(bs_labels(1), 1, 'bystats label1', failures)
  call check_int(bs_n(1), 2, 'bystats n group1', failures)
  call check_int(bs_n(3), 4, 'bystats n all', failures)
  call check_close(bs_mean(1,1), 3.0_dp, 1.0e-12_dp, 'bystats mean group1', failures)
  call check_close(bs_mean(2,2), 20.0_dp, 1.0e-12_dp, 'bystats mean group2 col2', failures)
  call check_close(bs_mean(3,1), 2.5_dp, 1.0e-12_dp, 'bystats overall', failures)

  call bystats2_mean(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[4,1]), &
                     [1,1,2,2], [1,2,1,2], b2_vlab, b2_hlab, b2_n, b2_miss, b2_mean, b2_status)
  call check_int(b2_status, 0, 'bystats2 status', failures)
  call check_int(b2_n(1,1), 1, 'bystats2 cell count', failures)
  call check_int(b2_n(3,3), 4, 'bystats2 all count', failures)
  call check_close(b2_mean(1,3,1), 1.5_dp, 1.0e-12_dp, 'bystats2 row margin', failures)
  call check_close(b2_mean(3,2,1), 3.0_dp, 1.0e-12_dp, 'bystats2 col margin', failures)
  call check_close(b2_mean(3,3,1), 2.5_dp, 1.0e-12_dp, 'bystats2 overall', failures)

  if (.not. all_digits('0123456789')) then
    print *, 'FAIL all.digits true'
    failures = failures + 1
  end if
  if (all_digits('12a3')) then
    print *, 'FAIL all.digits false'
    failures = failures + 1
  end if
  if (.not. all_is_numeric([character(len=4) :: ' 1.2', '3', 'NA'])) then
    print *, 'FAIL all.is.numeric true'
    failures = failures + 1
  end if
  if (all_is_numeric([character(len=3) :: '1', 'bad'])) then
    print *, 'FAIL all.is.numeric false'
    failures = failures + 1
  end if


  call find_matches(reshape([1.0_dp,1.0_dp],[1,2]), &
                    reshape([1.1_dp,1.4_dp,2.0_dp,1.2_dp,1.4_dp,2.0_dp],[3,2]), &
                    [0.5_dp,0.5_dp], [0.5_dp,0.5_dp], 2, fm_match, fm_dist, fm_n, fm_status)
  call check_int(fm_status, 0, 'find.matches status', failures)
  call check_int(fm_n(1), 2, 'find.matches count', failures)
  call check_int(fm_match(1,1), 1, 'find.matches first', failures)
  call check_int(fm_match(1,2), 2, 'find.matches second', failures)
  call check_close(fm_dist(1,1), 0.2_dp, 1.0e-12_dp, 'find.matches distance1', failures)
  call check_close(fm_dist(1,2), 1.28_dp, 1.0e-12_dp, 'find.matches distance2', failures)

  call seq_freq(reshape([1,1,0,0, 0,1,1,0, 0,0,1,0],[4,3]), &
                sf_assign, sf_order, sf_cond, sf_status)
  call check_int(sf_status, 0, 'seqFreq status', failures)
  call check_int(sf_order(1), 1, 'seqFreq order1', failures)
  call check_int(sf_order(2), 2, 'seqFreq order2', failures)
  call check_int(sf_assign(1), 1, 'seqFreq assignment1', failures)
  call check_int(sf_assign(3), 2, 'seqFreq assignment3', failures)
  call check_int(sf_assign(4), 0, 'seqFreq assignment none', failures)
  call check_int(sf_cond(1), 1, 'seqFreq zero positives', failures)
  call check_int(sf_cond(2), 1, 'seqFreq one positive', failures)
  call check_int(sf_cond(3), 2, 'seqFreq two positives', failures)

  call pair_up_diff([1.0_dp,4.0_dp,7.0_dp,2.0_dp,5.0_dp,3.0_dp,6.0_dp], &
                    [1,1,1,1,2,2,2], [1,1,2,2,3,3,4], [1,2,1,2,1,2,1], 1, &
                    pu_major, pu_minor, pu_diff, pu_mid, pu_sd, pu_lo, pu_hi, pu_lm, pu_um, pu_status, &
                    lower=[0.9_dp,3.9_dp,6.9_dp,1.9_dp,4.9_dp,2.9_dp,5.9_dp], &
                    upper=[1.1_dp,4.1_dp,7.1_dp,2.1_dp,5.1_dp,3.1_dp,6.1_dp])
  call check_int(pu_status, 0, 'pairUpDiff status', failures)
  call check_int(size(pu_diff), 4, 'pairUpDiff pair count', failures)
  call check_int(pu_major(1), 1, 'pairUpDiff major1', failures)
  call check_int(pu_minor(1), 1, 'pairUpDiff minor1', failures)
  call check_close(pu_diff(1), 3.0_dp, 1.0e-12_dp, 'pairUpDiff diff1', failures)
  call check_close(pu_diff(2), -5.0_dp, 1.0e-12_dp, 'pairUpDiff diff2', failures)
  call check_close(pu_mid(1), 2.5_dp, 1.0e-12_dp, 'pairUpDiff midpoint', failures)
  call check_close(pu_hi(1) - pu_diff(1), sqrt(0.02_dp), 2.0e-12_dp, &
                   'pairUpDiff CI halfwidth', failures)


  call smean_cl_boot([1.0_dp,2.0_dp,3.0_dp], reshape([1,1,1, 1,2,3, 3,3,3, 1,3,2],[3,4]), &
                     0.5_dp, boot_mean, boot_lo, boot_hi, boot_reps, boot_status)
  call check_int(boot_status, 0, 'smean.cl.boot status', failures)
  call check_close(boot_mean, 2.0_dp, 1.0e-12_dp, 'smean.cl.boot mean', failures)
  call check_close(boot_reps(1), 1.0_dp, 1.0e-12_dp, 'smean.cl.boot rep1', failures)
  call check_close(boot_reps(3), 3.0_dp, 1.0e-12_dp, 'smean.cl.boot rep3', failures)
  call check_close(boot_lo, 1.75_dp, 1.0e-12_dp, 'smean.cl.boot lower', failures)
  call check_close(boot_hi, 2.25_dp, 1.0e-12_dp, 'smean.cl.boot upper', failures)

  call smearing_est_tabulated([1.0_dp,2.0_dp], [0.0_dp,1.0_dp,2.0_dp,3.0_dp], &
                              [0.0_dp,10.0_dp,20.0_dp,30.0_dp], [-1.0_dp,1.0_dp], &
                              1, 0.5_dp, smear_out, smear_status)
  call check_int(smear_status, 0, 'smearingEst mean status', failures)
  call check_close(smear_out(1), 10.0_dp, 1.0e-12_dp, 'smearingEst mean1', failures)
  call check_close(smear_out(2), 20.0_dp, 1.0e-12_dp, 'smearingEst mean2', failures)
  call smearing_est_tabulated([1.0_dp,2.0_dp], [0.0_dp,1.0_dp,2.0_dp,3.0_dp], &
                              [0.0_dp,10.0_dp,20.0_dp,30.0_dp], [-1.0_dp,1.0_dp], &
                              2, 0.5_dp, smear_out, smear_status)
  call check_close(smear_out(1), 10.0_dp, 1.0e-12_dp, 'smearingEst quantile1', failures)


  call bpower_sim_counts([0,10,5,9], [10,0,5,1], 10, 10, 0.05_dp, &
                         bps_power, bps_lo, bps_hi, bps_status)
  call check_int(bps_status, 0, 'bpower.sim status', failures)
  call check_close(bps_power, 0.75_dp, 1.0e-12_dp, 'bpower.sim power', failures)
  call check_close(bps_lo, 0.75_dp - 1.96_dp*sqrt(0.75_dp*0.25_dp/4.0_dp), &
                   1.0e-12_dp, 'bpower.sim lower', failures)
  call check_close(bps_hi, 0.75_dp + 1.96_dp*sqrt(0.75_dp*0.25_dp/4.0_dp), &
                   1.0e-12_dp, 'bpower.sim upper', failures)

  call spower_simulated(reshape([1.0_dp,2.0_dp,3.0_dp],[3,1]), &
                        reshape([4.0_dp,5.0_dp,6.0_dp],[3,1]), &
                        reshape([6.0_dp,6.0_dp,6.0_dp,6.0_dp,6.0_dp,6.0_dp],[6,1]), &
                        0.05_dp, sps_power, sps_maxf, sps_maxc, sps_status)
  call check_int(sps_status, 0, 'spower status', failures)
  call check_close(sps_power, 1.0_dp, 1.0e-12_dp, 'spower power', failures)
  call check_close(sps_maxf, 6.0_dp, 1.0e-12_dp, 'spower max failure', failures)
  call check_close(sps_maxc, 6.0_dp, 1.0e-12_dp, 'spower max censor', failures)

  call inverse_function_all([-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp], &
                            [4.0_dp,1.0_dp,0.0_dp,1.0_dp,4.0_dp], [1.0_dp,2.0_dp], &
                            invf_roots, invf_nroots, invf_status)
  call check_int(invf_status, 0, 'inverseFunction status', failures)
  call check_int(invf_nroots(1), 2, 'inverseFunction roots1 count', failures)
  call check_int(invf_nroots(2), 2, 'inverseFunction roots2 count', failures)
  call check_close(invf_roots(1,1), -1.0_dp, 1.0e-12_dp, 'inverseFunction root1a', failures)
  call check_close(invf_roots(1,2), 1.0_dp, 1.0e-12_dp, 'inverseFunction root1b', failures)
  call check_close(invf_roots(2,1), -4.0_dp/3.0_dp, 1.0e-12_dp, 'inverseFunction root2a', failures)
  call check_close(invf_roots(2,2), 4.0_dp/3.0_dp, 1.0e-12_dp, 'inverseFunction root2b', failures)

  call bezier_curve([0.0_dp,1.0_dp,2.0_dp], [0.0_dp,2.0_dp,0.0_dp], 3, bz_x, bz_y, bz_status)
  call check_int(bz_status, 0, 'bezier status', failures)
  call check_close(bz_x(2), 1.0_dp, 1.0e-12_dp, 'bezier midpoint x', failures)
  call check_close(bz_y(2), 1.0_dp, 1.0e-12_dp, 'bezier midpoint y', failures)

  call cut2_explicit([-1.0_dp,0.0_dp,0.5_dp,1.0_dp,2.0_dp,3.0_dp], &
                     [0.0_dp,1.0_dp,2.0_dp], .true., cut_group, cut_eff, cut_status)
  call check_int(cut_status, 0, 'cut2 status', failures)
  call check_int(size(cut_eff), 5, 'cut2 effective cuts', failures)
  call check_int(cut_group(1), 1, 'cut2 extended low', failures)
  call check_int(cut_group(3), 2, 'cut2 interior', failures)
  call check_int(cut_group(6), 4, 'cut2 extended high', failures)

  call sim_po_cuts(8, [0.25_dp,0.5_dp,0.25_dp], 1.0_dp, &
                   reshape([0.1_dp,0.3_dp,0.6_dp,0.9_dp],[4,1]), &
                   reshape([0.1_dp,0.3_dp,0.6_dp,0.9_dp],[4,1]), &
                   poc_or, poc_status)
  call check_int(poc_status, 0, 'simPOcuts status', failures)
  call check_close(poc_or(1,1), 1.0_dp, 1.0e-12_dp, 'simPOcuts cut1', failures)
  call check_close(poc_or(1,2), 1.0_dp, 1.0e-12_dp, 'simPOcuts cut2', failures)

  call rcspline_function([0.0_dp,1.0_dp,2.0_dp,3.0_dp], [0.0_dp,1.0_dp,2.0_dp], &
                         [1.0_dp,2.0_dp,3.0_dp], 0, .false., rcf_y, rcf_status)
  call check_int(rcf_status, 0, 'rcsplineFunction status', failures)
  call check_close(rcf_y(1), 1.0_dp, 1.0e-12_dp, 'rcsplineFunction y0', failures)
  call check_close(rcf_y(2), 6.0_dp, 1.0e-12_dp, 'rcsplineFunction y1', failures)
  call check_close(rcf_y(4), 43.0_dp, 1.0e-12_dp, 'rcsplineFunction y3', failures)

  call rcspline_restate_coefficients([0.0_dp,1.0_dp,2.0_dp], [1.0_dp,2.0_dp,3.0_dp], &
                                     0, .false., rcr_coef, rcr_status)
  call check_int(rcr_status, 0, 'rcspline.restate status', failures)
  call check_int(size(rcr_coef), 5, 'rcspline.restate length', failures)
  call check_close(rcr_coef(1), 1.0_dp, 1.0e-12_dp, 'rcspline.restate intercept', failures)
  call check_close(rcr_coef(2), 2.0_dp, 1.0e-12_dp, 'rcspline.restate linear', failures)
  call check_close(rcr_coef(3), 3.0_dp, 1.0e-12_dp, 'rcspline.restate first cubic', failures)
  call check_close(rcr_coef(4), -6.0_dp, 1.0e-12_dp, 'rcspline.restate restricted1', failures)
  call check_close(rcr_coef(5), 3.0_dp, 1.0e-12_dp, 'rcspline.restate restricted2', failures)


  call spearman_test([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp], &
                     [1.0_dp,2.0_dp,4.0_dp,3.0_dp,5.0_dp], 1, &
                     st_r2, st_f, st_df1, st_df2, st_p, st_n, st_status)
  call check_int(st_status, 0, 'spearman.test status', failures)
  call check_int(st_n, 5, 'spearman.test n', failures)
  call check_close(st_r2, 0.81_dp, 1.0e-12_dp, 'spearman.test R2', failures)
  call check_close(st_df1, 1.0_dp, 1.0e-12_dp, 'spearman.test df1', failures)
  call check_close(st_df2, 3.0_dp, 1.0e-12_dp, 'spearman.test df2', failures)

  call chi_square_codes([1,1,2,2], [1,2,1,2], cs_stat, cs_df, cs_excess, cs_p, cs_status)
  call check_int(cs_status, 0, 'chiSquare status', failures)
  call check_close(cs_stat, 0.0_dp, 1.0e-12_dp, 'chiSquare statistic', failures)
  call check_close(cs_df, 1.0_dp, 1.0e-12_dp, 'chiSquare df', failures)
  call check_close(cs_excess, -1.0_dp, 1.0e-12_dp, 'chiSquare excess', failures)

  call combine_levels_codes([1,2,2,2,3,3,3,3,4,5], 3, .false., cl_out, cl_map, cl_status)
  call check_int(cl_status, 0, 'combine.levels unordered status', failures)
  call check_int(maxval(cl_out), 3, 'combine.levels unordered groups', failures)
  call check_int(cl_map(1), cl_map(4), 'combine.levels unordered rare A/D', failures)
  call check_int(cl_map(1), cl_map(5), 'combine.levels unordered rare A/E', failures)
  call combine_levels_codes([1,2,2,2,3,3,3,3,4,5], 3, .true., cl_out, cl_map, cl_status)
  call check_int(cl_status, 0, 'combine.levels ordered status', failures)
  call check_int(maxval(cl_out), 2, 'combine.levels ordered groups', failures)
  call check_int(cl_map(1), cl_map(2), 'combine.levels ordered first pair', failures)
  call check_int(cl_map(3), cl_map(5), 'combine.levels ordered tail', failures)

  call princmp_regular(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
                               2.0_dp,4.0_dp,6.0_dp,8.0_dp], [4,2]), &
                       .true., 2, pm_scores, pm_sload, pm_oload, pm_var, pm_scale, pm_status)
  call check_int(pm_status, 0, 'princmp status', failures)
  call check_close(pm_var(1), 2.0_dp, 1.0e-10_dp, 'princmp variance1', failures)
  call check_close(pm_var(2), 0.0_dp, 1.0e-10_dp, 'princmp variance2', failures)
  call check_close(abs(pm_sload(1,1)), 1.0_dp/sqrt(2.0_dp), 1.0e-10_dp, &
                   'princmp loading1', failures)


  call bootkm_resampled([1.0_dp,2.0_dp,3.0_dp], [.true.,.true.,.false.], &
                        reshape([1,2,3,1,1,3],[3,2]), .true., 0.5_dp, 2.0_dp, &
                        bkm_est, bkm_status)
  call check_int(bkm_status, 0, 'bootkm status', failures)
  call check_close(bkm_est(1), 1.0_dp/3.0_dp, 1.0e-12_dp, 'bootkm survival1', failures)
  call check_close(bkm_est(2), 1.0_dp/3.0_dp, 1.0e-12_dp, 'bootkm survival2', failures)

  allocate(mc_score(2,4))
  mc_score = reshape([0.4_dp,0.9_dp,0.2_dp,0.8_dp,0.3_dp,0.7_dp,0.1_dp,0.6_dp],[2,4])
  call match_cases_core([0.0_dp,10.0_dp], [0.1_dp,0.2_dp,9.0_dp,10.1_dp], &
                        0.5_dp, 2, .true., mc_score, mc_matches, mc_nmatch, mc_status)
  call check_int(mc_status, 0, 'matchCases status', failures)
  call check_int(mc_nmatch(1), 2, 'matchCases case1 count', failures)
  call check_int(mc_matches(1,1), 1, 'matchCases closest1', failures)
  call check_int(mc_matches(1,2), 2, 'matchCases closest2', failures)
  call check_int(mc_nmatch(2), 1, 'matchCases case2 count', failures)
  call check_int(mc_matches(2,1), 4, 'matchCases case2 closest', failures)

  call largest_empty_rexhaustive([0.5_dp], [0.5_dp], [0.0_dp,1.0_dp], [0.0_dp,1.0_dp], &
                                 0.0_dp, 0.0_dp, le_center, le_rect, le_area, le_status)
  call check_int(le_status, 0, 'largest.empty status', failures)
  call check_close(le_area, 0.5_dp, 1.0e-12_dp, 'largest.empty area', failures)
  call check_close(le_center(1), 0.25_dp, 1.0e-12_dp, 'largest.empty center x', failures)
  call check_close(le_center(2), 0.5_dp, 1.0e-12_dp, 'largest.empty center y', failures)

  call gbayes_update(0.0_dp, 4.0_dp, 2.0_dp, 1.0_dp, gb_mean_post, gb_var_post, gb_status, &
                     3.0_dp, gb_mean_pred, gb_var_pred)
  call check_int(gb_status, 0, 'gbayes status', failures)
  call check_close(gb_mean_post, 1.6_dp, 1.0e-12_dp, 'gbayes posterior mean', failures)
  call check_close(gb_var_post, 0.8_dp, 1.0e-12_dp, 'gbayes posterior variance', failures)
  call check_close(gb_mean_pred, 1.6_dp, 1.0e-12_dp, 'gbayes predictive mean', failures)
  call check_close(gb_var_pred, 3.8_dp, 1.0e-12_dp, 'gbayes predictive variance', failures)

  gb_pred_d = gbayes_mix_pred_density(0.5_dp, 1.0_dp, 0.25_dp, 0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp)
  gb_pred_c = gbayes_mix_pred_cdf(0.5_dp, 1.0_dp, 0.25_dp, 0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp)
  call check_close(gb_pred_d, 0.19182627573918287_dp, 1.0e-12_dp, 'gbayesMixPredNoData density', failures)
  call check_close(gb_pred_c, 0.22513717698504074_dp, 1.0e-12_dp, 'gbayesMixPredNoData cdf', failures)

  gb_post_d = gbayes_mix_post_density(0.5_dp, 1.0_dp, 1.0_dp, 0.25_dp, 0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp)
  gb_post_c = gbayes_mix_post_cdf(0.5_dp, 1.0_dp, 1.0_dp, 0.25_dp, 0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp)
  gb_post_m = gbayes_mix_post_mean(1.0_dp, 1.0_dp, 0.25_dp, 0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp)
  call check_close(gb_post_d, 0.40196128031168304_dp, 1.0e-12_dp, 'gbayesMixPost density', failures)
  call check_close(gb_post_c, 0.3052631029209445_dp, 1.0e-12_dp, 'gbayesMixPost cdf', failures)
  call check_close(gb_post_m, 0.9815521054074611_dp, 1.0e-12_dp, 'gbayesMixPost mean', failures)

  gb_power = gbayes1_power_np(0.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 0.0_dp, 0.05_dp)
  call check_close(gb_power, 0.01154906534430422_dp, 5.0e-9_dp, 'gbayes1PowerNP', failures)

  call gbayes_mix_power_np(0.5_dp, 1.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, &
                           [-10.0_dp,10.0_dp], 0.05_dp, gb_mix_critical, gb_mix_power, gb_status)
  call check_int(gb_status, 0, 'gbayesMixPowerNP status', failures)
  call check_close(gb_mix_critical, 1.9688692786833477_dp, 2.0e-10_dp, &
                   'gbayesMixPowerNP critical', failures)
  call check_close(gb_mix_power, 0.07093412437350288_dp, 2.0e-10_dp, &
                   'gbayesMixPowerNP power', failures)



  call gbayes2_tabulated([-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp], &
                         [0.25_dp,0.25_dp,0.25_dp,0.25_dp,0.25_dp], &
                         1.0_dp, 0.0_dp, 0.05_dp, gb2_prob, gb2_status)
  call check_int(gb2_status, 0, 'gbayes2 status', failures)
  call check_close(gb2_prob, 0.1097551418578421_dp, 5.0e-10_dp, 'gbayes2 tabulated', failures)

  call props_po_counts(reshape([20,10,30,20,50,70],[2,3]), [1.0_dp,2.0_dp], 1, &
                       ppo_obs, ppo_pred, ppo_status)
  call check_int(ppo_status, 0, 'propsPO status', failures)
  call check_close(ppo_obs(1,1), 0.2_dp, 1.0e-12_dp, 'propsPO observed', failures)
  call check_close(ppo_pred(1,2), 0.3_dp, 1.0e-12_dp, 'propsPO reference', failures)
  call check_close(sum(ppo_pred(2,:)), 1.0_dp, 1.0e-12_dp, 'propsPO predicted sum', failures)

  call props_trans_counts(reshape([1,1,2,2, 1,2,2,3, 2,2,3,3],[4,3]), 3, &
                          ptr_count, ptr_prop, ptr_status)
  call check_int(ptr_status, 0, 'propsTrans status', failures)
  call check_int(ptr_count(1,1,1), 1, 'propsTrans count11', failures)
  call check_int(ptr_count(1,2,1), 1, 'propsTrans count12', failures)
  call check_close(ptr_prop(1,1,1), 0.5_dp, 1.0e-12_dp, 'propsTrans prop11', failures)
  call check_close(ptr_prop(2,3,2), 0.5_dp, 1.0e-12_dp, 'propsTrans prop23', failures)

  call ynbind_numeric(reshape([0.0_dp,1.0_dp,0.0_dp,1.0_dp, &
                               1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
                               0.0_dp,0.0_dp,0.0_dp,1.0_dp],[4,3]), &
                      .true., yn_y, yn_order, yn_prop, yn_status)
  call check_int(yn_status, 0, 'ynbind status', failures)
  call check_int(yn_order(1), 3, 'ynbind order1', failures)
  call check_int(yn_order(2), 1, 'ynbind order2', failures)
  call check_int(yn_order(3), 2, 'ynbind order3', failures)
  call check_close(yn_prop(1), 0.25_dp, 1.0e-12_dp, 'ynbind prop1', failures)
  call check_close(yn_prop(3), 1.0_dp, 1.0e-12_dp, 'ynbind prop3', failures)

  rf_missing = .false.
  rf_missing(1,1) = .true.
  rf_missing(2,1) = .true.
  rf_missing(3,2) = .true.
  rf_missing(4,1) = .true.
  rf_missing(4,3) = .true.
  call reformm_missing_order(rf_missing, rf_order, rf_varmiss, rf_obsmiss, rf_nimp, rf_status)
  call check_int(rf_status, 0, 'reformM status', failures)
  call check_int(rf_varmiss(1), 3, 'reformM var missing1', failures)
  call check_int(rf_order(1), 1, 'reformM order1', failures)
  call check_int(rf_order(2), 2, 'reformM stable tie order', failures)
  call check_int(rf_nimp, 80, 'reformM recommended imputations', failures)

  ma_x = reshape([1.0_dp,3.0_dp,5.0_dp,7.0_dp, 2.0_dp,4.0_dp,6.0_dp,8.0_dp], [4,2])
  call mapply_group_mean(ma_x, [2,1,2,1], ma_mean, ma_levels, ma_count, ma_status)
  call check_int(ma_status, 0, 'mApply status', failures)
  call check_int(ma_levels(1), 1, 'mApply level1', failures)
  call check_int(ma_levels(2), 2, 'mApply level2', failures)
  call check_close(ma_mean(1,1), 5.0_dp, 1.0e-12_dp, 'mApply group1 mean1', failures)
  call check_close(ma_mean(2,2), 4.0_dp, 1.0e-12_dp, 'mApply group2 mean2', failures)
  call check_int(ma_count(1,1), 2, 'mApply count', failures)

  so_off = 0.0_dp
  so_absorb = [.false., .false., .true.]
  call soprob_markov_ord([1.0_dp,-1.0_dp], [0.0_dp,0.0_dp], so_off, 1, so_absorb, so_prob, so_status)
  call check_int(so_status, 0, 'soprobMarkovOrd status', failures)
  call check_close(so_prob(1,1), 0.2689414213699951_dp, 2.0e-12_dp, &
                   'soprobMarkovOrd initial state1', failures)
  call check_close(sum(so_prob(1,:)), 1.0_dp, 2.0e-12_dp, 'soprobMarkovOrd initial sum', failures)
  call check_close(sum(so_prob(3,:)), 1.0_dp, 2.0e-12_dp, 'soprobMarkovOrd final sum', failures)
  if (so_prob(2,3) <= so_prob(1,3)) then
    print *, 'FAIL soprobMarkovOrd absorbing accumulation'
    failures = failures + 1
  end if


  call ord_group_boot_mean([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp], [2,3], 0.0_dp, &
                           og_grouped, og_m, og_status)
  call check_int(og_status, 0, 'ordGroupBoot status', failures)
  call check_int(og_m, 2, 'ordGroupBoot m', failures)
  call check_close(og_grouped(1), 1.5_dp, 1.0e-12_dp, 'ordGroupBoot first mean', failures)
  call check_close(og_grouped(6), 5.5_dp, 1.0e-12_dp, 'ordGroupBoot last mean', failures)

  call clowess_smooth([3.0_dp,1.0_dp,2.0_dp,4.0_dp], [7.0_dp,3.0_dp,5.0_dp,9.0_dp], &
                      1.0_dp, 0, clo_x, clo_y, clo_status)
  call check_int(clo_status, 0, 'clowess status', failures)
  call check_close(clo_x(1), 1.0_dp, 1.0e-12_dp, 'clowess sorted x', failures)
  call check_close(clo_y(1), 3.0_dp, 1.0e-10_dp, 'clowess linear first', failures)
  call check_close(clo_y(4), 9.0_dp, 1.0e-10_dp, 'clowess linear last', failures)

  call curve_smooth_observed([1.0_dp,2.0_dp,3.0_dp,1.0_dp,2.0_dp], &
                             [3.0_dp,5.0_dp,7.0_dp,10.0_dp,20.0_dp], [1,1,1,2,2], &
                             1.0_dp, 0, cs_x, cs_y, cs_id, cvs_status)
  call check_int(cvs_status, 0, 'curveSmooth status', failures)
  call check_close(cs_y(2), 5.0_dp, 1.0e-10_dp, 'curveSmooth smoothed linear', failures)
  call check_close(cs_y(4), 10.0_dp, 1.0e-12_dp, 'curveSmooth short curve unchanged', failures)

  call wtd_loess_noiter([1.0_dp,2.0_dp,3.0_dp,4.0_dp], [3.0_dp,5.0_dp,7.0_dp,9.0_dp], &
                          [1.0_dp,2.0_dp,1.0_dp,3.0_dp], 1.0_dp, 1, wl_fit, wl_status)
  call check_int(wl_status, 0, 'wtd.loess.noiter status', failures)
  call check_close(wl_fit(1), 3.0_dp, 1.0e-10_dp, 'wtd.loess.noiter first', failures)
  call check_close(wl_fit(4), 9.0_dp, 1.0e-10_dp, 'wtd.loess.noiter last', failures)

  sm_off = 0.0_dp
  sm_u = reshape([0.9_dp,0.2_dp, 0.1_dp,0.5_dp], [2,2])
  sm_absorb = [.false., .false., .true.]
  call sim_markov_ord([1.0_dp,-1.0_dp], sm_off, [1,1], sm_absorb, sm_u, .true., sm_states, smk_status)
  call check_int(smk_status, 0, 'simMarkovOrd status', failures)
  call check_int(sm_states(1,1), 3, 'simMarkovOrd absorb hit', failures)
  call check_int(sm_states(1,2), 3, 'simMarkovOrd carry absorb', failures)
  call check_int(sm_states(2,1), 1, 'simMarkovOrd first draw', failures)
  call check_int(sm_states(2,2), 2, 'simMarkovOrd second draw', failures)

  qnan = ieee_value(0.0_dp, ieee_quiet_nan)
  call impute_median([1.0_dp,qnan,5.0_dp,3.0_dp,qnan], imp_y, imp_idx, imp_status)
  call check_int(imp_status, 0, 'impute median status', failures)
  call check_int(size(imp_idx), 2, 'impute median count', failures)
  call check_int(imp_idx(1), 2, 'impute median index1', failures)
  call check_close(imp_y(2), 3.0_dp, 1.0e-12_dp, 'impute median fill1', failures)
  call check_close(imp_y(5), 3.0_dp, 1.0e-12_dp, 'impute median fill2', failures)
  call impute_constant([qnan,2.0_dp,qnan], -1.0_dp, imp_y, imp_idx, imp_status)
  call check_close(imp_y(1), -1.0_dp, 1.0e-12_dp, 'impute constant fill', failures)

  call mov_stats_n_raw([5.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp], &
                       [50.0_dp,10.0_dp,20.0_dp,30.0_dp,40.0_dp], &
                       1, 2, 1, mov_x, mov_mean, mov_q1, mov_med, mov_q3, mov_n, mov_status)
  call check_int(mov_status, 0, 'movStats status', failures)
  call check_int(size(mov_n), 3, 'movStats windows', failures)
  call check_int(mov_n(1), 3, 'movStats first N', failures)
  call check_close(mov_x(1), 2.0_dp, 1.0e-12_dp, 'movStats center1', failures)
  call check_close(mov_mean(2), 30.0_dp, 1.0e-12_dp, 'movStats mean2', failures)
  call check_close(mov_med(3), 40.0_dp, 1.0e-12_dp, 'movStats median3', failures)
  call check_close(mov_q1(1), 15.0_dp, 1.0e-12_dp, 'movStats q1', failures)
  call check_close(mov_q3(1), 25.0_dp, 1.0e-12_dp, 'movStats q3', failures)

  call summarize_mean_codes(reshape([1.0_dp,3.0_dp,5.0_dp,7.0_dp, &
                                     2.0_dp,4.0_dp,6.0_dp,8.0_dp], [4,2]), &
                            reshape([2,1,2,1, 1,1,2,2], [4,2]), &
                            sum_levels, sum_means, sum_counts, sum_status)
  call check_int(sum_status, 0, 'summarize status', failures)
  call check_int(size(sum_levels,1), 4, 'summarize groups', failures)
  call check_int(sum_levels(1,1), 1, 'summarize lexicographic first key1', failures)
  call check_int(sum_levels(1,2), 1, 'summarize lexicographic first key2', failures)
  call check_close(sum_means(1,1), 3.0_dp, 1.0e-12_dp, 'summarize first mean', failures)
  call check_int(sum_counts(4,2), 1, 'summarize count', failures)


  call areg_tran_numeric([1.0_dp,2.0_dp,3.0_dp], 'c', [real(dp) ::], 3, at_out, at_status)
  call check_int(at_status, 0, 'aregTran categorical status', failures)
  call check_int(size(at_out,2), 2, 'aregTran categorical columns', failures)
  call check_close(at_out(1,1), 0.0_dp, 1.0e-12_dp, 'aregTran reference level', failures)
  call check_close(at_out(2,1), 1.0_dp, 1.0e-12_dp, 'aregTran level2', failures)
  call check_close(at_out(3,2), 1.0_dp, 1.0e-12_dp, 'aregTran level3', failures)
  call areg_tran_numeric([0.0_dp,0.5_dp,1.0_dp], 's', [0.0_dp,0.5_dp,1.0_dp], 0, at_out, at_status)
  call check_int(at_status, 0, 'aregTran spline status', failures)
  call check_int(size(at_out,2), 2, 'aregTran spline columns', failures)

  call areg_linear_fit(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp], [4,1]), &
                       [3.0_dp,5.0_dp,7.0_dp,9.0_dp], al_coef, al_fit, al_res, &
                       al_r2, al_mae, al_medae, al_status)
  call check_int(al_status, 0, 'areg linear status', failures)
  call check_close(al_coef(1), 1.0_dp, 1.0e-10_dp, 'areg linear intercept', failures)
  call check_close(al_coef(2), 2.0_dp, 1.0e-10_dp, 'areg linear slope', failures)
  call check_close(al_r2, 1.0_dp, 1.0e-12_dp, 'areg linear r2', failures)
  call check_close(al_mae, 0.0_dp, 1.0e-10_dp, 'areg linear mean abs error', failures)

  qnan = ieee_value(0.0_dp, ieee_quiet_nan)
  call dataframe_reduce_numeric(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
                                         qnan,qnan,qnan,1.0_dp, &
                                         0.0_dp,0.0_dp,0.0_dp,1.0_dp], [4,3]), &
                                0.5_dp, 0.3_dp, dfr_red, dfr_keep, dfr_reason, dfr_status)
  call check_int(dfr_status, 0, 'dataframeReduce status', failures)
  call check_int(size(dfr_keep), 1, 'dataframeReduce kept count', failures)
  call check_int(dfr_keep(1), 1, 'dataframeReduce kept first', failures)
  call check_int(dfr_reason(2), 1, 'dataframeReduce missing reason', failures)
  call check_int(dfr_reason(3), 2, 'dataframeReduce prevalence reason', failures)

  call data_rep_numeric(reshape([1.0_dp,1.0_dp,2.0_dp,qnan, &
                                 5.0_dp,5.0_dp,6.0_dp,7.0_dp], [4,2]), &
                        dr_unique, dr_counts, dr_n, dr_status)
  call check_int(dr_status, 0, 'dataRep status', failures)
  call check_int(dr_n, 3, 'dataRep complete rows', failures)
  call check_int(size(dr_counts), 2, 'dataRep unique rows', failures)
  call check_int(dr_counts(1), 2, 'dataRep duplicate count', failures)

  call varclus_numeric(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp, &
                                1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp, &
                                1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp], [5,3]), &
                       1, 1, 1, vc_sim, vc_np, vc_left, vc_right, vc_height, vc_status)
  call check_int(vc_status, 0, 'varclus status', failures)
  call check_close(vc_sim(1,2), 1.0_dp, 1.0e-12_dp, 'varclus identical similarity', failures)
  call check_close(vc_height(1), 0.0_dp, 1.0e-12_dp, 'varclus first merge height', failures)

  call redun_numeric(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp, &
                              2.0_dp,4.0_dp,6.0_dp,8.0_dp,10.0_dp, &
                              1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp], [5,3]), &
                     0.99_dp, .false., rd_keep, rd_remove, rd_r2, rd_status)
  call check_int(rd_status, 0, 'redun status', failures)
  call check_int(size(rd_remove), 1, 'redun removed count', failures)
  call check_close(rd_r2(1), 1.0_dp, 1.0e-10_dp, 'redun removed r2', failures)
  call check_int(size(rd_keep), 2, 'redun kept count', failures)

  call gbayes_seq([0.0_dp], [1.0_dp], [-1,0], [0.0_dp,-1.0_dp], [0.0_dp,1.0_dp], &
                  [0.0_dp,0.0_dp], [1.0_dp,1.0_dp], gb_prob, gb_pm, gb_ps, gbs_status)
  call check_int(gbs_status, 0, 'gbayesSeqSim status', failures)
  call check_close(gb_prob(1,1), 0.5_dp, 1.0e-12_dp, 'gbayesSeqSim tail', failures)
  call check_close(gb_pm(1,1), 0.0_dp, 1.0e-12_dp, 'gbayesSeqSim posterior mean', failures)
  call check_close(gb_ps(1,1), sqrt(0.5_dp), 1.0e-12_dp, 'gbayesSeqSim posterior sd', failures)
  call check_close(gb_prob(1,2), 0.8427007929497149_dp, 1.0e-10_dp, 'gbayesSeqSim interval', failures)

  block
    real(dp) :: tr(2,2,2)
    tr(:,:,1) = reshape([0.8_dp,0.1_dp,0.2_dp,0.9_dp], [2,2])
    tr(:,:,2) = tr(:,:,1)
    call soprob_markov_ordm_core([1.0_dp,0.0_dp], tr, mo_occ, mo_status)
  end block
  call check_int(mo_status, 0, 'soprobMarkovOrdm status', failures)
  call check_close(mo_occ(2,1), 0.8_dp, 1.0e-12_dp, 'soprobMarkovOrdm t2 state1', failures)
  call check_close(mo_occ(3,1), 0.66_dp, 1.0e-12_dp, 'soprobMarkovOrdm t3 state1', failures)
  call check_close(sum(mo_occ(3,:)), 1.0_dp, 1.0e-12_dp, 'soprobMarkovOrdm probability sum', failures)


  block
    real(dp), allocatable :: ab_coef(:), ab_boot(:,:), fmi_mean(:), fmi_cov(:,:), fmi_vi(:), fmi_mi(:), fmi_df(:)
    real(dp), allocatable :: imp_out(:), es_est(:,:,:), es_var(:,:,:), rmb_coef(:), rmb_boot(:,:), rmb_fit(:,:)
    real(dp) :: ab_r2, ab_r2v
    real(dp) :: xv(5,1), yv(5), vv(2,2,3), cc(2,3)
    integer :: bi(5,2), ab_nf, st2, pos(2), grp(6,1,1), looks2(2), ri(5,2)
    real(dp) :: resp(6,1,1)

    xv(:,1) = [0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    yv = 1.0_dp + 2.0_dp * xv(:,1)
    bi(:,1) = [1,2,3,4,5]
    bi(:,2) = [5,4,3,2,1]
    call areg_boot_linear(xv, yv, bi, ab_coef, ab_r2, ab_r2v, ab_boot, ab_nf, st2)
    call check_int(st2, 0, 'areg.boot linear status', failures)
    call check_int(ab_nf, 0, 'areg.boot linear failures', failures)
    call check_close(ab_coef(1), 1.0_dp, 1.0e-10_dp, 'areg.boot intercept', failures)
    call check_close(ab_coef(2), 2.0_dp, 1.0e-10_dp, 'areg.boot slope', failures)
    call check_close(ab_r2v, 1.0_dp, 1.0e-10_dp, 'areg.boot validated r2', failures)

    cc = reshape([1.0_dp,2.0_dp, 1.2_dp,2.2_dp, 0.8_dp,1.8_dp], [2,3])
    vv = 0.0_dp
    vv(1,1,:) = 0.04_dp
    vv(2,2,:) = 0.04_dp
    call fit_mult_impute_combine(cc, vv, fmi_mean, fmi_cov, fmi_vi, fmi_mi, fmi_df, st2)
    call check_int(st2, 0, 'fit.mult.impute status', failures)
    call check_close(fmi_mean(1), 1.0_dp, 1.0e-12_dp, 'fit.mult.impute mean 1', failures)
    call check_close(fmi_mean(2), 2.0_dp, 1.0e-12_dp, 'fit.mult.impute mean 2', failures)
    call check_close(fmi_cov(1,1), 0.0933333333333333_dp, 1.0e-12_dp, 'fit.mult.impute variance', failures)

    pos = [2,4]
    call impute_transcan_numeric([1.0_dp,qnan,3.0_dp,qnan], pos, &
                                 reshape([20.0_dp,40.0_dp,21.0_dp,41.0_dp],[2,2]), 2, imp_out, st2)
    call check_int(st2, 0, 'impute.transcan status', failures)
    call check_close(imp_out(2), 21.0_dp, 1.0e-12_dp, 'impute.transcan first', failures)
    call check_close(imp_out(4), 41.0_dp, 1.0e-12_dp, 'impute.transcan second', failures)

    grp(:,1,1) = [0,1,0,1,0,1]
    resp(:,1,1) = [1.0_dp,3.0_dp,2.0_dp,4.0_dp,3.0_dp,5.0_dp]
    looks2 = [4,6]
    call est_seq_sim_mean_difference(grp, resp, looks2, es_est, es_var, st2)
    call check_int(st2, 0, 'estSeqSim status', failures)
    call check_close(es_est(1,1,1), 2.0_dp, 1.0e-12_dp, 'estSeqSim look1 estimate', failures)
    call check_close(es_var(1,1,1), 0.5_dp, 1.0e-12_dp, 'estSeqSim look1 variance', failures)
    call check_close(es_var(2,1,1), 2.0_dp/3.0_dp, 1.0e-12_dp, 'estSeqSim look2 variance', failures)

    ri(:,1) = [1,2,3,4,5]
    ri(:,2) = [5,4,3,2,1]
    call rm_boot_linear(xv(:,1), yv, [0.5_dp,2.5_dp], ri, rmb_coef, rmb_boot, rmb_fit, st2)
    call check_int(st2, 0, 'rm.boot status', failures)
    call check_close(rmb_coef(1), 1.0_dp, 1.0e-10_dp, 'rm.boot intercept', failures)
    call check_close(rmb_coef(2), 2.0_dp, 1.0e-10_dp, 'rm.boot slope', failures)
    call check_close(rmb_fit(2,1), 6.0_dp, 1.0e-10_dp, 'rm.boot original fitted', failures)
  end block


  block
    real(dp) :: lr, pv, bpo, crit
    real(dp), allocatable :: probs(:,:), ints_fit(:)
    real(dp) :: grp_po(12), initoff(2), toff(1,3,2), targ(2,3)
    integer :: ypo(12), dfo, st3
    logical :: absorb3(3)

    grp_po = [0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
              1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp]
    ypo = [1,1,1,2,2,2, 2,2,3,3,3,3]
    call ord_test_po_numeric(grp_po, ypo, lr, dfo, pv, bpo, st3)
    call check_int(st3, 0, 'ordTestpo status', failures)
    call check_int(dfo, 1, 'ordTestpo df', failures)
    if (lr <= 1.0_dp .or. pv >= 0.5_dp .or. bpo <= 0.0_dp) then
      print *, 'FAIL ordTestpo association', lr, pv, bpo
      failures = failures + 1
    end if

    initoff = 0.0_dp
    toff = 0.0_dp
    absorb3 = .false.
    call soprob_markov_ord([1.0_dp,-1.0_dp], initoff, toff, 1, absorb3, probs, st3)
    call check_int(st3, 0, 'intMarkovOrd target generation', failures)
    targ(1,:) = probs(1,:)
    targ(2,:) = probs(2,:)
    call int_markov_ord_intercepts([0.4_dp,-0.4_dp], initoff, toff, 1, absorb3, [1,2], targ, ints_fit, crit, st3)
    call check_int(st3, 0, 'intMarkovOrd status', failures)
    if (crit > 1.0e-4_dp) then
      print *, 'FAIL intMarkovOrd criterion', crit
      failures = failures + 1
    end if
    call check_close(ints_fit(1), 1.0_dp, 2.0e-3_dp, 'intMarkovOrd intercept1', failures)
    call check_close(ints_fit(2), -1.0_dp, 2.0e-3_dp, 'intMarkovOrd intercept2', failures)
  end block


  block
    real(dp), allocatable :: sb(:), ss(:), sp(:)
    real(dp) :: tr(12), pw
    integer :: yo(12,3), st4
    tr = [0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
          1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp]
    yo(:,1) = [1,1,1,2,2,2, 2,2,3,3,3,3]
    yo(:,2) = [1,1,2,2,2,3, 1,2,2,2,3,3]
    yo(:,3) = [1,1,1,1,2,2, 2,3,3,3,3,3]
    call sim_reg_ord_supplied(tr, yo, 0.05_dp, sb, ss, sp, pw, st4)
    call check_int(st4, 0, 'simRegOrd status', failures)
    if (any(ieee_is_nan(sb)) .or. any(ieee_is_nan(sp)) .or. pw < 0.0_dp .or. pw > 1.0_dp) then
      print *, 'FAIL simRegOrd outputs', pw
      failures = failures + 1
    end if
    if (sb(1) <= 0.0_dp .or. sb(3) <= 0.0_dp) then
      print *, 'FAIL simRegOrd treatment slope', sb
      failures = failures + 1
    end if
  end block

  if (failures /= 0) error stop 1
  print *, 'All Hmisc tests passed.'
contains
  subroutine check_close(actual, expected, tol, label, failures)
    real(dp), intent(in) :: actual !! Computed scalar.
    real(dp), intent(in) :: expected !! Expected scalar.
    real(dp), intent(in) :: tol !! Absolute tolerance.
    character(*), intent(in) :: label !! Test label.
    integer, intent(inout) :: failures !! Running failure count.
    if (abs(actual-expected) > tol .or. ieee_is_nan(actual)) then
      print *, 'FAIL ', trim(label), actual, expected
      failures = failures + 1
    end if
  end subroutine check_close
  subroutine check_int(actual, expected, label, failures)
    integer, intent(in) :: actual !! Computed integer.
    integer, intent(in) :: expected !! Expected integer.
    character(*), intent(in) :: label !! Test label.
    integer, intent(inout) :: failures !! Running failure count.
    if (actual /= expected) then
      print *, 'FAIL ', trim(label), actual, expected
      failures = failures + 1
    end if
  end subroutine check_int
end program test_hmisc
