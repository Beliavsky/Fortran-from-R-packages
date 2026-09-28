program test_more
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use psych, only : dp, cta_dynamics, cta15_dynamics, stats_by_bootstrap, score_by_correlations, &
      mediate_from_covariance, factor_analysis_sapa, factor_analysis_pooled, parallel_analysis_poly, &
      simulate_omega_study, bass_ackward_factor, factor_analysis_by_groups, simulate_multi_subject, &
      best_scales_from_correlation, exploratory_sem, correlation_fiml, test_irt_benchmark, phi_demo_simulation, &
      glb_algebraic_admm
   implicit none
   real(dp), allocatable :: tt(:, :), at(:, :), frac(:), mt(:), ma(:), expf(:, :)
   real(dp), allocatable :: adj(:, :, :), al(:, :), avr(:, :), g6(:, :), svar(:, :)
   real(dp), allocatable :: a(:, :), b(:, :), ct(:, :), cd(:, :), ip(:, :), id(:, :)
   real(dp), allocatable :: ref(:, :), meanl(:, :), sdl(:, :), pooled(:, :, :)
   real(dp), allocatable :: opc(:), rpc(:), ofa(:), rfa(:), omres(:, :)
   real(dp), allocatable :: icc1(:, :), icc2(:, :), ftrial(:, :)
   real(dp), allocatable :: backl(:, :, :), backphi(:, :, :), trans(:, :, :)
   real(dp), allocatable :: gl(:, :, :), pld(:, :), gloads(:, :, :), mld(:, :), sld(:, :), mph(:, :)
   real(dp), allocatable :: multi(:, :, :), tm(:), bkeys(:, :), bvalid(:), bcut(:)
   integer, allocatable :: status(:), parents(:, :), backinfo(:), groupinfo(:)
   integer :: info, nc, nf, nseed
   integer, allocatable :: seed(:)
   integer :: groups(8)
   real(dp) :: x(8, 2), keys(4, 2), corrs(4, 4, 1), cov(3, 3), r(4, 4)
   real(dp) :: ordinal(12, 3), r5(5, 5)
   real(dp) :: r6(6, 6), l6(6, 2), ph0(2, 2), u6(6), nanv
   real(dp), allocatable :: esem_load(:, :), esem_structure(:, :), esem_phi(:, :), esem_resid(:, :)
   real(dp), allocatable :: fiml_mean(:), fiml_cov(:, :), fiml_cor(:, :)
   real(dp) :: esem_rms, esem_fit, fiml_obj, fiml_dat(4, 2), fiml_miss(6, 2)
   integer :: fiml_iters
   real(dp), allocatable :: irt_truth(:), irt_estimate(:), irt_params(:, :)
   real(dp), allocatable :: demo_tetra(:, :), demo_phis(:, :), demo_yule(:, :), demo_data(:, :)
   real(dp) :: irt_score_cor
   real(dp), allocatable :: glb_solution(:)
   real(dp) :: glb_value, glb_primal, glb_dual, glb_cov(3, 3)
   integer :: glb_iterations

   call cta15_dynamics([1.0_dp], 1, tt, at, frac, mt, ma, info=info)
   call check(info == 0, 'cta15 info')
   call close(tt(1, 1), 0.0475_dp, 1.0e-12_dp, 'cta15 tendency')
   call close(at(1, 1), 0.999875_dp, 1.0e-12_dp, 'cta15 action')
   call close(frac(1), 1.0_dp, 1.0e-12_dp, 'cta15 fraction')

   call cta_dynamics([1.0_dp], 1, tt, at, frac, expf, learn=.false., info=info)
   call check(info == 0, 'cta info')
   call close(tt(1, 1), 0.0485_dp, 1.0e-12_dp, 'cta tendency')
   call close(at(1, 1), 0.999925_dp, 1.0e-12_dp, 'cta action')

   keys = 0.0_dp
   keys(1:2, 1) = 1.0_dp
   keys(3:4, 2) = 1.0_dp
   corrs(:, :, 1) = reshape([ &
      1.0_dp, 0.5_dp, 0.1_dp, 0.1_dp, &
      0.5_dp, 1.0_dp, 0.1_dp, 0.1_dp, &
      0.1_dp, 0.1_dp, 1.0_dp, 0.5_dp, &
      0.1_dp, 0.1_dp, 0.5_dp, 1.0_dp], [4, 4])
   call score_by_correlations(keys, corrs, adj, al, avr, g6, svar, info=status)
   call check(all(status == 0), 'scoreBy info')
   call close(al(1, 1), 2.0_dp / 3.0_dp, 1.0e-10_dp, 'scoreBy alpha')
   call close(svar(1, 1), 3.0_dp, 1.0e-12_dp, 'scoreBy variance')
   call close(adj(1, 2, 1), 0.4_dp / 3.0_dp, 1.0e-10_dp, 'scoreBy correlation')

   cov = reshape([1.0_dp, 0.5_dp, 0.4_dp, 0.5_dp, 1.0_dp, 0.6_dp, 0.4_dp, 0.6_dp, 1.0_dp], [3, 3])
   call mediate_from_covariance(cov, [1], [2], [3], a, b, ct, cd, ip, id, info)
   call check(info == 0, 'mediate info')
   call close(a(1, 1), 0.5_dp, 1.0e-10_dp, 'mediate a')
   call close(b(1, 1), 0.5333333333333333_dp, 1.0e-8_dp, 'mediate b')
   call close(ct(1, 1), 0.4_dp, 1.0e-10_dp, 'mediate total')
   call close(cd(1, 1), 0.1333333333333333_dp, 1.0e-8_dp, 'mediate direct')
   call close(ip(1, 1), id(1, 1), 1.0e-8_dp, 'mediate indirect equality')

   r = reshape([ &
      1.0_dp, 0.56_dp, 0.48_dp, 0.40_dp, &
      0.56_dp, 1.0_dp, 0.42_dp, 0.35_dp, &
      0.48_dp, 0.42_dp, 1.0_dp, 0.30_dp, &
      0.40_dp, 0.35_dp, 0.30_dp, 1.0_dp], [4, 4])
   allocate(pooled(4, 4, 2))
   pooled(:, :, 1) = r
   pooled(:, :, 2) = r
   call factor_analysis_pooled(pooled, 1, ref, meanl, sdl, info)
   call check(info == 0, 'pooled info')
   call check(all(shape(meanl) == [4, 1]), 'pooled shape')
   call check(maxval(abs(sdl)) < 1.0e-10_dp, 'pooled identical sd')

   call random_seed(size=nseed)
   allocate(seed(nseed))
   seed = 24680
   call random_seed(put=seed)
   call factor_analysis_sapa(r, 1, 80, 3, ref, meanl, sdl, info)
   call check(info == 0, 'sapa info')
   call check(all(shape(meanl) == [4, 1]), 'sapa shape')
   call check(all(sdl >= 0.0_dp), 'sapa sd')

   ordinal = reshape([ &
      1.0_dp,1.0_dp,1.0_dp, 1.0_dp,1.0_dp,2.0_dp, 1.0_dp,2.0_dp,1.0_dp, 1.0_dp,2.0_dp,2.0_dp, &
      2.0_dp,1.0_dp,2.0_dp, 2.0_dp,2.0_dp,2.0_dp, 2.0_dp,2.0_dp,3.0_dp, 3.0_dp,2.0_dp,3.0_dp, &
      3.0_dp,3.0_dp,2.0_dp, 3.0_dp,3.0_dp,3.0_dp, 2.0_dp,3.0_dp,3.0_dp, 3.0_dp,2.0_dp,2.0_dp], [12, 3])
   seed = 13579
   call random_seed(put=seed)
   call parallel_analysis_poly(ordinal, 2, opc, rpc, ofa, rfa, nc, nf, info)
   call check(info == 0, 'parallel poly info')
   call check(size(opc) == 3 .and. size(rpc) == 3, 'parallel poly shape')
   call check(nc >= 0 .and. nf >= 0, 'parallel poly counts')

   x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
   x(:, 2) = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 6.0_dp, 5.0_dp, 8.0_dp, 7.0_dp]
   groups = [1, 1, 2, 2, 3, 3, 4, 4]
   seed = 999
   call random_seed(put=seed)
   call stats_by_bootstrap(x, groups, 3, icc1, icc2, ftrial, replace=.false., info=info)
   call check(info == 0, 'statsBy boot info')
   call check(all(shape(icc1) == [3, 2]), 'statsBy boot shape')

   seed = 54321
   call random_seed(put=seed)
   call simulate_omega_study(6, 2, 0, 2, 0.25_dp, omres, info)
   call check(info == 0, 'sim omega info')
   call check(all(shape(omres) == [2, 8]), 'sim omega shape')
   call check(all(omres(:, 2) >= 0.0_dp), 'sim omega nonnegative')

   call bass_ackward_factor(r, 2, backl, backphi, trans, parents, backinfo)
   call check(all(backinfo == 0), 'bassAckward info')
   call check(all(shape(backl) == [4, 2, 2]), 'bassAckward shape')
   call check(parents(1, 2) >= 1 .and. parents(1, 2) <= 2, 'bassAckward parent')

   allocate(gl(4, 4, 2))
   gl(:, :, 1) = r
   gl(:, :, 2) = r
   call factor_analysis_by_groups(r, gl, 1, pld, gloads, mld, sld, mph, groupinfo)
   call check(all(groupinfo == 0), 'faBy info')
   call check(all(shape(mld) == [4, 1]), 'faBy shape')
   call check(maxval(abs(sld)) < 1.0e-10_dp, 'faBy identical sd')

   seed = 31415
   call random_seed(put=seed)
   call simulate_multi_subject(2, 4, 2, 8, 2, multi, tm, info=info)
   call check(info == 0, 'sim.multi info')
   call check(all(shape(multi) == [8, 6, 2]), 'sim.multi shape')
   call close(tm(1), 6.0_dp, 1.0e-12_dp, 'sim.multi time')

   r5 = 0.0_dp
   r5(1, :) = [1.0_dp, -0.3_dp, 0.1_dp, 0.0_dp, 0.8_dp]
   r5(2, :) = [-0.3_dp, 1.0_dp, 0.0_dp, 0.1_dp, -0.6_dp]
   r5(3, :) = [0.1_dp, 0.0_dp, 1.0_dp, 0.2_dp, 0.2_dp]
   r5(4, :) = [0.0_dp, 0.1_dp, 0.2_dp, 1.0_dp, 0.05_dp]
   r5(5, :) = [0.8_dp, -0.6_dp, 0.2_dp, 0.05_dp, 1.0_dp]
   call best_scales_from_correlation(r5, [5], 2, 0.1_dp, bkeys, bvalid, bcut, info)
   call check(info == 0, 'bestScales info')
   call close(bcut(1), 0.6_dp, 1.0e-12_dp, 'bestScales cut')
   call close(bkeys(1, 1), 1.0_dp, 1.0e-12_dp, 'bestScales positive key')
   call close(bkeys(2, 1), -1.0_dp, 1.0e-12_dp, 'bestScales negative key')


   l6 = 0.0_dp
   l6(:, 1) = [0.8_dp, 0.7_dp, 0.75_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   l6(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp, 0.8_dp, 0.7_dp, 0.75_dp]
   ph0 = reshape([1.0_dp, 0.3_dp, 0.3_dp, 1.0_dp], [2, 2])
   r6 = matmul(l6, matmul(ph0, transpose(l6)))
   do nc = 1, 6
      u6(nc) = 1.0_dp - r6(nc, nc)
      r6(nc, nc) = r6(nc, nc) + u6(nc)
   end do
   call exploratory_sem(r6, [1, 2, 3], [4, 5, 6], 1, 1, esem_load, esem_structure, esem_phi, esem_resid, &
      esem_rms, esem_fit, info)
   call check(info == 0, 'esem info')
   call check(all(shape(esem_load) == [6, 2]), 'esem loading shape')
   call check(esem_rms < 0.08_dp .and. esem_fit > 0.90_dp, 'esem fit')
   call check(abs(esem_phi(1, 2)) > 0.05_dp .and. abs(esem_phi(1, 2)) < 0.8_dp, 'esem factor relation')

   fiml_dat = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 5.0_dp], [4, 2])
   call correlation_fiml(fiml_dat, fiml_mean, fiml_cov, fiml_cor, fiml_obj, fiml_iters, info)
   call check(info == 0, 'fiml complete info')
   call close(fiml_mean(1), 2.5_dp, 1.0e-10_dp, 'fiml mean')
   call close(fiml_cov(1, 1), 1.25_dp, 1.0e-10_dp, 'fiml variance 1')
   call close(fiml_cov(2, 2), 2.25_dp, 1.0e-10_dp, 'fiml variance 2')
   call close(fiml_cov(1, 2), 1.5_dp, 1.0e-10_dp, 'fiml covariance')
   call check(ieee_is_finite(fiml_obj), 'fiml complete objective')

   nanv = ieee_value(0.0_dp, ieee_quiet_nan)
   fiml_miss = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, nanv, 6.0_dp, &
      1.0_dp, 2.1_dp, nanv, 4.2_dp, 5.3_dp, 6.2_dp], [6, 2])
   call correlation_fiml(fiml_miss, fiml_mean, fiml_cov, fiml_cor, fiml_obj, fiml_iters, info, &
      max_iter=1000, tolerance=1.0e-7_dp)
   call check(info == 0, 'fiml missing info')
   call check(all([(ieee_is_finite(fiml_mean(nc)), nc = 1, size(fiml_mean))]), 'fiml missing means finite')
   call check(ieee_is_finite(fiml_obj), 'fiml missing objective')
   call check(fiml_cor(1, 2) > 0.9_dp .and. fiml_cor(1, 2) <= 1.0_dp, 'fiml missing correlation')


   seed = 24681357
   call random_seed(put=seed)
   call test_irt_benchmark(5, 800, irt_truth, irt_estimate, irt_params, irt_score_cor, info, &
      low=-1.5_dp, high=1.5_dp)
   call check(info == 0, 'test.irt info')
   call check(all(shape(irt_params) == [5, 4]), 'test.irt parameter shape')
   call check(irt_score_cor > 0.6_dp .and. irt_score_cor <= 1.0_dp, 'test.irt trait recovery')

   seed = 97531
   call random_seed(put=seed)
   call phi_demo_simulation(1500, 0.6_dp, [-1.0_dp, 0.0_dp, 1.0_dp], demo_tetra, demo_phis, demo_yule, &
      demo_data, info)
   call check(info == 0, 'phi.demo info')
   call check(all(shape(demo_tetra) == [3, 3]), 'phi.demo tetra shape')
   call check(all(shape(demo_phis) == [5, 5]), 'phi.demo phi shape')
   call check(all(shape(demo_yule) == [3, 3]), 'phi.demo yule shape')
   call check(abs(demo_phis(1, 2) - 0.6_dp) < 0.08_dp, 'phi.demo latent-observed correlation')


   glb_cov = 0.4_dp
   do nc = 1, 3
      glb_cov(nc, nc) = 1.0_dp
   end do
   call glb_algebraic_admm(glb_cov, glb_value, glb_solution, glb_iterations, glb_primal, glb_dual, info, &
      tolerance=1.0e-9_dp)
   call check(info == 0, 'glb algebraic info')
   call close(glb_value, 2.0_dp / 3.0_dp, 3.0e-7_dp, 'glb algebraic equicorrelation')
   call check(maxval(abs(glb_solution - 0.4_dp)) < 3.0e-7_dp, 'glb algebraic diagonal')

   print *, 'more extended focused tests passed'
contains
   subroutine check(ok, label)
      logical, intent(in) :: ok !! Condition expected to be true.
      character(len=*), intent(in) :: label !! Description printed if the assertion fails.
      if (.not. ok) then
         print *, 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check
   subroutine close(actual, expected, tol, label)
      real(dp), intent(in) :: actual !! Computed value to compare with the reference.
      real(dp), intent(in) :: expected !! Reference value for the comparison.
      real(dp), intent(in) :: tol !! Absolute tolerance accepted by the comparison.
      character(len=*), intent(in) :: label !! Description printed if the comparison fails.
      call check(abs(actual - expected) <= tol, label)
   end subroutine close
end program test_more
