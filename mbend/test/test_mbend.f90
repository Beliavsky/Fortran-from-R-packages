program test_mbend
   use mbend_api, only : bend, bend_result, dp, mbend_err_already_pd, mbend_err_method, mbend_success
   implicit none

   call test_hj_unweighted()
   call test_hj_weighted()
   call test_fixed_elements()
   call test_reciprocal_equivalence()
   call test_correlation_lrs()
   call test_errors()
   print '(a)', 'test_mbend: PASS'

contains

   subroutine test_hj_unweighted()
      real(dp) :: v(5, 5)
      type(bend_result) :: fit
      integer :: info

      call make_v(v)
      call bend(v, fit, info)
      call require(info == mbend_success, 'HJ unweighted status')
      call require(fit%converged, 'HJ unweighted convergence')
      call require(fit%iterations == 1, 'HJ unweighted iteration count')
      call require_close(fit%bent(1, 1), 103.16918_dp, 5.0e-5_dp, 'HJ bent(1,1)')
      call require_close(fit%bent(2, 4), 74.07039_dp, 5.0e-5_dp, 'HJ bent(2,4)')
      call require_close(fit%init_ev(5), -18.523500_dp, 5.0e-6_dp, 'HJ initial minimum eigenvalue')
      call require_close(fit%final_ev(5), 0.0001_dp, 5.0e-8_dp, 'HJ final minimum eigenvalue')
      call require_close(fit%min_dev, -5.929612_dp, 5.0e-6_dp, 'HJ minimum deviation')
      call require_close(fit%max_dev, 6.497343_dp, 5.0e-6_dp, 'HJ maximum deviation')
      call require_close(fit%bent(fit%loc_min_dev(1), fit%loc_min_dev(2)) - &
                         v(fit%loc_min_dev(1), fit%loc_min_dev(2)), fit%min_dev, 1.0e-12_dp, &
                         'HJ minimum location')
      call require_close(fit%bent(fit%loc_max_dev(1), fit%loc_max_dev(2)) - &
                         v(fit%loc_max_dev(1), fit%loc_max_dev(2)), fit%max_dev, 1.0e-12_dp, &
                         'HJ maximum location')
      call require_close(fit%aad, 3.372692_dp, 5.0e-6_dp, 'HJ AAD')
      call require_close(fit%cor, 0.985632_dp, 5.0e-6_dp, 'HJ correlation')
      call require_close(fit%rmsd, 3.927457_dp, 5.0e-6_dp, 'HJ RMSD')
   end subroutine test_hj_unweighted

   subroutine test_hj_weighted()
      real(dp) :: v(5, 5), w(5, 5)
      type(bend_result) :: fit
      integer :: info

      call make_v(v)
      call make_w(w)
      call bend(v, fit, info, wtmat=w, reciprocal=.true.)
      call require(info == mbend_success, 'HJ weighted status')
      call require(fit%converged, 'HJ weighted convergence')
      call require(fit%iterations == 428, 'HJ weighted iteration count')
      call require_close(fit%bent(1, 1), 100.16154_dp, 5.0e-5_dp, 'weighted bent(1,1)')
      call require_close(fit%bent(2, 4), 59.99921_dp, 5.0e-5_dp, 'weighted bent(2,4)')
      call require(fit%has_weighted_stats, 'weighted statistics flag')
      call require(fit%w_gt_0 == 15, 'weighted positive count')
      call require_close(fit%waad, 0.6100103_dp, 5.0e-6_dp, 'weighted AAD')
      call require_close(fit%wcor, 0.9955121_dp, 5.0e-6_dp, 'weighted correlation')
      call require_close(fit%wrmsd, 0.5326555_dp, 5.0e-6_dp, 'weighted RMSD')
   end subroutine test_hj_weighted

   subroutine test_fixed_elements()
      real(dp) :: v(5, 5), w(5, 5)
      type(bend_result) :: fit
      integer :: info

      call make_v(v)
      call make_w(w)
      w(1:2, 1:2) = 0.0_dp
      call bend(v, fit, info, wtmat=w, reciprocal=.true.)
      call require(info == mbend_success, 'fixed-element status')
      call require(fit%converged, 'fixed-element convergence')
      call require(maxval(abs(fit%bent(1:2, 1:2) - v(1:2, 1:2))) < 1.0e-12_dp, 'fixed block unchanged')
      call require(fit%w_gt_0 == 12, 'fixed-element positive weight count')
   end subroutine test_fixed_elements


   subroutine test_reciprocal_equivalence()
      real(dp) :: v(5, 5), w(5, 5), winv(5, 5)
      type(bend_result) :: fit_a, fit_b
      integer :: i, info_a, info_b, j

      call make_v(v)
      call make_w(w)
      winv = w
      do j = 1, 5
         do i = 1, 5
            if (winv(i, j) > 0.0_dp) winv(i, j) = 1.0_dp / winv(i, j)
         end do
      end do
      call bend(v, fit_a, info_a, wtmat=w, reciprocal=.true.)
      call bend(v, fit_b, info_b, wtmat=winv)
      call require(info_a == mbend_success .and. info_b == mbend_success, 'reciprocal equivalence status')
      call require(maxval(abs(fit_a%bent - fit_b%bent)) < 1.0e-10_dp, 'reciprocal equivalence matrix')
   end subroutine test_reciprocal_equivalence

   subroutine test_correlation_lrs()
      real(dp) :: d(5), v(5, 5), v2(5, 5), w(5, 5)
      type(bend_result) :: fit
      integer :: i, info, j

      call make_v(v)
      call make_w(w)
      do i = 1, 5
         d(i) = sqrt(v(i, i))
      end do
      do j = 1, 5
         do i = 1, 5
            v2(i, j) = v(i, j) / (d(i) * d(j))
         end do
      end do
      call bend(v2, fit, info, wtmat=w, reciprocal=.true., method='lrs')
      call require(info == mbend_success, 'LRS correlation status')
      call require(fit%converged, 'LRS correlation convergence')
      call require(fit%is_correlation, 'correlation detection')
      call require(fit%iterations > 0 .and. fit%iterations <= 10000, 'LRS correlation iteration count')
      call require_close(fit%bent(1, 2), 0.9447382578_dp, 5.0e-5_dp, 'LRS correlation bent(1,2)')
      call require(maxval(abs([(fit%bent(i, i) - 1.0_dp, i=1,5)])) < 1.0e-12_dp, 'correlation diagonal fixed')
      call require(minval(fit%final_ev) > 0.0_dp, 'LRS correlation positive definite')
   end subroutine test_correlation_lrs

   subroutine test_errors()
      real(dp) :: identity(2, 2), v(5, 5)
      type(bend_result) :: fit
      integer :: info

      identity = 0.0_dp
      identity(1, 1) = 1.0_dp
      identity(2, 2) = 1.0_dp
      call bend(identity, fit, info)
      call require(info == mbend_err_already_pd, 'already-PD status')

      call make_v(v)
      call bend(v, fit, info, method='bad')
      call require(info == mbend_err_method, 'invalid-method status')
   end subroutine test_errors

   subroutine make_v(v)
      real(dp), intent(out) :: v(5, 5) !! Canonical covariance matrix from the upstream mbend documentation.

      v = reshape([ &
         100.0_dp, 95.0_dp, 80.0_dp, 40.0_dp, 40.0_dp, &
          95.0_dp,100.0_dp, 95.0_dp, 80.0_dp, 40.0_dp, &
          80.0_dp, 95.0_dp,100.0_dp, 95.0_dp, 80.0_dp, &
          40.0_dp, 80.0_dp, 95.0_dp,100.0_dp, 95.0_dp, &
          40.0_dp, 40.0_dp, 80.0_dp, 95.0_dp,100.0_dp], [5, 5])
   end subroutine make_v

   subroutine make_w(w)
      real(dp), intent(out) :: w(5, 5) !! Canonical weight matrix from the upstream mbend documentation.

      w = reshape([ &
         1000.0_dp, 500.0_dp,  20.0_dp,  50.0_dp, 200.0_dp, &
          500.0_dp,1000.0_dp, 500.0_dp,   5.0_dp,  50.0_dp, &
           20.0_dp, 500.0_dp,1000.0_dp,  20.0_dp,  20.0_dp, &
           50.0_dp,   5.0_dp,  20.0_dp,1000.0_dp, 200.0_dp, &
          200.0_dp,  50.0_dp,  20.0_dp, 200.0_dp,1000.0_dp], [5, 5])
   end subroutine make_w

   subroutine require(condition, label)
      logical, intent(in) :: condition !! Condition that must be true for the regression test to pass.
      character(len=*), intent(in) :: label !! Human-readable name of the tested condition.

      if (.not. condition) then
         write (*, '(a)') 'FAIL: ' // trim(label)
         error stop 1
      end if
   end subroutine require

   subroutine require_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Value produced by the Fortran translation.
      real(dp), intent(in) :: expected !! Upstream reference value to compare against.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
      character(len=*), intent(in) :: label !! Human-readable name of the tested value.

      call require(abs(actual - expected) <= tolerance, label)
   end subroutine require_close

end program test_mbend
