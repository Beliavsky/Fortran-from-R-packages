! SPDX-License-Identifier: GPL-2.0-only
module marss_harvey_info
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_kf_result, marss_constraints
   use marss_types, only : marss_constraint_block
   use marss_kalman, only : marss_kfss
   use marss_parameters, only : marss_b_at, marss_u_at, marss_q_at
   use marss_parameters, only : marss_z_at, marss_a_at, marss_r_at, marss_v0_effective
   use marss_covariance, only : psd_inverse
   use marss_constraints_mod, only : marss_constraints_parameter_count, marss_apply_constraints
   use marss_utils, only : normal_quantile
   use r_linalg, only : inverse_matrix
   implicit none
   private
   public :: marss_fisher_i_harvey_linear
   public :: marss_param_cis_harvey_linear

contains

   pure subroutine marss_fisher_i_harvey_linear(template, constraints, theta, fisher, info)
      type(marss_model), intent(in) :: template !! Template model defining observations, fixed values, and constrained block shapes.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta constraints whose beta coordinates define derivatives.
      real(dp), intent(in) :: theta(:) !! Free beta coordinates at which Harvey observed information is evaluated.
      real(dp), allocatable, intent(out) :: fisher(:, :) !! Harvey (1989) observed Fisher information in beta coordinates.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid constraints, diffuse initialization, or filter failure.
      type(marss_model) :: model
      type(marss_kf_result) :: kf
      real(dp), allocatable :: finv(:, :)
      real(dp), allocatable :: dinnov(:, :)
      real(dp), allocatable :: df(:, :, :)
      real(dp), allocatable :: dx_prev(:, :)
      real(dp), allocatable :: dx_pred(:, :)
      real(dp), allocatable :: dx_filt(:, :)
      real(dp), allocatable :: dp_prev(:, :, :)
      real(dp), allocatable :: dp_pred(:, :, :)
      real(dp), allocatable :: dp_filt(:, :, :)
      real(dp) :: bt(size(template%b, 1), size(template%b, 2))
      real(dp) :: db(size(template%b, 1), size(template%b, 2))
      real(dp) :: du(size(template%u))
      real(dp) :: dq(size(template%q, 1), size(template%q, 2))
      real(dp) :: zt(size(template%z, 1), size(template%z, 2))
      real(dp) :: dz(size(template%z, 1), size(template%z, 2))
      real(dp) :: da(size(template%a))
      real(dp) :: rt(size(template%r, 1), size(template%r, 2))
      real(dp) :: dr(size(template%r, 1), size(template%r, 2))
      real(dp) :: dx0(size(template%x0))
      real(dp) :: dv0(size(template%v0, 1), size(template%v0, 2))
      real(dp) :: v0eff(size(template%v0, 1), size(template%v0, 2))
      real(dp) :: zwork(size(template%z, 1), size(template%z, 2))
      real(dp) :: dzwork(size(template%z, 1), size(template%z, 2))
      real(dp) :: rwork(size(template%r, 1), size(template%r, 2))
      real(dp) :: drwork(size(template%r, 1), size(template%r, 2))
      real(dp) :: dawork(size(template%a))
      real(dp) :: fmat(size(template%r, 1), size(template%r, 2))
      real(dp) :: pzt(size(template%b, 1), size(template%z, 1))
      real(dp) :: work_mn(size(template%b, 1), size(template%z, 1))
      real(dp) :: work_nm(size(template%z, 1), size(template%b, 1))
      real(dp) :: work_mm(size(template%b, 1), size(template%b, 2))
      real(dp) :: innov_work(size(template%z, 1))
      logical :: missing(size(template%z, 1))
      real(dp) :: value
      integer :: i
      integer :: j
      integer :: m
      integer :: n
      integer :: p
      integer :: rank
      integer :: stat
      integer :: t
      integer :: tt

      p = marss_constraints_parameter_count(constraints)
      if (size(theta) /= p) then
         allocate(fisher(0, 0))
         info = 1
         return
      end if
      call marss_apply_constraints(template, constraints, theta, model, stat)
      if (stat /= 0) then
         allocate(fisher(0, 0))
         info = 10 + stat
         return
      end if
      if (model%diffuse) then
         allocate(fisher(0, 0))
         info = 2
         return
      end if
      allocate(fisher(p, p))
      fisher = 0.0_dp
      if (p == 0) then
         info = 0
         return
      end if
      call marss_kfss(model, kf, smoother=.false.)
      if (.not. kf%ok) then
         info = 100 + kf%info
         return
      end if

      m = size(model%b, 1)
      n = size(model%z, 1)
      tt = size(model%y, 2)
      allocate(dinnov(n, p), df(n, n, p))
      allocate(dx_prev(m, p), dx_pred(m, p), dx_filt(m, p))
      allocate(dp_prev(m, m, p), dp_pred(m, m, p), dp_filt(m, m, p))
      dx_prev = 0.0_dp
      dp_prev = 0.0_dp
      v0eff = marss_v0_effective(model)

      do t = 1, tt
         bt = marss_b_at(model, t)
         zt = marss_z_at(model, t)
         rt = marss_r_at(model, t)
         zwork = zt
         rwork = rt
         innov_work = kf%innov(:, t)
         missing = ieee_is_nan(model%y(:, t))
         do i = 1, n
            if (missing(i)) then
               zwork(i, :) = 0.0_dp
               rwork(i, :) = 0.0_dp
               rwork(:, i) = 0.0_dp
               innov_work(i) = 0.0_dp
            end if
         end do
         fmat = matmul(matmul(zwork, kf%p_pred(:, :, t)), transpose(zwork)) + rwork
         fmat = 0.5_dp * (fmat + transpose(fmat))
         call psd_inverse(fmat, finv, rank, stat)
         if (stat /= 0) then
            info = 200 + t
            return
         end if
         pzt = matmul(kf%p_pred(:, :, t), transpose(zwork))

         do i = 1, p
            call constraint_derivatives_at(model, constraints, i, t, db, du, dq, dz, da, dr, dx0, dv0)
            dzwork = dz
            drwork = dr
            dawork = da
            do j = 1, n
               if (missing(j)) then
                  dzwork(j, :) = 0.0_dp
                  drwork(j, :) = 0.0_dp
                  drwork(:, j) = 0.0_dp
                  dawork(j) = 0.0_dp
               end if
            end do

            if (t == 1) then
               if (model%tinitx == 0) then
                  dx_pred(:, i) = matmul(db, model%x0) + matmul(bt, dx0) + du
                  dp_pred(:, :, i) = matmul(matmul(db, v0eff), transpose(bt)) + &
                     matmul(matmul(bt, dv0), transpose(bt)) + &
                     matmul(matmul(bt, v0eff), transpose(db)) + dq
               else
                  dx_pred(:, i) = dx0
                  dp_pred(:, :, i) = dv0
               end if
            else
               dx_pred(:, i) = matmul(db, kf%x_filt(:, t - 1)) + matmul(bt, dx_prev(:, i)) + du
               dp_pred(:, :, i) = matmul(matmul(db, kf%p_filt(:, :, t - 1)), transpose(bt)) + &
                  matmul(matmul(bt, dp_prev(:, :, i)), transpose(bt)) + &
                  matmul(matmul(bt, kf%p_filt(:, :, t - 1)), transpose(db)) + dq
            end if
            dp_pred(:, :, i) = 0.5_dp * (dp_pred(:, :, i) + transpose(dp_pred(:, :, i)))
            dinnov(:, i) = -matmul(zwork, dx_pred(:, i)) - matmul(dzwork, kf%x_pred(:, t)) - dawork
            df(:, :, i) = matmul(matmul(dzwork, kf%p_pred(:, :, t)), transpose(zwork)) + &
               matmul(matmul(zwork, dp_pred(:, :, i)), transpose(zwork)) + &
               matmul(matmul(zwork, kf%p_pred(:, :, t)), transpose(dzwork)) + drwork
            df(:, :, i) = 0.5_dp * (df(:, :, i) + transpose(df(:, :, i)))

            work_mn = matmul(dp_pred(:, :, i), transpose(zwork))
            dx_filt(:, i) = dx_pred(:, i) + matmul(matmul(work_mn, finv), innov_work)
            work_mn = matmul(kf%p_pred(:, :, t), transpose(dzwork))
            dx_filt(:, i) = dx_filt(:, i) + matmul(matmul(work_mn, finv), innov_work)
            dx_filt(:, i) = dx_filt(:, i) - &
               matmul(matmul(matmul(matmul(pzt, finv), df(:, :, i)), finv), innov_work)
            dx_filt(:, i) = dx_filt(:, i) + matmul(matmul(pzt, finv), dinnov(:, i))

            dp_filt(:, :, i) = dp_pred(:, :, i)
            work_mn = matmul(dp_pred(:, :, i), transpose(zwork))
            work_nm = matmul(zwork, kf%p_pred(:, :, t))
            dp_filt(:, :, i) = dp_filt(:, :, i) - matmul(matmul(work_mn, finv), work_nm)
            work_mn = matmul(kf%p_pred(:, :, t), transpose(dzwork))
            dp_filt(:, :, i) = dp_filt(:, :, i) - matmul(matmul(work_mn, finv), work_nm)
            work_mm = matmul(matmul(matmul(matmul(pzt, finv), df(:, :, i)), finv), work_nm)
            dp_filt(:, :, i) = dp_filt(:, :, i) + work_mm
            work_nm = matmul(dzwork, kf%p_pred(:, :, t))
            dp_filt(:, :, i) = dp_filt(:, :, i) - matmul(matmul(pzt, finv), work_nm)
            work_nm = matmul(zwork, dp_pred(:, :, i))
            dp_filt(:, :, i) = dp_filt(:, :, i) - matmul(matmul(pzt, finv), work_nm)
            dp_filt(:, :, i) = 0.5_dp * (dp_filt(:, :, i) + transpose(dp_filt(:, :, i)))
         end do

         do i = 1, p
            do j = i, p
               value = 0.5_dp * trace_product4(finv, df(:, :, i), finv, df(:, :, j)) + &
                  dot_product(dinnov(:, i), matmul(finv, dinnov(:, j)))
               fisher(i, j) = fisher(i, j) + value
               if (i /= j) fisher(j, i) = fisher(i, j)
            end do
         end do
         dx_prev = dx_filt
         dp_prev = dp_filt
      end do
      fisher = 0.5_dp * (fisher + transpose(fisher))
      info = 0
   end subroutine marss_fisher_i_harvey_linear

   pure subroutine marss_param_cis_harvey_linear(template, constraints, theta, alpha, se, lower, upper, info)
      type(marss_model), intent(in) :: template !! Template model defining fixed values and constrained block shapes.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta constraints for the fitted coordinates.
      real(dp), intent(in) :: theta(:) !! Fitted free beta coordinates.
      real(dp), intent(in) :: alpha !! Two-sided significance level, for example 0.05 for 95 percent intervals.
      real(dp), allocatable, intent(out) :: se(:) !! Harvey-information asymptotic standard errors for beta.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower normal-approximation confidence limits for beta.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper normal-approximation confidence limits for beta.
      integer, intent(out) :: info !! Zero on success; nonzero when Harvey information or its inverse fails.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fisher(:, :)
      real(dp) :: zcrit
      integer :: i
      integer :: stat

      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         allocate(se(0), lower(0), upper(0))
         info = 1
         return
      end if
      call marss_fisher_i_harvey_linear(template, constraints, theta, fisher, stat)
      if (stat /= 0) then
         allocate(se(0), lower(0), upper(0))
         info = 10 + stat
         return
      end if
      if (size(theta) == 0) then
         allocate(se(0), lower(0), upper(0))
         info = 0
         return
      end if
      call inverse_matrix(fisher, covariance, stat)
      if (stat /= 0) then
         allocate(se(0), lower(0), upper(0))
         info = 100 + stat
         return
      end if
      allocate(se(size(theta)), lower(size(theta)), upper(size(theta)))
      do i = 1, size(theta)
         se(i) = sqrt(max(covariance(i, i), 0.0_dp))
      end do
      zcrit = normal_quantile(1.0_dp - 0.5_dp * alpha)
      lower = theta - zcrit * se
      upper = theta + zcrit * se
      info = 0
   end subroutine marss_param_cis_harvey_linear

   pure subroutine constraint_derivatives_at(model, constraints, parameter, t, db, du, dq, dz, da, dr, dx0, dv0)
      type(marss_model), intent(in) :: model !! Applied model supplying active classic and marxss block shapes.
      type(marss_constraints), intent(in) :: constraints !! Affine constraint blocks defining exact beta derivatives.
      integer, intent(in) :: parameter !! One-based global beta-coordinate index in documented constraint order.
      integer, intent(in) :: t !! One-based observation time selecting derivatives of time-indexed blocks.
      real(dp), intent(out) :: db(:, :) !! Derivative of the effective B block at time t.
      real(dp), intent(out) :: du(:) !! Derivative of U+C*c(t) at time t.
      real(dp), intent(out) :: dq(:, :) !! Derivative of G*Q*G' or direct Q at time t.
      real(dp), intent(out) :: dz(:, :) !! Derivative of the effective Z block at time t.
      real(dp), intent(out) :: da(:) !! Derivative of A+D*d(t) at time t.
      real(dp), intent(out) :: dr(:, :) !! Derivative of H*R*H' or direct R at time t.
      real(dp), intent(out) :: dx0(:) !! Derivative of the initial-state mean.
      real(dp), intent(out) :: dv0(:, :) !! Derivative of L*V0*L' or direct V0.
      real(dp), allocatable :: dc(:, :)
      real(dp), allocatable :: dd(:, :)
      real(dp), allocatable :: dg(:, :)
      real(dp), allocatable :: dh(:, :)
      real(dp), allocatable :: dl(:, :)
      real(dp), allocatable :: dq_raw(:, :)
      real(dp), allocatable :: dr_raw(:, :)
      real(dp), allocatable :: dv0_raw(:, :)
      real(dp), allocatable :: gt(:, :)
      real(dp), allocatable :: ht(:, :)
      real(dp), allocatable :: lt(:, :)
      real(dp), allocatable :: qraw(:, :)
      real(dp), allocatable :: rraw(:, :)
      real(dp), allocatable :: v0raw(:, :)
      integer :: offset

      db = 0.0_dp
      du = 0.0_dp
      dq = 0.0_dp
      dz = 0.0_dp
      da = 0.0_dp
      dr = 0.0_dp
      dx0 = 0.0_dp
      dv0 = 0.0_dp
      offset = 0
      call matrix_derivative(constraints%b, parameter, offset, t, allocated(model%b_t), db)
      call vector_derivative(constraints%u, parameter, offset, t, allocated(model%u_t), du)

      if (allocated(model%q_noise)) then
         allocate(dq_raw(size(model%q_noise, 1), size(model%q_noise, 2)))
         call matrix_derivative(constraints%q, parameter, offset, t, allocated(model%q_noise_t), dq_raw)
      else
         call matrix_derivative(constraints%q, parameter, offset, t, allocated(model%q_t), dq)
      end if
      call matrix_derivative(constraints%z, parameter, offset, t, allocated(model%z_t), dz)
      call vector_derivative(constraints%a, parameter, offset, t, allocated(model%a_t), da)
      if (allocated(model%r_noise)) then
         allocate(dr_raw(size(model%r_noise, 1), size(model%r_noise, 2)))
         call matrix_derivative(constraints%r, parameter, offset, t, allocated(model%r_noise_t), dr_raw)
      else
         call matrix_derivative(constraints%r, parameter, offset, t, allocated(model%r_t), dr)
      end if
      call vector_derivative(constraints%x0, parameter, offset, 1, .false., dx0)
      if (allocated(model%v0_noise)) then
         allocate(dv0_raw(size(model%v0_noise, 1), size(model%v0_noise, 2)))
         call matrix_derivative(constraints%v0, parameter, offset, 1, .false., dv0_raw)
      else
         call matrix_derivative(constraints%v0, parameter, offset, 1, .false., dv0)
      end if

      if (allocated(model%c_coef)) then
         allocate(dc(size(model%c_coef, 1), size(model%c_coef, 2)))
         call matrix_derivative(constraints%c, parameter, offset, t, allocated(model%c_coef_t), dc)
         du = du + matmul(dc, model%state_covariates(:, t))
      else
         offset = offset + block_parameter_count_local(constraints%c)
      end if
      if (allocated(model%d_coef)) then
         allocate(dd(size(model%d_coef, 1), size(model%d_coef, 2)))
         call matrix_derivative(constraints%d, parameter, offset, t, allocated(model%d_coef_t), dd)
         da = da + matmul(dd, model%obs_covariates(:, t))
      else
         offset = offset + block_parameter_count_local(constraints%d)
      end if

      if (allocated(model%g)) then
         allocate(dg(size(model%g, 1), size(model%g, 2)))
         call matrix_derivative(constraints%g, parameter, offset, t, allocated(model%g_t), dg)
      else
         offset = offset + block_parameter_count_local(constraints%g)
      end if
      if (allocated(model%h)) then
         allocate(dh(size(model%h, 1), size(model%h, 2)))
         call matrix_derivative(constraints%h, parameter, offset, t, allocated(model%h_t), dh)
      else
         offset = offset + block_parameter_count_local(constraints%h)
      end if
      if (allocated(model%l)) then
         allocate(dl(size(model%l, 1), size(model%l, 2)))
         call matrix_derivative(constraints%l, parameter, offset, 1, .false., dl)
      else
         offset = offset + block_parameter_count_local(constraints%l)
      end if

      if (allocated(model%q_noise)) then
         allocate(gt(size(model%g, 1), size(model%g, 2)))
         allocate(qraw(size(model%q_noise, 1), size(model%q_noise, 2)))
         if (allocated(model%g_t)) then
            gt = model%g_t(:, :, t)
         else
            gt = model%g
         end if
         if (allocated(model%q_noise_t)) then
            qraw = model%q_noise_t(:, :, t)
         else
            qraw = model%q_noise
         end if
         dq = matmul(matmul(dg, qraw), transpose(gt)) + &
            matmul(matmul(gt, dq_raw), transpose(gt)) + &
            matmul(matmul(gt, qraw), transpose(dg))
      end if
      if (allocated(model%r_noise)) then
         allocate(ht(size(model%h, 1), size(model%h, 2)))
         allocate(rraw(size(model%r_noise, 1), size(model%r_noise, 2)))
         if (allocated(model%h_t)) then
            ht = model%h_t(:, :, t)
         else
            ht = model%h
         end if
         if (allocated(model%r_noise_t)) then
            rraw = model%r_noise_t(:, :, t)
         else
            rraw = model%r_noise
         end if
         dr = matmul(matmul(dh, rraw), transpose(ht)) + &
            matmul(matmul(ht, dr_raw), transpose(ht)) + &
            matmul(matmul(ht, rraw), transpose(dh))
      end if
      if (allocated(model%v0_noise)) then
         allocate(lt(size(model%l, 1), size(model%l, 2)))
         allocate(v0raw(size(model%v0_noise, 1), size(model%v0_noise, 2)))
         lt = model%l
         v0raw = model%v0_noise
         dv0 = matmul(matmul(dl, v0raw), transpose(lt)) + &
            matmul(matmul(lt, dv0_raw), transpose(lt)) + &
            matmul(matmul(lt, v0raw), transpose(dl))
      end if
   end subroutine constraint_derivatives_at

   pure subroutine matrix_derivative(block, parameter, offset, t, time_varying, derivative)
      type(marss_constraint_block), intent(in) :: block !! Affine matrix block whose design column is inspected.
      integer, intent(in) :: parameter !! One-based global beta-coordinate index.
      integer, intent(inout) :: offset !! Number of beta coordinates preceding this block, advanced on return.
      integer, intent(in) :: t !! One-based time slice used for a time-varying block.
      logical, intent(in) :: time_varying !! True when flattened design rows contain consecutive time slices.
      real(dp), intent(out) :: derivative(:, :) !! Selected matrix derivative, or zero when beta belongs elsewhere.
      integer :: first
      integer :: last
      integer :: local
      integer :: npar
      integer :: slice_size

      derivative = 0.0_dp
      npar = block_parameter_count_local(block)
      if (parameter > offset .and. parameter <= offset + npar) then
         local = parameter - offset
         if (time_varying) then
            slice_size = size(derivative)
            first = (t - 1) * slice_size + 1
            last = t * slice_size
            derivative = reshape(block%design(first:last, local), shape(derivative))
         else
            derivative = reshape(block%design(:, local), shape(derivative))
         end if
      end if
      offset = offset + npar
   end subroutine matrix_derivative

   pure subroutine vector_derivative(block, parameter, offset, t, time_varying, derivative)
      type(marss_constraint_block), intent(in) :: block !! Affine vector block whose design column is inspected.
      integer, intent(in) :: parameter !! One-based global beta-coordinate index.
      integer, intent(inout) :: offset !! Number of beta coordinates preceding this block, advanced on return.
      integer, intent(in) :: t !! One-based time slice used for a time-varying block.
      logical, intent(in) :: time_varying !! True when flattened design rows contain consecutive time slices.
      real(dp), intent(out) :: derivative(:) !! Selected vector derivative, or zero when beta belongs elsewhere.
      integer :: first
      integer :: last
      integer :: local
      integer :: npar
      integer :: slice_size

      derivative = 0.0_dp
      npar = block_parameter_count_local(block)
      if (parameter > offset .and. parameter <= offset + npar) then
         local = parameter - offset
         if (time_varying) then
            slice_size = size(derivative)
            first = (t - 1) * slice_size + 1
            last = t * slice_size
            derivative = block%design(first:last, local)
         else
            derivative = block%design(:, local)
         end if
      end if
      offset = offset + npar
   end subroutine vector_derivative

   pure integer function block_parameter_count_local(block) result(nparam)
      type(marss_constraint_block), intent(in) :: block !! Affine block whose design-column count is requested.

      nparam = 0
      if (allocated(block%design)) nparam = size(block%design, 2)
   end function block_parameter_count_local

   pure function trace_product4(a, b, c, d) result(value)
      real(dp), intent(in) :: a(:, :) !! First square matrix in trace(A B C D).
      real(dp), intent(in) :: b(:, :) !! Second square matrix in trace(A B C D).
      real(dp), intent(in) :: c(:, :) !! Third square matrix in trace(A B C D).
      real(dp), intent(in) :: d(:, :) !! Fourth square matrix in trace(A B C D).
      real(dp) :: value
      real(dp) :: product(size(a, 1), size(a, 2))
      integer :: i

      product = matmul(matmul(matmul(a, b), c), d)
      value = 0.0_dp
      do i = 1, size(product, 1)
         value = value + product(i, i)
      end do
   end function trace_product4

end module marss_harvey_info
