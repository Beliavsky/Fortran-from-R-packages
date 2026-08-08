! SPDX-License-Identifier: MPL-2.0
module trustoptim
   use trustoptim_kinds, only : dp
   use trustoptim_types
   use trustoptim_linalg, only : vecnorm, cholesky_lower, cholesky_solve, scaled_norm
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   public :: dp
   public :: trustoptim_control, trustoptim_result, sparse_symmetric_matrix
   public :: trust_method_sr1, trust_method_bfgs, trust_method_sparse
   public :: trust_success, trust_etolg, trust_emaxiter
   public :: trust_optim, trust_optim_sr1, trust_optim_bfgs, trust_optim_sparse
   public :: objective_fn, gradient_fn, sparse_hessian_fn

   interface trust_optim
      module procedure trust_optim_quasi_dispatch
      module procedure trust_optim_sparse
   end interface trust_optim

   abstract interface
      function objective_fn(x) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function objective_fn

      subroutine gradient_fn(x, g)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: g(:)
      end subroutine gradient_fn

      subroutine sparse_hessian_fn(x, h)
         import dp, sparse_symmetric_matrix
         real(dp), intent(in) :: x(:)
         type(sparse_symmetric_matrix), intent(inout) :: h
      end subroutine sparse_hessian_fn
   end interface

contains

   subroutine trust_optim_quasi_dispatch(x0, fn, gr, method, result, control)
      real(dp), intent(in) :: x0(:)
      procedure(objective_fn) :: fn
      procedure(gradient_fn) :: gr
      character(len=*), intent(in) :: method
      type(trustoptim_result), intent(out) :: result
      type(trustoptim_control), intent(in), optional :: control
      type(trustoptim_control) :: con

      con = trustoptim_control()
      if (present(control)) con = control
      select case (trim(method))
      case ('SR1', 'sr1')
         if (con%preconditioner == 1) con%preconditioner = 0
         call optimize_quasi(x0, fn, gr, trust_method_sr1, con, result)
      case ('BFGS', 'bfgs')
         call optimize_quasi(x0, fn, gr, trust_method_bfgs, con, result)
      case default
         error stop 'trustOptim: quasi method must be SR1 or BFGS'
      end select
   end subroutine trust_optim_quasi_dispatch

   subroutine trust_optim_sr1(x0, fn, gr, result, control)
      real(dp), intent(in) :: x0(:)
      procedure(objective_fn) :: fn
      procedure(gradient_fn) :: gr
      type(trustoptim_result), intent(out) :: result
      type(trustoptim_control), intent(in), optional :: control
      type(trustoptim_control) :: con

      con = trustoptim_control()
      if (present(control)) con = control
      if (con%preconditioner == 1) con%preconditioner = 0
      call optimize_quasi(x0, fn, gr, trust_method_sr1, con, result)
   end subroutine trust_optim_sr1

   subroutine trust_optim_bfgs(x0, fn, gr, result, control)
      real(dp), intent(in) :: x0(:)
      procedure(objective_fn) :: fn
      procedure(gradient_fn) :: gr
      type(trustoptim_result), intent(out) :: result
      type(trustoptim_control), intent(in), optional :: control
      type(trustoptim_control) :: con

      con = trustoptim_control()
      if (present(control)) con = control
      call optimize_quasi(x0, fn, gr, trust_method_bfgs, con, result)
   end subroutine trust_optim_bfgs

   subroutine trust_optim_sparse(x0, fn, gr, hs, result, control)
      real(dp), intent(in) :: x0(:)
      procedure(objective_fn) :: fn
      procedure(gradient_fn) :: gr
      procedure(sparse_hessian_fn) :: hs
      type(trustoptim_result), intent(out) :: result
      type(trustoptim_control), intent(in), optional :: control
      type(trustoptim_control) :: con

      con = trustoptim_control()
      if (present(control)) con = control
      call optimize_sparse(x0, fn, gr, hs, con, result)
   end subroutine trust_optim_sparse

   subroutine optimize_quasi(x0, fn, gr, method, con, result)
      real(dp), intent(in) :: x0(:)
      procedure(objective_fn) :: fn
      procedure(gradient_fn) :: gr
      integer, intent(in) :: method
      type(trustoptim_control), intent(in) :: con
      type(trustoptim_result), intent(out) :: result
      real(dp), allocatable :: x(:), g(:), b(:,:), l(:,:), s(:), y(:), gtry(:)
      real(dp), allocatable :: bs(:), w(:)
      real(dp) :: f, ftry, rad, nrm_g, ared, pred, ap, gs, sbs, denom, crit
      integer :: n, iter, step_status, cg_iters
      character(len=64) :: cg_reason
      logical :: accepted

      n = size(x0)
      call validate_control(con, n)
      allocate(x(n), g(n), b(n,n), l(n,n), s(n), y(n), gtry(n), bs(n), w(n))
      x = x0
      f = fn(x)
      call gr(x, g)
      call validate_start(f, g, n)
      f = con%function_scale_factor * f
      g = con%function_scale_factor * g
      b = 0.0_dp
      b = identity_matrix(n)
      rad = con%start_trust_radius
      call update_quasi_preconditioner(b, con%preconditioner, l)
      nrm_g = vecnorm(g)
      iter = 0
      result%status = trust_continue
      cg_reason = ''
      cg_iters = 0

      do while (result%status == trust_continue)
         iter = iter + 1
         call solve_trust_cg_dense(g, b, l, rad, con%cg_tol, con%trust_iter, &
                                   s, cg_iters, cg_reason)
         call assess_step_dense(x, f, g, b, l, s, fn, gr, con, rad, &
                                ftry, gtry, ared, pred, ap, step_status, accepted)

         if (accepted) then
            y = con%function_scale_factor * gtry - g
            x = x + s
            f = con%function_scale_factor * ftry
            g = con%function_scale_factor * gtry
            nrm_g = vecnorm(g)
         end if

         if (step_status == trust_contract .or. step_status == trust_failedcg .or. &
             step_status == trust_enegmove) rad = rad * con%contract_factor
         if (step_status == trust_expand) rad = rad * con%expand_factor

         result%status = step_status
         if (result%status == trust_failedcg .or. result%status == trust_enegmove) then
            result%status = trust_continue
         end if
         if (nrm_g / sqrt(real(n,dp)) <= con%prec) result%status = trust_success
         if (iter >= con%maxit) result%status = trust_emaxiter
         if (rad < con%stop_trust_radius) result%status = trust_etolg

         if (step_status == trust_moved .or. step_status == trust_expand) then
            if (result%status /= trust_success .and. result%status /= trust_emaxiter .and. &
                result%status /= trust_etolg) then
               if (method == trust_method_sr1) then
                  call update_sr1(b, s, y, w, denom, crit)
               else
                  call update_bfgs(b, s, y, bs, gs, sbs)
               end if
               if (mod(iter, max(1, con%precond_refresh_freq)) == 0) then
                  call update_quasi_preconditioner(b, con%preconditioner, l)
               end if
               result%status = trust_continue
            end if
         end if
         if (step_status == trust_contract .and. result%status == trust_contract) then
            result%status = trust_continue
         end if
      end do

      allocate(result%solution(n), result%gradient(n))
      result%solution = x
      result%fval = fn(x)
      call gr(x, result%gradient)
      result%iterations = iter
      result%trust_radius = rad
      result%method = method
      result%hessian_update_method = method
      result%last_cg_iterations = cg_iters
      result%last_cg_reason = cg_reason
   end subroutine optimize_quasi

   subroutine optimize_sparse(x0, fn, gr, hs, con, result)
      real(dp), intent(in) :: x0(:)
      procedure(objective_fn) :: fn
      procedure(gradient_fn) :: gr
      procedure(sparse_hessian_fn) :: hs
      type(trustoptim_control), intent(in) :: con
      type(trustoptim_result), intent(out) :: result
      real(dp), allocatable :: x(:), g(:), l(:,:), s(:), gtry(:)
      real(dp) :: f, ftry, rad, nrm_g, ared, pred, ap
      integer :: n, iter, step_status, cg_iters
      character(len=64) :: cg_reason
      logical :: accepted
      type(sparse_symmetric_matrix) :: h

      n = size(x0)
      call validate_control(con, n)
      allocate(x(n), g(n), l(n,n), s(n), gtry(n))
      x = x0
      f = fn(x)
      call gr(x, g)
      call validate_start(f, g, n)
      f = con%function_scale_factor * f
      g = con%function_scale_factor * g
      call hs(x, h)
      call validate_sparse_hessian(h, n)
      h%val = con%function_scale_factor * h%val
      rad = con%start_trust_radius
      call update_sparse_preconditioner(h, con%preconditioner, l)
      nrm_g = vecnorm(g)
      iter = 0
      result%status = trust_continue
      cg_reason = ''
      cg_iters = 0

      do while (result%status == trust_continue)
         iter = iter + 1
         call solve_trust_cg_sparse(g, h, l, rad, con%cg_tol, con%trust_iter, &
                                    s, cg_iters, cg_reason)
         call assess_step_sparse(x, f, g, h, l, s, fn, gr, con, rad, &
                                 ftry, gtry, ared, pred, ap, step_status, accepted)

         if (accepted) then
            x = x + s
            f = con%function_scale_factor * ftry
            g = con%function_scale_factor * gtry
            nrm_g = vecnorm(g)
         end if

         if (step_status == trust_contract .or. step_status == trust_failedcg .or. &
             step_status == trust_enegmove) rad = rad * con%contract_factor
         if (step_status == trust_expand) rad = rad * con%expand_factor

         result%status = step_status
         if (result%status == trust_failedcg .or. result%status == trust_enegmove) then
            result%status = trust_continue
         end if
         if (nrm_g / sqrt(real(n,dp)) <= con%prec) result%status = trust_success
         if (iter >= con%maxit) result%status = trust_emaxiter
         if (rad < con%stop_trust_radius) result%status = trust_etolg

         if (step_status == trust_moved .or. step_status == trust_expand) then
            if (result%status /= trust_success .and. result%status /= trust_emaxiter .and. &
                result%status /= trust_etolg) then
               call hs(x, h)
               call validate_sparse_hessian(h, n)
               h%val = con%function_scale_factor * h%val
               if (mod(iter, max(1, con%precond_refresh_freq)) == 0) then
                  call update_sparse_preconditioner(h, con%preconditioner, l)
               end if
               result%status = trust_continue
            end if
         end if
         if (step_status == trust_contract .and. result%status == trust_contract) then
            result%status = trust_continue
         end if
      end do

      allocate(result%solution(n), result%gradient(n))
      result%solution = x
      result%fval = fn(x)
      call gr(x, result%gradient)
      result%iterations = iter
      result%trust_radius = rad
      result%method = trust_method_sparse
      result%hessian_update_method = 0
      result%last_cg_iterations = cg_iters
      result%last_cg_reason = cg_reason
      call hs(x, result%hessian)
      result%nnz = result%hessian%nnz
   end subroutine optimize_sparse

   subroutine solve_trust_cg_dense(g, b, l, rad, tol, max_cg, p, num_iters, reason)
      real(dp), intent(in) :: g(:), b(:,:), l(:,:), rad, tol
      integer, intent(in) :: max_cg
      real(dp), intent(out) :: p(:)
      integer, intent(out) :: num_iters
      character(len=*), intent(out) :: reason
      real(dp), allocatable :: z(:), zold(:), r(:), y(:), d(:), bd(:)
      real(dp) :: norm_r, norm_g, dBd, a, beta, ry, ryold, tau
      integer :: j, n

      n = size(g)
      allocate(z(n), zold(n), r(n), y(n), d(n), bd(n))
      z = 0.0_dp
      r = -g
      norm_r = scaled_norm(l, r)
      norm_g = scaled_norm(l, g)
      call cholesky_solve(l, r, y)
      d = y
      p = 0.0_dp
      reason = 'Exceeded max CG iters'
      num_iters = max_cg

      do j = 1, max_cg
         bd = matmul(b, d)
         dBd = dot_product(d, bd)
         if (dBd <= 0.0_dp) then
            tau = find_tau(l, z, d, rad)
            p = z + tau * d
            num_iters = j
            reason = 'Negative curvature'
            return
         end if
         ry = dot_product(r, y)
         a = ry / dBd
         zold = z
         z = z + a * d
         if (scaled_norm(l, z) >= rad) then
            tau = find_tau(l, zold, d, rad)
            p = zold + tau * d
            num_iters = j
            reason = 'Intersect TR bound'
            return
         end if
         r = r - a * bd
         norm_r = scaled_norm(l, r)
         if (norm_g <= tiny(1.0_dp) .or. norm_r / norm_g < tol) then
            p = z
            num_iters = j
            reason = 'Reached tolerance'
            return
         end if
         ryold = ry
         call cholesky_solve(l, r, y)
         ry = dot_product(r, y)
         beta = ry / ryold
         d = y + beta * d
      end do
      p = z
   end subroutine solve_trust_cg_dense

   subroutine solve_trust_cg_sparse(g, h, l, rad, tol, max_cg, p, num_iters, reason)
      real(dp), intent(in) :: g(:), l(:,:), rad, tol
      type(sparse_symmetric_matrix), intent(in) :: h
      integer, intent(in) :: max_cg
      real(dp), intent(out) :: p(:)
      integer, intent(out) :: num_iters
      character(len=*), intent(out) :: reason
      real(dp), allocatable :: z(:), zold(:), r(:), y(:), d(:), hd(:)
      real(dp) :: norm_r, norm_g, dHd, a, beta, ry, ryold, tau
      integer :: j, n

      n = size(g)
      allocate(z(n), zold(n), r(n), y(n), d(n), hd(n))
      z = 0.0_dp
      r = -g
      norm_r = scaled_norm(l, r)
      norm_g = scaled_norm(l, g)
      call cholesky_solve(l, r, y)
      d = y
      p = 0.0_dp
      reason = 'Exceeded max CG iters'
      num_iters = max_cg

      do j = 1, max_cg
         call h%matvec(d, hd)
         dHd = dot_product(d, hd)
         if (dHd <= 0.0_dp) then
            tau = find_tau(l, z, d, rad)
            p = z + tau * d
            num_iters = j
            reason = 'Negative curvature'
            return
         end if
         ry = dot_product(r, y)
         a = ry / dHd
         zold = z
         z = z + a * d
         if (scaled_norm(l, z) >= rad) then
            tau = find_tau(l, zold, d, rad)
            p = zold + tau * d
            num_iters = j
            reason = 'Intersect TR bound'
            return
         end if
         r = r - a * hd
         norm_r = scaled_norm(l, r)
         if (norm_g <= tiny(1.0_dp) .or. norm_r / norm_g < tol) then
            p = z
            num_iters = j
            reason = 'Reached tolerance'
            return
         end if
         ryold = ry
         call cholesky_solve(l, r, y)
         ry = dot_product(r, y)
         beta = ry / ryold
         d = y + beta * d
      end do
      p = z
   end subroutine solve_trust_cg_sparse

   subroutine assess_step_dense(x, f, g, b, l, s, fn, gr, con, rad, &
                                ftry, gtry, ared, pred, ap, status, accepted)
      real(dp), intent(in) :: x(:), f, g(:), b(:,:), l(:,:), s(:), rad
      procedure(objective_fn) :: fn
      procedure(gradient_fn) :: gr
      type(trustoptim_control), intent(in) :: con
      real(dp), intent(out) :: ftry, gtry(:), ared, pred, ap
      integer, intent(out) :: status
      logical, intent(out) :: accepted
      real(dp) :: bs(size(s)), gs, sbs, stepnorm

      accepted = .false.
      status = trust_unknown
      stepnorm = scaled_norm(l, s)
      if (.not. ieee_is_finite(stepnorm)) then
         status = trust_failedcg
         return
      end if
      ftry = fn(x + s)
      if (.not. ieee_is_finite(ftry)) then
         status = trust_failedcg
         return
      end if
      ftry = con%function_scale_factor * ftry
      ared = f - ftry
      gs = dot_product(g, s)
      bs = matmul(b, s)
      sbs = dot_product(s, bs)
      pred = -(gs + 0.5_dp * sbs)
      if (pred < 0.0_dp) then
         status = trust_enegmove
         return
      end if
      if (abs(pred) <= tiny(1.0_dp)) then
         status = trust_failedcg
         return
      end if
      ap = ared / pred
      if (ap > con%contract_threshold) then
         call gr(x + s, gtry)
         if (.not. all(ieee_is_finite(gtry))) then
            status = trust_failedcg
            return
         end if
         ftry = ftry / con%function_scale_factor
         accepted = .true.
         if (ap > con%expand_threshold_ap .and. &
             stepnorm >= con%expand_threshold_radius * rad) then
            status = trust_expand
         else
            status = trust_moved
         end if
      else
         ftry = ftry / con%function_scale_factor
         status = trust_contract
      end if
   end subroutine assess_step_dense

   subroutine assess_step_sparse(x, f, g, h, l, s, fn, gr, con, rad, &
                                 ftry, gtry, ared, pred, ap, status, accepted)
      real(dp), intent(in) :: x(:), f, g(:), l(:,:), s(:), rad
      type(sparse_symmetric_matrix), intent(in) :: h
      procedure(objective_fn) :: fn
      procedure(gradient_fn) :: gr
      type(trustoptim_control), intent(in) :: con
      real(dp), intent(out) :: ftry, gtry(:), ared, pred, ap
      integer, intent(out) :: status
      logical, intent(out) :: accepted
      real(dp) :: hs(size(s)), gs, shs, stepnorm

      accepted = .false.
      status = trust_unknown
      stepnorm = scaled_norm(l, s)
      if (.not. ieee_is_finite(stepnorm)) then
         status = trust_failedcg
         return
      end if
      ftry = fn(x + s)
      if (.not. ieee_is_finite(ftry)) then
         status = trust_failedcg
         return
      end if
      ftry = con%function_scale_factor * ftry
      ared = f - ftry
      gs = dot_product(g, s)
      call h%matvec(s, hs)
      shs = dot_product(s, hs)
      pred = -(gs + 0.5_dp * shs)
      if (pred < 0.0_dp) then
         status = trust_enegmove
         return
      end if
      if (abs(pred) <= tiny(1.0_dp)) then
         status = trust_failedcg
         return
      end if
      ap = ared / pred
      if (ap > con%contract_threshold) then
         call gr(x + s, gtry)
         if (.not. all(ieee_is_finite(gtry))) then
            status = trust_failedcg
            return
         end if
         ftry = ftry / con%function_scale_factor
         accepted = .true.
         if (ap > con%expand_threshold_ap .and. &
             stepnorm >= con%expand_threshold_radius * rad) then
            status = trust_expand
         else
            status = trust_moved
         end if
      else
         ftry = ftry / con%function_scale_factor
         status = trust_contract
      end if
   end subroutine assess_step_sparse

   subroutine update_sr1(b, s, y, w, denom, crit)
      real(dp), intent(inout) :: b(:,:)
      real(dp), intent(in) :: s(:), y(:)
      real(dp), intent(out) :: w(:), denom, crit
      integer :: i, j, n

      w = y - matmul(b, s)
      denom = dot_product(w, s)
      crit = 1.0e-7_dp * vecnorm(s) * vecnorm(w)
      if (abs(denom) > crit) then
         n = size(s)
         do j = 1, n
            do i = 1, n
               b(i,j) = b(i,j) + w(i) * w(j) / denom
            end do
         end do
      end if
   end subroutine update_sr1

   subroutine update_bfgs(b, s, y, bs, ys, sbs)
      real(dp), intent(inout) :: b(:,:)
      real(dp), intent(in) :: s(:), y(:)
      real(dp), intent(out) :: bs(:), ys, sbs
      integer :: i, j, n

      ys = dot_product(y, s)
      bs = matmul(b, s)
      sbs = dot_product(s, bs)
      if (abs(ys) <= tiny(1.0_dp) .or. abs(sbs) <= tiny(1.0_dp)) return
      n = size(s)
      do j = 1, n
         do i = 1, n
            b(i,j) = b(i,j) - bs(i) * bs(j) / sbs + y(i) * y(j) / ys
         end do
      end do
   end subroutine update_bfgs

   subroutine update_quasi_preconditioner(b, precond_id, l)
      real(dp), intent(in) :: b(:,:)
      integer, intent(in) :: precond_id
      real(dp), intent(out) :: l(:,:)
      logical :: ok
      integer :: n, i

      n = size(b,1)
      if (precond_id == 1) then
         call cholesky_lower(b, l, ok)
         if (ok) return
      end if
      l = 0.0_dp
      do i = 1, n
         l(i,i) = 1.0_dp
      end do
   end subroutine update_quasi_preconditioner

   subroutine update_sparse_preconditioner(h, precond_id, l)
      type(sparse_symmetric_matrix), intent(in) :: h
      integer, intent(in) :: precond_id
      real(dp), intent(out) :: l(:,:)
      real(dp), allocatable :: a(:,:), t(:)
      real(dp) :: beta, bmin, alpha, inc
      logical :: ok
      integer :: n, i, j, attempt

      n = h%n
      if (precond_id /= 1) then
         l = 0.0_dp
         do i = 1, n
            l(i,i) = 1.0_dp
         end do
         return
      end if

      allocate(a(n,n), t(n))
      call h%to_dense(a)
      do j = 1, n
         t(j) = sqrt(sum(a(:,j) * a(:,j)))
         if (t(j) <= tiny(1.0_dp)) t(j) = 1.0_dp
      end do
      do j = 1, n
         do i = 1, n
            a(i,j) = a(i,j) / sqrt(t(i) * t(j))
         end do
      end do
      beta = sqrt(sum(a * a))
      bmin = minval([(a(i,i), i=1,n)])
      if (bmin > 0.0_dp) then
         alpha = 0.0_dp
      else
         alpha = 0.5_dp * beta
      end if
      do attempt = 1, 60
         call cholesky_lower(a, l, ok)
         if (ok) return
         inc = max(2.0_dp * alpha, 0.5_dp * beta) - alpha
         if (inc <= 0.0_dp) inc = max(sqrt(epsilon(1.0_dp)), 0.5_dp * beta)
         alpha = alpha + inc
         do i = 1, n
            a(i,i) = a(i,i) + inc
         end do
      end do
      l = 0.0_dp
      do i = 1, n
         l(i,i) = 1.0_dp
      end do
   end subroutine update_sparse_preconditioner

   pure function find_tau(l, z, d, rad) result(tau)
      real(dp), intent(in) :: l(:,:), z(:), d(:), rad
      real(dp) :: tau
      real(dp) :: wd(size(d)), wz(size(z)), d2, z2, zd, root
      integer :: i, j, n

      n = size(d)
      wd = 0.0_dp
      wz = 0.0_dp
      do i = 1, n
         do j = i, n
            wd(i) = wd(i) + l(j,i) * d(j)
            wz(i) = wz(i) + l(j,i) * z(j)
         end do
      end do
      d2 = dot_product(wd, wd)
      z2 = dot_product(wz, wz)
      zd = dot_product(wd, wz)
      root = max(0.0_dp, zd * zd - d2 * (z2 - rad * rad))
      if (d2 <= tiny(1.0_dp)) then
         tau = 0.0_dp
      else
         tau = (sqrt(root) - zd) / d2
      end if
   end function find_tau

   function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   subroutine validate_start(f, g, n)
      real(dp), intent(in) :: f, g(:)
      integer, intent(in) :: n
      if (n <= 0) error stop 'trustOptim: number of variables must be positive'
      if (.not. ieee_is_finite(f)) error stop 'trustOptim: initial objective is not finite'
      if (size(g) /= n) error stop 'trustOptim: gradient length mismatch'
      if (.not. all(ieee_is_finite(g))) error stop 'trustOptim: initial gradient is not finite'
   end subroutine validate_start

   subroutine validate_control(con, n)
      type(trustoptim_control), intent(in) :: con
      integer, intent(in) :: n
      if (n <= 0) error stop 'trustOptim: number of variables must be positive'
      if (con%start_trust_radius <= 0.0_dp) error stop 'trustOptim: bad start trust radius'
      if (con%stop_trust_radius <= 0.0_dp) error stop 'trustOptim: bad stop trust radius'
      if (con%cg_tol <= 0.0_dp) error stop 'trustOptim: bad CG tolerance'
      if (con%prec <= 0.0_dp) error stop 'trustOptim: bad precision'
      if (con%maxit < 0) error stop 'trustOptim: maxit must be nonnegative'
      if (con%contract_threshold < 0.0_dp) error stop 'trustOptim: bad contract threshold'
      if (con%expand_threshold_ap < 0.0_dp) error stop 'trustOptim: bad expansion threshold'
      if (con%expand_threshold_radius < 0.0_dp) error stop 'trustOptim: bad radius expansion threshold'
      if (con%expand_factor < 0.0_dp) error stop 'trustOptim: bad expansion factor'
      if (abs(con%function_scale_factor) <= tiny(1.0_dp) .or. &
          .not. ieee_is_finite(con%function_scale_factor)) then
         error stop 'trustOptim: invalid function scale factor'
      end if
      if (con%trust_iter <= 0) error stop 'trustOptim: trust_iter must be positive'
   end subroutine validate_control

   subroutine validate_sparse_hessian(h, n)
      type(sparse_symmetric_matrix), intent(in) :: h
      integer, intent(in) :: n
      integer :: k
      if (h%n /= n) error stop 'trustOptim: sparse Hessian dimension mismatch'
      if (h%nnz < 0) error stop 'trustOptim: invalid sparse Hessian nnz'
      if (h%nnz == 0) return
      if (.not. allocated(h%row) .or. .not. allocated(h%col) .or. .not. allocated(h%val)) then
         error stop 'trustOptim: incomplete sparse Hessian storage'
      end if
      if (size(h%row) /= h%nnz .or. size(h%col) /= h%nnz .or. size(h%val) /= h%nnz) then
         error stop 'trustOptim: sparse Hessian storage length mismatch'
      end if
      do k = 1, h%nnz
         if (h%row(k) < h%col(k)) error stop 'trustOptim: sparse Hessian must store lower triangle'
         if (h%row(k) < 1 .or. h%row(k) > n .or. h%col(k) < 1 .or. h%col(k) > n) then
            error stop 'trustOptim: sparse Hessian index out of range'
         end if
      end do
      if (.not. all(ieee_is_finite(h%val))) error stop 'trustOptim: sparse Hessian is not finite'
   end subroutine validate_sparse_hessian

end module trustoptim
