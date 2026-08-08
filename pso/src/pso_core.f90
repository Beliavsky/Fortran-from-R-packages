! SPDX-License-Identifier: LGPL-3.0-only
module pso_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use pso_kinds, only : dp
   use pso_types, only : pso_control, pso_result, pso_objective, pso_gradient, &
      pso_spso2007, pso_spso2011, pso_hybrid_on, pso_hybrid_improved
   use pso_random, only : random_box_matrix, random_uniform_vector, &
      random_r_sphere_like, random_r_spheres_like, random_permutation
   use pso_lbfgsb, only : bounded_lbfgs_grad, bounded_lbfgs_numgrad
   implicit none
   private
   public :: psoptim

   interface psoptim
      module procedure psoptim_vector_bounds
      module procedure psoptim_scalar_bounds
   end interface psoptim

contains

   subroutine psoptim_vector_bounds(par, fn, lower, upper, result, control, gr)
      real(dp), intent(in) :: par(:)
      procedure(pso_objective) :: fn
      real(dp), intent(in) :: lower(:), upper(:)
      type(pso_result), intent(out) :: result
      type(pso_control), intent(in), optional :: control
      procedure(pso_gradient), optional :: gr

      type(pso_control) :: con
      integer :: n, swarm, iter, feval, restarts, stagnate, ibest, i, j
      integer :: idx, local_fev, count_other
      integer, allocatable :: order(:), informant(:)
      real(dp) :: prob, diameter, vmax, reltol, error, old_error
      real(dp) :: inertia, progress, nrm, local_f
      real(dp) :: cp2, cp3, cg3, cpg3
      real(dp), allocatable :: x(:,:), v(:,:), pbest(:,:), fx(:), fp(:)
      real(dp), allocatable :: temp(:), sphere(:), rvec(:), xlocal(:)
      real(dp), allocatable :: tmat(:,:), smat(:,:), rmat(:,:)
      logical, allocatable :: links(:,:), improved(:)
      logical :: init_links, any_start, improved_global

      con = pso_control()
      if (present(control)) con = control
      call validate_inputs(par, lower, upper, con)

      n = size(par)
      if (con%swarm_size > 0) then
         swarm = con%swarm_size
      else if (con%pso_type == pso_spso2007) then
         swarm = floor(10.0_dp + 2.0_dp * sqrt(real(n, dp)))
      else
         swarm = 40
      end if
      if (swarm < 1) error stop "psoptim: swarm_size must be positive"

      if (con%informant_p >= 0.0_dp) then
         prob = con%informant_p
      else
         prob = 1.0_dp - (1.0_dp - 1.0_dp / real(swarm, dp)) ** con%k
      end if
      if (prob < 0.0_dp .or. prob > 1.0_dp) &
         error stop "psoptim: informant_p must lie in [0,1]"

      if (con%diameter >= 0.0_dp) then
         diameter = con%diameter
      else
         diameter = sqrt(sum((upper - lower) ** 2))
      end if
      if (con%v_max >= 0.0_dp) then
         vmax = con%v_max * diameter
      else
         vmax = -1.0_dp
      end if
      reltol = con%reltol * diameter

      cp2 = con%c_p / 2.0_dp
      cp3 = con%c_p / 3.0_dp
      cg3 = con%c_g / 3.0_dp
      cpg3 = cp3 + cg3

      allocate(x(n,swarm), v(n,swarm), pbest(n,swarm), fx(swarm), fp(swarm))
      allocate(temp(n), sphere(n), rvec(n), xlocal(n), order(swarm))
      allocate(informant(swarm), links(swarm,swarm), improved(swarm))

      call random_box_matrix(x, lower, upper)
      any_start = all(.not. ieee_is_nan(par)) .and. all(par >= lower) .and. all(par <= upper)
      if (any_start) x(:,1) = par

      if (con%pso_type == pso_spso2007) then
         allocate(tmat(n,swarm))
         call random_box_matrix(tmat, lower, upper)
         v = (tmat - x) / 2.0_dp
         deallocate(tmat)
      else
         call random_number(v)
         do j = 1, swarm
            v(:,j) = lower - x(:,j) + v(:,j) * (upper - lower)
         end do
      end if
      call clamp_velocity_matrix(v, vmax)

      do i = 1, swarm
         fx(i) = fn(x(:,i)) / con%fnscale
      end do
      feval = swarm
      pbest = x
      fp = fx
      ibest = minloc(fp, dim=1)
      error = fp(ibest)
      init_links = .true.
      iter = 1
      restarts = 0
      stagnate = 0

      if (con%trace) then
         call print_settings(con, swarm, prob, diameter, vmax)
         if (con%report == 1) then
            call print_progress(iter, error, -1.0_dp, .false.)
            if (con%trace_stats) call append_trace(result, iter, error, fx, x)
         end if
      end if

      do while (iter < con%maxit .and. feval < con%maxf .and. &
                error > con%abstol .and. restarts < con%max_restart .and. &
                stagnate < con%maxit_stagnate)
         iter = iter + 1

         if (prob < 1.0_dp .and. init_links) call make_links(links, prob)

         if (.not. con%vectorize) then
            if (con%rand_order) then
               call random_permutation(order)
            else
               do i = 1, swarm
                  order(i) = i
               end do
            end if

            do idx = 1, swarm
               i = order(idx)
               if (prob >= 1.0_dp) then
                  j = ibest
               else
                  j = best_informant(links(:,i), fp)
               end if

               progress = max(safe_ratio(iter, con%maxit), safe_ratio(feval, con%maxf))
               inertia = con%w0 + (con%w1 - con%w0) * progress
               v(:,i) = inertia * v(:,i)

               if (con%pso_type == pso_spso2007) then
                  call random_uniform_vector(rvec, 0.0_dp, con%c_p)
                  v(:,i) = v(:,i) + rvec * (pbest(:,i) - x(:,i))
                  if (i /= j) then
                     call random_uniform_vector(rvec, 0.0_dp, con%c_g)
                     v(:,i) = v(:,i) + rvec * (pbest(:,j) - x(:,i))
                  end if
               else
                  if (i /= j) then
                     temp = cp3 * pbest(:,i) + cg3 * pbest(:,j) - cpg3 * x(:,i)
                  else
                     temp = cp2 * (pbest(:,i) - x(:,i))
                  end if
                  call random_r_sphere_like(sphere, sqrt(sum(temp * temp)))
                  v(:,i) = v(:,i) + temp + sphere
               end if

               call clamp_velocity(v(:,i), vmax)
               x(:,i) = x(:,i) + v(:,i)
               call enforce_bounds(x(:,i), v(:,i), lower, upper)

               if (con%hybrid == pso_hybrid_on) then
                  xlocal = x(:,i)
                  call local_refine(fn, gr, xlocal, lower, upper, con, local_f, local_fev)
                  v(:,i) = v(:,i) + xlocal - x(:,i)
                  x(:,i) = xlocal
                  fx(i) = local_f
                  feval = feval + local_fev
               else
                  fx(i) = fn(x(:,i)) / con%fnscale
                  feval = feval + 1
               end if

               if (fx(i) < fp(i)) then
                  pbest(:,i) = x(:,i)
                  fp(i) = fx(i)
                  if (fp(i) < fp(ibest)) then
                     ibest = i
                     if (con%hybrid == pso_hybrid_improved) then
                        xlocal = x(:,i)
                        call local_refine(fn, gr, xlocal, lower, upper, con, local_f, local_fev)
                        v(:,i) = v(:,i) + xlocal - x(:,i)
                        x(:,i) = xlocal
                        pbest(:,i) = xlocal
                        fx(i) = local_f
                        fp(i) = local_f
                        feval = feval + local_fev
                     end if
                  end if
               end if
               if (feval >= con%maxf) exit
            end do
         else
            if (prob >= 1.0_dp) then
               informant = ibest
            else
               do i = 1, swarm
                  informant(i) = best_informant(links(:,i), fp)
               end do
            end if

            progress = max(safe_ratio(iter, con%maxit), safe_ratio(feval, con%maxf))
            inertia = con%w0 + (con%w1 - con%w0) * progress
            v = inertia * v

            if (con%pso_type == pso_spso2007) then
               allocate(rmat(n,swarm))
               call random_number(rmat)
               rmat = con%c_p * rmat
               v = v + rmat * (pbest - x)
               deallocate(rmat)

               count_other = count(informant /= [(i, i=1,swarm)])
               if (count_other > 0) then
                  allocate(rmat(n,count_other))
                  call random_number(rmat)
                  ! The upstream vectorized R branch uses c.p here, not c.g.
                  rmat = con%c_p * rmat
                  j = 0
                  do i = 1, swarm
                     if (informant(i) /= i) then
                        j = j + 1
                        v(:,i) = v(:,i) + rmat(:,j) * &
                           (pbest(:,informant(i)) - x(:,i))
                     end if
                  end do
                  deallocate(rmat)
               end if
            else
               allocate(tmat(n,swarm), smat(n,swarm))
               do i = 1, swarm
                  if (informant(i) == i) then
                     tmat(:,i) = cp2 * (pbest(:,i) - x(:,i))
                  else
                     tmat(:,i) = cp3 * pbest(:,i) + &
                        cg3 * pbest(:,informant(i)) - cpg3 * x(:,i)
                  end if
               end do
               do i = 1, swarm
                  fx(i) = sqrt(sum(tmat(:,i) * tmat(:,i)))
               end do
               call random_r_spheres_like(smat, fx)
               v = v + tmat + smat
               deallocate(tmat, smat)
            end if

            call clamp_velocity_matrix(v, vmax)
            x = x + v
            do i = 1, swarm
               call enforce_bounds(x(:,i), v(:,i), lower, upper)
            end do

            if (con%hybrid == pso_hybrid_on) then
               do i = 1, swarm
                  xlocal = x(:,i)
                  call local_refine(fn, gr, xlocal, lower, upper, con, local_f, local_fev)
                  v(:,i) = v(:,i) + xlocal - x(:,i)
                  x(:,i) = xlocal
                  fx(i) = local_f
                  feval = feval + local_fev
               end do
            else
               do i = 1, swarm
                  fx(i) = fn(x(:,i)) / con%fnscale
               end do
               feval = feval + swarm
            end if

            improved = fx < fp
            if (any(improved)) then
               do i = 1, swarm
                  if (improved(i)) then
                     pbest(:,i) = x(:,i)
                     fp(i) = fx(i)
                  end if
               end do
               ibest = minloc(fp, dim=1)
               improved_global = improved(ibest)
               if (improved_global .and. con%hybrid == pso_hybrid_improved) then
                  i = ibest
                  xlocal = x(:,i)
                  call local_refine(fn, gr, xlocal, lower, upper, con, local_f, local_fev)
                  v(:,i) = v(:,i) + xlocal - x(:,i)
                  x(:,i) = xlocal
                  pbest(:,i) = xlocal
                  fx(i) = local_f
                  fp(i) = local_f
                  feval = feval + local_fev
               end if
            end if
         end if

         old_error = error
         if (abs(reltol) > tiny(1.0_dp)) then
            nrm = swarm_diameter_from_best(x, pbest(:,ibest))
            if (nrm < reltol) then
               call random_box_matrix(x, lower, upper)
               allocate(tmat(n,swarm))
               call random_box_matrix(tmat, lower, upper)
               v = (tmat - x) / 2.0_dp
               deallocate(tmat)
               call clamp_velocity_matrix(v, vmax)
               restarts = restarts + 1
               if (con%trace) write(*,'(a,i0,a)') "It ", iter, ": restarting"
            end if
         else
            nrm = -1.0_dp
         end if

         init_links = .not. (fp(ibest) < old_error .or. fp(ibest) > old_error)
         if (init_links) then
            stagnate = stagnate + 1
         else
            stagnate = 0
         end if
         error = fp(ibest)

         if (con%trace .and. modulo(iter, con%report) == 0) then
            call print_progress(iter, error, nrm, abs(reltol) > tiny(1.0_dp))
            if (con%trace_stats) call append_trace(result, iter, error, fx, x)
         end if
      end do

      allocate(result%par(n))
      result%par = pbest(:,ibest)
      result%value = fp(ibest)
      result%function_evaluations = feval
      result%iterations = iter
      result%restarts = restarts

      if (error <= con%abstol) then
         result%convergence = 0
         result%message = "Converged"
      else if (feval >= con%maxf) then
         result%convergence = 1
         result%message = "Maximal number of function evaluations reached"
      else if (iter >= con%maxit) then
         result%convergence = 2
         result%message = "Maximal number of iterations reached"
      else if (restarts >= con%max_restart) then
         result%convergence = 3
         result%message = "Maximal number of restarts reached"
      else
         result%convergence = 4
         result%message = "Maximal number of iterations without improvement reached"
      end if
      if (con%trace) write(*,'(a)') result%message
   end subroutine psoptim_vector_bounds

   subroutine psoptim_scalar_bounds(par, fn, lower, upper, result, control, gr)
      real(dp), intent(in) :: par(:)
      procedure(pso_objective) :: fn
      real(dp), intent(in) :: lower, upper
      type(pso_result), intent(out) :: result
      type(pso_control), intent(in), optional :: control
      procedure(pso_gradient), optional :: gr
      real(dp), allocatable :: lower_v(:), upper_v(:)

      allocate(lower_v(size(par)), upper_v(size(par)))
      lower_v = lower
      upper_v = upper
      if (present(control)) then
         if (present(gr)) then
            call psoptim_vector_bounds(par, fn, lower_v, upper_v, result, control, gr)
         else
            call psoptim_vector_bounds(par, fn, lower_v, upper_v, result, control)
         end if
      else
         if (present(gr)) then
            call psoptim_vector_bounds(par, fn, lower_v, upper_v, result, gr=gr)
         else
            call psoptim_vector_bounds(par, fn, lower_v, upper_v, result)
         end if
      end if
   end subroutine psoptim_scalar_bounds

   subroutine validate_inputs(par, lower, upper, con)
      real(dp), intent(in) :: par(:), lower(:), upper(:)
      type(pso_control), intent(in) :: con

      if (size(par) < 1) error stop "psoptim: par must be nonempty"
      if (size(lower) /= size(par) .or. size(upper) /= size(par)) &
         error stop "psoptim: bounds must have the same size as par"
      if (any(.not. ieee_is_finite(lower)) .or. any(.not. ieee_is_finite(upper))) &
         error stop "psoptim: fixed finite bounds must be provided"
      if (any(lower > upper)) error stop "psoptim: lower exceeds upper"
      if (abs(con%fnscale) <= tiny(1.0_dp)) error stop "psoptim: fnscale must be nonzero"
      if (con%maxit < 1 .or. con%maxf < 1) error stop "psoptim: invalid iteration/evaluation limit"
      if (con%report < 1) error stop "psoptim: report must be positive"
      if (con%pso_type /= pso_spso2007 .and. con%pso_type /= pso_spso2011) &
         error stop "psoptim: invalid pso_type"
      if (con%hybrid < 0 .or. con%hybrid > 2) error stop "psoptim: invalid hybrid mode"
   end subroutine validate_inputs

   subroutine local_refine(fn, gr, x, lower, upper, con, f, nfev)
      procedure(pso_objective) :: fn
      procedure(pso_gradient), optional :: gr
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: lower(:), upper(:)
      type(pso_control), intent(in) :: con
      real(dp), intent(out) :: f
      integer, intent(out) :: nfev

      if (present(gr)) then
         call bounded_lbfgs_grad(fn, gr, x, lower, upper, con%fnscale, &
            con%hybrid_maxit, con%hybrid_memory, con%hybrid_reltol, f, nfev)
      else
         call bounded_lbfgs_numgrad(fn, x, lower, upper, con%fnscale, &
            con%hybrid_maxit, con%hybrid_memory, con%hybrid_reltol, f, nfev)
      end if
   end subroutine local_refine

   subroutine make_links(links, prob)
      logical, intent(out) :: links(:,:)
      real(dp), intent(in) :: prob
      real(dp), allocatable :: u(:,:)
      integer :: i

      allocate(u(size(links,1), size(links,2)))
      call random_number(u)
      links = u <= prob
      do i = 1, size(links,1)
         links(i,i) = .true.
      end do
   end subroutine make_links

   integer function best_informant(mask, fp) result(idx)
      logical, intent(in) :: mask(:)
      real(dp), intent(in) :: fp(:)
      integer :: i
      real(dp) :: best

      idx = 1
      best = huge(1.0_dp)
      do i = 1, size(fp)
         if (mask(i) .and. fp(i) < best) then
            best = fp(i)
            idx = i
         end if
      end do
   end function best_informant

   subroutine clamp_velocity(v, vmax)
      real(dp), intent(inout) :: v(:)
      real(dp), intent(in) :: vmax
      real(dp) :: nrm

      if (vmax < 0.0_dp) return
      nrm = sqrt(sum(v * v))
      if (nrm > vmax .and. nrm > 0.0_dp) v = (vmax / nrm) * v
   end subroutine clamp_velocity

   subroutine clamp_velocity_matrix(v, vmax)
      real(dp), intent(inout) :: v(:,:)
      real(dp), intent(in) :: vmax
      integer :: j
      if (vmax < 0.0_dp) return
      do j = 1, size(v,2)
         call clamp_velocity(v(:,j), vmax)
      end do
   end subroutine clamp_velocity_matrix

   subroutine enforce_bounds(x, v, lower, upper)
      real(dp), intent(inout) :: x(:), v(:)
      real(dp), intent(in) :: lower(:), upper(:)
      integer :: i

      do i = 1, size(x)
         if (x(i) < lower(i)) then
            x(i) = lower(i)
            v(i) = 0.0_dp
         else if (x(i) > upper(i)) then
            x(i) = upper(i)
            v(i) = 0.0_dp
         end if
      end do
   end subroutine enforce_bounds

   real(dp) function swarm_diameter_from_best(x, best) result(d)
      real(dp), intent(in) :: x(:,:), best(:)
      integer :: j
      real(dp) :: d2

      d2 = 0.0_dp
      do j = 1, size(x,2)
         d2 = max(d2, sum((x(:,j) - best) ** 2))
      end do
      d = sqrt(d2)
   end function swarm_diameter_from_best

   real(dp) function safe_ratio(a, b) result(r)
      integer, intent(in) :: a, b
      r = real(a, dp) / real(max(1,b), dp)
   end function safe_ratio

   subroutine append_trace(result, it, error, fx, x)
      type(pso_result), intent(inout) :: result
      integer, intent(in) :: it
      real(dp), intent(in) :: error, fx(:), x(:,:)
      integer, allocatable :: ti(:)
      real(dp), allocatable :: te(:), tf(:,:), tx(:,:,:)
      integer :: old, new

      old = result%ntrace
      new = old + 1
      allocate(ti(new), te(new), tf(size(fx),new), tx(size(x,1),size(x,2),new))
      if (old > 0) then
         ti(1:old) = result%trace_it
         te(1:old) = result%trace_error
         tf(:,1:old) = result%trace_f
         tx(:,:,1:old) = result%trace_x
      end if
      ti(new) = it
      te(new) = error
      tf(:,new) = fx
      tx(:,:,new) = x
      call move_alloc(ti, result%trace_it)
      call move_alloc(te, result%trace_error)
      call move_alloc(tf, result%trace_f)
      call move_alloc(tx, result%trace_x)
      result%ntrace = new
   end subroutine append_trace

   subroutine print_settings(con, swarm, prob, diameter, vmax)
      type(pso_control), intent(in) :: con
      integer, intent(in) :: swarm
      real(dp), intent(in) :: prob, diameter, vmax

      write(*,'(a,i0,a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
         "S=", swarm, ", K=", con%k, ", p=", prob, ", w0=", con%w0, &
         ", w1=", con%w1, ", c.p=", con%c_p, ", c.g=", con%c_g
      write(*,'(a,es12.4,a,es12.4,a,l1,a,i0)') &
         "v.max=", vmax, ", d=", diameter, ", vectorize=", con%vectorize, &
         ", hybrid=", con%hybrid
   end subroutine print_settings

   subroutine print_progress(it, error, diameter, show_diameter)
      integer, intent(in) :: it
      real(dp), intent(in) :: error, diameter
      logical, intent(in) :: show_diameter

      if (show_diameter) then
         write(*,'(a,i0,a,es14.6,a,es14.6)') &
            "It ", it, ": fitness=", error, ", swarm diam.=", diameter
      else
         write(*,'(a,i0,a,es14.6)') "It ", it, ": fitness=", error
      end if
   end subroutine print_progress

end module pso_core
