! SPDX-License-Identifier: GPL-3.0-only
module garchito_optimizer
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchito_kinds, only : dp
   use garchito_callbacks, only : objective_callback, projection_callback
   use garchito_types, only : garchito_control, garchito_success, &
      garchito_max_iterations, garchito_numerical_failure
   implicit none
   private

   public :: bounded_nelder_mead

contains

   subroutine bounded_nelder_mead(fn, project, data, x0, lower, upper, control, &
                                  xbest, fbest, status, iterations, evaluations)
      procedure(objective_callback) :: fn
      procedure(projection_callback) :: project
      class(*), intent(in) :: data
      real(dp), intent(in) :: x0(:), lower(:), upper(:)
      type(garchito_control), intent(in) :: control
      real(dp), intent(out) :: xbest(:)
      real(dp), intent(out) :: fbest
      integer, intent(out) :: status, iterations, evaluations

      real(dp), allocatable :: simplex(:, :), values(:), centroid(:)
      real(dp), allocatable :: xr(:), xe(:), xc(:), trial(:), step(:)
      real(dp) :: fr, fe, fc, ftrial, fspread, xspread, scale
      integer :: n, i, j
      logical :: improved

      n = size(x0)
      allocate(simplex(n, n + 1), values(n + 1), centroid(n))
      allocate(xr(n), xe(n), xc(n), trial(n), step(n))
      evaluations = 0
      iterations = 0
      status = garchito_max_iterations

      simplex(:, 1) = max(lower, min(upper, x0))
      call project(simplex(:, 1), data)
      call evaluate(simplex(:, 1), values(1))

      do j = 1, n
         simplex(:, j + 1) = simplex(:, 1)
         scale = max(control%simplex_scale * max(abs(x0(j)), 1.0e-8_dp), &
                     1.0e-6_dp * max(upper(j) - lower(j), 1.0_dp))
         simplex(j, j + 1) = min(upper(j), simplex(j, j + 1) + scale)
         call project(simplex(:, j + 1), data)
         if (maxval(abs(simplex(:, j + 1) - simplex(:, 1))) <= epsilon(1.0_dp)) then
            simplex(j, j + 1) = max(lower(j), simplex(j, j + 1) - 2.0_dp*scale)
            call project(simplex(:, j + 1), data)
         end if
         call evaluate(simplex(:, j + 1), values(j + 1))
      end do

      do iterations = 1, control%max_iterations
         call sort_simplex(simplex, values)
         fspread = maxval(abs(values - values(1)))
         xspread = 0.0_dp
         do j = 2, n + 1
            xspread = max(xspread, maxval(abs(simplex(:, j) - simplex(:, 1))))
         end do
         if (fspread <= control%tolerance * (1.0_dp + abs(values(1))) .and. &
             xspread <= sqrt(control%tolerance) * &
             (1.0_dp + maxval(abs(simplex(:, 1))))) then
            status = garchito_success
            exit
         end if
         if (evaluations >= control%max_evaluations) exit

         centroid = sum(simplex(:, 1:n), dim=2) / real(n, dp)
         xr = centroid + (centroid - simplex(:, n + 1))
         xr = max(lower, min(upper, xr))
         call project(xr, data)
         call evaluate(xr, fr)

         if (fr < values(1)) then
            xe = centroid + 2.0_dp * (xr - centroid)
            xe = max(lower, min(upper, xe))
            call project(xe, data)
            call evaluate(xe, fe)
            if (fe < fr) then
               simplex(:, n + 1) = xe
               values(n + 1) = fe
            else
               simplex(:, n + 1) = xr
               values(n + 1) = fr
            end if
         else if (fr < values(n)) then
            simplex(:, n + 1) = xr
            values(n + 1) = fr
         else
            if (fr < values(n + 1)) then
               xc = centroid + 0.5_dp * (xr - centroid)
            else
               xc = centroid + 0.5_dp * (simplex(:, n + 1) - centroid)
            end if
            xc = max(lower, min(upper, xc))
            call project(xc, data)
            call evaluate(xc, fc)
            if (fc < min(fr, values(n + 1))) then
               simplex(:, n + 1) = xc
               values(n + 1) = fc
            else
               do j = 2, n + 1
                  simplex(:, j) = simplex(:, 1) + &
                     0.5_dp * (simplex(:, j) - simplex(:, 1))
                  simplex(:, j) = max(lower, min(upper, simplex(:, j)))
                  call project(simplex(:, j), data)
                  call evaluate(simplex(:, j), values(j))
               end do
            end if
         end if

         if (control%trace > 0) then
            if (mod(iterations, control%trace) == 0) then
               call sort_simplex(simplex, values)
               write(*, '(a,i0,a,es16.8)') 'iteration ', iterations, &
                  ': objective = ', values(1)
            end if
         end if
      end do

      iterations = min(iterations, control%max_iterations)
      call sort_simplex(simplex, values)
      xbest = simplex(:, 1)
      fbest = values(1)

      ! A bounded coordinate polish is useful when projected simplex vertices
      ! become nearly collinear at the stationarity boundary.
      step = max(1.0e-7_dp * max(upper - lower, 1.0_dp), &
                 1.0e-5_dp * max(abs(xbest), 1.0_dp))
      do i = 1, 80
         improved = .false.
         do j = 1, n
            trial = xbest
            trial(j) = min(upper(j), trial(j) + step(j))
            call project(trial, data)
            call evaluate(trial, ftrial)
            if (ftrial < fbest) then
               xbest = trial
               fbest = ftrial
               improved = .true.
               cycle
            end if
            trial = xbest
            trial(j) = max(lower(j), trial(j) - step(j))
            call project(trial, data)
            call evaluate(trial, ftrial)
            if (ftrial < fbest) then
               xbest = trial
               fbest = ftrial
               improved = .true.
            end if
         end do
         if (.not. improved) step = 0.5_dp * step
         if (maxval(step) <= control%tolerance * &
             (1.0_dp + maxval(abs(xbest)))) exit
         if (evaluations >= control%max_evaluations) exit
      end do

      if (.not. ieee_is_finite(fbest)) status = garchito_numerical_failure

   contains

      subroutine evaluate(x, value)
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: value
         call fn(x, value, data)
         evaluations = evaluations + 1
         if (.not. ieee_is_finite(value)) value = huge(1.0_dp) / 16.0_dp
      end subroutine evaluate

   end subroutine bounded_nelder_mead

   subroutine sort_simplex(simplex, values)
      real(dp), intent(inout) :: simplex(:, :), values(:)
      real(dp), allocatable :: column(:)
      real(dp) :: key
      integer :: i, j

      allocate(column(size(simplex, 1)))
      do i = 2, size(values)
         key = values(i)
         column = simplex(:, i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j + 1) = values(j)
            simplex(:, j + 1) = simplex(:, j)
            j = j - 1
         end do
         values(j + 1) = key
         simplex(:, j + 1) = column
      end do
   end subroutine sort_simplex

end module garchito_optimizer
