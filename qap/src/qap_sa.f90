module qap_sa_mod
   use qap_kinds, only : dp, i64
   use qap_rng, only : qap_rng_t
   use qap_types, only : qap_control_t, qap_result_t
   use qap_core, only : qap_obj, qap_swap_delta, qap_is_permutation, qap_validate
   implicit none
   private
   public :: qap_sa, qap_solve

contains

   subroutine qap_solve(A, B, result, control, initial_perm)
      real(dp), intent(in) :: A(:,:), B(:,:)
      type(qap_result_t), intent(out) :: result
      type(qap_control_t), intent(in), optional :: control
      integer, intent(in), optional :: initial_perm(:)

      type(qap_control_t) :: ctl
      type(qap_rng_t) :: rng
      type(qap_result_t) :: one
      integer, allocatable :: perm(:)
      integer :: r, n, i

      call qap_validate(A, B)
      n = size(A,1)
      ctl = qap_control_t()
      if (present(control)) ctl = control
      if (ctl%rep < 1) error stop "qap_solve: rep must be positive"
      if (ctl%fiter <= 0.0_dp) error stop "qap_solve: fiter must be positive"
      if (ctl%ft < 0.0_dp .or. ctl%ft >= 1.0_dp) then
         error stop "qap_solve: ft must be in [0,1)"
      end if
      if (ctl%maxsteps < 1) error stop "qap_solve: maxsteps must be positive"
      if (ctl%miter < 0) ctl%miter = 2*n

      allocate(result%permutation(n), perm(n))
      result%permutation = [(i, i=1,n)]
      result%objective = huge(1.0_dp)
      result%best_rep = 0

      call rng%seed(ctl%seed)
      if (ctl%verbose) then
         write(*,'(a)') "Simulated Annealing Heuristic by Burkard and Rendl."
         write(*,'(5x,a,2(1x,a12))') "rep", "best_obj", "current_obj"
      end if

      do r = 1, ctl%rep
         if (present(initial_perm) .and. r == 1) then
            if (size(initial_perm) /= n) error stop "qap_solve: initial_perm has wrong size"
            if (.not. qap_is_permutation(initial_perm)) then
               error stop "qap_solve: initial_perm is not a permutation"
            end if
            perm = initial_perm
         else
            perm = [(i, i=1,size(perm))]
            call rng%shuffle(perm)
         end if

         call qap_sa(A, B, perm, one, ctl%miter, ctl%fiter, ctl%ft, ctl%maxsteps, rng)
         result%attempted_swaps = result%attempted_swaps + one%attempted_swaps
         result%accepted_swaps = result%accepted_swaps + one%accepted_swaps
         result%duplicate_trials = result%duplicate_trials + one%duplicate_trials
         result%cooling_steps = result%cooling_steps + one%cooling_steps
         if (one%objective < result%objective) then
            result%objective = one%objective
            result%permutation = one%permutation
            result%best_rep = r
         end if
         if (ctl%verbose) then
            write(*,'(i8,2(1x,f12.0))') r, result%objective, one%objective
         end if
      end do
   end subroutine qap_solve

   subroutine qap_sa(A, B, perm, result, miter, fiter, ft, maxsteps, rng)
      real(dp), intent(in) :: A(:,:), B(:,:)
      integer, intent(inout) :: perm(:)
      type(qap_result_t), intent(out) :: result
      integer, intent(in) :: miter, maxsteps
      real(dp), intent(in) :: fiter, ft
      type(qap_rng_t), intent(inout) :: rng

      integer :: n, i, i1, i2, m1, step, tmp
      real(dp) :: ia, ib, t1, ol, min_obj, max_obj, delta, p, u
      real(dp) :: best_obj, initial_obj
      logical :: found_best
      integer, allocatable :: best_perm(:)

      n = size(perm)
      allocate(best_perm(n), result%permutation(n))

      ia = sum(A)
      ib = sum(B)
      t1 = ia / real(n*n-n, dp)
      t1 = t1 * ib
      m1 = miter
      ol = qap_obj(A, B, perm)
      initial_obj = ol

      ! qapsim.f initializes its incumbent objective to the rough mean t1 and
      ! does not store the starting permutation immediately. Preserve that
      ! behavior. The initialized best_perm is only a safe fallback for the
      ! pathological case in which no observed state ever reaches t1.
      best_obj = t1
      best_perm = perm
      found_best = .false.

      step = 0
      do
         min_obj = ol
         max_obj = -huge(1.0_dp)

         do i = 1, m1
            result%attempted_swaps = result%attempted_swaps + 1_i64

            u = rng%uniform()
            i1 = int(u * real(n,dp) + 0.5_dp)
            if (i1 < 1) i1 = 1
            if (i1 > n) i1 = n

            u = rng%uniform()
            i2 = int(u * real(n,dp) + 0.5_dp)
            if (i2 < 1) i2 = 1
            if (i2 > n) i2 = n

            if (i1 == i2) then
               result%duplicate_trials = result%duplicate_trials + 1_i64
               call update_observed()
               cycle
            end if

            delta = qap_swap_delta(A, B, perm, i1, i2)
            if (delta > 0.0_dp) then
               if (t1 <= 0.0_dp) then
                  p = 0.0_dp
               else if (delta / t1 > 10.0_dp) then
                  p = 0.0_dp
               else
                  p = exp(-delta / t1)
               end if
               if (rng%uniform() > p) cycle
            end if

            tmp = perm(i1)
            perm(i1) = perm(i2)
            perm(i2) = tmp
            ol = ol + delta
            result%accepted_swaps = result%accepted_swaps + 1_i64
            call update_observed()
         end do

         t1 = t1 * ft
         m1 = int(real(m1,dp) * fiter)
         step = step + 1
         result%cooling_steps = result%cooling_steps + 1_i64
         if (.not. (max_obj > min_obj .and. step < maxsteps)) exit
      end do

      if (.not. found_best) best_obj = initial_obj
      result%permutation = best_perm
      result%objective = best_obj
      result%best_rep = 1

   contains

      subroutine update_observed()
         min_obj = min(min_obj, ol)
         max_obj = max(max_obj, ol)
         if (ol <= best_obj) then
            best_obj = ol
            best_perm = perm
            found_best = .true.
         end if
      end subroutine update_observed

   end subroutine qap_sa

end module qap_sa_mod
