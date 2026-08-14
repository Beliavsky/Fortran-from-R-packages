module rebayes_kw
   use rebayes_kinds, only : dp
   use rebayes_math, only : normalize_prob, safe_log
   implicit none
   private
   public :: kw_control, kw_result, kw_fit, kw_primal, kw_dual

   type :: kw_control
      integer :: max_iter = 20000
      real(dp) :: tol = 1.0e-9_dp
      real(dp) :: weight_floor = 1.0e-15_dp
      integer :: vertex_every = 25
   end type kw_control

   type :: kw_result
      real(dp), allocatable :: f(:)
      real(dp), allocatable :: g(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: kkt_gap = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 1
   end type kw_result
contains
   subroutine kw_fit(a, d, weights, result, control, f0)
      real(dp), intent(in) :: a(:,:), d(:), weights(:)
      type(kw_result), intent(out) :: result
      type(kw_control), intent(in), optional :: control
      real(dp), intent(in), optional :: f0(:)
      type(kw_control) :: ctl
      real(dp), allocatable :: b(:,:), w(:), f(:), g(:), score(:), fnew(:)
      real(dp) :: ll, llnew, gap, alpha
      integer :: n, m, it, jbest
      n = size(a,1); m = size(a,2)
      if (size(d) /= m .or. size(weights) /= n) error stop "kw_fit: dimension mismatch"
      ctl = kw_control(); if (present(control)) ctl = control
      allocate(b(n,m),w(n),f(m),g(n),score(m),fnew(m))
      b = a*spread(d,1,n)
      if (any(b < 0.0_dp)) error stop "kw_fit: A*d must be nonnegative"
      w = max(weights,0.0_dp); call normalize_prob(w)
      if (present(f0)) then
         if (size(f0) /= m) error stop "kw_fit: bad f0 size"
         f = max(f0,ctl%weight_floor); call normalize_prob(f)
      else
         f = 1.0_dp/real(m,dp)
      end if
      g = matmul(b,f)
      if (any(g <= 0.0_dp)) then
         result%status = 2
         allocate(result%f(m),result%g(n)); result%f=f; result%g=g
         return
      end if
      ll = sum(w*safe_log(g))
      do it = 1, ctl%max_iter
         score = matmul(transpose(b),w/g)
         gap = maxval(score)-1.0_dp
         if (gap <= ctl%tol) exit
         fnew = max(f*score,ctl%weight_floor)
         call normalize_prob(fnew)
         g = matmul(b,fnew)
         llnew = sum(w*safe_log(g))
         if (llnew + 100.0_dp*epsilon(1.0_dp) < ll) then
            fnew = 0.5_dp*(f+fnew); call normalize_prob(fnew)
            g = matmul(b,fnew); llnew = sum(w*safe_log(g))
         end if
         f = fnew; ll = llnew
         if (ctl%vertex_every > 0 .and. mod(it,ctl%vertex_every) == 0) then
            score = matmul(transpose(b),w/g)
            jbest = maxloc(score,dim=1)
            if (score(jbest)-1.0_dp > 10.0_dp*ctl%tol) then
               call best_vertex_step(g,b(:,jbest),w,alpha)
               if (alpha > 0.0_dp) then
                  f = (1.0_dp-alpha)*f
                  f(jbest) = f(jbest)+alpha
                  g = matmul(b,f)
                  ll = sum(w*safe_log(g))
               end if
            end if
         end if
      end do
      score = matmul(transpose(b),w/g)
      gap = maxval(score)-1.0_dp
      allocate(result%f(m),result%g(n))
      result%f = f; result%g = g
      result%loglik = sum(w*safe_log(g))
      result%kkt_gap = gap
      result%iterations = min(it,ctl%max_iter)
      result%status = merge(0,1,gap <= max(10.0_dp*ctl%tol,1.0e-7_dp))
   end subroutine kw_fit

   subroutine best_vertex_step(g, bj, w, alpha)
      real(dp), intent(in) :: g(:), bj(:), w(:)
      real(dp), intent(out) :: alpha
      real(dp) :: lo, hi, mid, der, delta(size(g))
      integer :: k
      delta = bj-g
      der = sum(w*delta/g)
      if (der <= 0.0_dp) then
         alpha = 0.0_dp; return
      end if
      lo = 0.0_dp; hi = 1.0_dp-1.0e-12_dp
      der = sum(w*delta/(g+hi*delta))
      if (der >= 0.0_dp) then
         alpha = hi; return
      end if
      do k = 1, 60
         mid = 0.5_dp*(lo+hi)
         der = sum(w*delta/(g+mid*delta))
         if (der > 0.0_dp) then
            lo = mid
         else
            hi = mid
         end if
      end do
      alpha = 0.5_dp*(lo+hi)
   end subroutine best_vertex_step

   subroutine kw_dual(a,d,w,result,control)
      real(dp), intent(in) :: a(:,:), d(:), w(:)
      type(kw_result), intent(out) :: result
      type(kw_control), intent(in), optional :: control
      call kw_fit(a,d,w,result,control)
   end subroutine kw_dual

   subroutine kw_primal(a,d,w,result,control)
      real(dp), intent(in) :: a(:,:), d(:), w(:)
      type(kw_result), intent(out) :: result
      type(kw_control), intent(in), optional :: control
      call kw_fit(a,d,w,result,control)
   end subroutine kw_primal
end module rebayes_kw
