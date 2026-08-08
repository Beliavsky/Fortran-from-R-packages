! Upstream OOR license declaration: LGPL (version unspecified).
module oor_stosoo
   use oor_kinds, only : dp
   use oor_interfaces, only : vector_objective
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
   implicit none
   private
   public :: stosoo, stosoo_options, stosoo_result, soo_node, soo_level

   type :: stosoo_options
      logical :: stochastic = .true.
      logical :: maximize = .false.
      logical :: keep_tree = .false.
      integer :: verbose = 0
      integer :: k_max = 0
      integer :: h_max = 0
      real(dp) :: delta = -1.0_dp
   end type stosoo_options

   type :: soo_node
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: x_min(:)
      real(dp), allocatable :: x_max(:)
      logical :: leaf = .true.
      logical :: is_new = .false.
      real(dp) :: sums = 0.0_dp
      real(dp) :: b = 0.0_dp
      integer :: k = 0
      real(dp), allocatable :: values(:)
   end type soo_node

   type :: soo_level
      type(soo_node), allocatable :: nodes(:)
   end type soo_level

   type :: stosoo_result
      real(dp), allocatable :: par(:)
      real(dp) :: value = 0.0_dp
      integer :: evaluations = 0
      real(dp), allocatable :: xs(:,:)
      real(dp), allocatable :: ys(:)
      type(soo_level), allocatable :: tree(:)
   end type stosoo_result

contains
   subroutine stosoo(fn, lower, upper, nb_iter, result, options)
      procedure(vector_objective) :: fn
      real(dp), intent(in) :: lower(:), upper(:)
      integer, intent(in) :: nb_iter
      type(stosoo_result), intent(out) :: result
      type(stosoo_options), intent(in), optional :: options
      type(stosoo_options) :: opt
      type(soo_level), allocatable :: tree(:)
      integer :: d, kmax, hmax, h, i, imax, splitd, n, cap, j, final_idx
      logical :: at_least_one, have_final
      real(dp) :: delta, ucbk, fnscale, sampled, b_hi, b_hi_max, v_max, pinf, ninf
      real(dp), allocatable :: xx(:), xg(:), xd(:), newmin(:), newmax(:), finalx(:)
      real(dp), allocatable :: hx(:,:), hy(:)
      real(dp) :: finaly, best_sum
      type(soo_node) :: child

      if (size(lower) /= size(upper)) error stop "stosoo: lower and upper must have the same size"
      d = size(lower)
      if (d < 1) error stop "stosoo: dimension must be positive"
      if (nb_iter < 2) error stop "stosoo: nb_iter must be at least 2"
      if (any(upper <= lower)) error stop "stosoo: every upper bound must exceed its lower bound"

      opt = stosoo_options()
      if (present(options)) opt = options
      if (opt%delta > 0.0_dp) then
         delta = opt%delta
      else
         delta = 1.0_dp / sqrt(real(nb_iter, dp))
      end if
      if (delta <= 0.0_dp) error stop "stosoo: delta must be positive"

      if (opt%stochastic) then
         if (opt%k_max > 0) then
            kmax = opt%k_max
         else
            kmax = max(1, ceiling(real(nb_iter, dp) / log(real(nb_iter, dp))**3))
         end if
         if (opt%h_max > 0) then
            hmax = opt%h_max
         else
            hmax = max(1, ceiling(sqrt(real(nb_iter, dp) / real(kmax, dp))))
         end if
      else
         kmax = 1
         if (opt%h_max > 0) then
            hmax = opt%h_max
         else
            hmax = max(1, ceiling(sqrt(real(nb_iter, dp))))
         end if
      end if
      if (kmax < 1 .or. hmax < 1) error stop "stosoo: k_max and h_max must be positive"

      if (opt%maximize) then
         fnscale = 1.0_dp
      else
         fnscale = -1.0_dp
      end if
      ucbk = log(real(nb_iter, dp)**2 / delta) / 2.0_dp
      pinf = ieee_value(0.0_dp, ieee_positive_inf)
      ninf = ieee_value(0.0_dp, ieee_negative_inf)

      allocate(tree(hmax))
      do h = 1, hmax
         allocate(tree(h)%nodes(0))
      end do

      child = make_node(fill_vector(d, 0.5_dp), fill_vector(d, 0.0_dp), &
                        fill_vector(d, 1.0_dp), .true., .false.)
      sampled = eval_scaled(fn, child%x, lower, upper, fnscale)
      child%k = 1
      child%sums = sampled
      child%b = sampled + sqrt(ucbk)
      allocate(child%values(1))
      child%values(1) = sampled
      call append_node(tree(1), child)

      cap = nb_iter + 4
      allocate(hx(d, cap), hy(cap))
      hx(:,1) = child%x
      hy(1) = sampled
      n = 1
      finaly = ninf
      allocate(finalx(d))
      finalx = child%x
      have_final = .false.
      at_least_one = .true.

      allocate(xx(d), xg(d), xd(d), newmin(d), newmax(d))

      do while (n < nb_iter)
         if (.not. at_least_one) exit
         v_max = ninf
         at_least_one = .false.

         do h = 1, hmax
            if (n > nb_iter) exit
            imax = 0
            b_hi_max = ninf
            do i = 1, size(tree(h)%nodes)
               if (tree(h)%nodes(i)%leaf .and. .not. tree(h)%nodes(i)%is_new) then
                  if (opt%stochastic) then
                     b_hi = tree(h)%nodes(i)%b
                  else
                     b_hi = tree(h)%nodes(i)%sums / real(tree(h)%nodes(i)%k, dp)
                  end if
                  if (b_hi > b_hi_max) then
                     b_hi_max = b_hi
                     imax = i
                  end if
               end if
            end do

            if (imax > 0) then
               if (h + 1 <= hmax) then
                  if (b_hi_max >= v_max) then
                     at_least_one = .true.
                     xx = tree(h)%nodes(imax)%x
                     if (tree(h)%nodes(imax)%k < kmax) then
                        sampled = eval_scaled(fn, xx, lower, upper, fnscale)
                        call store_history(hx, hy, n, xx, sampled)
                        if (sampled > finaly) then
                           finalx = xx
                           finaly = sampled
                           have_final = .true.
                        end if
                        call append_value(tree(h)%nodes(imax)%values, sampled)
                        tree(h)%nodes(imax)%sums = tree(h)%nodes(imax)%sums + sampled
                        tree(h)%nodes(imax)%k = tree(h)%nodes(imax)%k + 1
                        tree(h)%nodes(imax)%b = tree(h)%nodes(imax)%sums / &
                           real(tree(h)%nodes(imax)%k, dp) + &
                           sqrt(ucbk / real(tree(h)%nodes(imax)%k, dp))
                     else
                        tree(h)%nodes(imax)%leaf = .false.
                        splitd = max_range_dimension(tree(h)%nodes(imax)%x_min, &
                                                     tree(h)%nodes(imax)%x_max)
                        xg = xx
                        xg(splitd) = (5.0_dp * tree(h)%nodes(imax)%x_min(splitd) + &
                                      tree(h)%nodes(imax)%x_max(splitd)) / 6.0_dp
                        xd = xx
                        xd(splitd) = (tree(h)%nodes(imax)%x_min(splitd) + &
                                      5.0_dp * tree(h)%nodes(imax)%x_max(splitd)) / 6.0_dp

                        newmin = tree(h)%nodes(imax)%x_min
                        newmax = tree(h)%nodes(imax)%x_max
                        newmax(splitd) = (2.0_dp * tree(h)%nodes(imax)%x_min(splitd) + &
                                          tree(h)%nodes(imax)%x_max(splitd)) / 3.0_dp
                        child = make_node(xg, newmin, newmax, .true., .true.)
                        sampled = eval_scaled(fn, xg, lower, upper, fnscale)
                        call store_history(hx, hy, n, xg, sampled)
                        if (sampled > finaly) then
                           finalx = xg
                           finaly = sampled
                           have_final = .true.
                        end if
                        child%k = 1
                        child%sums = sampled
                        child%b = sampled + sqrt(ucbk)
                        allocate(child%values(1))
                        child%values(1) = sampled
                        call append_node(tree(h + 1), child)

                        newmin = tree(h)%nodes(imax)%x_min
                        newmax = tree(h)%nodes(imax)%x_max
                        newmin(splitd) = (tree(h)%nodes(imax)%x_min(splitd) + &
                                          2.0_dp * tree(h)%nodes(imax)%x_max(splitd)) / 3.0_dp
                        child = make_node(xd, newmin, newmax, .true., .true.)
                        sampled = eval_scaled(fn, xd, lower, upper, fnscale)
                        call store_history(hx, hy, n, xd, sampled)
                        if (sampled > finaly) then
                           finalx = xd
                           finaly = sampled
                           have_final = .true.
                        end if
                        child%k = 1
                        child%sums = sampled
                        child%b = sampled + sqrt(ucbk)
                        allocate(child%values(1))
                        child%values(1) = sampled
                        call append_node(tree(h + 1), child)

                        newmin = tree(h)%nodes(imax)%x_min
                        newmax = tree(h)%nodes(imax)%x_max
                        newmin(splitd) = (2.0_dp * tree(h)%nodes(imax)%x_min(splitd) + &
                                          tree(h)%nodes(imax)%x_max(splitd)) / 3.0_dp
                        newmax(splitd) = (tree(h)%nodes(imax)%x_min(splitd) + &
                                          2.0_dp * tree(h)%nodes(imax)%x_max(splitd)) / 3.0_dp
                        child = make_node(xx, newmin, newmax, .true., .true.)
                        child%k = tree(h)%nodes(imax)%k
                        child%sums = tree(h)%nodes(imax)%sums
                        child%b = tree(h)%nodes(imax)%b
                        if (allocated(tree(h)%nodes(imax)%values)) then
                           allocate(child%values(size(tree(h)%nodes(imax)%values)))
                           child%values = tree(h)%nodes(imax)%values
                        end if
                        call append_node(tree(h + 1), child)
                        v_max = b_hi_max
                     end if
                  end if
               end if
            end if
         end do

         do h = 1, hmax
            do i = 1, size(tree(h)%nodes)
               tree(h)%nodes(i)%is_new = .false.
            end do
         end do
      end do

      if (opt%stochastic) then
         do h = hmax, 1, -1
            final_idx = 0
            best_sum = ninf
            do i = 1, size(tree(h)%nodes)
               if (.not. tree(h)%nodes(i)%leaf) then
                  if (tree(h)%nodes(i)%sums > best_sum) then
                     best_sum = tree(h)%nodes(i)%sums
                     final_idx = i
                  end if
               end if
            end do
            if (final_idx > 0) then
               finalx = tree(h)%nodes(final_idx)%x
               finaly = tree(h)%nodes(final_idx)%sums / real(tree(h)%nodes(final_idx)%k, dp)
               have_final = .true.
               exit
            end if
         end do
      end if

      if (.not. have_final) then
         finalx = tree(1)%nodes(1)%x
         finaly = tree(1)%nodes(1)%sums / real(tree(1)%nodes(1)%k, dp)
      end if

      allocate(result%par(d))
      result%par = finalx * (upper - lower) + lower
      result%value = fnscale * finaly
      result%evaluations = n
      allocate(result%xs(d, n), result%ys(n))
      do j = 1, n
         result%xs(:,j) = hx(:,j) * (upper - lower) + lower
         result%ys(j) = fnscale * hy(j)
      end do
      if (opt%keep_tree) then
         allocate(result%tree(hmax))
         result%tree = tree
      end if
   end subroutine stosoo

   function make_node(x, xmin, xmax, leaf, is_new) result(node)
      real(dp), intent(in) :: x(:), xmin(:), xmax(:)
      logical, intent(in) :: leaf, is_new
      type(soo_node) :: node
      allocate(node%x(size(x)), node%x_min(size(xmin)), node%x_max(size(xmax)))
      node%x = x
      node%x_min = xmin
      node%x_max = xmax
      node%leaf = leaf
      node%is_new = is_new
   end function make_node

   subroutine append_node(level, node)
      type(soo_level), intent(inout) :: level
      type(soo_node), intent(in) :: node
      type(soo_node), allocatable :: tmp(:)
      integer :: nold
      nold = size(level%nodes)
      allocate(tmp(nold + 1))
      if (nold > 0) tmp(1:nold) = level%nodes
      tmp(nold + 1) = node
      call move_alloc(tmp, level%nodes)
   end subroutine append_node

   subroutine append_value(a, value)
      real(dp), allocatable, intent(inout) :: a(:)
      real(dp), intent(in) :: value
      real(dp), allocatable :: tmp(:)
      integer :: nold
      if (.not. allocated(a)) then
         allocate(a(1))
         a(1) = value
         return
      end if
      nold = size(a)
      allocate(tmp(nold + 1))
      if (nold > 0) tmp(1:nold) = a
      tmp(nold + 1) = value
      call move_alloc(tmp, a)
   end subroutine append_value

   subroutine store_history(hx, hy, n, x, y)
      real(dp), intent(inout) :: hx(:,:), hy(:)
      integer, intent(inout) :: n
      real(dp), intent(in) :: x(:), y
      if (n >= size(hy)) error stop "stosoo: internal history capacity exceeded"
      n = n + 1
      hx(:,n) = x
      hy(n) = y
   end subroutine store_history

   function eval_scaled(fn, xnorm, lower, upper, fnscale) result(y)
      procedure(vector_objective) :: fn
      real(dp), intent(in) :: xnorm(:), lower(:), upper(:), fnscale
      real(dp) :: y
      real(dp), allocatable :: x(:)
      allocate(x(size(xnorm)))
      x = xnorm * (upper - lower) + lower
      y = fn(x) / fnscale
   end function eval_scaled

   function fill_vector(n, value) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: value
      real(dp), allocatable :: x(:)
      allocate(x(n))
      x = value
   end function fill_vector

   function max_range_dimension(xmin, xmax) result(idx)
      real(dp), intent(in) :: xmin(:), xmax(:)
      integer :: idx, i
      real(dp) :: best, width
      idx = 1
      best = xmax(1) - xmin(1)
      do i = 2, size(xmin)
         width = xmax(i) - xmin(i)
         if (width > best) then
            best = width
            idx = i
         end if
      end do
   end function max_range_dimension
end module oor_stosoo
