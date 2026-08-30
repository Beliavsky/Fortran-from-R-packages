module polyclip_boolean
   use polyclip_kinds, only: dp, i8
   use polyclip_types, only: clip_intersection, clip_union, clip_difference, clip_xor, &
      fill_evenodd
   use polyclip_integer, only: int_path, int_set, append_int_path
   use polyclip_geometry, only: fill_contains, segment_intersections, simplify_collinear_path, int_path_area
   implicit none
   private

   type :: split_edge
      integer(i8) :: ax = 0_i8, ay = 0_i8, bx = 0_i8, by = 0_i8
      integer(i8), allocatable :: px(:), py(:)
      integer :: npt = 0
   end type split_edge

   type :: boundary_seg
      integer(i8) :: ax = 0_i8, ay = 0_i8, bx = 0_i8, by = 0_i8
   end type boundary_seg

   type :: seg_list
      type(boundary_seg), allocatable :: v(:)
      integer :: n = 0
   end type seg_list

   public :: clip_closed_int, clip_open_int, simplify_int_set

contains

   subroutine clip_closed_int(a, b, op, fill_a, fill_b, result, ierr)
      type(int_set), intent(in) :: a, b
      integer, intent(in) :: op, fill_a, fill_b
      type(int_set), intent(out) :: result
      integer, intent(out), optional :: ierr
      type(split_edge), allocatable :: edges(:)
      type(seg_list) :: boundary
      integer :: nedge, i, j, k, nout
      integer(i8) :: ix(2), iy(2)
      integer(i8), allocatable :: sx(:), sy(:)
      integer :: ns
      real(dp) :: mx, my, dx, dy, len, delta, lx, ly, rx, ry, den, t1, t2
      logical :: al, ar, bl, br, rl, rr

      if (present(ierr)) ierr = 0
      if (op < clip_intersection .or. op > clip_xor) then
         if (present(ierr)) then
            ierr = 1; allocate(result%path(0)); return
         else
            error stop "polyclip: invalid clipping operation"
         end if
      end if
      call build_closed_edges(a, b, edges, nedge)
      if (nedge == 0) then
         allocate(result%path(0)); return
      end if

      do i = 1, nedge - 1
         do j = i + 1, nedge
            call segment_intersections(edges(i)%ax, edges(i)%ay, edges(i)%bx, edges(i)%by, &
               edges(j)%ax, edges(j)%ay, edges(j)%bx, edges(j)%by, nout, ix, iy)
            do k = 1, nout
               call add_split_point(edges(i), ix(k), iy(k))
               call add_split_point(edges(j), ix(k), iy(k))
            end do
         end do
      end do

      do i = 1, nedge
         call sorted_split_points(edges(i), sx, sy, ns)
         if (ns < 2) cycle
         do j = 1, ns - 1
            if (sx(j) == sx(j + 1) .and. sy(j) == sy(j + 1)) cycle
            dx = real(edges(i)%bx - edges(i)%ax, dp)
            dy = real(edges(i)%by - edges(i)%ay, dp)
            len = hypot(dx, dy)
            if (len <= 0.0_dp) cycle
            den = dx * dx + dy * dy
            t1 = (real(sx(j) - edges(i)%ax, dp) * dx + &
                  real(sy(j) - edges(i)%ay, dp) * dy) / den
            t2 = (real(sx(j + 1) - edges(i)%ax, dp) * dx + &
                  real(sy(j + 1) - edges(i)%ay, dp) * dy) / den
            mx = real(edges(i)%ax, dp) + 0.5_dp * (t1 + t2) * dx
            my = real(edges(i)%ay, dp) + 0.5_dp * (t1 + t2) * dy
            delta = min(1.0e-4_dp, 0.1_dp * len)
            lx = mx - delta * dy / len
            ly = my + delta * dx / len
            rx = mx + delta * dy / len
            ry = my - delta * dx / len
            al = fill_contains(a, lx, ly, fill_a)
            ar = fill_contains(a, rx, ry, fill_a)
            bl = fill_contains(b, lx, ly, fill_b)
            br = fill_contains(b, rx, ry, fill_b)
            rl = bool_value(al, bl, op)
            rr = bool_value(ar, br, op)
            if (rl .eqv. rr) cycle
            if (rl) then
               call append_unique_seg(boundary, sx(j), sy(j), sx(j + 1), sy(j + 1))
            else
               call append_unique_seg(boundary, sx(j + 1), sy(j + 1), sx(j), sy(j))
            end if
         end do
      end do
      call trace_boundaries(boundary, result)
   end subroutine clip_closed_int

   subroutine simplify_int_set(a, fill_type, result, ierr)
      type(int_set), intent(in) :: a
      integer, intent(in) :: fill_type
      type(int_set), intent(out) :: result
      integer, intent(out), optional :: ierr
      type(int_set) :: empty
      integer :: stat
      allocate(empty%path(0))
      call clip_closed_int(a, empty, clip_union, fill_type, fill_evenodd, result, stat)
      if (present(ierr)) ierr = stat
   end subroutine simplify_int_set

   subroutine clip_open_int(a, b, op, fill_b, result, ierr)
      type(int_set), intent(in) :: a, b
      integer, intent(in) :: op, fill_b
      type(int_set), intent(out) :: result
      integer, intent(out), optional :: ierr
      type(split_edge) :: edge
      integer :: ip, iseg, ib, jb, nb, nout, k, ns
      integer(i8) :: ix(2), iy(2)
      integer(i8), allocatable :: sx(:), sy(:), runx(:), runy(:)
      integer :: nrun
      real(dp) :: mx, my, dx, dy, den, t1, t2
      logical :: inside, keep

      if (present(ierr)) ierr = 0
      allocate(result%path(0))
      if (.not. allocated(a%path)) return
      do ip = 1, size(a%path)
         if (a%path(ip)%size() < 2) cycle
         nrun = 0
         allocate(runx(max(2, 4 * a%path(ip)%size() + 16)), runy(max(2, 4 * a%path(ip)%size() + 16)))
         do iseg = 1, a%path(ip)%size() - 1
            call init_split_edge(edge, a%path(ip)%x(iseg), a%path(ip)%y(iseg), &
               a%path(ip)%x(iseg + 1), a%path(ip)%y(iseg + 1), max(4, 2 * count_edges(b) + 4))
            if (allocated(b%path)) then
               do ib = 1, size(b%path)
                  nb = b%path(ib)%size()
                  if (nb < 2) cycle
                  do jb = 1, nb
                     call segment_intersections(edge%ax, edge%ay, edge%bx, edge%by, &
                        b%path(ib)%x(jb), b%path(ib)%y(jb), &
                        b%path(ib)%x(merge(jb + 1, 1, jb < nb)), &
                        b%path(ib)%y(merge(jb + 1, 1, jb < nb)), nout, ix, iy)
                     do k = 1, nout
                        call add_split_point(edge, ix(k), iy(k))
                     end do
                  end do
               end do
            end if
            call sorted_split_points(edge, sx, sy, ns)
            do k = 1, ns - 1
               if (sx(k) == sx(k + 1) .and. sy(k) == sy(k + 1)) cycle
               dx = real(edge%bx - edge%ax, dp)
               dy = real(edge%by - edge%ay, dp)
               den = dx * dx + dy * dy
               if (den <= 0.0_dp) cycle
               t1 = (real(sx(k) - edge%ax, dp) * dx + real(sy(k) - edge%ay, dp) * dy) / den
               t2 = (real(sx(k + 1) - edge%ax, dp) * dx + real(sy(k + 1) - edge%ay, dp) * dy) / den
               mx = real(edge%ax, dp) + 0.5_dp * (t1 + t2) * dx
               my = real(edge%ay, dp) + 0.5_dp * (t1 + t2) * dy
               inside = fill_contains(b, mx, my, fill_b)
               select case (op)
               case (clip_intersection)
                  keep = inside
               case (clip_union, clip_difference, clip_xor)
                  keep = .not. inside
               case default
                  keep = .false.
               end select
               if (keep) then
                  call append_run_point(runx, runy, nrun, sx(k), sy(k))
                  call append_run_point(runx, runy, nrun, sx(k + 1), sy(k + 1))
               else
                  call flush_run(result, runx, runy, nrun)
               end if
            end do
         end do
         call flush_run(result, runx, runy, nrun)
         deallocate(runx, runy)
      end do
   end subroutine clip_open_int

   pure logical function bool_value(a, b, op) result(v)
      logical, intent(in) :: a, b
      integer, intent(in) :: op
      select case (op)
      case (clip_intersection); v = a .and. b
      case (clip_union);        v = a .or. b
      case (clip_difference);   v = a .and. (.not. b)
      case (clip_xor);          v = a .neqv. b
      case default;             v = .false.
      end select
   end function bool_value

   subroutine build_closed_edges(a, b, edges, nedge)
      type(int_set), intent(in) :: a, b
      type(split_edge), allocatable, intent(out) :: edges(:)
      integer, intent(out) :: nedge
      integer :: total, i, j, n, cap
      total = count_edges(a) + count_edges(b)
      nedge = 0
      allocate(edges(total))
      cap = max(4, 2 * total + 4)
      if (allocated(a%path)) then
         do i = 1, size(a%path)
            n = a%path(i)%size(); if (n < 2) cycle
            do j = 1, n
               nedge = nedge + 1
               call init_split_edge(edges(nedge), a%path(i)%x(j), a%path(i)%y(j), &
                  a%path(i)%x(merge(j + 1, 1, j < n)), a%path(i)%y(merge(j + 1, 1, j < n)), cap)
            end do
         end do
      end if
      if (allocated(b%path)) then
         do i = 1, size(b%path)
            n = b%path(i)%size(); if (n < 2) cycle
            do j = 1, n
               nedge = nedge + 1
               call init_split_edge(edges(nedge), b%path(i)%x(j), b%path(i)%y(j), &
                  b%path(i)%x(merge(j + 1, 1, j < n)), b%path(i)%y(merge(j + 1, 1, j < n)), cap)
            end do
         end do
      end if
   end subroutine build_closed_edges

   pure integer function count_edges(a) result(nedge)
      type(int_set), intent(in) :: a
      integer :: i
      nedge = 0
      if (.not. allocated(a%path)) return
      do i = 1, size(a%path)
         if (a%path(i)%size() >= 2) nedge = nedge + a%path(i)%size()
      end do
   end function count_edges

   subroutine init_split_edge(e, ax, ay, bx, by, cap)
      type(split_edge), intent(out) :: e
      integer(i8), intent(in) :: ax, ay, bx, by
      integer, intent(in) :: cap
      e%ax = ax; e%ay = ay; e%bx = bx; e%by = by; e%npt = 2
      allocate(e%px(cap), e%py(cap))
      e%px(1) = ax; e%py(1) = ay
      e%px(2) = bx; e%py(2) = by
   end subroutine init_split_edge

   subroutine add_split_point(e, x, y)
      type(split_edge), intent(inout) :: e
      integer(i8), intent(in) :: x, y
      integer :: i, oldcap, newcap
      integer(i8), allocatable :: tx(:), ty(:)
      do i = 1, e%npt
         if (e%px(i) == x .and. e%py(i) == y) return
      end do
      if (e%npt >= size(e%px)) then
         oldcap = size(e%px); newcap = max(oldcap + 1, 2 * oldcap)
         allocate(tx(newcap), ty(newcap))
         tx(1:oldcap) = e%px; ty(1:oldcap) = e%py
         call move_alloc(tx, e%px); call move_alloc(ty, e%py)
      end if
      e%npt = e%npt + 1
      e%px(e%npt) = x; e%py(e%npt) = y
   end subroutine add_split_point

   subroutine sorted_split_points(e, x, y, n)
      type(split_edge), intent(in) :: e
      integer(i8), allocatable, intent(out) :: x(:), y(:)
      integer, intent(out) :: n
      real(dp), allocatable :: t(:)
      real(dp) :: dx, dy, den, tv
      integer(i8) :: xv, yv
      integer :: i, j
      n = e%npt
      allocate(x(n), y(n), t(n))
      dx = real(e%bx - e%ax, dp); dy = real(e%by - e%ay, dp)
      den = dx * dx + dy * dy
      do i = 1, n
         x(i) = e%px(i); y(i) = e%py(i)
         if (den > 0.0_dp) then
            t(i) = (real(x(i) - e%ax, dp) * dx + real(y(i) - e%ay, dp) * dy) / den
         else
            t(i) = 0.0_dp
         end if
      end do
      do i = 2, n
         tv = t(i); xv = x(i); yv = y(i); j = i - 1
         do while (j >= 1)
            if (t(j) <= tv) exit
            t(j + 1) = t(j); x(j + 1) = x(j); y(j + 1) = y(j); j = j - 1
         end do
         t(j + 1) = tv; x(j + 1) = xv; y(j + 1) = yv
      end do
      deallocate(t)
   end subroutine sorted_split_points

   subroutine append_unique_seg(list, ax, ay, bx, by)
      type(seg_list), intent(inout) :: list
      integer(i8), intent(in) :: ax, ay, bx, by
      type(boundary_seg), allocatable :: tmp(:)
      integer :: i, cap
      if (ax == bx .and. ay == by) return
      do i = 1, list%n
         if (list%v(i)%ax == ax .and. list%v(i)%ay == ay .and. &
             list%v(i)%bx == bx .and. list%v(i)%by == by) return
         if (list%v(i)%ax == bx .and. list%v(i)%ay == by .and. &
             list%v(i)%bx == ax .and. list%v(i)%by == ay) then
            list%v(i) = list%v(list%n)
            list%n = list%n - 1
            return
         end if
      end do
      if (.not. allocated(list%v)) allocate(list%v(64))
      if (list%n == size(list%v)) then
         cap = 2 * size(list%v); allocate(tmp(cap)); tmp(1:list%n) = list%v(1:list%n)
         call move_alloc(tmp, list%v)
      end if
      list%n = list%n + 1
      list%v(list%n)%ax = ax; list%v(list%n)%ay = ay
      list%v(list%n)%bx = bx; list%v(list%n)%by = by
   end subroutine append_unique_seg

   subroutine trace_boundaries(list, result)
      type(seg_list), intent(in) :: list
      type(int_set), intent(out) :: result
      logical, allocatable :: used(:)
      integer(i8), allocatable :: xx(:), yy(:)
      integer :: i, npt, current, nxt, steps
      integer(i8) :: sx, sy, cx, cy
      type(int_path) :: p
      if (list%n == 0) then
         allocate(result%path(0)); return
      end if
      allocate(result%path(0), used(list%n), xx(list%n + 2), yy(list%n + 2))
      used = .false.
      do i = 1, list%n
         if (used(i)) cycle
         npt = 2
         xx(1) = list%v(i)%ax; yy(1) = list%v(i)%ay
         xx(2) = list%v(i)%bx; yy(2) = list%v(i)%by
         sx = xx(1); sy = yy(1); cx = xx(2); cy = yy(2)
         used(i) = .true.; current = i; steps = 0
         do while (.not. (cx == sx .and. cy == sy))
            nxt = choose_next(list, used, current, cx, cy)
            if (nxt == 0) exit
            used(nxt) = .true.
            npt = npt + 1
            if (npt > size(xx)) exit
            xx(npt) = list%v(nxt)%bx; yy(npt) = list%v(nxt)%by
            cx = xx(npt); cy = yy(npt); current = nxt
            steps = steps + 1
            if (steps > list%n + 1) exit
         end do
         if (npt >= 4 .and. xx(npt) == sx .and. yy(npt) == sy) then
            if (allocated(p%x)) deallocate(p%x, p%y)
            allocate(p%x(npt - 1), p%y(npt - 1))
            p%x = xx(1:npt - 1); p%y = yy(1:npt - 1)
            call simplify_collinear_path(p, .true.)
            if (p%size() >= 3 .and. abs(int_path_area(p)) > 0.0_dp) call append_int_path(result, p)
         end if
      end do
   end subroutine trace_boundaries

   integer function choose_next(list, used, current, x, y) result(best)
      type(seg_list), intent(in) :: list
      logical, intent(in) :: used(:)
      integer, intent(in) :: current
      integer(i8), intent(in) :: x, y
      integer :: j
      real(dp) :: dinx, diny, backx, backy, doutx, douty, crossv, dotv, cw, bestcw
      real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)
      best = 0; bestcw = huge(1.0_dp)
      dinx = real(list%v(current)%bx - list%v(current)%ax, dp)
      diny = real(list%v(current)%by - list%v(current)%ay, dp)
      backx = -dinx; backy = -diny
      do j = 1, list%n
         if (used(j)) cycle
         if (list%v(j)%ax /= x .or. list%v(j)%ay /= y) cycle
         doutx = real(list%v(j)%bx - list%v(j)%ax, dp)
         douty = real(list%v(j)%by - list%v(j)%ay, dp)
         crossv = doutx * backy - douty * backx
         dotv = doutx * backx + douty * backy
         cw = atan2(crossv, dotv)
         if (cw < 0.0_dp) cw = cw + twopi
         if (cw < bestcw) then
            bestcw = cw; best = j
         end if
      end do
   end function choose_next

   subroutine append_run_point(x, y, n, px, py)
      integer(i8), allocatable, intent(inout) :: x(:), y(:)
      integer, intent(inout) :: n
      integer(i8), intent(in) :: px, py
      integer(i8), allocatable :: tx(:), ty(:)
      integer :: cap
      if (n > 0) then
         if (x(n) == px .and. y(n) == py) return
      end if
      if (n == size(x)) then
         cap = 2 * size(x); allocate(tx(cap), ty(cap)); tx(1:n) = x; ty(1:n) = y
         call move_alloc(tx, x); call move_alloc(ty, y)
      end if
      n = n + 1; x(n) = px; y(n) = py
   end subroutine append_run_point

   subroutine flush_run(result, x, y, n)
      type(int_set), intent(inout) :: result
      integer(i8), intent(in) :: x(:), y(:)
      integer, intent(inout) :: n
      type(int_path) :: p
      if (n >= 2) then
         allocate(p%x(n), p%y(n))
         p%x = x(n:1:-1)
         p%y = y(n:1:-1)
         call append_int_path(result, p)
      end if
      n = 0
   end subroutine flush_run

end module polyclip_boolean
