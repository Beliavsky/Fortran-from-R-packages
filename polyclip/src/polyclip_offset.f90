module polyclip_offset
   use polyclip_kinds, only: dp, i8
   use polyclip_types, only: join_square, join_round, join_miter, &
      end_closed_polygon, end_closed_line, end_open_butt, end_open_square, end_open_round, &
      fill_positive, fill_negative
   use polyclip_integer, only: int_path, int_set, append_int_path, round_away_i8
   use polyclip_geometry, only: int_orientation, int_path_area, simplify_collinear_path
   use polyclip_boolean, only: simplify_int_set
   implicit none
   private

   type :: path_builder
      integer(i8), allocatable :: x(:), y(:)
      integer :: n = 0
   end type path_builder

   public :: offset_int_set

contains

   subroutine offset_int_set(a, delta, join_type, end_type, miter_limit, arc_tolerance, result, ierr)
      type(int_set), intent(in) :: a
      real(dp), intent(in) :: delta, miter_limit, arc_tolerance
      integer, intent(in) :: join_type, end_type
      type(int_set), intent(out) :: result
      integer, intent(out), optional :: ierr
      type(int_set) :: src, raw, tmp
      integer :: stat, i, largest
      real(dp) :: best

      if (present(ierr)) ierr = 0
      if (join_type < join_square .or. join_type > join_miter) then
         if (present(ierr)) then
            ierr = 1; allocate(result%path(0)); return
         else
            error stop "polyclip: invalid join type"
         end if
      end if
      if (end_type < end_closed_polygon .or. end_type > end_open_round) then
         if (present(ierr)) then
            ierr = 2; allocate(result%path(0)); return
         else
            error stop "polyclip: invalid end type"
         end if
      end if
      src = a
      call fix_orientations(src, end_type)
      call do_offset(src, delta, join_type, end_type, miter_limit, arc_tolerance, raw)
      if (abs(delta) <= tiny(1.0_dp)) then
         result = raw
         return
      end if
      if (delta > 0.0_dp) then
         call simplify_int_set(raw, fill_positive, result, stat)
      else
         if (end_type /= end_closed_polygon) then
            allocate(result%path(0)); return
         end if
         call add_negative_outer(raw, tmp)
         call simplify_int_set(tmp, fill_negative, result, stat)
         if (result%size() > 0) then
            largest = 1; best = abs(int_path_area(result%path(1)))
            do i = 2, result%size()
               if (abs(int_path_area(result%path(i))) > best) then
                  largest = i; best = abs(int_path_area(result%path(i)))
               end if
            end do
            call remove_path(result, largest)
            do i = 1, result%size()
               call reverse_path(result%path(i))
            end do
         end if
      end if
      if (present(ierr)) ierr = stat
   end subroutine offset_int_set

   subroutine fix_orientations(a, end_type)
      type(int_set), intent(inout) :: a
      integer, intent(in) :: end_type
      integer :: i, j, low_path, low_vertex, n
      integer(i8) :: lx, ly
      logical :: have_low
      if (.not. allocated(a%path)) return
      if (end_type == end_closed_polygon) then
         have_low = .false.; low_path = 0; low_vertex = 0
         do i = 1, size(a%path)
            n = a%path(i)%size(); if (n < 3) cycle
            do j = 1, n
               if (.not. have_low) then
                  have_low = .true.; low_path = i; low_vertex = j
                  lx = a%path(i)%x(j); ly = a%path(i)%y(j)
               else if (a%path(i)%y(j) > ly .or. &
                   (a%path(i)%y(j) == ly .and. a%path(i)%x(j) < lx)) then
                  low_path = i; low_vertex = j
                  lx = a%path(i)%x(j); ly = a%path(i)%y(j)
               end if
            end do
         end do
         if (have_low) then
            if (.not. int_orientation(a%path(low_path))) then
               do i = 1, size(a%path)
                  call reverse_path(a%path(i))
               end do
            end if
         end if
      else if (end_type == end_closed_line) then
         do i = 1, size(a%path)
            if (.not. int_orientation(a%path(i))) call reverse_path(a%path(i))
         end do
      end if
   end subroutine fix_orientations

   subroutine do_offset(a, delta, join_type, end_type, miter_limit, arc_tolerance, dest)
      type(int_set), intent(in) :: a
      real(dp), intent(in) :: delta, miter_limit, arc_tolerance
      integer, intent(in) :: join_type, end_type
      type(int_set), intent(out) :: dest
      integer :: ip, n, j, k, steps
      real(dp), allocatable :: nx(:), ny(:)
      real(dp) :: miter_lim2, arc_y, step_count, sin_step, cos_step, steps_per_rad
      real(dp) :: x, y, x2
      type(path_builder) :: out
      type(int_path) :: p

      allocate(dest%path(0))
      if (.not. allocated(a%path)) return
      if (abs(delta) <= tiny(1.0_dp)) then
         if (end_type == end_closed_polygon) dest = a
         return
      end if
      if (delta < 0.0_dp .and. end_type /= end_closed_polygon) return

      if (miter_limit > 2.0_dp) then
         miter_lim2 = 2.0_dp / (miter_limit * miter_limit)
      else
         miter_lim2 = 0.5_dp
      end if
      if (arc_tolerance <= 0.0_dp) then
         arc_y = 0.25_dp
      else if (arc_tolerance > abs(delta) * 0.25_dp) then
         arc_y = abs(delta) * 0.25_dp
      else
         arc_y = arc_tolerance
      end if
      if (arc_y <= 0.0_dp .or. abs(delta) <= tiny(1.0_dp)) then
         step_count = 4.0_dp
      else
         step_count = acos(-1.0_dp) / acos(max(-1.0_dp, min(1.0_dp, 1.0_dp - arc_y / abs(delta))))
      end if
      if (step_count > abs(delta) * acos(-1.0_dp)) step_count = abs(delta) * acos(-1.0_dp)
      step_count = max(1.0_dp, step_count)
      sin_step = sin(2.0_dp * acos(-1.0_dp) / step_count)
      cos_step = cos(2.0_dp * acos(-1.0_dp) / step_count)
      steps_per_rad = step_count / (2.0_dp * acos(-1.0_dp))
      if (delta < 0.0_dp) sin_step = -sin_step

      do ip = 1, size(a%path)
         p = a%path(ip)
         n = p%size()
         if (n == 0) cycle
         call builder_init(out, max(16, 2 * n + 16))
         if (n == 1) then
            if (join_type == join_round) then
               x = 1.0_dp; y = 0.0_dp
               steps = max(1, round_away_int(step_count))
               do j = 1, steps
                  call builder_add(out, round_away_i8(real(p%x(1), dp) + x * delta), &
                     round_away_i8(real(p%y(1), dp) + y * delta))
                  x2 = x
                  x = x * cos_step - sin_step * y
                  y = x2 * sin_step + y * cos_step
               end do
            else
               call builder_add(out, round_away_i8(real(p%x(1), dp) - delta), round_away_i8(real(p%y(1), dp) - delta))
               call builder_add(out, round_away_i8(real(p%x(1), dp) + delta), round_away_i8(real(p%y(1), dp) - delta))
               call builder_add(out, round_away_i8(real(p%x(1), dp) + delta), round_away_i8(real(p%y(1), dp) + delta))
               call builder_add(out, round_away_i8(real(p%x(1), dp) - delta), round_away_i8(real(p%y(1), dp) + delta))
            end if
            call builder_to_path(out, p)
            if (p%size() >= 3) call append_int_path(dest, p)
            cycle
         end if

         allocate(nx(n), ny(n))
         do j = 1, n - 1
            call unit_normal(p%x(j), p%y(j), p%x(j + 1), p%y(j + 1), nx(j), ny(j))
         end do
         if (end_type == end_closed_line .or. end_type == end_closed_polygon) then
            call unit_normal(p%x(n), p%y(n), p%x(1), p%y(1), nx(n), ny(n))
         else
            nx(n) = nx(n - 1); ny(n) = ny(n - 1)
         end if

         select case (end_type)
         case (end_closed_polygon)
            k = n
            do j = 1, n
               call offset_point(out, p, nx, ny, j, k, delta, join_type, miter_lim2, &
                  sin_step, cos_step, steps_per_rad)
            end do
            call builder_to_path(out, p)
            call simplify_collinear_path(p, .true.)
            if (p%size() >= 3) call append_int_path(dest, p)

         case (end_closed_line)
            k = n
            do j = 1, n
               call offset_point(out, p, nx, ny, j, k, delta, join_type, miter_lim2, &
                  sin_step, cos_step, steps_per_rad)
            end do
            call builder_to_path(out, p)
            if (p%size() >= 3) call append_int_path(dest, p)
            call builder_init(out, max(16, 2 * n + 16))
            call reverse_normals_for_closed(nx, ny)
            k = 1
            do j = n, 1, -1
               call offset_point(out, a%path(ip), nx, ny, j, k, delta, join_type, miter_lim2, &
                  sin_step, cos_step, steps_per_rad)
            end do
            call builder_to_path(out, p)
            if (p%size() >= 3) call append_int_path(dest, p)

         case default
            k = 1
            do j = 2, n - 1
               call offset_point(out, p, nx, ny, j, k, delta, join_type, miter_lim2, &
                  sin_step, cos_step, steps_per_rad)
            end do
            call add_end_cap(out, p, nx, ny, n, n - 1, delta, end_type, sin_step, cos_step, steps_per_rad, .true.)
            call reverse_normals_for_open(nx, ny)
            k = n
            do j = n - 1, 2, -1
               call offset_point(out, p, nx, ny, j, k, delta, join_type, miter_lim2, &
                  sin_step, cos_step, steps_per_rad)
            end do
            call add_end_cap(out, p, nx, ny, 1, 2, delta, end_type, sin_step, cos_step, steps_per_rad, .false.)
            call builder_to_path(out, p)
            if (p%size() >= 3) call append_int_path(dest, p)
         end select
         deallocate(nx, ny)
      end do
   end subroutine do_offset

   subroutine offset_point(out, p, nx, ny, j, k, delta, join_type, miter_lim2, sin_step, cos_step, steps_per_rad)
      type(path_builder), intent(inout) :: out
      type(int_path), intent(in) :: p
      real(dp), intent(in) :: nx(:), ny(:), delta, miter_lim2, sin_step, cos_step, steps_per_rad
      integer, intent(in) :: j, join_type
      integer, intent(inout) :: k
      real(dp) :: sin_a, cos_a, r
      sin_a = nx(k) * ny(j) - nx(j) * ny(k)
      if (abs(sin_a * delta) < 1.0_dp) then
         cos_a = nx(k) * nx(j) + ny(j) * ny(k)
         if (cos_a > 0.0_dp) then
            call builder_add(out, round_away_i8(real(p%x(j), dp) + nx(k) * delta), &
               round_away_i8(real(p%y(j), dp) + ny(k) * delta))
            return
         end if
      else
         sin_a = max(-1.0_dp, min(1.0_dp, sin_a))
      end if
      if (sin_a * delta < 0.0_dp) then
         call builder_add(out, round_away_i8(real(p%x(j), dp) + nx(k) * delta), &
            round_away_i8(real(p%y(j), dp) + ny(k) * delta))
         call builder_add(out, p%x(j), p%y(j))
         call builder_add(out, round_away_i8(real(p%x(j), dp) + nx(j) * delta), &
            round_away_i8(real(p%y(j), dp) + ny(j) * delta))
      else
         select case (join_type)
         case (join_miter)
            r = 1.0_dp + nx(j) * nx(k) + ny(j) * ny(k)
            if (r >= miter_lim2) then
               call do_miter(out, p%x(j), p%y(j), nx(k), ny(k), nx(j), ny(j), delta, r)
            else
               call do_square(out, p%x(j), p%y(j), nx(k), ny(k), nx(j), ny(j), delta, sin_a)
            end if
         case (join_square)
            call do_square(out, p%x(j), p%y(j), nx(k), ny(k), nx(j), ny(j), delta, sin_a)
         case (join_round)
            call do_round(out, p%x(j), p%y(j), nx(k), ny(k), nx(j), ny(j), delta, sin_a, &
               sin_step, cos_step, steps_per_rad)
         end select
      end if
      k = j
   end subroutine offset_point

   subroutine do_square(out, px, py, nkx, nky, njx, njy, delta, sin_a)
      type(path_builder), intent(inout) :: out
      integer(i8), intent(in) :: px, py
      real(dp), intent(in) :: nkx, nky, njx, njy, delta, sin_a
      real(dp) :: dx, dotv
      dotv = nkx * njx + nky * njy
      dx = tan(atan2(sin_a, dotv) / 4.0_dp)
      call builder_add(out, round_away_i8(real(px, dp) + delta * (nkx - nky * dx)), &
         round_away_i8(real(py, dp) + delta * (nky + nkx * dx)))
      call builder_add(out, round_away_i8(real(px, dp) + delta * (njx + njy * dx)), &
         round_away_i8(real(py, dp) + delta * (njy - njx * dx)))
   end subroutine do_square

   subroutine do_miter(out, px, py, nkx, nky, njx, njy, delta, r)
      type(path_builder), intent(inout) :: out
      integer(i8), intent(in) :: px, py
      real(dp), intent(in) :: nkx, nky, njx, njy, delta, r
      real(dp) :: q
      q = delta / r
      call builder_add(out, round_away_i8(real(px, dp) + (nkx + njx) * q), &
         round_away_i8(real(py, dp) + (nky + njy) * q))
   end subroutine do_miter

   subroutine do_round(out, px, py, nkx, nky, njx, njy, delta, sin_a, sin_step, cos_step, steps_per_rad)
      type(path_builder), intent(inout) :: out
      integer(i8), intent(in) :: px, py
      real(dp), intent(in) :: nkx, nky, njx, njy, delta, sin_a, sin_step, cos_step, steps_per_rad
      real(dp) :: a, x, y, x2, dotv
      integer :: i, steps
      dotv = nkx * njx + nky * njy
      a = atan2(sin_a, dotv)
      steps = max(round_away_int(steps_per_rad * abs(a)), 1)
      x = nkx; y = nky
      do i = 1, steps
         call builder_add(out, round_away_i8(real(px, dp) + x * delta), &
            round_away_i8(real(py, dp) + y * delta))
         x2 = x
         x = x * cos_step - sin_step * y
         y = x2 * sin_step + y * cos_step
      end do
      call builder_add(out, round_away_i8(real(px, dp) + njx * delta), &
         round_away_i8(real(py, dp) + njy * delta))
   end subroutine do_round

   subroutine add_end_cap(out, p, nx, ny, j, k, delta, end_type, sin_step, cos_step, steps_per_rad, at_end)
      type(path_builder), intent(inout) :: out
      type(int_path), intent(in) :: p
      real(dp), intent(inout) :: nx(:), ny(:)
      integer, intent(in) :: j, k, end_type
      real(dp), intent(in) :: delta, sin_step, cos_step, steps_per_rad
      logical, intent(in) :: at_end
      real(dp) :: sin_a
      if (end_type == end_open_butt) then
         if (at_end) then
            call builder_add(out, round_away_i8(real(p%x(j), dp) + nx(j) * delta), &
               round_away_i8(real(p%y(j), dp) + ny(j) * delta))
            call builder_add(out, round_away_i8(real(p%x(j), dp) - nx(j) * delta), &
               round_away_i8(real(p%y(j), dp) - ny(j) * delta))
         else
            call builder_add(out, round_away_i8(real(p%x(j), dp) - nx(j) * delta), &
               round_away_i8(real(p%y(j), dp) - ny(j) * delta))
            call builder_add(out, round_away_i8(real(p%x(j), dp) + nx(j) * delta), &
               round_away_i8(real(p%y(j), dp) + ny(j) * delta))
         end if
      else
         sin_a = 0.0_dp
         if (at_end) then
            nx(j) = -nx(j); ny(j) = -ny(j)
         end if
         if (end_type == end_open_square) then
            call do_square(out, p%x(j), p%y(j), nx(k), ny(k), nx(j), ny(j), delta, sin_a)
         else
            call do_round(out, p%x(j), p%y(j), nx(k), ny(k), nx(j), ny(j), delta, sin_a, &
               sin_step, cos_step, steps_per_rad)
         end if
      end if
   end subroutine add_end_cap

   subroutine unit_normal(x1, y1, x2, y2, nx, ny)
      integer(i8), intent(in) :: x1, y1, x2, y2
      real(dp), intent(out) :: nx, ny
      real(dp) :: dx, dy, f
      if (x1 == x2 .and. y1 == y2) then
         nx = 0.0_dp; ny = 0.0_dp; return
      end if
      dx = real(x2 - x1, dp); dy = real(y2 - y1, dp)
      f = 1.0_dp / hypot(dx, dy)
      nx = dy * f; ny = -dx * f
   end subroutine unit_normal

   subroutine reverse_normals_for_closed(nx, ny)
      real(dp), intent(inout) :: nx(:), ny(:)
      real(dp) :: oldx, oldy, lastx, lasty
      integer :: j, n
      n = size(nx); lastx = nx(n); lasty = ny(n)
      do j = n, 2, -1
         oldx = nx(j - 1); oldy = ny(j - 1)
         nx(j) = -oldx; ny(j) = -oldy
      end do
      nx(1) = -lastx; ny(1) = -lasty
   end subroutine reverse_normals_for_closed

   subroutine reverse_normals_for_open(nx, ny)
      real(dp), intent(inout) :: nx(:), ny(:)
      integer :: j, n
      n = size(nx)
      do j = n, 2, -1
         nx(j) = -nx(j - 1); ny(j) = -ny(j - 1)
      end do
      nx(1) = -nx(2); ny(1) = -ny(2)
   end subroutine reverse_normals_for_open

   subroutine builder_init(b, cap)
      type(path_builder), intent(inout) :: b
      integer, intent(in) :: cap
      if (allocated(b%x)) deallocate(b%x, b%y)
      allocate(b%x(max(2, cap)), b%y(max(2, cap)))
      b%n = 0
   end subroutine builder_init

   subroutine builder_add(b, x, y)
      type(path_builder), intent(inout) :: b
      integer(i8), intent(in) :: x, y
      integer(i8), allocatable :: tx(:), ty(:)
      integer :: cap
      if (b%n > 0) then
         if (b%x(b%n) == x .and. b%y(b%n) == y) return
      end if
      if (b%n == size(b%x)) then
         cap = 2 * size(b%x); allocate(tx(cap), ty(cap)); tx(1:b%n) = b%x; ty(1:b%n) = b%y
         call move_alloc(tx, b%x); call move_alloc(ty, b%y)
      end if
      b%n = b%n + 1; b%x(b%n) = x; b%y(b%n) = y
   end subroutine builder_add

   subroutine builder_to_path(b, p)
      type(path_builder), intent(in) :: b
      type(int_path), intent(out) :: p
      integer :: n
      n = b%n
      if (n > 1) then
         if (b%x(n) == b%x(1) .and. b%y(n) == b%y(1)) n = n - 1
      end if
      allocate(p%x(n), p%y(n))
      if (n > 0) then
         p%x = b%x(1:n); p%y = b%y(1:n)
      end if
   end subroutine builder_to_path

   pure integer function round_away_int(x) result(v)
      real(dp), intent(in) :: x
      if (x < 0.0_dp) then
         v = int(x - 0.5_dp)
      else
         v = int(x + 0.5_dp)
      end if
   end function round_away_int

   subroutine add_negative_outer(raw, tmp)
      type(int_set), intent(in) :: raw
      type(int_set), intent(out) :: tmp
      type(int_path) :: outer
      integer(i8) :: xmin, xmax, ymin, ymax
      integer :: i
      if (raw%size() == 0) then
         allocate(tmp%path(0)); return
      end if
      xmin = minval(raw%path(1)%x); xmax = maxval(raw%path(1)%x)
      ymin = minval(raw%path(1)%y); ymax = maxval(raw%path(1)%y)
      do i = 2, raw%size()
         xmin = min(xmin, minval(raw%path(i)%x)); xmax = max(xmax, maxval(raw%path(i)%x))
         ymin = min(ymin, minval(raw%path(i)%y)); ymax = max(ymax, maxval(raw%path(i)%y))
      end do
      tmp = raw
      allocate(outer%x(4), outer%y(4))
      ! Clockwise enclosing rectangle, matching ClipperOffset's negative-delta cleanup.
      outer%x = [xmin - 10_i8, xmax + 10_i8, xmax + 10_i8, xmin - 10_i8]
      outer%y = [ymax + 10_i8, ymax + 10_i8, ymin - 10_i8, ymin - 10_i8]
      call append_int_path(tmp, outer)
   end subroutine add_negative_outer

   subroutine remove_path(s, idx)
      type(int_set), intent(inout) :: s
      integer, intent(in) :: idx
      type(int_path), allocatable :: tmp(:)
      integer :: n
      n = s%size()
      if (idx < 1 .or. idx > n) return
      allocate(tmp(n - 1))
      if (idx > 1) tmp(1:idx - 1) = s%path(1:idx - 1)
      if (idx < n) tmp(idx:n - 1) = s%path(idx + 1:n)
      call move_alloc(tmp, s%path)
   end subroutine remove_path

   subroutine reverse_path(p)
      type(int_path), intent(inout) :: p
      integer(i8), allocatable :: x(:), y(:)
      integer :: i, n
      n = p%size(); allocate(x(n), y(n))
      do i = 1, n
         x(i) = p%x(n + 1 - i); y(i) = p%y(n + 1 - i)
      end do
      p%x = x; p%y = y
   end subroutine reverse_path

end module polyclip_offset
