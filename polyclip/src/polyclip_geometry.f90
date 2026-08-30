module polyclip_geometry
   use polyclip_kinds, only: dp, i8
   use polyclip_types, only: fill_evenodd, fill_nonzero, fill_positive, fill_negative
   use polyclip_integer, only: int_path, int_set, round_away_i8
   implicit none
   private

   public :: int_path_area, int_orientation, point_in_int_path
   public :: fill_contains, segment_intersections, point_on_segment_real
   public :: line_intersection, normalize_path_orientation
   public :: simplify_collinear_path

contains

   pure real(dp) function int_path_area(p) result(a)
      type(int_path), intent(in) :: p
      integer :: i, j, n
      a = 0.0_dp
      n = p%size()
      if (n < 3) return
      j = n
      do i = 1, n
         a = a + real(p%x(j) + p%x(i), dp) * real(p%y(j) - p%y(i), dp)
         j = i
      end do
      a = -0.5_dp * a
   end function int_path_area

   pure logical function int_orientation(p) result(ccw)
      type(int_path), intent(in) :: p
      ccw = int_path_area(p) >= 0.0_dp
   end function int_orientation

   pure integer function point_in_int_path(px, py, path) result(ans)
      integer(i8), intent(in) :: px, py
      type(int_path), intent(in) :: path
      integer :: i, j, n
      real(dp) :: xint
      logical :: inside
      n = path%size()
      if (n < 3) then
         ans = 0
         return
      end if
      inside = .false.
      j = n
      do i = 1, n
         if (point_on_segment_i8(px, py, path%x(j), path%y(j), path%x(i), path%y(i))) then
            ans = -1
            return
         end if
         if ((path%y(i) > py) .neqv. (path%y(j) > py)) then
            xint = real(path%x(j), dp) + real(py - path%y(j), dp) * &
               real(path%x(i) - path%x(j), dp) / real(path%y(i) - path%y(j), dp)
            if (xint > real(px, dp)) inside = .not. inside
         end if
         j = i
      end do
      if (inside) then
         ans = 1
      else
         ans = 0
      end if
   end function point_in_int_path

   pure logical function point_on_segment_i8(px, py, ax, ay, bx, by) result(on)
      integer(i8), intent(in) :: px, py, ax, ay, bx, by
      real(dp) :: crossv, scale
      crossv = real(bx - ax, dp) * real(py - ay, dp) - &
               real(by - ay, dp) * real(px - ax, dp)
      scale = abs(real(bx - ax, dp) * real(py - ay, dp)) + &
              abs(real(by - ay, dp) * real(px - ax, dp))
      on = abs(crossv) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, scale)
      if (.not. on) return
      on = real(px, dp) >= real(min(ax, bx), dp) - 0.25_dp .and. &
           real(px, dp) <= real(max(ax, bx), dp) + 0.25_dp .and. &
           real(py, dp) >= real(min(ay, by), dp) - 0.25_dp .and. &
           real(py, dp) <= real(max(ay, by), dp) + 0.25_dp
   end function point_on_segment_i8

   pure logical function point_on_segment_real(px, py, ax, ay, bx, by, tol) result(on)
      real(dp), intent(in) :: px, py, ax, ay, bx, by, tol
      real(dp) :: crossv, scale
      crossv = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
      scale = abs((bx - ax) * (py - ay)) + abs((by - ay) * (px - ax))
      on = abs(crossv) <= tol * max(1.0_dp, scale)
      if (.not. on) return
      on = px >= min(ax, bx) - tol .and. px <= max(ax, bx) + tol .and. &
           py >= min(ay, by) - tol .and. py <= max(ay, by) + tol
   end function point_on_segment_real

   pure logical function fill_contains(paths, px, py, fill_type) result(inside)
      type(int_set), intent(in) :: paths
      real(dp), intent(in) :: px, py
      integer, intent(in) :: fill_type
      integer :: i, j, k, n, winding, crossings
      real(dp) :: ax, ay, bx, by, crossv
      winding = 0
      crossings = 0
      if (.not. allocated(paths%path)) then
         inside = .false.
         return
      end if
      do k = 1, size(paths%path)
         n = paths%path(k)%size()
         if (n < 3) cycle
         j = n
         do i = 1, n
            ax = real(paths%path(k)%x(j), dp)
            ay = real(paths%path(k)%y(j), dp)
            bx = real(paths%path(k)%x(i), dp)
            by = real(paths%path(k)%y(i), dp)
            if ((ay > py) .neqv. (by > py)) then
               if (ax + (py - ay) * (bx - ax) / (by - ay) > px) crossings = crossings + 1
            end if
            crossv = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
            if (ay <= py) then
               if (by > py .and. crossv > 0.0_dp) winding = winding + 1
            else
               if (by <= py .and. crossv < 0.0_dp) winding = winding - 1
            end if
            j = i
         end do
      end do
      select case (fill_type)
      case (fill_evenodd)
         inside = mod(crossings, 2) /= 0
      case (fill_nonzero)
         inside = winding /= 0
      case (fill_positive)
         inside = winding > 0
      case (fill_negative)
         inside = winding < 0
      case default
         inside = .false.
      end select
   end function fill_contains

   subroutine segment_intersections(ax, ay, bx, by, cx, cy, dx, dy, nout, ix, iy)
      integer(i8), intent(in) :: ax, ay, bx, by, cx, cy, dx, dy
      integer, intent(out) :: nout
      integer(i8), intent(out) :: ix(2), iy(2)
      real(dp) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(dp) :: rx, ry, sx, sy, qpx, qpy, den, numt, numu, t, u
      real(dp) :: scale, tol, px, py
      integer(i8) :: tx(4), ty(4)
      integer :: i, m
      logical :: collinear
      x1 = real(ax, dp); y1 = real(ay, dp)
      x2 = real(bx, dp); y2 = real(by, dp)
      x3 = real(cx, dp); y3 = real(cy, dp)
      x4 = real(dx, dp); y4 = real(dy, dp)
      rx = x2 - x1; ry = y2 - y1
      sx = x4 - x3; sy = y4 - y3
      qpx = x3 - x1; qpy = y3 - y1
      den = rx * sy - ry * sx
      scale = abs(rx * sy) + abs(ry * sx)
      tol = 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, scale)
      nout = 0
      ix = 0_i8; iy = 0_i8
      if (abs(den) > tol) then
         numt = qpx * sy - qpy * sx
         numu = qpx * ry - qpy * rx
         t = numt / den
         u = numu / den
         if (t >= -1.0e-12_dp .and. t <= 1.0_dp + 1.0e-12_dp .and. &
             u >= -1.0e-12_dp .and. u <= 1.0_dp + 1.0e-12_dp) then
            t = max(0.0_dp, min(1.0_dp, t))
            px = x1 + t * rx
            py = y1 + t * ry
            nout = 1
            ix(1) = round_away_i8(px)
            iy(1) = round_away_i8(py)
         end if
         return
      end if

      collinear = abs(qpx * ry - qpy * rx) <= &
         32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(qpx * ry) + abs(qpy * rx))
      if (.not. collinear) return

      m = 0
      call maybe_add(ax, ay, cx, cy, dx, dy, tx, ty, m)
      call maybe_add(bx, by, cx, cy, dx, dy, tx, ty, m)
      call maybe_add(cx, cy, ax, ay, bx, by, tx, ty, m)
      call maybe_add(dx, dy, ax, ay, bx, by, tx, ty, m)
      do i = 1, m
         if (nout == 0) then
            nout = 1; ix(1) = tx(i); iy(1) = ty(i)
         else if (tx(i) /= ix(1) .or. ty(i) /= iy(1)) then
            nout = 2; ix(2) = tx(i); iy(2) = ty(i)
            exit
         end if
      end do
   contains
      subroutine maybe_add(px0, py0, ux, uy, vx, vy, ox, oy, nm)
         integer(i8), intent(in) :: px0, py0, ux, uy, vx, vy
         integer(i8), intent(inout) :: ox(4), oy(4)
         integer, intent(inout) :: nm
         integer :: z
         if (.not. point_on_segment_i8(px0, py0, ux, uy, vx, vy)) return
         do z = 1, nm
            if (ox(z) == px0 .and. oy(z) == py0) return
         end do
         nm = nm + 1
         ox(nm) = px0; oy(nm) = py0
      end subroutine maybe_add
   end subroutine segment_intersections

   subroutine line_intersection(ax, ay, bx, by, cx, cy, dx, dy, ok, px, py)
      real(dp), intent(in) :: ax, ay, bx, by, cx, cy, dx, dy
      logical, intent(out) :: ok
      real(dp), intent(out) :: px, py
      real(dp) :: rx, ry, sx, sy, den, t
      rx = bx - ax; ry = by - ay
      sx = dx - cx; sy = dy - cy
      den = rx * sy - ry * sx
      if (abs(den) <= 64.0_dp * epsilon(1.0_dp) * &
          max(1.0_dp, abs(rx * sy) + abs(ry * sx))) then
         ok = .false.; px = 0.0_dp; py = 0.0_dp; return
      end if
      t = ((cx - ax) * sy - (cy - ay) * sx) / den
      px = ax + t * rx
      py = ay + t * ry
      ok = .true.
   end subroutine line_intersection

   subroutine normalize_path_orientation(p, ccw)
      type(int_path), intent(inout) :: p
      logical, intent(in) :: ccw
      integer(i8), allocatable :: tx(:), ty(:)
      integer :: n, i
      if (p%size() < 3) return
      if (int_orientation(p) .eqv. ccw) return
      n = p%size()
      allocate(tx(n), ty(n))
      do i = 1, n
         tx(i) = p%x(n + 1 - i)
         ty(i) = p%y(n + 1 - i)
      end do
      p%x = tx; p%y = ty
   end subroutine normalize_path_orientation

   subroutine simplify_collinear_path(p, closed)
      type(int_path), intent(inout) :: p
      logical, intent(in) :: closed
      integer(i8), allocatable :: xx(:), yy(:)
      integer :: n, i, prev, nxt, k
      real(dp) :: crossv, scale
      logical :: keep
      n = p%size()
      if (n <= merge(3, 2, closed)) return
      allocate(xx(n), yy(n))
      k = 0
      do i = 1, n
         if (closed) then
            prev = i - 1; if (prev < 1) prev = n
            nxt = i + 1; if (nxt > n) nxt = 1
         else
            if (i == 1 .or. i == n) then
               k = k + 1; xx(k) = p%x(i); yy(k) = p%y(i); cycle
            end if
            prev = i - 1; nxt = i + 1
         end if
         crossv = real(p%x(i) - p%x(prev), dp) * real(p%y(nxt) - p%y(i), dp) - &
                  real(p%y(i) - p%y(prev), dp) * real(p%x(nxt) - p%x(i), dp)
         scale = abs(real(p%x(i) - p%x(prev), dp) * real(p%y(nxt) - p%y(i), dp)) + &
                 abs(real(p%y(i) - p%y(prev), dp) * real(p%x(nxt) - p%x(i), dp))
         keep = abs(crossv) > 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, scale)
         if (keep) then
            k = k + 1; xx(k) = p%x(i); yy(k) = p%y(i)
         end if
      end do
      if (k >= merge(3, 2, closed)) then
         deallocate(p%x, p%y)
         allocate(p%x(k), p%y(k))
         p%x = xx(1:k); p%y = yy(1:k)
      end if
   end subroutine simplify_collinear_path

end module polyclip_geometry
