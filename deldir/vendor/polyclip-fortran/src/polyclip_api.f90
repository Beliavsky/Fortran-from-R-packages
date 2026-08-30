module polyclip_api
   use polyclip_kinds, only: dp, i8
   use polyclip_types, only: poly_path, poly_set, clip_intersection, clip_union, clip_difference, clip_xor, &
      fill_evenodd, join_square, end_closed_polygon
   use polyclip_integer, only: int_set, scale_parameters_one, scale_parameters_two, scale_to_int, scale_from_int
   use polyclip_geometry, only: point_in_int_path
   use polyclip_boolean, only: clip_closed_int, clip_open_int, simplify_int_set
   use polyclip_offset, only: offset_int_set
   use polyclip_minkowski, only: minkowski_sum_int
   implicit none
   private

   public :: polyclip_apply, polysimplify, polyoffset, polylineoffset, polyminkowski, pointinpolygon

contains

   subroutine polyclip_apply(a, b, result, op, fill_a, fill_b, closed, eps, x0, y0, ierr)
      type(poly_set), intent(in) :: a, b
      type(poly_set), intent(out) :: result
      integer, intent(in), optional :: op, fill_a, fill_b
      logical, intent(in), optional :: closed
      real(dp), intent(in), optional :: eps, x0, y0
      integer, intent(out), optional :: ierr
      type(int_set) :: ia, ib, ir
      real(dp) :: de, dx0, dy0
      integer :: iop, ifa, ifb, stat
      logical :: is_closed
      call scale_parameters_two(a, b, de, dx0, dy0)
      if (present(eps)) de = eps
      if (present(x0)) dx0 = x0
      if (present(y0)) dy0 = y0
      iop = clip_intersection; if (present(op)) iop = op
      ifa = fill_evenodd; if (present(fill_a)) ifa = fill_a
      ifb = fill_evenodd; if (present(fill_b)) ifb = fill_b
      is_closed = .true.; if (present(closed)) is_closed = closed
      call scale_to_int(a, dx0, dy0, de, ia, stat)
      if (stat /= 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = stat; return
      end if
      call scale_to_int(b, dx0, dy0, de, ib, stat)
      if (stat /= 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = stat; return
      end if
      if (is_closed) then
         call clip_closed_int(ia, ib, iop, ifa, ifb, ir, stat)
      else
         call clip_open_int(ia, ib, iop, ifb, ir, stat)
      end if
      call scale_from_int(ir, dx0, dy0, de, result)
      if (present(ierr)) ierr = stat
   end subroutine polyclip_apply

   subroutine polysimplify(a, result, fill_type, eps, x0, y0, ierr)
      type(poly_set), intent(in) :: a
      type(poly_set), intent(out) :: result
      integer, intent(in), optional :: fill_type
      real(dp), intent(in), optional :: eps, x0, y0
      integer, intent(out), optional :: ierr
      type(int_set) :: ia, ir
      real(dp) :: de, dx0, dy0
      integer :: ft, stat
      call scale_parameters_one(a, de, dx0, dy0)
      if (present(eps)) de = eps
      if (present(x0)) dx0 = x0
      if (present(y0)) dy0 = y0
      ft = fill_evenodd; if (present(fill_type)) ft = fill_type
      call scale_to_int(a, dx0, dy0, de, ia, stat)
      if (stat /= 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = stat; return
      end if
      call simplify_int_set(ia, ft, ir, stat)
      call scale_from_int(ir, dx0, dy0, de, result)
      if (present(ierr)) ierr = stat
   end subroutine polysimplify

   subroutine polyoffset(a, delta, result, join_type, miter_limit, arc_tolerance, eps, x0, y0, ierr)
      type(poly_set), intent(in) :: a
      real(dp), intent(in) :: delta
      type(poly_set), intent(out) :: result
      integer, intent(in), optional :: join_type
      real(dp), intent(in), optional :: miter_limit, arc_tolerance, eps, x0, y0
      integer, intent(out), optional :: ierr
      type(int_set) :: ia, ir
      real(dp) :: de, dx0, dy0, ml, at
      integer :: jt, stat
      call scale_parameters_one(a, de, dx0, dy0)
      if (present(eps)) de = eps
      if (present(x0)) dx0 = x0
      if (present(y0)) dy0 = y0
      jt = join_square; if (present(join_type)) jt = join_type
      ml = 2.0_dp; if (present(miter_limit)) ml = miter_limit
      at = abs(delta) / 100.0_dp; if (present(arc_tolerance)) at = arc_tolerance
      at = max(de / 4.0_dp, at)
      call scale_to_int(a, dx0, dy0, de, ia, stat)
      if (stat /= 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = stat; return
      end if
      call offset_int_set(ia, delta / de, jt, end_closed_polygon, ml, at / de, ir, stat)
      call scale_from_int(ir, dx0, dy0, de, result)
      if (present(ierr)) ierr = stat
   end subroutine polyoffset

   subroutine polylineoffset(a, delta, result, join_type, end_type, miter_limit, arc_tolerance, eps, x0, y0, ierr)
      type(poly_set), intent(in) :: a
      real(dp), intent(in) :: delta
      type(poly_set), intent(out) :: result
      integer, intent(in), optional :: join_type, end_type
      real(dp), intent(in), optional :: miter_limit, arc_tolerance, eps, x0, y0
      integer, intent(out), optional :: ierr
      type(int_set) :: ia, ir
      real(dp) :: de, dx0, dy0, ml, at
      integer :: jt, et, stat
      call scale_parameters_one(a, de, dx0, dy0)
      if (present(eps)) de = eps
      if (present(x0)) dx0 = x0
      if (present(y0)) dy0 = y0
      jt = join_square; if (present(join_type)) jt = join_type
      et = end_closed_polygon; if (present(end_type)) et = end_type
      ml = 2.0_dp; if (present(miter_limit)) ml = miter_limit
      at = abs(delta) / 100.0_dp; if (present(arc_tolerance)) at = arc_tolerance
      at = max(de / 4.0_dp, at)
      call scale_to_int(a, dx0, dy0, de, ia, stat)
      if (stat /= 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = stat; return
      end if
      call offset_int_set(ia, delta / de, jt, et, ml, at / de, ir, stat)
      call scale_from_int(ir, dx0, dy0, de, result)
      if (present(ierr)) ierr = stat
   end subroutine polylineoffset

   subroutine polyminkowski(a, b, result, closed, eps, x0, y0, ierr)
      type(poly_set), intent(in) :: a, b
      type(poly_set), intent(out) :: result
      logical, intent(in), optional :: closed
      real(dp), intent(in), optional :: eps, x0, y0
      integer, intent(out), optional :: ierr
      type(int_set) :: ia, ib, ir
      real(dp) :: de, dx0, dy0, xmin, xmax, ymin, ymax
      integer :: stat
      logical :: cl
      if (a%size() /= 1) then
         allocate(result%path(0)); if (present(ierr)) ierr = 10; return
      end if
      call combined_minkowski_defaults(a, b, de, dx0, dy0, xmin, xmax, ymin, ymax)
      if (present(eps)) de = eps
      if (present(x0)) dx0 = x0
      if (present(y0)) dy0 = y0
      cl = .true.; if (present(closed)) cl = closed
      call scale_to_int(a, dx0, dy0, de, ia, stat)
      if (stat /= 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = stat; return
      end if
      call scale_to_int(b, dx0, dy0, de, ib, stat)
      if (stat /= 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = stat; return
      end if
      call minkowski_sum_int(ia%path(1), ib, cl, ir, stat)
      call scale_from_int(ir, dx0, dy0, de, result)
      if (present(ierr)) ierr = stat
   end subroutine polyminkowski

   subroutine pointinpolygon(px, py, a, answer, eps, x0, y0, ierr)
      real(dp), intent(in) :: px(:), py(:)
      type(poly_path), intent(in) :: a
      integer, allocatable, intent(out) :: answer(:)
      real(dp), intent(in), optional :: eps, x0, y0
      integer, intent(out), optional :: ierr
      type(poly_set) :: aset
      type(int_set) :: ia
      real(dp) :: de, dx0, dy0, sx, sy, lim
      integer :: i, stat
      if (size(px) /= size(py)) then
         allocate(answer(0)); if (present(ierr)) ierr = 20; return
      end if
      allocate(aset%path(1)); aset%path(1) = a
      call scale_parameters_one(aset, de, dx0, dy0)
      if (present(eps)) de = eps
      if (present(x0)) dx0 = x0
      if (present(y0)) dy0 = y0
      call scale_to_int(aset, dx0, dy0, de, ia, stat)
      if (stat /= 0) then
         allocate(answer(0)); if (present(ierr)) ierr = stat; return
      end if
      allocate(answer(size(px)))
      lim = 0.9_dp * real(huge(0_i8), dp)
      do i = 1, size(px)
         sx = (px(i) - dx0) / de; sy = (py(i) - dy0) / de
         if (abs(sx) > lim .or. abs(sy) > lim) then
            answer(i) = 0
         else
            answer(i) = point_in_int_path(int(sx, i8), int(sy, i8), ia%path(1))
         end if
      end do
      if (present(ierr)) ierr = 0
   end subroutine pointinpolygon

   subroutine combined_minkowski_defaults(a, b, eps, x0, y0, xmin, xmax, ymin, ymax)
      type(poly_set), intent(in) :: a, b
      real(dp), intent(out) :: eps, x0, y0, xmin, xmax, ymin, ymax
      integer :: i
      logical :: first
      first = .true.; xmin = 0.0_dp; xmax = 0.0_dp; ymin = 0.0_dp; ymax = 0.0_dp
      call absorb(a); call absorb(b)
      if (first) then
         xmin = 0.0_dp; xmax = 1.0_dp; ymin = 0.0_dp; ymax = 1.0_dp
      end if
      eps = max(xmax - xmin, ymax - ymin) / 1.0e9_dp
      if (eps <= 0.0_dp) eps = 1.0e-9_dp
      x0 = xmin
      y0 = ymin - (ymax - ymin) / 16.0_dp
   contains
      subroutine absorb(s)
         type(poly_set), intent(in) :: s
         if (.not. allocated(s%path)) return
         do i = 1, size(s%path)
            if (s%path(i)%size() == 0) cycle
            if (first) then
               xmin = minval(s%path(i)%x); xmax = maxval(s%path(i)%x)
               ymin = minval(s%path(i)%y); ymax = maxval(s%path(i)%y); first = .false.
            else
               xmin = min(xmin, minval(s%path(i)%x)); xmax = max(xmax, maxval(s%path(i)%x))
               ymin = min(ymin, minval(s%path(i)%y)); ymax = max(ymax, maxval(s%path(i)%y))
            end if
         end do
      end subroutine absorb
   end subroutine combined_minkowski_defaults

end module polyclip_api
