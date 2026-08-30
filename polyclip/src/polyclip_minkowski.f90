module polyclip_minkowski
   use polyclip_kinds, only: i8
   use polyclip_types, only: fill_nonzero
   use polyclip_integer, only: int_path, int_set, append_int_path
   use polyclip_geometry, only: int_orientation
   use polyclip_boolean, only: simplify_int_set
   implicit none
   private
   public :: minkowski_sum_int, minkowski_diff_int

contains

   subroutine minkowski_sum_int(pattern, paths, closed, result, ierr)
      type(int_path), intent(in) :: pattern
      type(int_set), intent(in) :: paths
      logical, intent(in) :: closed
      type(int_set), intent(out) :: result
      integer, intent(out), optional :: ierr
      type(int_set) :: raw
      type(int_path) :: q, translated
      integer(i8), allocatable :: ppx(:, :), ppy(:, :)
      integer :: ip, i, j, ni, np, delta, stat
      allocate(raw%path(0))
      np = pattern%size()
      if (np == 0 .or. paths%size() == 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = 0; return
      end if
      do ip = 1, paths%size()
         ni = paths%path(ip)%size(); if (ni == 0) cycle
         allocate(ppx(ni, np), ppy(ni, np))
         do i = 1, ni
            do j = 1, np
               ppx(i, j) = paths%path(ip)%x(i) + pattern%x(j)
               ppy(i, j) = paths%path(ip)%y(i) + pattern%y(j)
            end do
         end do
         delta = merge(1, 0, closed)
         do i = 1, ni - 1 + delta
            do j = 1, np
               allocate(q%x(4), q%y(4))
               q%x = [ppx(mod(i - 1, ni) + 1, mod(j - 1, np) + 1), &
                      ppx(mod(i, ni) + 1, mod(j - 1, np) + 1), &
                      ppx(mod(i, ni) + 1, mod(j, np) + 1), &
                      ppx(mod(i - 1, ni) + 1, mod(j, np) + 1)]
               q%y = [ppy(mod(i - 1, ni) + 1, mod(j - 1, np) + 1), &
                      ppy(mod(i, ni) + 1, mod(j - 1, np) + 1), &
                      ppy(mod(i, ni) + 1, mod(j, np) + 1), &
                      ppy(mod(i - 1, ni) + 1, mod(j, np) + 1)]
               if (.not. int_orientation(q)) call reverse_path(q)
               call append_int_path(raw, q)
               deallocate(q%x, q%y)
            end do
         end do
         if (closed) then
            allocate(translated%x(ni), translated%y(ni))
            translated%x = paths%path(ip)%x + pattern%x(1)
            translated%y = paths%path(ip)%y + pattern%y(1)
            call append_int_path(raw, translated)
            deallocate(translated%x, translated%y)
         end if
         deallocate(ppx, ppy)
      end do
      call simplify_int_set(raw, fill_nonzero, result, stat)
      if (present(ierr)) ierr = stat
   end subroutine minkowski_sum_int

   subroutine minkowski_diff_int(a, b, result, ierr)
      type(int_path), intent(in) :: a, b
      type(int_set), intent(out) :: result
      integer, intent(out), optional :: ierr
      type(int_set) :: raw, paths
      type(int_path) :: q
      integer(i8), allocatable :: ppx(:, :), ppy(:, :)
      integer :: i, j, na, nb, stat
      allocate(raw%path(0), paths%path(1))
      paths%path(1) = b
      na = a%size(); nb = b%size()
      if (na == 0 .or. nb == 0) then
         allocate(result%path(0)); if (present(ierr)) ierr = 0; return
      end if
      allocate(ppx(nb, na), ppy(nb, na))
      do i = 1, nb
         do j = 1, na
            ppx(i, j) = b%x(i) - a%x(j)
            ppy(i, j) = b%y(i) - a%y(j)
         end do
      end do
      do i = 1, nb
         do j = 1, na
            allocate(q%x(4), q%y(4))
            q%x = [ppx(i, j), ppx(mod(i, nb) + 1, j), ppx(mod(i, nb) + 1, mod(j, na) + 1), &
                   ppx(i, mod(j, na) + 1)]
            q%y = [ppy(i, j), ppy(mod(i, nb) + 1, j), ppy(mod(i, nb) + 1, mod(j, na) + 1), &
                   ppy(i, mod(j, na) + 1)]
            if (.not. int_orientation(q)) call reverse_path(q)
            call append_int_path(raw, q)
            deallocate(q%x, q%y)
         end do
      end do
      call simplify_int_set(raw, fill_nonzero, result, stat)
      if (present(ierr)) ierr = stat
   end subroutine minkowski_diff_int

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

end module polyclip_minkowski
