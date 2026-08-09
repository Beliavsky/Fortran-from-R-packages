! SPDX-License-Identifier: GPL-2.0-or-later
module coneproj_core
   use coneproj_kinds, only : dp
   use coneproj_types, only : cone_result, qprog_result, coneproj_success, coneproj_max_iter, &
      coneproj_invalid_input, coneproj_singular
   use coneproj_linalg, only : solve_spd, least_squares, cholesky_upper, inverse_upper, &
      solve_upper, solve_upper_transpose
   implicit none
   private
   public :: cone_a, cone_b, qprog

contains

   subroutine cone_a(y, amat, result, weights, start_face, max_iter)
      real(dp), intent(in) :: y(:), amat(:,:)
      type(cone_result), intent(out) :: result
      real(dp), intent(in), optional :: weights(:)
      integer, intent(in), optional :: start_face(:)
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: yt(:), at(:,:), sw(:)
      integer :: n, m, imax

      n = size(y)
      m = size(amat,1)
      if (size(amat,2) /= n .or. n < 1 .or. m < 1) then
         result%status = coneproj_invalid_input
         return
      end if
      allocate(yt(n), at(m,n))
      yt = y
      at = amat
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            result%status = coneproj_invalid_input
            return
         end if
         allocate(sw(n))
         sw = sqrt(weights)
         yt = sw * y
         at = amat
         at = at / spread(sw, dim=1, ncopies=m)
      end if
      imax = 10000
      if (present(max_iter)) imax = max_iter
      call cone_a_core(yt, at, result, start_face, imax)
      if (present(weights) .and. allocated(result%fit)) result%fit = result%fit / sw
   end subroutine cone_a

   subroutine cone_a_core(y, amat, result, start_face, max_iter)
      real(dp), intent(in) :: y(:), amat(:,:)
      type(cone_result), intent(out) :: result
      integer, intent(in), optional :: start_face(:)
      integer, intent(in) :: max_iter
      real(dp), allocatable :: delta(:,:), b2(:), theta(:), xmat(:,:), a(:), gram(:,:), rhs(:), avec(:)
      logical, allocatable :: active(:)
      integer, allocatable :: idx(:)
      real(dp) :: rownorm, vmax, amin
      real(dp), parameter :: sm = 1.0e-11_dp
      integer :: n, m, i, k, na, info, add_i, rem_i
      logical :: check

      n = size(y)
      m = size(amat,1)
      allocate(delta(m,n), b2(m), theta(n), active(m), avec(m))
      do i = 1, m
         rownorm = sqrt(dot_product(amat(i,:), amat(i,:)))
         if (rownorm <= tiny(1.0_dp)) then
            result%status = coneproj_invalid_input
            return
         end if
         delta(i,:) = -amat(i,:) / rownorm
      end do
      b2 = matmul(delta, y)
      active = .false.
      theta = 0.0_dp
      check = .false.
      vmax = maxval(b2)
      if (vmax > 2.0_dp * sm) then
         add_i = first_max_index(b2)
         active(add_i) = .true.
         if (present(start_face)) then
            do k = 1, size(start_face)
               if (start_face(k) >= 1 .and. start_face(k) <= m) active(start_face(k)) = .true.
            end do
         end if
      else
         check = .true.
      end if
      result%steps = 0
      allocate(xmat(0,n))

      do while (.not. check .and. result%steps < max_iter)
         result%steps = result%steps + 1
         call active_indices(active, idx)
         na = size(idx)
         deallocate(xmat)
         allocate(xmat(na,n))
         do k = 1, na
            xmat(k,:) = delta(idx(k),:)
         end do
         allocate(gram(na,na), rhs(na))
         gram = matmul(xmat, transpose(xmat))
         rhs = matmul(xmat, y)
         call solve_spd(gram, rhs, a, info)
         if (info /= 0) then
            call least_squares(transpose(xmat), y, a, info)
         end if
         deallocate(gram, rhs)
         if (info /= 0) then
            result%status = coneproj_singular
            return
         end if
         amin = minval(a)
         if (amin < -sm) then
            avec = huge(1.0_dp)
            do k = 1, na
               avec(idx(k)) = a(k)
            end do
            rem_i = first_min_active(avec, active)
            active(rem_i) = .false.
            check = .false.
         else
            theta = matmul(transpose(xmat), a)
            b2 = matmul(delta, y - theta) / real(n, dp)
            vmax = maxval(b2)
            if (vmax > 2.0_dp * sm) then
               add_i = first_max_index(b2)
               active(add_i) = .true.
               check = .false.
            else
               check = .true.
            end if
         end if
         if (allocated(a)) deallocate(a)
         if (allocated(idx)) deallocate(idx)
      end do

      result%status = coneproj_success
      if (.not. check) result%status = coneproj_max_iter
      result%df = n - count(active)
      allocate(result%fit(n), result%coefs(m), result%xmat(size(xmat,1),n))
      result%fit = y - theta
      result%coefs = 0.0_dp
      result%xmat = xmat
      call active_indices(active, result%face)
   end subroutine cone_a_core

   subroutine cone_b(y, delta, result, vmat, weights, start_face, max_iter)
      real(dp), intent(in) :: y(:), delta(:,:)
      type(cone_result), intent(out) :: result
      real(dp), intent(in), optional :: vmat(:,:)
      real(dp), intent(in), optional :: weights(:)
      integer, intent(in), optional :: start_face(:)
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: yt(:), dt(:,:), vt(:,:), sw(:)
      integer :: n, m, p, imax

      n = size(y)
      m = size(delta,2)
      p = 0
      if (present(vmat)) p = size(vmat,2)
      if (size(delta,1) /= n .or. m < 1) then
         result%status = coneproj_invalid_input
         return
      end if
      if (present(vmat)) then
         if (size(vmat,1) /= n) then
            result%status = coneproj_invalid_input
            return
         end if
      end if
      allocate(yt(n), dt(n,m), vt(n,p))
      yt = y
      dt = delta
      if (p > 0) vt = vmat
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            result%status = coneproj_invalid_input
            return
         end if
         allocate(sw(n))
         sw = sqrt(weights)
         yt = sw * y
         dt = spread(sw, dim=2, ncopies=m) * delta
         if (p > 0) vt = spread(sw, dim=2, ncopies=p) * vmat
      end if
      imax = n * n
      if (present(max_iter)) imax = max_iter
      call cone_b_core(yt, dt, vt, result, start_face, imax)
      if (present(weights) .and. allocated(result%fit)) result%fit = result%fit / sw
   end subroutine cone_b

   subroutine cone_b_core(y, delta, vmat, result, start_face, max_iter)
      real(dp), intent(in) :: y(:), delta(:,:), vmat(:,:)
      type(cone_result), intent(out) :: result
      integer, intent(in), optional :: start_face(:)
      integer, intent(in) :: max_iter
      real(dp), allocatable :: sigma(:,:), scalars(:), theta(:), b2(:), a(:), xmat(:,:), gram(:,:), rhs(:), avec(:)
      logical, allocatable :: active(:)
      integer, allocatable :: idx(:)
      real(dp), parameter :: sm = 1.0e-5_dp
      real(dp) :: rownorm, vmax, amin
      integer :: n, m, p, nt, i, k, na, info, add_i, rem_i
      logical :: check

      n = size(y)
      m = size(delta,2)
      p = size(vmat,2)
      nt = p + m
      allocate(sigma(nt,n), scalars(m), active(nt), theta(n), b2(nt), avec(nt))
      active = .false.
      theta = 0.0_dp
      if (p > 0) then
         sigma(1:p,:) = transpose(vmat)
         active(1:p) = .true.
         call least_squares(vmat, y, a, info)
         if (info /= 0) then
            result%status = coneproj_singular
            return
         end if
         theta = matmul(vmat, a)
         deallocate(a)
      end if
      do i = 1, m
         rownorm = sqrt(dot_product(delta(:,i), delta(:,i)))
         if (rownorm <= tiny(1.0_dp)) then
            result%status = coneproj_invalid_input
            return
         end if
         scalars(i) = rownorm
         sigma(p+i,:) = delta(:,i) / rownorm
      end do
      b2 = matmul(sigma, y - theta) / real(n, dp)
      check = .false.
      vmax = maxval(b2)
      if (vmax > 2.0_dp * sm) then
         add_i = first_max_index(b2)
         active(add_i) = .true.
         if (present(start_face)) then
            do k = 1, size(start_face)
               if (start_face(k) >= 1 .and. start_face(k) <= nt) active(start_face(k)) = .true.
            end do
         end if
      else
         check = .true.
      end if
      result%steps = 0
      allocate(xmat(0,n))

      if (check) then
         allocate(result%fit(n), result%coefs(nt), result%xmat(0,n))
         result%fit = theta
         result%coefs = 0.0_dp
         if (p > 0) then
            call least_squares(vmat, y, a, info)
            if (info == 0) result%coefs(1:p) = a
         end if
         result%df = count(active)
         call positive_indices(result%coefs, result%face)
         result%status = coneproj_success
         return
      end if

      do while (.not. check .and. result%steps < max_iter)
         result%steps = result%steps + 1
         call active_indices(active, idx)
         na = size(idx)
         deallocate(xmat)
         allocate(xmat(na,n))
         do k = 1, na
            xmat(k,:) = sigma(idx(k),:)
         end do
         allocate(gram(na,na), rhs(na))
         gram = matmul(xmat, transpose(xmat))
         rhs = matmul(xmat, y)
         call solve_spd(gram, rhs, a, info)
         if (info /= 0) call least_squares(transpose(xmat), y, a, info)
         deallocate(gram, rhs)
         if (info /= 0) then
            result%status = coneproj_singular
            return
         end if
         if (na > p) then
            amin = huge(1.0_dp)
            do k = 1, na
               if (idx(k) > p) amin = min(amin, a(k))
            end do
            if (amin < -sm) then
               avec = huge(1.0_dp)
               do k = 1, na
                  if (idx(k) > p) avec(idx(k)) = a(k)
               end do
               rem_i = first_min_edge(avec, active, p)
               active(rem_i) = .false.
               check = .false.
            else
               theta = matmul(transpose(xmat), a)
               b2 = matmul(sigma, y - theta) / real(n,dp)
               vmax = maxval(b2)
               if (vmax > 2.0_dp * sm) then
                  add_i = first_max_index(b2)
                  active(add_i) = .true.
                  check = .false.
               else
                  check = .true.
               end if
            end if
         else
            theta = matmul(transpose(xmat), a)
            check = .true.
         end if
         if (allocated(a)) deallocate(a)
         if (allocated(idx)) deallocate(idx)
      end do

      call active_indices(active, idx)
      na = size(idx)
      if (na > 0) then
         deallocate(xmat)
         allocate(xmat(na,n), gram(na,na), rhs(na))
         do k = 1, na
            xmat(k,:) = sigma(idx(k),:)
         end do
         gram = matmul(xmat, transpose(xmat))
         rhs = matmul(xmat, y)
         call solve_spd(gram, rhs, a, info)
         if (info /= 0) call least_squares(transpose(xmat), y, a, info)
         theta = matmul(transpose(xmat), a)
      end if
      allocate(result%fit(n), result%coefs(nt), result%xmat(size(xmat,1),n))
      result%fit = theta
      result%coefs = 0.0_dp
      if (na > 0) then
         do k = 1, na
            result%coefs(idx(k)) = a(k)
         end do
      end if
      do i = 1, m
         result%coefs(p+i) = result%coefs(p+i) / scalars(i)
      end do
      result%xmat = xmat
      result%df = count(active)
      call positive_indices(result%coefs, result%face)
      result%status = coneproj_success
      if (.not. check) result%status = coneproj_max_iter
   end subroutine cone_b_core

   subroutine qprog(q, c, amat, b, result, start_face, max_iter)
      real(dp), intent(in) :: q(:,:), c(:), amat(:,:), b(:)
      type(qprog_result), intent(out) :: result
      integer, intent(in), optional :: start_face(:)
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: u(:,:), uinv(:,:), theta0(:), nnc(:), z(:), atil(:,:), beta(:)
      type(cone_result) :: cres
      integer :: n, m, info, imax
      logical :: shifted

      n = size(c)
      m = size(amat,1)
      if (size(q,1) /= n .or. size(q,2) /= n .or. size(amat,2) /= n .or. size(b) /= m) then
         result%status = coneproj_invalid_input
         return
      end if
      call cholesky_upper(q, u, info)
      if (info /= 0) then
         result%status = coneproj_singular
         return
      end if
      call inverse_upper(u, uinv, info)
      if (info /= 0) then
         result%status = coneproj_singular
         return
      end if
      allocate(theta0(n), nnc(n), atil(m,n))
      theta0 = 0.0_dp
      shifted = any(abs(b) > 0.0_dp)
      if (shifted) then
         call least_squares(amat, b, theta0, info)
         if (info /= 0) then
            result%status = coneproj_singular
            return
         end if
         nnc = c - matmul(q, theta0)
      else
         nnc = c
      end if
      call solve_upper_transpose(u, nnc, z, info)
      if (info /= 0) then
         result%status = coneproj_singular
         return
      end if
      atil = matmul(amat, uinv)
      imax = 10000
      if (present(max_iter)) imax = max_iter
      call cone_a_core(z, atil, cres, start_face, imax)
      call solve_upper(u, cres%fit, beta, info)
      if (info /= 0) then
         result%status = coneproj_singular
         return
      end if
      beta = beta + theta0
      allocate(result%theta(n), result%xmat(size(cres%xmat,1),n), result%face(size(cres%face)))
      result%theta = beta
      result%xmat = cres%xmat
      result%face = cres%face
      result%df = cres%df
      result%steps = cres%steps
      result%status = cres%status
      result%objective = 0.5_dp * dot_product(beta, matmul(q,beta)) - dot_product(c,beta)
   end subroutine qprog

   subroutine active_indices(active, idx)
      logical, intent(in) :: active(:)
      integer, allocatable, intent(out) :: idx(:)
      integer :: i, k
      allocate(idx(count(active)))
      k = 0
      do i = 1, size(active)
         if (active(i)) then
            k = k + 1
            idx(k) = i
         end if
      end do
   end subroutine active_indices

   subroutine positive_indices(x, idx)
      real(dp), intent(in) :: x(:)
      integer, allocatable, intent(out) :: idx(:)
      integer :: i, k
      allocate(idx(count(x > 0.0_dp)))
      k = 0
      do i = 1, size(x)
         if (x(i) > 0.0_dp) then
            k = k + 1
            idx(k) = i
         end if
      end do
   end subroutine positive_indices

   integer function first_max_index(x) result(idx)
      real(dp), intent(in) :: x(:)
      real(dp) :: vmax
      integer :: i
      vmax = maxval(x)
      idx = 1
      do i = 1, size(x)
         if (x(i) >= vmax) then
            idx = i
            exit
         end if
      end do
   end function first_max_index

   integer function first_min_active(x, active) result(idx)
      real(dp), intent(in) :: x(:)
      logical, intent(in) :: active(:)
      real(dp) :: vmin
      integer :: i
      vmin = huge(1.0_dp)
      idx = 1
      do i = 1, size(x)
         if (active(i) .and. x(i) < vmin) then
            vmin = x(i)
            idx = i
         end if
      end do
   end function first_min_active

   integer function first_min_edge(x, active, p) result(idx)
      real(dp), intent(in) :: x(:)
      logical, intent(in) :: active(:)
      integer, intent(in) :: p
      real(dp) :: vmin
      integer :: i
      vmin = huge(1.0_dp)
      idx = p + 1
      do i = p + 1, size(x)
         if (active(i) .and. x(i) < vmin) then
            vmin = x(i)
            idx = i
         end if
      end do
   end function first_min_edge

end module coneproj_core
