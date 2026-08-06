! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_smooths
   use mgcv_kinds, only : dp, pi_dp
   use mgcv_utils, only : center_columns, difference_matrix, &
      tensor_product_model_matrix, tensor_product_penalties, standardize_vector
   use splines, only : bs_basis, natural_spline_basis, type7_quantile
   implicit none
   private

   integer, parameter, public :: smooth_cr = 1
   integer, parameter, public :: smooth_ps = 2
   integer, parameter, public :: smooth_cyclic = 3
   integer, parameter, public :: smooth_tp1 = 4
   integer, parameter, public :: smooth_tp2 = 5
   integer, parameter, public :: smooth_random = 6

   type, public :: smooth_spec_t
      integer :: kind = 0
      integer :: degree = 3
      integer :: penalty_order = 2
      integer :: ncoef = 0
      real(dp) :: bounds(2) = 0.0_dp
      real(dp) :: means(2) = 0.0_dp
      real(dp) :: scales(2) = 1.0_dp
      real(dp) :: period = 0.0_dp
      real(dp), allocatable :: knots(:)
      real(dp), allocatable :: column_means(:)
      real(dp), allocatable :: centers(:, :)
      integer, allocatable :: levels(:)
      real(dp), allocatable :: penalties(:, :, :)
   end type smooth_spec_t

   public :: construct_cr_smooth, construct_ps_smooth
   public :: construct_cyclic_smooth, construct_tp_smooth_1d
   public :: construct_tp_smooth_2d, construct_random_effect
   public :: predict_smooth, tensor_smooth, append_columns
   public :: embed_penalty, identity_penalty, place_knots

contains

   subroutine construct_cr_smooth(x, k, basis, spec, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: basis(:, :)
      type(smooth_spec_t), intent(out) :: spec
      integer, intent(out) :: status
      real(dp), allocatable :: raw(:, :), grid(:), bg(:, :), d2(:, :), s(:, :)
      integer :: i, ng
      real(dp) :: h

      status = 0
      if (size(x) < 4 .or. k < 3 .or. maxval(x) <= minval(x)) then
         allocate(basis(0, 0)); status = 1; return
      end if
      spec%kind = smooth_cr; spec%ncoef = k; spec%bounds = [minval(x), maxval(x)]
      call natural_spline_basis(x, raw, df=k, intercept=.false., &
                                boundary_knots=spec%bounds, interior_knots=spec%knots, status=status)
      if (status /= 0) then
         allocate(basis(0, 0)); return
      end if
      basis = raw
      call center_columns(basis, spec%column_means)
      ng = max(100, 10 * k)
      allocate(grid(ng))
      do i = 1, ng
         grid(i) = spec%bounds(1) + real(i - 1, dp) * (spec%bounds(2) - spec%bounds(1)) / real(ng - 1, dp)
      end do
      h = grid(2) - grid(1)
      call natural_spline_basis(grid, bg, knots=spec%knots, intercept=.false., &
                                boundary_knots=spec%bounds, status=status)
      if (status /= 0) return
      allocate(d2(ng - 2, size(bg, 2)))
      d2 = (bg(3:ng, :) - 2.0_dp * bg(2:ng - 1, :) + bg(1:ng - 2, :)) / (h * h)
      allocate(s(k, k)); s = matmul(transpose(d2), d2) * h
      allocate(spec%penalties(k, k, 1)); spec%penalties(:, :, 1) = s
   end subroutine construct_cr_smooth

   subroutine construct_ps_smooth(x, k, basis, spec, status, degree, penalty_order)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: basis(:, :)
      type(smooth_spec_t), intent(out) :: spec
      integer, intent(out) :: status
      integer, intent(in), optional :: degree, penalty_order
      real(dp), allocatable :: raw(:, :), d(:, :)
      integer :: deg, po

      deg = 3; if (present(degree)) deg = degree
      po = 2; if (present(penalty_order)) po = penalty_order
      status = 0
      if (size(x) < 2 .or. k <= deg .or. po < 1 .or. po >= k .or. maxval(x) <= minval(x)) then
         allocate(basis(0, 0)); status = 1; return
      end if
      spec%kind = smooth_ps; spec%degree = deg; spec%penalty_order = po
      spec%ncoef = k; spec%bounds = [minval(x), maxval(x)]
      call bs_basis(x, raw, degree=deg, df=k, intercept=.true., &
                    boundary_knots=spec%bounds, interior_knots=spec%knots, status=status)
      if (status /= 0) then
         allocate(basis(0, 0)); return
      end if
      basis = raw
      call center_columns(basis, spec%column_means)
      d = difference_matrix(k, po)
      allocate(spec%penalties(k, k, 1))
      spec%penalties(:, :, 1) = matmul(transpose(d), d)
   end subroutine construct_ps_smooth

   subroutine construct_cyclic_smooth(x, k, basis, spec, status, period, origin)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: basis(:, :)
      type(smooth_spec_t), intent(out) :: spec
      integer, intent(out) :: status
      real(dp), intent(in), optional :: period, origin
      real(dp) :: p, x0, angle, freq
      integer :: i, j, harmonic

      status = 0
      if (size(x) < 2 .or. k < 2) then
         allocate(basis(0, 0)); status = 1; return
      end if
      x0 = minval(x); if (present(origin)) x0 = origin
      p = maxval(x) - minval(x); if (present(period)) p = period
      if (p <= 0.0_dp) then
         allocate(basis(0, 0)); status = 2; return
      end if
      spec%kind = smooth_cyclic; spec%ncoef = k; spec%period = p
      spec%bounds = [x0, x0 + p]
      allocate(basis(size(x), k), spec%column_means(k), spec%penalties(k, k, 1))
      basis = 0.0_dp; spec%penalties = 0.0_dp
      do j = 1, k
         harmonic = (j + 1) / 2
         freq = 2.0_dp * pi_dp * real(harmonic, dp) / p
         do i = 1, size(x)
            angle = freq * (x(i) - x0)
            if (mod(j, 2) == 1) then
               basis(i, j) = cos(angle)
            else
               basis(i, j) = sin(angle)
            end if
         end do
         spec%penalties(j, j, 1) = freq**4
      end do
      spec%column_means = 0.0_dp
   end subroutine construct_cyclic_smooth

   subroutine construct_tp_smooth_1d(x, k, basis, spec, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: basis(:, :)
      type(smooth_spec_t), intent(out) :: spec
      integer, intent(out) :: status
      real(dp), allocatable :: z(:)
      integer :: i, j, nc
      real(dp) :: p

      status = 0
      if (size(x) < 4 .or. k < 3 .or. maxval(x) <= minval(x)) then
         allocate(basis(0, 0)); status = 1; return
      end if
      spec%kind = smooth_tp1; spec%ncoef = k
      call standardize_vector(x, spec%means(1), spec%scales(1), z)
      nc = k - 1
      allocate(spec%centers(1, nc), basis(size(x), k))
      basis(:, 1) = z
      do j = 1, nc
         p = real(j, dp) / real(nc + 1, dp)
         spec%centers(1, j) = type7_quantile(z, p)
         do i = 1, size(x)
            basis(i, j + 1) = abs(z(i) - spec%centers(1, j))**3
         end do
      end do
      call center_columns(basis, spec%column_means)
      allocate(spec%penalties(k, k, 1)); spec%penalties = 0.0_dp
      do j = 2, k
         spec%penalties(j, j, 1) = 1.0_dp
      end do
   end subroutine construct_tp_smooth_1d

   subroutine construct_tp_smooth_2d(x, y, k, basis, spec, status)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: basis(:, :)
      type(smooth_spec_t), intent(out) :: spec
      integer, intent(out) :: status
      real(dp), allocatable :: zx(:), zy(:)
      real(dp) :: p, r2
      integer :: i, j, nc, idx

      status = 0
      if (size(x) /= size(y) .or. size(x) < 6 .or. k < 5) then
         allocate(basis(0, 0)); status = 1; return
      end if
      spec%kind = smooth_tp2; spec%ncoef = k
      call standardize_vector(x, spec%means(1), spec%scales(1), zx)
      call standardize_vector(y, spec%means(2), spec%scales(2), zy)
      nc = k - 2
      allocate(spec%centers(2, nc), basis(size(x), k))
      basis(:, 1) = zx; basis(:, 2) = zy
      do j = 1, nc
         p = real(j, dp) / real(nc + 1, dp)
         spec%centers(1, j) = type7_quantile(zx, p)
         idx = 1 + modulo(nint(real(j, dp) * 0.6180339887498949_dp * real(size(x), dp)), size(x))
         spec%centers(2, j) = zy(idx)
         do i = 1, size(x)
            r2 = (zx(i) - spec%centers(1, j))**2 + (zy(i) - spec%centers(2, j))**2
            if (r2 <= 1.0e-28_dp) then
               basis(i, j + 2) = 0.0_dp
            else
               basis(i, j + 2) = 0.5_dp * r2 * log(r2)
            end if
         end do
      end do
      call center_columns(basis, spec%column_means)
      allocate(spec%penalties(k, k, 1)); spec%penalties = 0.0_dp
      do j = 3, k
         spec%penalties(j, j, 1) = 1.0_dp
      end do
   end subroutine construct_tp_smooth_2d

   subroutine construct_random_effect(group, basis, spec, status)
      integer, intent(in) :: group(:)
      real(dp), allocatable, intent(out) :: basis(:, :)
      type(smooth_spec_t), intent(out) :: spec
      integer, intent(out) :: status
      integer, allocatable :: temp(:)
      integer :: i, j, nlev, value
      logical :: found

      status = 0
      if (size(group) == 0) then
         allocate(basis(0, 0)); status = 1; return
      end if
      allocate(temp(size(group))); nlev = 0
      do i = 1, size(group)
         value = group(i); found = .false.
         do j = 1, nlev
            if (temp(j) == value) then; found = .true.; exit; end if
         end do
         if (.not. found) then; nlev = nlev + 1; temp(nlev) = value; end if
      end do
      allocate(spec%levels(nlev)); spec%levels = temp(1:nlev)
      allocate(basis(size(group), nlev)); basis = 0.0_dp
      do i = 1, size(group)
         do j = 1, nlev
            if (group(i) == spec%levels(j)) then; basis(i, j) = 1.0_dp; exit; end if
         end do
      end do
      spec%kind = smooth_random; spec%ncoef = nlev
      allocate(spec%column_means(nlev)); spec%column_means = 0.0_dp
      allocate(spec%penalties(nlev, nlev, 1)); spec%penalties = 0.0_dp
      do i = 1, nlev; spec%penalties(i, i, 1) = 1.0_dp; end do
   end subroutine construct_random_effect

   subroutine predict_smooth(spec, x, basis, status, y, group)
      type(smooth_spec_t), intent(in) :: spec
      real(dp), intent(in), optional :: x(:), y(:)
      integer, intent(in), optional :: group(:)
      real(dp), allocatable, intent(out) :: basis(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: raw(:, :), z(:), zy(:)
      real(dp) :: angle, freq, r2
      integer :: i, j, harmonic

      status = 0
      select case (spec%kind)
      case (smooth_cr)
         if (.not. present(x)) then; allocate(basis(0, 0)); status = 1; return; end if
         call natural_spline_basis(x, raw, knots=spec%knots, intercept=.false., &
                                   boundary_knots=spec%bounds, status=status)
         if (status /= 0) then; allocate(basis(0, 0)); return; end if
         allocate(basis(size(raw, 1), size(raw, 2)))
         basis = raw - spread(spec%column_means, 1, size(raw, 1))
      case (smooth_ps)
         if (.not. present(x)) then; allocate(basis(0, 0)); status = 1; return; end if
         call bs_basis(x, raw, degree=spec%degree, knots=spec%knots, intercept=.true., &
                       boundary_knots=spec%bounds, status=status)
         if (status /= 0) then; allocate(basis(0, 0)); return; end if
         allocate(basis(size(raw, 1), size(raw, 2)))
         basis = raw - spread(spec%column_means, 1, size(raw, 1))
      case (smooth_cyclic)
         if (.not. present(x)) then; allocate(basis(0, 0)); status = 1; return; end if
         allocate(basis(size(x), spec%ncoef))
         do j = 1, spec%ncoef
            harmonic = (j + 1) / 2
            freq = 2.0_dp * pi_dp * real(harmonic, dp) / spec%period
            do i = 1, size(x)
               angle = freq * (x(i) - spec%bounds(1))
               if (mod(j, 2) == 1) then
                  basis(i, j) = cos(angle)
               else
                  basis(i, j) = sin(angle)
               end if
            end do
         end do
      case (smooth_tp1)
         if (.not. present(x)) then; allocate(basis(0, 0)); status = 1; return; end if
         allocate(z(size(x))); z = (x - spec%means(1)) / spec%scales(1)
         allocate(raw(size(x), spec%ncoef)); raw(:, 1) = z
         do j = 1, size(spec%centers, 2)
            raw(:, j + 1) = abs(z - spec%centers(1, j))**3
         end do
         allocate(basis(size(raw, 1), size(raw, 2)))
         basis = raw - spread(spec%column_means, 1, size(raw, 1))
      case (smooth_tp2)
         if (.not. present(x) .or. .not. present(y)) then
            allocate(basis(0, 0)); status = 1; return
         end if
         if (size(x) /= size(y)) then; allocate(basis(0, 0)); status = 2; return; end if
         allocate(z(size(x)), zy(size(y)))
         z = (x - spec%means(1)) / spec%scales(1)
         zy = (y - spec%means(2)) / spec%scales(2)
         allocate(raw(size(x), spec%ncoef)); raw(:, 1) = z; raw(:, 2) = zy
         do j = 1, size(spec%centers, 2)
            do i = 1, size(x)
               r2 = (z(i) - spec%centers(1, j))**2 + (zy(i) - spec%centers(2, j))**2
               raw(i, j + 2) = merge(0.5_dp * r2 * log(r2), 0.0_dp, r2 > 1.0e-28_dp)
            end do
         end do
         allocate(basis(size(raw, 1), size(raw, 2)))
         basis = raw - spread(spec%column_means, 1, size(raw, 1))
      case (smooth_random)
         if (.not. present(group)) then; allocate(basis(0, 0)); status = 1; return; end if
         allocate(basis(size(group), spec%ncoef)); basis = 0.0_dp
         do i = 1, size(group)
            do j = 1, spec%ncoef
               if (group(i) == spec%levels(j)) then; basis(i, j) = 1.0_dp; exit; end if
            end do
         end do
      case default
         allocate(basis(0, 0)); status = 3
      end select
   end subroutine predict_smooth

   subroutine tensor_smooth(basis1, s1, basis2, s2, basis, penalties, status)
      real(dp), intent(in) :: basis1(:, :), s1(:, :), basis2(:, :), s2(:, :)
      real(dp), allocatable, intent(out) :: basis(:, :), penalties(:, :, :)
      integer, intent(out) :: status
      if (size(basis1, 1) /= size(basis2, 1) .or. &
          size(s1, 1) /= size(basis1, 2) .or. size(s2, 1) /= size(basis2, 2)) then
         allocate(basis(0, 0), penalties(0, 0, 0)); status = 1; return
      end if
      basis = tensor_product_model_matrix(basis1, basis2)
      call tensor_product_penalties(s1, s2, penalties)
      status = 0
   end subroutine tensor_smooth

   function append_columns(a, b) result(c)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable :: c(:, :)
      if (size(a, 1) /= size(b, 1)) then
         allocate(c(0, 0)); return
      end if
      allocate(c(size(a, 1), size(a, 2) + size(b, 2)))
      if (size(a, 2) > 0) c(:, 1:size(a, 2)) = a
      if (size(b, 2) > 0) c(:, size(a, 2) + 1:) = b
   end function append_columns

   function embed_penalty(s, first, total) result(full)
      real(dp), intent(in) :: s(:, :)
      integer, intent(in) :: first, total
      real(dp), allocatable :: full(:, :)
      integer :: last
      allocate(full(total, total)); full = 0.0_dp
      last = first + size(s, 1) - 1
      if (first >= 1 .and. last <= total .and. size(s, 1) == size(s, 2)) then
         full(first:last, first:last) = s
      end if
   end function embed_penalty


   subroutine place_knots(x, n_knots, knots)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: n_knots
      real(dp), allocatable, intent(out) :: knots(:)
      integer :: i
      if (n_knots <= 0 .or. size(x) == 0) then
         allocate(knots(0)); return
      end if
      allocate(knots(n_knots))
      if (n_knots == 1) then
         knots(1) = type7_quantile(x, 0.5_dp)
      else
         do i = 1, n_knots
            knots(i) = type7_quantile(x, real(i - 1, dp) / real(n_knots - 1, dp))
         end do
      end if
   end subroutine place_knots

   function identity_penalty(n) result(s)
      integer, intent(in) :: n
      real(dp), allocatable :: s(:, :)
      integer :: i
      allocate(s(n, n)); s = 0.0_dp
      do i = 1, n; s(i, i) = 1.0_dp; end do
   end function identity_penalty

end module mgcv_smooths
