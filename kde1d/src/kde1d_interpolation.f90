! SPDX-License-Identifier: MIT
! Derived from kde1d-cpp interpolation code by Thomas Nagler and Thibault Vatter.
! Upstream code and this translation are distributed under the MIT license.
module kde1d_interpolation
   use r_kinds, only : dp
   implicit none
   private

   type, public :: interpolation_grid
      real(dp), allocatable :: grid_points(:)
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: cell_coefs(:, :)
      real(dp), allocatable :: cumulative_integrals(:)
   contains
      procedure :: initialize => interpolation_initialize
      procedure :: interpolate => interpolation_interpolate
      procedure :: integrate => interpolation_integrate
      procedure :: quantile => interpolation_quantile
      procedure :: normalize => interpolation_normalize
   end type interpolation_grid

contains

   pure subroutine interpolation_initialize(self, grid_points, values, norm_times, ierr)
      class(interpolation_grid), intent(inout) :: self !! Interpolation object to initialize.
      real(dp), intent(in) :: grid_points(:) !! Strictly increasing interpolation abscissae.
      real(dp), intent(in) :: values(:) !! Density ordinates corresponding to grid_points.
      integer, intent(in) :: norm_times !! Number of spline-mass normalization passes, normally 0 or 3.
      integer, intent(out) :: ierr !! Zero on success; nonzero for invalid grid input.

      ierr = 0
      if (size(grid_points) /= size(values) .or. size(grid_points) < 2) then
         ierr = 1
         return
      end if
      if (any(grid_points(2:) <= grid_points(:size(grid_points) - 1))) then
         ierr = 2
         return
      end if

      self%grid_points = grid_points
      self%values = values
      call self%normalize(norm_times, ierr)
   end subroutine interpolation_initialize

   pure subroutine interpolation_normalize(self, times, ierr)
      class(interpolation_grid), intent(inout) :: self !! Interpolation object whose ordinates are normalized.
      integer, intent(in) :: times !! Number of normalization iterations to apply.
      integer, intent(out) :: ierr !! Zero on success; nonzero when total spline mass is invalid.

      integer :: iteration
      integer :: k
      real(dp) :: integral

      ierr = 0
      if (.not. allocated(self%grid_points) .or. .not. allocated(self%values)) then
         ierr = 1
         return
      end if

      do iteration = 1, max(0, times)
         call update_cell_coefs(self)
         integral = 0.0_dp
         do k = 1, size(self%grid_points) - 1
            integral = integral + cubic_integral(0.0_dp, 1.0_dp, self%cell_coefs(:, k)) &
               * (self%grid_points(k + 1) - self%grid_points(k))
         end do
         if (.not. (integral > 0.0_dp)) then
            ierr = 2
            return
         end if
         self%values = self%values / integral
      end do

      call update_cell_coefs(self)
      call update_cumulative_integrals(self)
   end subroutine interpolation_normalize

   pure function interpolation_interpolate(self, x) result(y)
      class(interpolation_grid), intent(in) :: self !! Initialized spline interpolation object.
      real(dp), intent(in) :: x(:) !! Evaluation points; NaNs are propagated by arithmetic comparisons.
      real(dp) :: y(size(x))

      integer :: i
      integer :: k
      real(dp) :: xev
      real(dp) :: width

      do i = 1, size(x)
         k = find_cell(self, x(i))
         width = self%grid_points(k + 1) - self%grid_points(k)
         xev = (x(i) - self%grid_points(k)) / width
         if (xev <= 0.0_dp) then
            y(i) = self%values(k) * exp(-0.5_dp * xev * xev)
         else if (xev >= 1.0_dp) then
            y(i) = self%values(k + 1) * exp(-0.5_dp * (xev - 1.0_dp)**2)
         else
            y(i) = cubic_poly(xev, self%cell_coefs(:, k))
         end if
      end do
   end function interpolation_interpolate

   pure function interpolation_integrate(self, x, normalize) result(y)
      class(interpolation_grid), intent(in) :: self !! Initialized spline interpolation object.
      real(dp), intent(in) :: x(:) !! Upper limits of integration on the original grid scale.
      logical, intent(in), optional :: normalize !! If true, divide by the spline's total grid mass.
      real(dp) :: y(size(x))

      integer :: i
      integer :: k
      logical :: do_normalize
      real(dp) :: position
      real(dp) :: total_mass
      real(dp) :: width

      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize
      total_mass = self%cumulative_integrals(size(self%cumulative_integrals))

      do i = 1, size(x)
         if (x(i) <= self%grid_points(1)) then
            y(i) = 0.0_dp
         else if (x(i) >= self%grid_points(size(self%grid_points))) then
            y(i) = total_mass
         else
            k = find_cell(self, x(i))
            width = self%grid_points(k + 1) - self%grid_points(k)
            position = (x(i) - self%grid_points(k)) / width
            y(i) = self%cumulative_integrals(k) &
               + cubic_integral(0.0_dp, position, self%cell_coefs(:, k)) * width
         end if
      end do

      if (do_normalize .and. total_mass > 0.0_dp) y = y / total_mass
   end function interpolation_integrate

   pure function interpolation_quantile(self, probabilities) result(q)
      class(interpolation_grid), intent(in) :: self !! Initialized spline interpolation object.
      real(dp), intent(in) :: probabilities(:) !! Probabilities in [0,1]; endpoints map to grid endpoints.
      real(dp) :: q(size(probabilities))

      integer :: i
      real(dp) :: total_mass

      total_mass = self%cumulative_integrals(size(self%cumulative_integrals))
      do i = 1, size(probabilities)
         q(i) = invert_integral(self, probabilities(i), total_mass)
      end do
   end function interpolation_quantile

   pure real(dp) function cubic_poly(x, a) result(value)
      real(dp), intent(in) :: x !! Unit-cell coordinate at which the cubic is evaluated.
      real(dp), intent(in) :: a(4) !! Cubic coefficients ordered from constant through cubic terms.

      value = a(1) + a(2) * x + a(3) * x * x + a(4) * x * x * x
   end function cubic_poly

   pure real(dp) function cubic_indef_integral(x, a) result(value)
      real(dp), intent(in) :: x !! Unit-cell coordinate at which the antiderivative is evaluated.
      real(dp), intent(in) :: a(4) !! Cubic coefficients ordered from constant through cubic terms.

      value = a(1) * x + 0.5_dp * a(2) * x**2 + a(3) * x**3 / 3.0_dp &
         + 0.25_dp * a(4) * x**4
   end function cubic_indef_integral

   pure real(dp) function cubic_integral(lower, upper, a) result(value)
      real(dp), intent(in) :: lower !! Lower unit-cell integration limit.
      real(dp), intent(in) :: upper !! Upper unit-cell integration limit.
      real(dp), intent(in) :: a(4) !! Cubic coefficients ordered from constant through cubic terms.

      value = cubic_indef_integral(upper, a) - cubic_indef_integral(lower, a)
   end function cubic_integral

   pure real(dp) function invert_integral(self, probability, total_mass) result(value)
      class(interpolation_grid), intent(in) :: self !! Initialized spline interpolation object.
      real(dp), intent(in) :: probability !! Probability in [0,1] to invert.
      real(dp), intent(in) :: total_mass !! Positive total spline mass over the interpolation grid.

      integer :: cell
      integer :: iteration
      real(dp) :: cell_mass
      real(dp) :: derivative
      real(dp) :: lower_position
      real(dp) :: newton_position
      real(dp) :: position
      real(dp) :: residual
      real(dp) :: target_mass
      real(dp) :: upper_position
      real(dp) :: width

      if (probability <= 0.0_dp) then
         value = self%grid_points(1)
         return
      end if
      if (probability >= 1.0_dp) then
         value = self%grid_points(size(self%grid_points))
         return
      end if

      target_mass = probability * total_mass
      cell = 1
      do while (cell < size(self%grid_points) - 1)
         if (self%cumulative_integrals(cell + 1) >= target_mass) exit
         cell = cell + 1
      end do

      width = self%grid_points(cell + 1) - self%grid_points(cell)
      cell_mass = self%cumulative_integrals(cell + 1) - self%cumulative_integrals(cell)
      if (.not. (cell_mass > 0.0_dp)) then
         value = self%grid_points(cell + 1)
         return
      end if

      position = max(0.0_dp, min(1.0_dp, &
         (target_mass - self%cumulative_integrals(cell)) / cell_mass))
      lower_position = 0.0_dp
      upper_position = 1.0_dp

      do iteration = 1, 35
         residual = cubic_indef_integral(position, self%cell_coefs(:, cell)) * width &
            - (target_mass - self%cumulative_integrals(cell))
         if (abs(residual) <= 8.0_dp * epsilon(1.0_dp) * total_mass) exit

         if (residual < 0.0_dp) then
            lower_position = position
         else
            upper_position = position
         end if

         derivative = cubic_poly(position, self%cell_coefs(:, cell)) * width
         if (derivative > 0.0_dp) then
            newton_position = position - residual / derivative
         else
            newton_position = lower_position
         end if
         if (derivative > 0.0_dp .and. newton_position > lower_position &
            .and. newton_position < upper_position) then
            position = newton_position
         else
            position = 0.5_dp * (lower_position + upper_position)
         end if
      end do

      value = self%grid_points(cell) + position * width
   end function invert_integral

   pure integer function find_cell(self, x) result(cell)
      class(interpolation_grid), intent(in) :: self !! Initialized spline interpolation object.
      real(dp), intent(in) :: x !! Evaluation point whose containing grid cell is requested.

      integer :: low
      integer :: high
      integer :: mid

      if (x <= self%grid_points(1)) then
         cell = 1
         return
      end if
      if (x >= self%grid_points(size(self%grid_points))) then
         cell = size(self%grid_points) - 1
         return
      end if

      low = 1
      high = size(self%grid_points)
      do while (low < high - 1)
         mid = low + (high - low) / 2
         if (x < self%grid_points(mid)) then
            high = mid
         else
            low = mid
         end if
      end do
      cell = low
   end function find_cell

   pure subroutine update_cell_coefs(self)
      class(interpolation_grid), intent(inout) :: self !! Interpolation object whose cubic cell coefficients are rebuilt.

      integer :: k

      if (allocated(self%cell_coefs)) deallocate(self%cell_coefs)
      allocate(self%cell_coefs(4, size(self%grid_points) - 1))
      do k = 1, size(self%grid_points) - 1
         self%cell_coefs(:, k) = find_cell_coefs(self, k)
      end do
   end subroutine update_cell_coefs

   pure function find_cell_coefs(self, k) result(a)
      class(interpolation_grid), intent(in) :: self !! Interpolation object supplying neighboring knots and ordinates.
      integer, intent(in) :: k !! One-based cell index between adjacent grid points.
      real(dp) :: a(4)

      integer :: k0
      integer :: k2
      integer :: k3
      real(dp) :: dt0
      real(dp) :: dt1
      real(dp) :: dt2
      real(dp) :: dx1
      real(dp) :: dx2

      k0 = max(k - 1, 1)
      k2 = k + 1
      k3 = min(k + 2, size(self%grid_points))
      dt0 = self%grid_points(k) - self%grid_points(k0)
      dt1 = self%grid_points(k2) - self%grid_points(k)
      dt2 = self%grid_points(k3) - self%grid_points(k2)

      dx1 = 0.0_dp
      dx2 = 0.0_dp
      if (dt0 > 0.0_dp) then
         dx1 = (self%values(k) - self%values(k0)) / dt0
         dx1 = dx1 - (self%values(k2) - self%values(k0)) / (dt0 + dt1)
         dx1 = dx1 + (self%values(k2) - self%values(k)) / dt1
      end if
      if (dt2 > 0.0_dp) then
         dx2 = (self%values(k2) - self%values(k)) / dt1
         dx2 = dx2 - (self%values(k3) - self%values(k)) / (dt1 + dt2)
         dx2 = dx2 + (self%values(k3) - self%values(k2)) / dt2
      end if

      dx1 = dx1 * dt1
      dx2 = dx2 * dt1
      dx1 = max(dx1, -3.0_dp * self%values(k))
      dx2 = min(dx2, 3.0_dp * self%values(k2))

      a(1) = self%values(k)
      a(2) = dx1
      a(3) = -3.0_dp * (self%values(k) - self%values(k2)) - 2.0_dp * dx1 - dx2
      a(4) = 2.0_dp * (self%values(k) - self%values(k2)) + dx1 + dx2
   end function find_cell_coefs

   pure subroutine update_cumulative_integrals(self)
      class(interpolation_grid), intent(inout) :: self !! Interpolation object whose cumulative spline masses are rebuilt.

      integer :: k

      if (allocated(self%cumulative_integrals)) deallocate(self%cumulative_integrals)
      allocate(self%cumulative_integrals(size(self%grid_points)))
      self%cumulative_integrals = 0.0_dp
      do k = 1, size(self%grid_points) - 1
         self%cumulative_integrals(k + 1) = self%cumulative_integrals(k) &
            + cubic_integral(0.0_dp, 1.0_dp, self%cell_coefs(:, k)) &
            * (self%grid_points(k + 1) - self%grid_points(k))
      end do
   end subroutine update_cumulative_integrals

end module kde1d_interpolation
