! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_simulation
   use mgcv_kinds, only : dp, pi_dp
   use mgcv_distributions, only : random_normal, random_poisson
   implicit none
   private

   integer, parameter, public :: sim_normal = 1
   integer, parameter, public :: sim_poisson = 2
   integer, parameter, public :: sim_binary = 3

   type, public :: gam_sim_data_t
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: f(:)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: components(:, :)
   end type gam_sim_data_t

   public :: gam_sim

contains

   subroutine gam_sim(n, example, distribution, scale, data, status, seed)
      integer, intent(in) :: n
      integer, intent(in), optional :: example, distribution, seed
      real(dp), intent(in), optional :: scale
      type(gam_sim_data_t), intent(out) :: data
      integer, intent(out) :: status
      integer :: eg, dist, i, nseed, j
      integer, allocatable :: seed_values(:)
      real(dp) :: sigma, u, probability

      status = 0; eg = 1; dist = sim_normal; sigma = 2.0_dp
      if (present(example)) eg = example
      if (present(distribution)) dist = distribution
      if (present(scale)) sigma = scale
      if (n < 1 .or. (eg /= 1 .and. eg /= 2)) then; status = 1; return; end if
      if (present(seed)) then
         call random_seed(size=nseed); allocate(seed_values(nseed))
         do j = 1, nseed
            seed_values(j) = modulo(seed + 104729 * j, huge(1) - 1)
            if (seed_values(j) <= 0) seed_values(j) = j
         end do
         call random_seed(put=seed_values)
      end if

      select case (eg)
      case (1)
         allocate(data%x(n, 4), data%components(n, 4), data%f(n), data%y(n))
         call random_number(data%x)
         data%components(:, 1) = 2.0_dp * sin(pi_dp * data%x(:, 1))
         data%components(:, 2) = exp(2.0_dp * data%x(:, 2))
         data%components(:, 3) = 0.2_dp * data%x(:, 3)**11 * (10.0_dp * (1.0_dp - data%x(:, 3)))**6 + &
             10.0_dp * (10.0_dp * data%x(:, 3))**3 * (1.0_dp - data%x(:, 3))**10
         data%components(:, 4) = 0.0_dp
         data%f = sum(data%components, dim=2)
      case (2)
         allocate(data%x(n, 2), data%components(n, 1), data%f(n), data%y(n))
         call random_number(data%x)
         data%f = pi_dp * 0.3_dp * 0.4_dp * (1.2_dp * exp(-(data%x(:,1)-0.2_dp)**2/0.3_dp**2 - &
            (data%x(:,2)-0.3_dp)**2/0.4_dp**2) + 0.8_dp * exp(-(data%x(:,1)-0.7_dp)**2/0.3_dp**2 - &
            (data%x(:,2)-0.8_dp)**2/0.4_dp**2))
         data%components(:, 1) = data%f
      end select

      select case (dist)
      case (sim_normal)
         do i = 1, n
            data%y(i) = data%f(i) + sigma * random_normal()
         end do
      case (sim_poisson)
         data%f = sigma * data%f
         do i = 1, n
            data%y(i) = real(random_poisson(exp(min(50.0_dp, data%f(i)))), dp)
         end do
      case (sim_binary)
         data%f = sigma * (data%f - 5.0_dp)
         do i = 1, n
            probability = 1.0_dp / (1.0_dp + exp(-max(-50.0_dp, min(50.0_dp, data%f(i)))))
            call random_number(u)
            data%y(i) = merge(1.0_dp, 0.0_dp, u < probability)
         end do
      case default
         status = 2
      end select
   end subroutine gam_sim

end module mgcv_simulation
