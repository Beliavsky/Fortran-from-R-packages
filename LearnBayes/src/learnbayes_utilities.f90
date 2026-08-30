module learnbayes_utilities
   use learnbayes_distributions, only: beta_logpdf
   use learnbayes_kinds, only: dp
   use learnbayes_sampling, only: eval_callback
   use learnbayes_types, only: bayes_discrete_result, bayes_grid_result, likelihood_callback, log_density_callback
   implicit none
   private

   public :: careertraj_setup
   public :: contour_grid
   public :: discrete_bayes
   public :: discrete_bayes_2
   public :: regroup
   public :: triplot_data

contains

   subroutine discrete_bayes(callback, parameter, prior, observations, result)
      type(likelihood_callback), intent(in) :: callback !! Bound likelihood callback plus optional numeric context.
      real(dp), intent(in) :: parameter(:) !! Discrete one-parameter support values.
      real(dp), intent(in) :: prior(:) !! Prior masses corresponding one-for-one with parameter.
      real(dp), intent(in) :: observations(:, :) !! Observation matrix; each row contributes independently to the likelihood.
      type(bayes_discrete_result), intent(out) :: result !! Posterior masses and prior predictive probability of the observations.
      real(dp), allocatable :: likelihood(:)
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: params(:)
      integer :: i
      integer :: j

      allocate(result%prob(size(parameter)), likelihood(size(parameter)))
      call callback_context(callback, data, params)
      likelihood = 1.0_dp
      do j = 1, size(parameter)
         do i = 1, size(observations, 1)
            likelihood(j) = likelihood(j)*callback%eval(observations(i, :), parameter(j), 0.0_dp, data, params)
         end do
      end do
      result%predictive = sum(prior*likelihood)
      result%prob = prior*likelihood/result%predictive
   end subroutine discrete_bayes

   subroutine discrete_bayes_2(callback, parameter1, parameter2, prior, observations, result)
      type(likelihood_callback), intent(in) :: callback !! Bound two-parameter likelihood callback plus optional numeric context.
      real(dp), intent(in) :: parameter1(:) !! First discrete support grid, indexing rows of prior and posterior matrices.
      real(dp), intent(in) :: parameter2(:) !! Second discrete support grid, indexing columns of prior and posterior matrices.
      real(dp), intent(in) :: prior(:, :) !! Joint prior masses on the Cartesian product of parameter1 and parameter2.
      real(dp), intent(in) :: observations(:, :) !! Observation matrix; zero rows represent an empty dataset and likelihood one.
      type(bayes_grid_result), intent(out) :: result !! Joint posterior masses and prior predictive probability.
      real(dp), allocatable :: likelihood(:, :)
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: params(:)
      integer :: i
      integer :: j
      integer :: k

      allocate(result%prob(size(parameter1), size(parameter2)), likelihood(size(parameter1), size(parameter2)))
      call callback_context(callback, data, params)
      likelihood = 1.0_dp
      do j = 1, size(parameter2)
         do i = 1, size(parameter1)
            do k = 1, size(observations, 1)
               likelihood(i, j) = likelihood(i, j)* &
                  callback%eval(observations(k, :), parameter1(i), parameter2(j), data, params)
            end do
         end do
      end do
      result%predictive = sum(prior*likelihood)
      result%prob = prior*likelihood/result%predictive
   end subroutine discrete_bayes_2

   subroutine callback_context(callback, data, params)
      type(likelihood_callback), intent(in) :: callback !! Likelihood callback whose stored numeric context is copied.
      real(dp), allocatable, intent(out) :: data(:, :) !! Copy of callback data, or a zero-by-zero matrix if absent.
      real(dp), allocatable, intent(out) :: params(:) !! Copy of callback constants, or a zero-length vector if absent.

      if (allocated(callback%data)) then
         allocate(data(size(callback%data, 1), size(callback%data, 2)))
         data = callback%data
      else
         allocate(data(0, 0))
      end if
      if (allocated(callback%params)) then
         allocate(params(size(callback%params)))
         params = callback%params
      else
         allocate(params(0))
      end if
   end subroutine callback_context

   subroutine regroup(data, group_size, grouped)
      real(dp), intent(in) :: data(:, :) !! Input matrix whose consecutive rows are accumulated into groups.
      integer, intent(in) :: group_size !! Positive number of source rows per group, matching LearnBayes regroup().
      real(dp), allocatable, intent(out) :: grouped(:, :) !! Row-summed grouped matrix; final remainders join the last full group.
      integer :: n
      integer :: n_group
      integer :: i
      integer :: j
      integer :: g

      n = size(data, 1)
      if (group_size <= 0 .or. n < group_size) then
         allocate(grouped(0, size(data, 2)))
         return
      end if
      n_group = n/group_size
      allocate(grouped(n_group, size(data, 2)))
      grouped = 0.0_dp
      do g = 1, n_group
         do i = 1, group_size
            j = (g - 1)*group_size + i
            grouped(g, :) = grouped(g, :) + data(j, :)
         end do
      end do
      if (n > n_group*group_size) then
         do j = n_group*group_size + 1, n
            grouped(n_group, :) = grouped(n_group, :) + data(j, :)
         end do
      end if
   end subroutine regroup

   subroutine careertraj_setup(player, data, player_ids, y, n_ab, age, seasons, info)
      integer, intent(in) :: player(:) !! Integer player identifier for each input row, replacing R character factor labels.
      real(dp), intent(in) :: data(:, :) !! Career data matrix using LearnBayes columns 3, 5, 10, and 13.
      integer, allocatable, intent(out) :: player_ids(:) !! Sorted-by-first-occurrence unique player identifiers.
      real(dp), allocatable, intent(out) :: y(:, :) !! Career outcome from source column 10, padded with zeros.
      real(dp), allocatable, intent(out) :: n_ab(:, :) !! At-bat quantity source column 5 minus source column 13, padded with zeros.
      real(dp), allocatable, intent(out) :: age(:, :) !! Source column 3 for each player's career rows, padded with zeros.
      integer, allocatable, intent(out) :: seasons(:) !! Count of rows with positive n_ab for each unique player.
      integer, intent(out) :: info !! Zero on success; nonzero if dimensions or required source columns are missing.
      integer, allocatable :: count_by_player(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: nplayer
      integer :: max_rows
      logical :: seen

      info = 0
      if (size(player) /= size(data, 1) .or. size(data, 2) < 13) then
         info = 1
         allocate(player_ids(0), y(0, 0), n_ab(0, 0), age(0, 0), seasons(0))
         return
      end if
      allocate(player_ids(size(player)))
      nplayer = 0
      do i = 1, size(player)
         seen = .false.
         do j = 1, nplayer
            if (player_ids(j) == player(i)) then
               seen = .true.
               exit
            end if
         end do
         if (.not. seen) then
            nplayer = nplayer + 1
            player_ids(nplayer) = player(i)
         end if
      end do
      player_ids = player_ids(1:nplayer)
      allocate(count_by_player(nplayer))
      count_by_player = 0
      do i = 1, size(player)
         do j = 1, nplayer
            if (player(i) == player_ids(j)) then
               count_by_player(j) = count_by_player(j) + 1
               exit
            end if
         end do
      end do
      max_rows = maxval(count_by_player)
      allocate(y(nplayer, max_rows), n_ab(nplayer, max_rows), age(nplayer, max_rows), seasons(nplayer))
      y = 0.0_dp
      n_ab = 0.0_dp
      age = 0.0_dp
      seasons = 0
      count_by_player = 0
      do i = 1, size(player)
         do j = 1, nplayer
            if (player(i) == player_ids(j)) then
               count_by_player(j) = count_by_player(j) + 1
               k = count_by_player(j)
               y(j, k) = data(i, 10)
               n_ab(j, k) = data(i, 5) - data(i, 13)
               age(j, k) = data(i, 3)
               if (n_ab(j, k) > 0.0_dp) seasons(j) = seasons(j) + 1
               exit
            end if
         end do
      end do
   end subroutine careertraj_setup

   subroutine contour_grid(logf, limits, ng, x, y, z)
      type(log_density_callback), intent(in) :: logf !! Two-parameter log density evaluated on a rectangular grid.
      real(dp), intent(in) :: limits(4) !! Rectangle xmin, xmax, ymin, ymax.
      integer, intent(in) :: ng !! Number of equally spaced grid points along each coordinate.
      real(dp), intent(out) :: x(:) !! First-coordinate grid of length ng.
      real(dp), intent(out) :: y(:) !! Second-coordinate grid of length ng.
      real(dp), intent(out) :: z(:, :) !! Shifted log-density grid with maximum equal to zero.
      real(dp) :: theta(2)
      real(dp) :: mx
      integer :: i
      integer :: j

      do i = 1, ng
         x(i) = limits(1) + real(i - 1, dp)*(limits(2) - limits(1))/real(ng - 1, dp)
         y(i) = limits(3) + real(i - 1, dp)*(limits(4) - limits(3))/real(ng - 1, dp)
      end do
      do j = 1, ng
         theta(2) = y(j)
         do i = 1, ng
            theta(1) = x(i)
            z(i, j) = eval_callback(logf, theta)
         end do
      end do
      mx = maxval(z)
      z = z - mx
   end subroutine contour_grid

   subroutine triplot_data(a, b, successes, failures, p, prior_density, likelihood_density, posterior_density)
      real(dp), intent(in) :: a !! Positive first beta prior shape parameter.
      real(dp), intent(in) :: b !! Positive second beta prior shape parameter.
      integer, intent(in) :: successes !! Observed binomial successes.
      integer, intent(in) :: failures !! Observed binomial failures.
      real(dp), intent(out) :: p(:) !! Probability grid from 0.005 to 0.995.
      real(dp), intent(out) :: prior_density(:) !! Beta(a,b) prior density on p.
      real(dp), intent(out) :: likelihood_density(:) !! Normalized beta(successes+1,failures+1) likelihood shape.
      real(dp), intent(out) :: posterior_density(:) !! Beta(a+successes,b+failures) posterior density on p.
      integer :: i
      integer :: n

      n = size(p)
      do i = 1, n
         p(i) = 0.005_dp + 0.99_dp*real(i - 1, dp)/real(n - 1, dp)
         prior_density(i) = exp(beta_logpdf(p(i), a, b))
         likelihood_density(i) = exp(beta_logpdf(p(i), real(successes + 1, dp), real(failures + 1, dp)))
         posterior_density(i) = exp(beta_logpdf(p(i), a + real(successes, dp), b + real(failures, dp)))
      end do
   end subroutine triplot_data

end module learnbayes_utilities
