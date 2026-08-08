! SPDX-License-Identifier: GPL-2.0-only
module soma_types
    use soma_kinds, only : dp
    implicit none
    private

    character(len=*), parameter, public :: strategy_all2one = "all2one"
    character(len=*), parameter, public :: strategy_t3a = "t3a"
    character(len=*), parameter, public :: strategy_pareto = "pareto"

    type, public :: soma_bounds
        real(dp), allocatable :: lower(:)
        real(dp), allocatable :: upper(:)
    end type soma_bounds

    type, public :: soma_options
        character(len=8) :: strategy = strategy_all2one
        integer :: population_size = 10
        integer :: n_migrations = 20
        real(dp) :: path_length = 3.0_dp
        real(dp) :: step_length = 0.11_dp
        real(dp) :: perturbation_chance = 0.1_dp
        real(dp) :: min_absolute_sep = 0.0_dp
        real(dp) :: min_relative_sep = 1.0e-3_dp
        integer :: n_steps = 0
        integer :: migrant_pool_size = 0
        integer :: leader_pool_size = 0
        integer :: n_migrants = 0
        real(dp) :: perturbation_frequency = 1.0_dp
        real(dp) :: step_frequency = 1.0_dp
    end type soma_options

    type, public :: soma_result
        integer :: leader = 0
        real(dp), allocatable :: population(:,:)
        real(dp), allocatable :: cost(:)
        real(dp), allocatable :: history(:)
        integer, allocatable :: evaluations(:)
        integer :: migrations = 0
        integer :: status = 0
        character(len=160) :: message = ""
    end type soma_result

    abstract interface
        function soma_cost_function(x) result(value)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function soma_cost_function
    end interface

    public :: soma_cost_function
    public :: make_bounds, bounds, all2one, t3a, pareto

    interface bounds
        module procedure make_bounds
    end interface bounds

contains

    function make_bounds(lower, upper) result(bnds)
        real(dp), intent(in) :: lower(:), upper(:)
        type(soma_bounds) :: bnds

        allocate(bnds%lower(size(lower)), bnds%upper(size(upper)))
        bnds%lower = lower
        bnds%upper = upper
    end function make_bounds

    function all2one(population_size, n_migrations, path_length, step_length, &
                     perturbation_chance, min_absolute_sep, min_relative_sep) result(options)
        integer, intent(in), optional :: population_size, n_migrations
        real(dp), intent(in), optional :: path_length, step_length, perturbation_chance
        real(dp), intent(in), optional :: min_absolute_sep, min_relative_sep
        type(soma_options) :: options

        options%strategy = strategy_all2one
        options%population_size = 10
        options%n_migrations = 20
        options%path_length = 3.0_dp
        options%step_length = 0.11_dp
        options%perturbation_chance = 0.1_dp
        options%min_absolute_sep = 0.0_dp
        options%min_relative_sep = 1.0e-3_dp
        options%n_steps = 0
        options%migrant_pool_size = 0
        options%leader_pool_size = 0
        options%n_migrants = 0
        if (present(population_size)) options%population_size = population_size
        if (present(n_migrations)) options%n_migrations = n_migrations
        if (present(path_length)) options%path_length = path_length
        if (present(step_length)) options%step_length = step_length
        if (present(perturbation_chance)) options%perturbation_chance = perturbation_chance
        if (present(min_absolute_sep)) options%min_absolute_sep = min_absolute_sep
        if (present(min_relative_sep)) options%min_relative_sep = min_relative_sep
    end function all2one

    function t3a(population_size, n_migrations, n_steps, migrant_pool_size, &
                 leader_pool_size, n_migrants, min_absolute_sep, min_relative_sep) result(options)
        integer, intent(in), optional :: population_size, n_migrations, n_steps
        integer, intent(in), optional :: migrant_pool_size, leader_pool_size, n_migrants
        real(dp), intent(in), optional :: min_absolute_sep, min_relative_sep
        type(soma_options) :: options

        options%strategy = strategy_t3a
        options%population_size = 30
        options%n_migrations = 20
        options%n_steps = 45
        options%migrant_pool_size = 10
        options%leader_pool_size = 10
        options%n_migrants = 4
        options%min_absolute_sep = 0.0_dp
        options%min_relative_sep = 1.0e-3_dp
        if (present(population_size)) options%population_size = population_size
        if (present(n_migrations)) options%n_migrations = n_migrations
        if (present(n_steps)) options%n_steps = n_steps
        if (present(migrant_pool_size)) options%migrant_pool_size = migrant_pool_size
        if (present(leader_pool_size)) options%leader_pool_size = leader_pool_size
        if (present(n_migrants)) options%n_migrants = n_migrants
        if (present(min_absolute_sep)) options%min_absolute_sep = min_absolute_sep
        if (present(min_relative_sep)) options%min_relative_sep = min_relative_sep
    end function t3a

    function pareto(population_size, n_migrations, n_steps, perturbation_frequency, &
                    step_frequency, min_absolute_sep, min_relative_sep) result(options)
        integer, intent(in), optional :: population_size, n_migrations, n_steps
        real(dp), intent(in), optional :: perturbation_frequency, step_frequency
        real(dp), intent(in), optional :: min_absolute_sep, min_relative_sep
        type(soma_options) :: options

        options%strategy = strategy_pareto
        options%population_size = 100
        options%n_migrations = 20
        options%n_steps = 10
        options%perturbation_frequency = 1.0_dp
        options%step_frequency = 1.0_dp
        options%min_absolute_sep = 0.0_dp
        options%min_relative_sep = 1.0e-3_dp
        if (present(population_size)) options%population_size = population_size
        if (present(n_migrations)) options%n_migrations = n_migrations
        if (present(n_steps)) options%n_steps = n_steps
        if (present(perturbation_frequency)) options%perturbation_frequency = perturbation_frequency
        if (present(step_frequency)) options%step_frequency = step_frequency
        if (present(min_absolute_sep)) options%min_absolute_sep = min_absolute_sep
        if (present(min_relative_sep)) options%min_relative_sep = min_relative_sep
    end function pareto

end module soma_types
