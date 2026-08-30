module e1071_discrete
    use e1071_kinds, only: dp
    use e1071_rng, only: rng_state, rng_uniform
    implicit none
    private

    public :: ddiscrete, pdiscrete, qdiscrete, rdiscrete

contains

    function ddiscrete(x, probs, values) result(density)
        real(dp), intent(in) :: x(:) !! Points at which the discrete probability mass function is evaluated.
        real(dp), intent(in) :: probs(:) !! Nonnegative unnormalized masses, one for each entry of values.
        real(dp), intent(in) :: values(:) !! Support values associated one-to-one with probs.
        real(dp), allocatable :: density(:)
        real(dp) :: total
        integer :: i
        integer :: j

        call validate_distribution(probs, values)
        total = sum(probs)
        allocate(density(size(x)))
        density = 0.0_dp
        do i = 1, size(x)
            do j = 1, size(values)
                if (abs(x(i) - values(j)) <= 0.0_dp) density(i) = density(i) + probs(j) / total
            end do
        end do
    end function ddiscrete

    function pdiscrete(q, probs, values) result(probability)
        real(dp), intent(in) :: q(:) !! Quantiles at which the cumulative discrete probability is evaluated.
        real(dp), intent(in) :: probs(:) !! Nonnegative unnormalized masses, one for each entry of values.
        real(dp), intent(in) :: values(:) !! Support values used for the less-than-or-equal cumulative comparison.
        real(dp), allocatable :: probability(:)
        real(dp) :: total
        integer :: i

        call validate_distribution(probs, values)
        total = sum(probs)
        allocate(probability(size(q)))
        do i = 1, size(q)
            probability(i) = sum(probs, mask=values <= q(i)) / total
        end do
    end function pdiscrete

    function qdiscrete(p, probs, values) result(quantile)
        real(dp), intent(in) :: p(:) !! Probabilities to invert; values outside [0,1] follow the finite-support edge convention.
        real(dp), intent(in) :: probs(:) !! Nonnegative unnormalized masses, one for each entry of values.
        real(dp), intent(in) :: values(:) !! Support values returned according to cumulative probability order.
        real(dp), allocatable :: quantile(:)
        real(dp), allocatable :: cumulative(:)
        integer :: i
        integer :: j

        call validate_distribution(probs, values)
        allocate(cumulative(size(probs)), quantile(size(p)))
        cumulative = probs
        do j = 2, size(cumulative)
            cumulative(j) = cumulative(j) + cumulative(j - 1)
        end do
        cumulative = cumulative / cumulative(size(cumulative))
        do i = 1, size(p)
            j = 1
            do while (j <= size(cumulative))
                if (p(i) <= cumulative(j)) exit
                j = j + 1
            end do
            if (j > size(values)) j = size(values)
            quantile(i) = values(j)
        end do
    end function qdiscrete

    function rdiscrete(n, probs, values, rng) result(sample)
        integer, intent(in) :: n !! Number of independent draws; must be nonnegative.
        real(dp), intent(in) :: probs(:) !! Nonnegative unnormalized sampling masses, one for each support value.
        real(dp), intent(in) :: values(:) !! Support values returned by the generated draws.
        type(rng_state), intent(inout) :: rng !! Mutable pseudo-random-number generator state used for sampling.
        real(dp), allocatable :: sample(:)
        real(dp), allocatable :: cumulative(:)
        real(dp) :: u
        integer :: i
        integer :: j

        if (n < 0) error stop "rdiscrete: n must be nonnegative"
        call validate_distribution(probs, values)
        allocate(cumulative(size(probs)), sample(n))
        cumulative = probs
        do j = 2, size(cumulative)
            cumulative(j) = cumulative(j) + cumulative(j - 1)
        end do
        cumulative = cumulative / cumulative(size(cumulative))
        do i = 1, n
            u = rng_uniform(rng)
            j = 1
            do while (j < size(cumulative))
                if (u <= cumulative(j)) exit
                j = j + 1
            end do
            sample(i) = values(j)
        end do
    end function rdiscrete

    subroutine validate_distribution(probs, values)
        real(dp), intent(in) :: probs(:) !! Candidate nonnegative masses to validate.
        real(dp), intent(in) :: values(:) !! Support vector whose length must match probs.

        if (size(probs) /= size(values)) error stop "discrete distribution: probs and values have different lengths"
        if (size(probs) == 0) error stop "discrete distribution: empty support"
        if (any(probs < 0.0_dp)) error stop "discrete distribution: negative probability"
        if (sum(probs) <= 0.0_dp) error stop "discrete distribution: total probability must be positive"
    end subroutine validate_distribution

end module e1071_discrete
