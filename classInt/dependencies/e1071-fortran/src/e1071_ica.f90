module e1071_ica
    use e1071_kinds, only: dp
    use e1071_rng, only: rng_state, rng_normal
    implicit none
    private

    integer, parameter, public :: ica_negative_kurtosis = 1
    integer, parameter, public :: ica_positive_kurtosis = 2
    integer, parameter, public :: ica_fourth_moment = 3

    type, public :: ica_model
        real(dp), allocatable :: weights(:, :)
        real(dp), allocatable :: projection(:, :)
        real(dp), allocatable :: initial_weights(:, :)
        integer :: epochs = 0
        real(dp) :: learning_rate = 0.0_dp
        integer :: nonlinearity = ica_negative_kurtosis
    end type ica_model

    public :: ica_fit, ica_transform

contains

    subroutine ica_fit(x, learning_rate, rng, model, epochs, ncomp, nonlinearity, initial_weights)
        real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric matrix processed in its supplied row order.
        real(dp), intent(in) :: learning_rate !! Fixed stochastic ICA learning rate used for every observation and epoch.
        type(rng_state), intent(inout) :: rng !! Mutable random generator used only when initial_weights is not supplied.
        type(ica_model), intent(out) :: model !! Fitted weight matrix, projected observations, initial weights, and training
        !! metadata.
        integer, intent(in), optional :: epochs !! Number of full stochastic training passes; defaults to 100 and must be positive.
        integer, intent(in), optional :: ncomp !! Number of independent components; defaults to the number of input variables.
        integer, intent(in), optional :: nonlinearity !! Built-in score code: negative kurtosis, positive kurtosis, or fourth
        !! moment.
        real(dp), intent(in), optional :: initial_weights(:, :) !! Optional ncomp-by-nvar starting matrix replacing
        !! random-normal initialization.
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: gy(:)
        real(dp), allocatable :: correction(:)
        integer :: use_epochs
        integer :: components
        integer :: fun
        integer :: i
        integer :: j
        integer :: c
        integer :: epoch

        use_epochs = 100
        if (present(epochs)) use_epochs = epochs
        components = size(x, 2)
        if (present(ncomp)) components = ncomp
        fun = ica_negative_kurtosis
        if (present(nonlinearity)) fun = nonlinearity
        if (use_epochs < 1 .or. components < 1) error stop "ica_fit: invalid epochs or ncomp"
        if (fun < 1 .or. fun > 3) error stop "ica_fit: invalid nonlinearity"
        if (learning_rate <= 0.0_dp) error stop "ica_fit: learning_rate must be positive"

        allocate(model%weights(components, size(x, 2)))
        if (present(initial_weights)) then
            if (size(initial_weights, 1) /= components .or. size(initial_weights, 2) /= size(x, 2)) then
                error stop "ica_fit: initial_weights shape mismatch"
            end if
            model%weights = initial_weights
        else
            do c = 1, components
                do j = 1, size(x, 2)
                    model%weights(c, j) = rng_normal(rng)
                end do
            end do
        end if
        model%initial_weights = model%weights
        allocate(y(components), gy(components), correction(size(x, 2)))
        do epoch = 1, use_epochs
            do i = 1, size(x, 1)
                y = matmul(model%weights, x(i, :))
                select case (fun)
                case (ica_negative_kurtosis)
                    gy = tanh(y)
                case (ica_positive_kurtosis)
                    gy = y - tanh(y)
                case (ica_fourth_moment)
                    gy = sign(1.0_dp, y) * y**2
                    where (abs(y) <= 0.0_dp) gy = 0.0_dp
                end select
                correction = x(i, :) - matmul(gy, model%weights)
                do c = 1, components
                    model%weights(c, :) = model%weights(c, :) + learning_rate * gy(c) * correction
                end do
            end do
        end do
        allocate(model%projection(size(x, 1), components))
        model%projection = matmul(x, transpose(model%weights))
        model%epochs = use_epochs
        model%learning_rate = learning_rate
        model%nonlinearity = fun
    end subroutine ica_fit

    function ica_transform(model, x) result(projection)
        type(ica_model), intent(in) :: model !! Fitted ICA model supplying the learned component weight matrix.
        real(dp), intent(in) :: x(:, :) !! New observation-by-variable matrix with the same variable count as the training data.
        real(dp), allocatable :: projection(:, :)

        if (size(x, 2) /= size(model%weights, 2)) error stop "ica_transform: variable count mismatch"
        allocate(projection(size(x, 1), size(model%weights, 1)))
        projection = matmul(x, transpose(model%weights))
    end function ica_transform

end module e1071_ica
