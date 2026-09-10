! SPDX-License-Identifier: GPL-2.0-or-later
module dlm_types
    use r_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: dlm_success = 0
    integer, parameter, public :: dlm_invalid_shape = 1
    integer, parameter, public :: dlm_invalid_model = 2
    integer, parameter, public :: dlm_linalg_failure = 3
    integer, parameter, public :: dlm_invalid_argument = 4
    integer, parameter, public :: dlm_sampling_failure = 5

    type, public :: dlm_model
        real(dp), allocatable :: m0(:)
        real(dp), allocatable :: c0(:, :)
        real(dp), allocatable :: ff(:, :)
        real(dp), allocatable :: v(:, :)
        real(dp), allocatable :: gg(:, :)
        real(dp), allocatable :: w(:, :)
        integer, allocatable :: jff(:, :)
        integer, allocatable :: jv(:, :)
        integer, allocatable :: jgg(:, :)
        integer, allocatable :: jw(:, :)
        real(dp), allocatable :: x(:, :)
    end type dlm_model

    type, public :: dlm_filter_result
        type(dlm_model) :: model
        real(dp), allocatable :: y(:, :)
        real(dp), allocatable :: m(:, :)
        real(dp), allocatable :: c(:, :, :)
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: f(:, :)
        real(dp), allocatable :: q(:, :, :)
    end type dlm_filter_result

    type, public :: dlm_smooth_result
        real(dp), allocatable :: s(:, :)
        real(dp), allocatable :: covariance(:, :, :)
    end type dlm_smooth_result

    type, public :: dlm_forecast_result
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: f(:, :)
        real(dp), allocatable :: q(:, :, :)
    end type dlm_forecast_result

    type, public :: dlm_gibbs_dig_result
        real(dp), allocatable :: d_v(:)
        real(dp), allocatable :: d_w(:, :)
        real(dp), allocatable :: theta(:, :, :)
    end type dlm_gibbs_dig_result

    public :: dp
    public :: dlm_model_is_valid
    public :: dlm_model_is_time_varying
    public :: dlm_matrices_at

contains

    pure logical function dlm_model_is_valid(model) result(valid)
        type(dlm_model), intent(in) :: model !! Dynamic linear model whose allocated components and shapes are checked.

        integer :: m_obs, p_state

        valid = .false.
        if (.not. allocated(model%m0)) return
        if (.not. allocated(model%c0)) return
        if (.not. allocated(model%ff)) return
        if (.not. allocated(model%v)) return
        if (.not. allocated(model%gg)) return
        if (.not. allocated(model%w)) return

        p_state = size(model%m0)
        m_obs = size(model%ff, 1)
        if (p_state < 1 .or. m_obs < 1) return
        if (size(model%ff, 2) /= p_state) return
        if (size(model%c0, 1) /= p_state .or. size(model%c0, 2) /= p_state) return
        if (size(model%v, 1) /= m_obs .or. size(model%v, 2) /= m_obs) return
        if (size(model%gg, 1) /= p_state .or. size(model%gg, 2) /= p_state) return
        if (size(model%w, 1) /= p_state .or. size(model%w, 2) /= p_state) return

        if (allocated(model%jff)) then
            if (any(shape(model%jff) /= shape(model%ff))) return
            if (.not. allocated(model%x)) return
        end if
        if (allocated(model%jv)) then
            if (any(shape(model%jv) /= shape(model%v))) return
            if (.not. allocated(model%x)) return
        end if
        if (allocated(model%jgg)) then
            if (any(shape(model%jgg) /= shape(model%gg))) return
            if (.not. allocated(model%x)) return
        end if
        if (allocated(model%jw)) then
            if (any(shape(model%jw) /= shape(model%w))) return
            if (.not. allocated(model%x)) return
        end if
        valid = .true.
    end function dlm_model_is_valid

    pure logical function dlm_model_is_time_varying(model) result(time_varying)
        type(dlm_model), intent(in) :: model !! Dynamic linear model tested for any allocated time-varying index matrix.

        time_varying = allocated(model%jff) .or. allocated(model%jv) .or. &
                       allocated(model%jgg) .or. allocated(model%jw)
    end function dlm_model_is_time_varying

    pure subroutine dlm_matrices_at(model, time_index, ff, v, gg, w, info)
        type(dlm_model), intent(in) :: model !! Dynamic linear model supplying baseline matrices and optional index maps.
        integer, intent(in) :: time_index !! One-based row of `model%x` used for time-varying substitutions.
        real(dp), allocatable, intent(out) :: ff(:, :) !! Observation matrix after substitutions at `time_index`.
        real(dp), allocatable, intent(out) :: v(:, :) !! Observation covariance after substitutions at `time_index`.
        real(dp), allocatable, intent(out) :: gg(:, :) !! State transition matrix after substitutions at `time_index`.
        real(dp), allocatable, intent(out) :: w(:, :) !! State innovation covariance after substitutions at `time_index`.
        integer, intent(out) :: info !! Zero on success or a `dlm_*` status code on invalid model/index data.

        integer :: i, j, k

        info = dlm_success
        if (.not. dlm_model_is_valid(model)) then
            info = dlm_invalid_model
            allocate(ff(0, 0), v(0, 0), gg(0, 0), w(0, 0))
            return
        end if

        ff = model%ff
        v = model%v
        gg = model%gg
        w = model%w

        if (.not. dlm_model_is_time_varying(model)) return
        if (.not. allocated(model%x)) then
            info = dlm_invalid_model
            return
        end if
        if (time_index < 1 .or. time_index > size(model%x, 1)) then
            info = dlm_invalid_argument
            return
        end if

        if (allocated(model%jff)) then
            do j = 1, size(model%jff, 2)
                do i = 1, size(model%jff, 1)
                    k = model%jff(i, j)
                    if (k > 0) then
                        if (k > size(model%x, 2)) then
                            info = dlm_invalid_model
                            return
                        end if
                        ff(i, j) = model%x(time_index, k)
                    end if
                end do
            end do
        end if

        if (allocated(model%jv)) then
            do j = 1, size(model%jv, 2)
                do i = 1, size(model%jv, 1)
                    k = model%jv(i, j)
                    if (k > 0) then
                        if (k > size(model%x, 2)) then
                            info = dlm_invalid_model
                            return
                        end if
                        v(i, j) = model%x(time_index, k)
                    end if
                end do
            end do
        end if

        if (allocated(model%jgg)) then
            do j = 1, size(model%jgg, 2)
                do i = 1, size(model%jgg, 1)
                    k = model%jgg(i, j)
                    if (k > 0) then
                        if (k > size(model%x, 2)) then
                            info = dlm_invalid_model
                            return
                        end if
                        gg(i, j) = model%x(time_index, k)
                    end if
                end do
            end do
        end if

        if (allocated(model%jw)) then
            do j = 1, size(model%jw, 2)
                do i = 1, size(model%jw, 1)
                    k = model%jw(i, j)
                    if (k > 0) then
                        if (k > size(model%x, 2)) then
                            info = dlm_invalid_model
                            return
                        end if
                        w(i, j) = model%x(time_index, k)
                    end if
                end do
            end do
        end if
    end subroutine dlm_matrices_at

end module dlm_types
