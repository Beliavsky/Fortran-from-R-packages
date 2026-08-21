! SPDX-License-Identifier: GPL-2.0-or-later
module tensora_types
    use tensora_kinds, only : dp, axis_name_len
    implicit none
    private

    type, public :: tensor_t
        complex(dp), allocatable :: data(:)
        integer, allocatable :: shape(:)
        character(len=axis_name_len), allocatable :: axis(:)
    contains
        procedure :: rank => tensor_rank
        procedure :: nelem => tensor_nelem
        procedure :: valid => tensor_valid
        procedure :: axis_pos => tensor_axis_pos
    end type tensor_t

    interface tensor
        module procedure tensor_from_real
        module procedure tensor_from_complex
    end interface tensor

    interface to_tensor
        module procedure tensor_from_real
        module procedure tensor_from_complex
    end interface to_tensor

    interface scalar_tensor
        module procedure scalar_tensor_real
        module procedure scalar_tensor_complex
    end interface scalar_tensor

    public :: tensor, to_tensor, tensor_from_real, tensor_from_complex
    public :: tensor_zeros, tensor_ones, scalar_tensor
    public :: tensor_is_real, real_data

contains

    function tensor_from_real(data, shape, axis) result(t)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: shape(:)
        character(len=*), intent(in), optional :: axis(:)
        type(tensor_t) :: t

        call init_metadata(t, size(data), shape, axis)
        allocate(t%data(size(data)))
        t%data = cmplx(data, 0.0_dp, dp)
    end function tensor_from_real

    function tensor_from_complex(data, shape, axis) result(t)
        complex(dp), intent(in) :: data(:)
        integer, intent(in) :: shape(:)
        character(len=*), intent(in), optional :: axis(:)
        type(tensor_t) :: t

        call init_metadata(t, size(data), shape, axis)
        allocate(t%data(size(data)))
        t%data = data
    end function tensor_from_complex

    function tensor_zeros(shape, axis) result(t)
        integer, intent(in) :: shape(:)
        character(len=*), intent(in), optional :: axis(:)
        type(tensor_t) :: t
        integer :: n

        n = product_or_one(shape)
        call init_metadata(t, n, shape, axis)
        allocate(t%data(n))
        t%data = cmplx(0.0_dp, 0.0_dp, dp)
    end function tensor_zeros

    function tensor_ones(shape, axis) result(t)
        integer, intent(in) :: shape(:)
        character(len=*), intent(in), optional :: axis(:)
        type(tensor_t) :: t
        integer :: n

        n = product_or_one(shape)
        call init_metadata(t, n, shape, axis)
        allocate(t%data(n))
        t%data = cmplx(1.0_dp, 0.0_dp, dp)
    end function tensor_ones

    function scalar_tensor_real(x) result(t)
        real(dp), intent(in) :: x
        type(tensor_t) :: t

        allocate(t%shape(0), t%axis(0), t%data(1))
        t%data(1) = cmplx(x, 0.0_dp, dp)
    end function scalar_tensor_real

    function scalar_tensor_complex(x) result(t)
        complex(dp), intent(in) :: x
        type(tensor_t) :: t

        allocate(t%shape(0), t%axis(0), t%data(1))
        t%data(1) = x
    end function scalar_tensor_complex

    pure integer function tensor_rank(self) result(r)
        class(tensor_t), intent(in) :: self

        if (allocated(self%shape)) then
            r = size(self%shape)
        else
            r = 0
        end if
    end function tensor_rank

    pure integer function tensor_nelem(self) result(n)
        class(tensor_t), intent(in) :: self

        if (allocated(self%data)) then
            n = size(self%data)
        else
            n = 0
        end if
    end function tensor_nelem

    pure logical function tensor_valid(self) result(ok)
        class(tensor_t), intent(in) :: self
        integer :: n

        ok = allocated(self%data) .and. allocated(self%shape) .and. allocated(self%axis)
        if (.not. ok) return
        if (size(self%shape) /= size(self%axis)) then
            ok = .false.
            return
        end if
        n = product_or_one(self%shape)
        ok = n == size(self%data) .and. all(self%shape > 0)
    end function tensor_valid

    pure integer function tensor_axis_pos(self, name) result(pos)
        class(tensor_t), intent(in) :: self
        character(len=*), intent(in) :: name
        integer :: k

        pos = 0
        if (.not. allocated(self%axis)) return
        do k = 1, size(self%axis)
            if (trim(self%axis(k)) == trim(name)) then
                if (pos /= 0) then
                    pos = -1
                    return
                end if
                pos = k
            end if
        end do
    end function tensor_axis_pos

    pure logical function tensor_is_real(t, tol) result(ans)
        type(tensor_t), intent(in) :: t
        real(dp), intent(in), optional :: tol
        real(dp) :: eps

        eps = 0.0_dp
        if (present(tol)) eps = tol
        ans = all(abs(aimag(t%data)) <= eps)
    end function tensor_is_real

    function real_data(t, tol) result(x)
        type(tensor_t), intent(in) :: t
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: x(:)
        real(dp) :: eps

        eps = 100.0_dp * epsilon(1.0_dp)
        if (present(tol)) eps = tol
        if (.not. tensor_is_real(t, eps)) error stop "real_data: tensor is complex"
        allocate(x(size(t%data)))
        x = real(t%data, dp)
    end function real_data

    subroutine init_metadata(t, ndata, shape, axis)
        type(tensor_t), intent(out) :: t
        integer, intent(in) :: ndata
        integer, intent(in) :: shape(:)
        character(len=*), intent(in), optional :: axis(:)
        integer :: k, n
        character(len=axis_name_len) :: tmp

        if (any(shape <= 0)) error stop "tensor: all dimensions must be positive"
        n = product_or_one(shape)
        if (n /= ndata) error stop "tensor: data size does not match shape"
        allocate(t%shape(size(shape)), t%axis(size(shape)))
        t%shape = shape
        if (present(axis)) then
            if (size(axis) /= size(shape)) error stop "tensor: axis-name count mismatch"
            do k = 1, size(axis)
                t%axis(k) = trim(axis(k))
            end do
        else
            do k = 1, size(shape)
                write(tmp, '(a,i0)') 'I', k
                t%axis(k) = trim(tmp)
            end do
        end if
    end subroutine init_metadata

    pure integer function product_or_one(x) result(p)
        integer, intent(in) :: x(:)

        if (size(x) == 0) then
            p = 1
        else
            p = product(x)
        end if
    end function product_or_one

end module tensora_types
