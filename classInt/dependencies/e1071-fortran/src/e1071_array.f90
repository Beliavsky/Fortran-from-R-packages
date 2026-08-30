module e1071_array
    implicit none
    private

    public :: array_linear_index

contains

    pure function array_linear_index(dimensions, subscripts) result(index_value)
        integer, intent(in) :: dimensions(:) !! Positive extent of each array dimension in Fortran/R column-major order.
        integer, intent(in) :: subscripts(:) !! One-based subscript in every dimension; length must match `dimensions`.
        integer :: index_value
        integer :: stride
        integer :: k

        if (size(subscripts) /= size(dimensions)) error stop "array_linear_index: wrong number of subscripts"
        if (any(dimensions < 1)) error stop "array_linear_index: dimensions must be positive"
        if (any(subscripts < 1) .or. any(subscripts > dimensions)) error stop "array_linear_index: subscript out of range"
        index_value = 1
        stride = 1
        do k = 1, size(dimensions)
            index_value = index_value + (subscripts(k) - 1) * stride
            stride = stride * dimensions(k)
        end do
    end function array_linear_index

end module e1071_array
