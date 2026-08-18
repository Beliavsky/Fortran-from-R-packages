! SPDX-License-Identifier: GPL-2.0-or-later
module tensora
    use tensora_kinds, only : dp, i8, axis_name_len
    use tensora_types
    use tensora_index, only : contraname, is_covariate_name, is_contravariate_name, &
        as_covariate_name, as_contravariate_name, positions_by_name, product_shape
    use tensora_core
    use tensora_linalg
    use tensora_stats
    implicit none
    public

    interface operator(+)
        module procedure tensor_plus
    end interface
    interface operator(-)
        module procedure tensor_minus
        module procedure tensor_negate
    end interface
    interface operator(*)
        module procedure tensor_times
        module procedure tensor_times_scalar
        module procedure scalar_times_tensor
        module procedure tensor_times_real
        module procedure real_times_tensor
    end interface
    interface operator(/)
        module procedure tensor_divide
        module procedure tensor_divide_scalar
        module procedure tensor_divide_real
    end interface

contains

    function tensor_plus(a, b) result(c)
        type(tensor_t), intent(in) :: a, b
        type(tensor_t) :: c
        c = add_tensor(a,b)
    end function tensor_plus

    function tensor_minus(a, b) result(c)
        type(tensor_t), intent(in) :: a, b
        type(tensor_t) :: c
        c = sub_tensor(a,b)
    end function tensor_minus

    function tensor_negate(a) result(c)
        type(tensor_t), intent(in) :: a
        type(tensor_t) :: c
        c = scale_tensor(a, cmplx(-1.0_dp,0.0_dp,dp))
    end function tensor_negate

    function tensor_times(a, b) result(c)
        type(tensor_t), intent(in) :: a, b
        type(tensor_t) :: c
        c = elem_mul_tensor(a,b)
    end function tensor_times

    function tensor_times_scalar(a, b) result(c)
        type(tensor_t), intent(in) :: a
        complex(dp), intent(in) :: b
        type(tensor_t) :: c
        c = scale_tensor(a,b)
    end function tensor_times_scalar

    function scalar_times_tensor(a, b) result(c)
        complex(dp), intent(in) :: a
        type(tensor_t), intent(in) :: b
        type(tensor_t) :: c
        c = scale_tensor(b,a)
    end function scalar_times_tensor

    function tensor_times_real(a, b) result(c)
        type(tensor_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(tensor_t) :: c
        c = scale_tensor(a, cmplx(b,0.0_dp,dp))
    end function tensor_times_real

    function real_times_tensor(a, b) result(c)
        real(dp), intent(in) :: a
        type(tensor_t), intent(in) :: b
        type(tensor_t) :: c
        c = scale_tensor(b, cmplx(a,0.0_dp,dp))
    end function real_times_tensor

    function tensor_divide(a, b) result(c)
        type(tensor_t), intent(in) :: a, b
        type(tensor_t) :: c
        c = elem_div_tensor(a,b)
    end function tensor_divide

    function tensor_divide_real(a, b) result(c)
        type(tensor_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(tensor_t) :: c
        c = scale_tensor(a, cmplx(1.0_dp/b,0.0_dp,dp))
    end function tensor_divide_real

    function tensor_divide_scalar(a, b) result(c)
        type(tensor_t), intent(in) :: a
        complex(dp), intent(in) :: b
        type(tensor_t) :: c
        c = scale_tensor(a,1.0_dp/b)
    end function tensor_divide_scalar

end module tensora
