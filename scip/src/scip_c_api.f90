! Apache-2.0
module scip_c_api
    use, intrinsic :: iso_c_binding
    implicit none
    private

    public :: c_scipf_model_create, c_scipf_model_free, c_scipf_add_var
    public :: c_scipf_add_linear_cons, c_scipf_add_quadratic_cons
    public :: c_scipf_add_sos1_cons, c_scipf_add_sos2_cons
    public :: c_scipf_add_indicator_cons
    public :: c_scipf_set_param_bool, c_scipf_set_param_int
    public :: c_scipf_set_param_long, c_scipf_set_param_real
    public :: c_scipf_set_param_char, c_scipf_set_param_string
    public :: c_scipf_set_heuristics, c_scipf_set_objective_sense
    public :: c_scipf_set_param_integer_auto, c_scipf_set_param_real_auto
    public :: c_scipf_set_param_text_auto
    public :: c_scipf_optimize, c_scipf_get_status, c_scipf_get_nsols
    public :: c_scipf_get_nvars, c_scipf_get_nconss
    public :: c_scipf_get_best_solution, c_scipf_get_solution_k, c_scipf_get_info

    interface
        integer(c_int) function c_scipf_model_create(name, handle) &
            bind(C, name="scipf_model_create")
            import :: c_int, c_char, c_ptr
            character(c_char), intent(in) :: name(*)
            type(c_ptr), intent(out) :: handle
        end function

        integer(c_int) function c_scipf_model_free(handle) &
            bind(C, name="scipf_model_free")
            import :: c_int, c_ptr
            type(c_ptr), value :: handle
        end function

        integer(c_int) function c_scipf_add_var(handle, obj, lb, ub, vtype, name, index1) &
            bind(C, name="scipf_add_var")
            import :: c_int, c_ptr, c_double, c_char
            type(c_ptr), value :: handle
            real(c_double), value :: obj, lb, ub
            integer(c_int), value :: vtype
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: index1
        end function

        integer(c_int) function c_scipf_add_linear_cons(handle, nv, vars1, coefs, &
            lhs, rhs, name, index1) bind(C, name="scipf_add_linear_cons")
            import :: c_int, c_ptr, c_double, c_char
            type(c_ptr), value :: handle
            integer(c_int), value :: nv
            integer(c_int), intent(in) :: vars1(*)
            real(c_double), intent(in) :: coefs(*)
            real(c_double), value :: lhs, rhs
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: index1
        end function

        integer(c_int) function c_scipf_add_quadratic_cons(handle, nlin, linvars1, &
            lincoefs, nquad, qvars1, qvars2, qcoefs, lhs, rhs, name, index1) &
            bind(C, name="scipf_add_quadratic_cons")
            import :: c_int, c_ptr, c_double, c_char
            type(c_ptr), value :: handle
            integer(c_int), value :: nlin, nquad
            integer(c_int), intent(in) :: linvars1(*), qvars1(*), qvars2(*)
            real(c_double), intent(in) :: lincoefs(*), qcoefs(*)
            real(c_double), value :: lhs, rhs
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: index1
        end function

        integer(c_int) function c_scipf_add_sos1_cons(handle, nv, vars1, weights, &
            name, index1) bind(C, name="scipf_add_sos1_cons")
            import :: c_int, c_ptr, c_double, c_char
            type(c_ptr), value :: handle
            integer(c_int), value :: nv
            integer(c_int), intent(in) :: vars1(*)
            real(c_double), intent(in) :: weights(*)
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: index1
        end function

        integer(c_int) function c_scipf_add_sos2_cons(handle, nv, vars1, weights, &
            name, index1) bind(C, name="scipf_add_sos2_cons")
            import :: c_int, c_ptr, c_double, c_char
            type(c_ptr), value :: handle
            integer(c_int), value :: nv
            integer(c_int), intent(in) :: vars1(*)
            real(c_double), intent(in) :: weights(*)
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: index1
        end function

        integer(c_int) function c_scipf_add_indicator_cons(handle, binvar1, nv, &
            vars1, coefs, rhs, name, index1) bind(C, name="scipf_add_indicator_cons")
            import :: c_int, c_ptr, c_double, c_char
            type(c_ptr), value :: handle
            integer(c_int), value :: binvar1, nv
            integer(c_int), intent(in) :: vars1(*)
            real(c_double), intent(in) :: coefs(*)
            real(c_double), value :: rhs
            character(c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: index1
        end function

        integer(c_int) function c_scipf_set_param_bool(handle, name, value) &
            bind(C, name="scipf_set_param_bool")
            import :: c_int, c_ptr, c_char
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
            integer(c_int), value :: value
        end function

        integer(c_int) function c_scipf_set_param_int(handle, name, value) &
            bind(C, name="scipf_set_param_int")
            import :: c_int, c_ptr, c_char
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
            integer(c_int), value :: value
        end function

        integer(c_int) function c_scipf_set_param_long(handle, name, value) &
            bind(C, name="scipf_set_param_long")
            import :: c_int, c_ptr, c_char, c_long_long
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
            integer(c_long_long), value :: value
        end function

        integer(c_int) function c_scipf_set_param_real(handle, name, value) &
            bind(C, name="scipf_set_param_real")
            import :: c_int, c_ptr, c_char, c_double
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
            real(c_double), value :: value
        end function

        integer(c_int) function c_scipf_set_param_char(handle, name, value) &
            bind(C, name="scipf_set_param_char")
            import :: c_int, c_ptr, c_char
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
            character(c_char), value :: value
        end function

        integer(c_int) function c_scipf_set_param_string(handle, name, value) &
            bind(C, name="scipf_set_param_string")
            import :: c_int, c_ptr, c_char
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*), value(*)
        end function

        integer(c_int) function c_scipf_set_param_integer_auto(handle, name, value) &
            bind(C, name="scipf_set_param_integer_auto")
            import :: c_int, c_ptr, c_char, c_long_long
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
            integer(c_long_long), value :: value
        end function

        integer(c_int) function c_scipf_set_param_real_auto(handle, name, value) &
            bind(C, name="scipf_set_param_real_auto")
            import :: c_int, c_ptr, c_char, c_double
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*)
            real(c_double), value :: value
        end function

        integer(c_int) function c_scipf_set_param_text_auto(handle, name, value) &
            bind(C, name="scipf_set_param_text_auto")
            import :: c_int, c_ptr, c_char
            type(c_ptr), value :: handle
            character(c_char), intent(in) :: name(*), value(*)
        end function

        integer(c_int) function c_scipf_set_heuristics(handle, setting) &
            bind(C, name="scipf_set_heuristics")
            import :: c_int, c_ptr
            type(c_ptr), value :: handle
            integer(c_int), value :: setting
        end function

        integer(c_int) function c_scipf_set_objective_sense(handle, maximize) &
            bind(C, name="scipf_set_objective_sense")
            import :: c_int, c_ptr
            type(c_ptr), value :: handle
            integer(c_int), value :: maximize
        end function

        integer(c_int) function c_scipf_optimize(handle) bind(C, name="scipf_optimize")
            import :: c_int, c_ptr
            type(c_ptr), value :: handle
        end function

        integer(c_int) function c_scipf_get_status(handle, buf, nbuf) &
            bind(C, name="scipf_get_status")
            import :: c_int, c_ptr, c_char
            type(c_ptr), value :: handle
            character(c_char), intent(out) :: buf(*)
            integer(c_int), value :: nbuf
        end function

        integer(c_int) function c_scipf_get_nsols(handle) bind(C, name="scipf_get_nsols")
            import :: c_int, c_ptr
            type(c_ptr), value :: handle
        end function

        integer(c_int) function c_scipf_get_nvars(handle) bind(C, name="scipf_get_nvars")
            import :: c_int, c_ptr
            type(c_ptr), value :: handle
        end function

        integer(c_int) function c_scipf_get_nconss(handle) bind(C, name="scipf_get_nconss")
            import :: c_int, c_ptr
            type(c_ptr), value :: handle
        end function

        integer(c_int) function c_scipf_get_best_solution(handle, obj, x, n) &
            bind(C, name="scipf_get_best_solution")
            import :: c_int, c_ptr, c_double
            type(c_ptr), value :: handle
            real(c_double), intent(out) :: obj
            real(c_double), intent(out) :: x(*)
            integer(c_int), value :: n
        end function

        integer(c_int) function c_scipf_get_solution_k(handle, k1, obj, x, n) &
            bind(C, name="scipf_get_solution_k")
            import :: c_int, c_ptr, c_double
            type(c_ptr), value :: handle
            integer(c_int), value :: k1
            real(c_double), intent(out) :: obj
            real(c_double), intent(out) :: x(*)
            integer(c_int), value :: n
        end function

        integer(c_int) function c_scipf_get_info(handle, solve_time, nodes, &
            iterations, gap, sol_count) bind(C, name="scipf_get_info")
            import :: c_int, c_ptr, c_double, c_long_long
            type(c_ptr), value :: handle
            real(c_double), intent(out) :: solve_time, gap
            integer(c_long_long), intent(out) :: nodes, iterations
            integer(c_int), intent(out) :: sol_count
        end function
    end interface
end module scip_c_api
