module deoptim
    use deoptim_kinds, only : dp, i8
    use deoptim_types, only : de_control, de_result, de_objective, de_map, &
        de_success, de_invalid_input, de_unsupported, de_objective_nan, de_map_error
    use deoptim_core, only : deoptim_solve
    implicit none
    public
end module deoptim
