module ppso
    use ppso_kinds, only : dp
    use ppso_types
    use ppso_pso, only : optim_pso, init_pso_state, run_pso_state
    use ppso_dds, only : optim_dds, init_dds_state, run_dds_state
    use ppso_checkpoint, only : save_pso_state, load_pso_state, save_dds_state, load_dds_state
    use ppso_benchmarks, only : rastrigin_function, ackley_function, griewank_function, &
        sample_function, sample_function2
    implicit none
    public
end module ppso
