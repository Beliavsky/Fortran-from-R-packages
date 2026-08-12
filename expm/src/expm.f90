! Public convenience module for expm-fortran.
module expm_module
    use expm_kinds, only : dp
    use expm_linalg, only : balance_real_result, balance_complex_result, balance_real, balance_complex
    use expm_matrix_functions, only : expm, expm_higham08, expm_pade, expm_taylor, &
        expm_almohy09, expm_rbs, expm_ward77, matrix_power
    use expm_frechet, only : expm_frechet_sps, expm_frechet_block
    use expm_condition, only : expm_cond_exact, expm_cond_1_est, expm_cond_f_est
    use expm_action, only : exp_action_result, exp_at_v
    use expm_logsqrt, only : sqrtm, logm
    use expm_eigen_methods, only : expm_eigen, logm_eigen, expm_hybrid_eigen_ward
    implicit none
    public
end module expm_module
