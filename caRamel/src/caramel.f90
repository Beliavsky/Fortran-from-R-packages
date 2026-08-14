module caramel
    use caramel_kinds, only: dp
    use caramel_pareto, only: pareto, pareto_2d, pareto_3d, dominate, dominated
    use caramel_utils, only: val2rank, boxes, rselect, matvcov, vol_splx, dimprove, downsize, decrease_pop
    use caramel_generation, only: index_block, cinterp, cextrap, crecombination, cusecovar, new_xval
    use caramel_optimizer, only: caramel_options, caramel_result, objective_function, caramel_optimize
    use caramel_random, only: seed_random
    implicit none
    public
end module caramel
