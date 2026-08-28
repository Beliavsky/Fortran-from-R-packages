module actuar
    use actuar_kinds
    use actuar_special
    use actuar_special_v02
    use actuar_continuous
    use actuar_continuous_v02
    use actuar_supplement
    use actuar_supplement_v02
    use actuar_discrete
    use actuar_phase_type
    use actuar_risk
    use actuar_aggregate_v02
    use actuar_ruin_v02
    use actuar_credibility
    use actuar_credibility_v02
    use actuar_hachemeister_v02
    use actuar_grouped_v02
    use actuar_mde_v03
    use actuar_coverage_v03
    use actuar_hache_bary_v03
    use actuar_hierarc_exact_v03
    use actuar_hier_sim_v03
    implicit none
    public
contains
    pure real(dp) function mbeta(order,shape1,shape2) result(moment)
        real(dp), intent(in) :: order,shape1,shape2
        moment=mbeta_act(order,shape1,shape2)
    end function mbeta

    pure real(dp) function levbeta(limit,shape1,shape2,order) result(moment)
        real(dp), intent(in) :: limit,shape1,shape2,order
        moment=levbeta_act(limit,shape1,shape2,order)
    end function levbeta

    pure real(dp) function munif(order,xmin,xmax) result(moment)
        real(dp), intent(in) :: order,xmin,xmax
        moment=munif_act(order,xmin,xmax)
    end function munif

    pure real(dp) function levunif(limit,xmin,xmax,order) result(moment)
        real(dp), intent(in) :: limit,xmin,xmax,order
        moment=levunif_act(limit,xmin,xmax,order)
    end function levunif
end module actuar
