module proxy
    use proxy_kinds, only: dp
    use proxy_ieee, only: proxy_nan, proxy_is_missing
    use proxy_numeric_measures
    use proxy_binary_measures
    use proxy_nominal_measures, only: nominal_similarity
    use proxy_gower
    use proxy_string_measures
    use proxy_utils, only: proxy_convert_default, proxy_convert_one_minus, proxy_convert_none, &
                           proxy_simil_to_dist, proxy_dist_to_simil
    use proxy_registry
    use proxy_engine
    use proxy_dist_utils
    implicit none
    public
end module proxy
