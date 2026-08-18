! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

module partitions
    use partitions_kinds, only : i8, dp
    use partitions_counts
    use partitions_enumeration, only : triangular_index, first_part, next_part, is_last_part, &
        distinct_parts, first_distinct_part, next_distinct_part, is_last_distinct_part, &
        restricted_parts, first_restricted_part, next_restricted_part, is_last_restricted_part, &
        block_parts, first_block_part, next_block_part, is_last_block_part, first_composition, &
        next_composition, is_last_composition, to_binary, to_decimal, composition_to_binary, &
        binary_to_composition, conjugate, durfee
    use partitions_permutations, only : next_permutation, permutations, plain_permutations, &
        multiset_permutations
    use partitions_sets, only : set_partitions, restricted_set_partitions, multiset_sequences, &
        multinomial_permutations, all_binomial, generalized_riffles, riffles
    use partitions_compat
    implicit none
    public
end module partitions
