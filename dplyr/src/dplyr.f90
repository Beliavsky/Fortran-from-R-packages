! SPDX-License-Identifier: MIT
module dplyr
   !! Public umbrella module for the restricted modern-Fortran dplyr API.
   use dplyr_expression, only: evaluate_expression
   use dplyr_groups, only: count_rows, group_by, group_indices, group_keys, grouped_df, &
      group_size, n_groups, summarise, summary, summary_spec, ungroup
   use dplyr_helpers, only: between, coalesce, consecutive_id, first, lag, last, lead, &
      n_distinct, na_if, near, nth
   use dplyr_joins, only: anti_join, cross_join, full_join, inner_join, left_join, &
      right_join, semi_join
   use dplyr_sets, only: intersect, setdiff, setequal, symdiff, union, union_all
   use dplyr_vector_ops
   use dplyr_verbs, only: arrange, bind_cols, bind_rows, distinct, filter, mutate, pull, &
      relocate, rename, select, slice, slice_head, slice_tail
   use tibble
   use vctrs
   implicit none
   public
end module dplyr
