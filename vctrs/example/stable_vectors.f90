program stable_vectors
   use vctrs, only: dp, new_vctr, vctr_type, vctrs_real, vec_c, vec_cast, &
      vec_group_id, vec_sort
   implicit none

   type(vctr_type) :: inputs(2), combined, sorted
   integer, allocatable :: groups(:)

   inputs(1) = new_vctr("value", [1, 2, 1])
   inputs(2) = new_vctr("value", [3.5_dp, 2.0_dp])
   combined = vec_c(inputs)
   sorted = vec_sort(combined)
   groups = vec_group_id(combined)

   if (combined%type_code /= vctrs_real) error stop "unexpected common type"
   write (*, '(a,*(g0,1x))') "combined: ", combined%real_values
   write (*, '(a,*(g0,1x))') "sorted:   ", sorted%real_values
   write (*, '(a,*(i0,1x))') "groups:   ", groups

   combined = vec_cast(inputs(1), vctrs_real)
   write (*, '(a,*(g0,1x))') "cast:     ", combined%real_values
end program stable_vectors
