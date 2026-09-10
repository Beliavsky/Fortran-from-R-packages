program test_slice_combine
   use vctrs, only: dp, new_vctr, vctr_type, vctrs_real, vec_c, vec_copy_element, &
      vec_init, vec_insert, vec_order, vec_rep, vec_rep_each, vec_slice, vec_sort, vec_unique, &
      vec_unique_count, vec_unique_loc
   implicit none

   type(vctr_type) :: character_vector, sliced, combined, sorted, unique
   type(vctr_type) :: numeric_vectors(2), character_vectors(2)
   integer, allocatable :: locations(:), order(:)

   character_vector = new_vctr("word", [character(len=5) :: "beta", "alpha", "beta"])
   sliced = vec_slice(character_vector, [3, 2, 2])
   if (any(sliced%character_values /= [character(len=5) :: "beta", "alpha", "alpha"])) then
      error stop "character slicing failed"
   end if

   numeric_vectors(1) = new_vctr("x", [1, 2])
   numeric_vectors(2) = new_vctr("x", [3.5_dp])
   combined = vec_c(numeric_vectors)
   if (combined%type_code /= vctrs_real) error stop "concatenation promotion failed"
   if (any(abs(combined%real_values - [1.0_dp, 2.0_dp, 3.5_dp]) > 0.0_dp)) then
      error stop "numeric concatenation failed"
   end if

   character_vectors(1) = new_vctr("x", [character(len=2) :: "a", "bb"])
   character_vectors(2) = new_vctr("x", [character(len=8) :: "longword"])
   combined = vec_c(character_vectors)
   if (combined%character_values(3) /= "longword") error stop "character width lost"
   character_vectors(1) = new_vctr("x", [character(len=1) ::])
   character_vectors(2) = new_vctr("x", [character(len=3) :: "abc"])
   combined = vec_c(character_vectors)
   if (combined%size() /= 1 .or. combined%character_values(1) /= "abc") then
      error stop "empty-vector concatenation failed"
   end if

   character_vectors(1) = new_vctr("x", [character(len=2) :: "a", "bb"])
   character_vectors(2) = new_vctr("x", [character(len=8) :: "longword"])

   unique = vec_unique(character_vector)
   locations = vec_unique_loc(character_vector)
   if (vec_unique_count(character_vector) /= 2) error stop "unique count failed"
   if (any(locations /= [1, 2])) error stop "unique locations failed"
   if (any(unique%character_values /= [character(len=5) :: "beta", "alpha"])) then
      error stop "unique values failed"
   end if

   character_vector%missing = [.false., .true., .false.]
   order = vec_order(character_vector)
   if (any(order /= [1, 3, 2])) error stop "stable missing-last order failed"
   sorted = vec_sort(character_vector, descending=.true.)
   if (sorted%missing(3) .neqv. .true.) error stop "sorted missing position failed"

   combined = vec_rep(new_vctr("x", [1, 2]), 2)
   if (any(combined%integer_values /= [1, 2, 1, 2])) error stop "whole-vector repeat failed"
   combined = vec_rep_each(new_vctr("x", [1, 2]), 2)
   if (any(combined%integer_values /= [1, 1, 2, 2])) error stop "element repeat failed"

   combined = vec_insert(new_vctr("x", [1, 2]), new_vctr("x", [8.5_dp]), 2)
   if (combined%type_code /= vctrs_real) error stop "insertion promotion failed"
   if (any(abs(combined%real_values - [1.0_dp, 8.5_dp, 2.0_dp]) > 0.0_dp)) then
      error stop "middle insertion failed"
   end if
   combined = vec_insert(character_vectors(1), character_vectors(2), 1)
   if (combined%character_values(1) /= "longword") error stop "front insertion failed"
   combined = vec_insert(character_vectors(1), character_vectors(2), 3)
   if (combined%character_values(3) /= "longword") error stop "append insertion failed"

   sliced = vec_init(character_vector, 2, missing=.true.)
   call vec_copy_element(character_vector, 1, sliced, 2)
   if (.not. sliced%missing(1) .or. sliced%missing(2)) error stop "initialized mask copy failed"
   if (sliced%character_values(2) /= character_vector%character_values(1)) then
      error stop "single-element copy failed"
   end if
end program test_slice_combine
