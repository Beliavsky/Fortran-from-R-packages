program test_compare_match
   use vctrs, only: new_vctr, vctr_type, vec_compare, vec_duplicate_any, &
      vec_duplicate_detect, vec_equal, vec_group_id, vec_identify_runs, &
      vec_if_else, vec_in, vec_match
   implicit none

   type(vctr_type) :: x, y, result
   integer, allocatable :: groups(:), runs(:)
   logical, allocatable :: duplicates(:)

   x = new_vctr("x", [1, 2, 0], [.false., .false., .true.])
   y = new_vctr("y", [1, 3, 0], [.false., .false., .true.])
   result = vec_equal(x, y)
   if (any(result%logical_values .neqv. [.true., .false., .false.])) then
      error stop "elementwise equality values failed"
   end if
   if (any(result%missing .neqv. [.false., .false., .true.])) then
      error stop "elementwise equality missingness failed"
   end if
   result = vec_equal(x, y, na_equal=.true.)
   if (any(result%logical_values .neqv. [.true., .false., .true.])) then
      error stop "missing-equal comparison failed"
   end if

   result = vec_compare(new_vctr("", [1, 4]), new_vctr("", [2, 4]))
   if (any(result%integer_values /= [-1, 0])) error stop "three-way comparison failed"

   result = vec_match(new_vctr("", [3, 1, 9]), new_vctr("", [1, 2, 3]))
   if (any(result%integer_values /= [3, 1, 0])) error stop "match locations failed"
   if (any(result%missing .neqv. [.false., .false., .true.])) then
      error stop "unmatched missingness failed"
   end if
   result = vec_in(new_vctr("", [3, 9]), new_vctr("", [1, 2, 3]))
   if (any(result%logical_values .neqv. [.true., .false.])) error stop "membership failed"

   x = new_vctr("", [2, 1, 2, 3, 1])
   duplicates = vec_duplicate_detect(x)
   if (any(duplicates .neqv. [.true., .true., .true., .false., .true.])) then
      error stop "duplicate detection failed"
   end if
   if (.not. vec_duplicate_any(x)) error stop "duplicate-any failed"
   groups = vec_group_id(x)
   if (any(groups /= [1, 2, 1, 3, 2])) error stop "group identifiers failed"
   runs = vec_identify_runs(new_vctr("", [1, 1, 2, 2, 1]))
   if (any(runs /= [1, 1, 2, 2, 3])) error stop "run identifiers failed"

   result = vec_if_else(new_vctr("", [.true., .false., .true.]), &
      new_vctr("", 10, 1), new_vctr("", [1, 2, 3]))
   if (any(result%integer_values /= [10, 2, 10])) error stop "type-stable if-else failed"
end program test_compare_match
