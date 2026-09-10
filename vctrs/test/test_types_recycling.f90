program test_types_recycling
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   use vctrs, only: dp, new_vctr, vctr_type, vctrs_integer, vctrs_real, &
      vec_any_missing, vec_cast, vec_detect_missing, vec_is, vec_ptype2, &
      vec_ptype_common, vec_recycle, vec_recycle_common, vec_size_common
   implicit none

   type(vctr_type) :: logical_vector, integer_vector, real_vector, recycled
   type(vctr_type) :: inputs(3)
   type(vctr_type), allocatable :: common(:)

   logical_vector = new_vctr("flag", [.true., .false.], [.false., .true.])
   integer_vector = vec_cast(logical_vector, vctrs_integer)
   real_vector = vec_cast(integer_vector, vctrs_real)
   if (.not. vec_is(logical_vector)) error stop "valid logical vector rejected"
   if (.not. vec_any_missing(logical_vector)) error stop "missing value not detected"
   if (any(vec_detect_missing(logical_vector) .neqv. [.false., .true.])) then
      error stop "missing mask changed"
   end if
   if (any(integer_vector%integer_values /= [1, 0])) error stop "logical-to-integer cast failed"
   if (any(abs(real_vector%real_values - [1.0_dp, 0.0_dp]) > 0.0_dp)) then
      error stop "integer-to-real cast failed"
   end if
   if (vec_ptype2(logical_vector, integer_vector) /= vctrs_integer) then
      error stop "logical/integer common type failed"
   end if
   if (vec_ptype2(integer_vector, real_vector) /= vctrs_real) then
      error stop "integer/real common type failed"
   end if

   inputs(1) = new_vctr("a", 2, 1)
   inputs(2) = new_vctr("b", [1, 2, 3])
   inputs(3) = new_vctr("c", 4.0_dp, 1)
   if (vec_size_common(inputs) /= 3) error stop "common size failed"
   if (vec_ptype_common(inputs) /= vctrs_real) error stop "common type failed"
   common = vec_recycle_common(inputs)
   if (any(common(1)%integer_values /= 2)) error stop "common recycling failed"
   if (common(3)%size() /= 3) error stop "real scalar recycling size failed"
   recycled = vec_recycle(new_vctr("empty", 7, 1), 0)
   if (recycled%size() /= 0) error stop "size-one to empty recycling failed"

   real_vector = new_vctr("masked", [ieee_value(0.0_dp, ieee_quiet_nan), 2.0_dp], &
      [.true., .false.])
   integer_vector = vec_cast(real_vector, vctrs_integer)
   if (any(integer_vector%integer_values /= [0, 2])) then
      error stop "masked real-to-integer cast failed"
   end if
end program test_types_recycling
