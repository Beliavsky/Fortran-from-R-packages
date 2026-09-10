program test_construction
   use tibble, only: dp, get_character, get_integer, get_logical, get_missing, get_real, &
      make_column, new_tibble, tibble_character, tibble_column, tibble_integer, &
      tibble_logical, tibble_real, tibble_type, validate_tibble
   implicit none

   type(tibble_column) :: columns(4)
   type(tibble_type) :: table
   character(len=:), allocatable :: labels(:), message
   integer, allocatable :: identifiers(:)
   real(dp), allocatable :: scores(:)
   logical, allocatable :: active(:), missing(:)

   columns(1) = make_column("id", [1, 2, 3])
   columns(2) = make_column("score", [2.5_dp, 4.0_dp, 8.5_dp], [.false., .true., .false.])
   columns(3) = make_column("active", [.true., .false., .true.])
   columns(4) = make_column("label", [character(len=5) :: "alpha", "beta", "gamma"])
   table = new_tibble(columns)

   if (table%nrow() /= 3 .or. table%ncol() /= 4) error stop "wrong table dimensions"
   if (.not. table%has_name("score")) error stop "exact name lookup failed"
   if (table%has_name("sco")) error stop "partial name matching must not occur"
   if (table%columns(1)%type_code /= tibble_integer) error stop "integer type tag failed"
   if (table%columns(2)%type_code /= tibble_real) error stop "real type tag failed"
   if (table%columns(3)%type_code /= tibble_logical) error stop "logical type tag failed"
   if (table%columns(4)%type_code /= tibble_character) error stop "character type tag failed"

   identifiers = get_integer(table, "id")
   scores = get_real(table, "score")
   missing = get_missing(table, "score")
   active = get_logical(table, "active")
   labels = get_character(table, "label")
   if (any(identifiers /= [1, 2, 3])) error stop "integer extraction failed"
   if (any(abs(scores - [2.5_dp, 4.0_dp, 8.5_dp]) > 1.0e-12_dp)) error stop "real extraction failed"
   if (any(missing .neqv. [.false., .true., .false.])) error stop "missing mask failed"
   if (any(active .neqv. [.true., .false., .true.])) error stop "logical extraction failed"
   if (any(labels /= [character(len=5) :: "alpha", "beta", "gamma"])) then
      error stop "character extraction failed"
   end if
   if (.not. validate_tibble(table, message)) error stop "valid table rejected: " // message

   table%columns(4)%missing = [.false., .true.]
   if (validate_tibble(table, message)) error stop "invalid missing mask accepted"
   if (index(message, "missing-value mask") == 0) error stop "unexpected validation message"
end program test_construction
