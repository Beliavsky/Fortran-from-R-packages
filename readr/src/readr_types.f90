! SPDX-License-Identifier: MIT
! SPDX-FileComment: Public data types for the Fortran readr translation.
module readr_types
   !! Defines collectors, token tables, diagnostics, and parse results.
   use forcats, only: factor_type
   use tibble, only: dp, tibble_column, tibble_type
   implicit none
   private

   integer, parameter, public :: col_guess_code = 0, col_logical_code = 1
   integer, parameter, public :: col_integer_code = 2, col_double_code = 3
   integer, parameter, public :: col_character_code = 4, col_factor_code = 5, col_skip_code = 6, col_number_code = 7
   public :: dp

   type, public :: locale_type
      !! Stores the portable numeric subset of readr locale settings.
      character(len=1) :: decimal_mark = "."
      character(len=1) :: grouping_mark = ","
      logical :: trim_whitespace = .true.
   end type locale_type

   type, public :: collector_type
      !! Specifies how one input column is parsed.
      integer :: type_code = col_guess_code
      character(len=:), allocatable :: levels(:)
      logical :: ordered = .false.
   end type collector_type

   type, public :: col_spec_type
      !! Associates optional column names with collectors and a default collector.
      character(len=:), allocatable :: names(:)
      type(collector_type), allocatable :: collectors(:)
      type(collector_type) :: default
   end type col_spec_type

   type, public :: parse_problem_type
      !! Describes one field that could not be tokenized or converted as requested.
      integer :: row = 0
      integer :: column = 0
      character(len=:), allocatable :: expected
      character(len=:), allocatable :: actual
      character(len=:), allocatable :: message
   end type parse_problem_type

   type, public :: token_table_type
      !! Stores rectangular text fields after quote-aware tokenization.
      character(len=:), allocatable :: field(:,:)
      logical, allocatable :: quoted(:,:)
      integer, allocatable :: fields_per_row(:)
      logical :: valid = .true.
      character(len=:), allocatable :: message
   contains
      procedure :: nrow => token_nrow
      procedure :: ncol => token_ncol
   end type token_table_type

   type, public :: parse_result_type
      !! Holds one parsed tibble column and any conversion problems.
      type(tibble_column) :: values
      type(parse_problem_type), allocatable :: problem(:)
   contains
      procedure :: ok => parse_result_ok
   end type parse_result_type

   type, public :: formatted_lines_type
      !! Owns formatted text records without compiler-dependent deferred-length results.
      character(len=:), allocatable :: line(:)
   end type formatted_lines_type

   type, public :: factor_parse_result_type
      !! Holds one parsed factor and any unknown-level problems.
      type(factor_type) :: values
      type(parse_problem_type), allocatable :: problem(:)
   contains
      procedure :: ok => factor_parse_result_ok
   end type factor_parse_result_type

   type, public :: read_result_type
      !! Holds a parsed rectangular table, its effective specification, and diagnostics.
      type(tibble_type) :: data
      type(col_spec_type) :: specification
      type(parse_problem_type), allocatable :: problem(:)
   contains
      procedure :: ok => read_result_ok
   end type read_result_type

   public :: col_character, col_double, col_factor, col_guess, col_integer
   public :: col_logical, col_number, col_skip, cols, default_locale, locale

contains

   pure elemental function col_guess() result(collector)
      !! Constructs a type-guessing collector.
      type(collector_type) :: collector

      collector%type_code = col_guess_code
   end function col_guess

   pure elemental function col_logical() result(collector)
      !! Constructs a logical collector.
      type(collector_type) :: collector

      collector%type_code = col_logical_code
   end function col_logical

   pure elemental function col_integer() result(collector)
      !! Constructs an integer collector.
      type(collector_type) :: collector

      collector%type_code = col_integer_code
   end function col_integer

   pure elemental function col_double() result(collector)
      !! Constructs a double-precision real collector.
      type(collector_type) :: collector

      collector%type_code = col_double_code
   end function col_double

   pure elemental function col_character() result(collector)
      !! Constructs a character collector.
      type(collector_type) :: collector

      collector%type_code = col_character_code
   end function col_character

   pure elemental function col_number() result(collector)
      !! Constructs a flexible decorated-number collector.
      type(collector_type) :: collector

      collector%type_code = col_number_code
   end function col_number

   pure function col_factor(levels, ordered) result(collector)
      !! Constructs a categorical collector with optional explicit levels.
      character(len=*), intent(in), optional :: levels(:) !! Permitted level labels.
      logical, intent(in), optional :: ordered             !! Whether the factor is ordinal.
      type(collector_type) :: collector

      collector%type_code = col_factor_code
      if (present(levels)) collector%levels = levels
      if (present(ordered)) collector%ordered = ordered
   end function col_factor

   pure elemental function col_skip() result(collector)
      !! Constructs a collector that omits its column.
      type(collector_type) :: collector

      collector%type_code = col_skip_code
   end function col_skip

   pure function cols(names, collectors, default) result(specification)
      !! Constructs a column specification from parallel name and collector arrays.
      character(len=*), intent(in), optional :: names(:) !! Names receiving explicit collectors.
      type(collector_type), intent(in), optional :: collectors(:) !! Explicit collectors.
      type(collector_type), intent(in), optional :: default !! Collector for unspecified columns.
      type(col_spec_type) :: specification

      specification%default = col_guess()
      if (present(default)) specification%default = default
      if (present(names) .neqv. present(collectors)) then
         error stop "cols: names and collectors must be supplied together"
      end if
      if (present(names)) then
         if (size(names) /= size(collectors)) error stop "cols: names and collectors differ in size"
         specification%names = names
         specification%collectors = collectors
      end if
   end function cols

   pure elemental function default_locale() result(settings)
      !! Returns the default decimal, grouping, and whitespace settings.
      type(locale_type) :: settings
   end function default_locale

   pure elemental function locale(decimal_mark, grouping_mark, trim_whitespace) result(settings)
      !! Constructs portable locale settings for numeric parsing.
      character(len=1), intent(in), optional :: decimal_mark  !! Decimal separator.
      character(len=1), intent(in), optional :: grouping_mark !! Thousands separator.
      logical, intent(in), optional :: trim_whitespace        !! Trim unquoted fields.
      type(locale_type) :: settings

      if (present(decimal_mark)) settings%decimal_mark = decimal_mark
      if (present(grouping_mark)) settings%grouping_mark = grouping_mark
      if (present(trim_whitespace)) settings%trim_whitespace = trim_whitespace
   end function locale

   pure elemental integer function token_nrow(self) result(n)
      !! Returns the token row count.
      class(token_table_type), intent(in) :: self !! Token table to inspect.

      n = 0
      if (allocated(self%field)) n = size(self%field, 1)
   end function token_nrow

   pure elemental integer function token_ncol(self) result(n)
      !! Returns the token column count.
      class(token_table_type), intent(in) :: self !! Token table to inspect.

      n = 0
      if (allocated(self%field)) n = size(self%field, 2)
   end function token_ncol

   pure elemental logical function parse_result_ok(self) result(ok)
      !! Reports whether scalar/vector parsing produced no diagnostics.
      class(parse_result_type), intent(in) :: self !! Parse result to inspect.

      ok = .not. allocated(self%problem)
      if (allocated(self%problem)) ok = size(self%problem) == 0
   end function parse_result_ok

   pure elemental logical function factor_parse_result_ok(self) result(ok)
      !! Reports whether factor parsing produced no diagnostics.
      class(factor_parse_result_type), intent(in) :: self !! Parse result to inspect.

      ok = .not. allocated(self%problem)
      if (allocated(self%problem)) ok = size(self%problem) == 0
   end function factor_parse_result_ok

   pure elemental logical function read_result_ok(self) result(ok)
      !! Reports whether rectangular parsing produced no diagnostics.
      class(read_result_type), intent(in) :: self !! Read result to inspect.

      ok = .not. allocated(self%problem)
      if (allocated(self%problem)) ok = size(self%problem) == 0
   end function read_result_ok

end module readr_types
