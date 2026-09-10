! SPDX-License-Identifier: MIT
! SPDX-FileComment: Public facade for the Fortran readr translation.
module readr
   !! Re-exports the supported eager readr API.
   use readr_parsers, only: guess_parser, parse_character, parse_double, parse_factor, &
      parse_guess, parse_integer, parse_logical, parse_number, parse_vector
   use readr_reader, only: problems, read_csv, read_csv2, read_delim, read_table, read_tsv, &
      spec, spec_csv, spec_csv2, spec_delim, spec_table, spec_tsv, stop_for_problems
   use readr_tokenizer, only: count_fields, read_file, read_lines, tokenize, write_file, write_lines
   use readr_types
   use readr_writer, only: format_csv, format_csv2, format_delim, format_tsv, &
      write_csv, write_csv2, write_delim, write_tsv
   implicit none
   public
end module readr
