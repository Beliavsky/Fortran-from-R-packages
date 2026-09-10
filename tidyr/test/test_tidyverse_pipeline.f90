! SPDX-License-Identifier: MIT
! SPDX-FileComment: Cross-package integration test for translated tidyverse components.
program test_tidyverse_pipeline
   !! Exercises readr, tibble, tidyselect, dplyr, tidyr, and stringr together.
   use dplyr, only: select
   use readr, only: read_csv, read_result_type
   use stringr, only: str_detect, str_to_upper
   use tibble, only: tibble_type
   use tidyselect, only: all_of
   use tidyr, only: drop_na, pivot_longer
   implicit none
   character(len=16) :: names(4)
   integer, allocatable :: selected(:)
   type(read_result_type) :: parsed
   type(tibble_type) :: data, long

   parsed = read_csv('../readr/data/sample.csv')
   if (.not. parsed%ok()) error stop 'pipeline: readr reported a problem'
   data = parsed%data
   names = [character(len=16) :: 'id', 'active', 'score', 'name']
   selected = all_of(names, [character(len=5) :: 'id', 'score'])
   data = select(data, selected)
   data = drop_na(data, [character(len=5) :: 'score'])
   long = pivot_longer(data, all_of([character(len=5) :: 'id', 'score'], ['score']), &
                       names_to='measure', values_to='value')
   if (long%nrow() /= 2) error stop 'pipeline: unexpected row count'
   if (.not. all(str_detect(str_to_upper(long%columns(2)%character_values), 'SCORE'))) then
      error stop 'pipeline: string transformation failed'
   end if
   print '(a)', 'Translated tidyverse integration test passed.'
end program test_tidyverse_pipeline
