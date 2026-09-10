! SPDX-License-Identifier: MIT
! SPDX-FileComment: Pivot example for the Fortran tidyr translation.
program reshape
   !! Converts a small wide tibble to tidy long form.
   use tidyr
   implicit none
   type(tibble_type) :: data, tidy

   data = new_tibble([make_column('station', [1, 2]), &
                      make_column('morning', [12.0d0, 14.0d0]), &
                      make_column('evening', [17.0d0, 18.0d0])])
   tidy = pivot_longer(data, ['morning', 'evening'], 'time', 'temperature')
   call print_tibble(tidy)
end program reshape
