module abind
   use abind_types, only : dp, array_value, int_vector, string_vector, init_array, set_dimnames, &
      set_dimname_names, clear_dimnames, array_rank, array_size, has_dimnames, make_int_vector
   use abind_core, only : abind_arrays, asub_array, afill_array, adrop_array, acorn_array, &
      indices_from_names, indices_from_mask
   implicit none
   private

   public :: dp
   public :: array_value
   public :: int_vector
   public :: string_vector
   public :: init_array
   public :: set_dimnames
   public :: set_dimname_names
   public :: clear_dimnames
   public :: array_rank
   public :: array_size
   public :: has_dimnames
   public :: make_int_vector
   public :: abind_arrays
   public :: asub_array
   public :: afill_array
   public :: adrop_array
   public :: acorn_array
   public :: indices_from_names
   public :: indices_from_mask

end module abind
