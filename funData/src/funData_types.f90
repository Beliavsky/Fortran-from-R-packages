module funData_types
   use funData_kinds, only : dp
   implicit none
   private

   type, public :: real_vector
      real(dp), allocatable :: values(:)
   end type real_vector

   type, public :: int_vector
      integer, allocatable :: values(:)
   end type int_vector

   type, public :: fun_data
      type(real_vector), allocatable :: argvals(:)
      integer, allocatable :: dims(:)
      real(dp), allocatable :: x(:)
   end type fun_data

   type, public :: irreg_curve
      real(dp), allocatable :: argvals(:)
      real(dp), allocatable :: x(:)
   end type irreg_curve

   type, public :: irreg_fun_data
      type(irreg_curve), allocatable :: curves(:)
   end type irreg_fun_data

   type, public :: multi_fun_data
      type(fun_data), allocatable :: components(:)
   end type multi_fun_data


   type, public :: basis_spec
      type(real_vector), allocatable :: argvals(:)
      integer, allocatable :: m(:)
      character(len=16), allocatable :: basis_type(:)
      type(int_vector), allocatable :: ignore_deg(:)
   end type basis_spec

   type, public :: sim_fun_result
      type(fun_data) :: sim_data
      type(fun_data) :: true_funs
      real(dp), allocatable :: true_vals(:)
   end type sim_fun_result

   type, public :: sim_multi_result
      type(multi_fun_data) :: sim_data
      type(multi_fun_data) :: true_funs
      real(dp), allocatable :: true_vals(:)
   end type sim_multi_result

end module funData_types
