module varmapack
   use varmapack_kinds, only : dp
   use varmapack_analysis, only : varmapack_autocov, varmapack_cov2corr
   use varmapack_model_mod, only : make_varmapack_model, varmapack_model_type
   use varmapack_testcases_mod, only : varmapack_testcase, varmapack_testcases
   implicit none
   private
   public :: dp
   public :: varmapack_autocov
   public :: varmapack_cov2corr
   public :: varmapack_model
   public :: varmapack_model_type
   public :: varmapack_testcase
   public :: varmapack_testcases

   interface varmapack_model
      procedure :: make_varmapack_model
   end interface varmapack_model
end module varmapack
