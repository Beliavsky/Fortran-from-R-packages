module nspmix_types
   use nspmix_kinds, only : dp
   implicit none
   private
   integer, parameter, public :: NSP_NORMAL=1, NSP_POISSON=2, NSP_GEOM=3
   integer, parameter, public :: NSP_NBINOM=4, NSP_CVPS=5, NSP_MLOGIT=6

   type, public :: disc_dist
      real(dp), allocatable :: pt(:)
      real(dp), allocatable :: pr(:)
   end type disc_dist

   type, public :: nsp_data
      integer :: family = 0
      real(dp), allocatable :: v(:), w(:)
      real(dp) :: size = 1.0_dp
      real(dp), allocatable :: ni(:), mi(:), ri(:)
      integer, allocatable :: group(:)
      real(dp), allocatable :: y(:), trials(:), xmat(:,:)
   end type nsp_data

   type, public :: nspmix_result
      type(disc_dist) :: mix
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: gradient(:)
      real(dp) :: ll = -huge(1.0_dp)
      real(dp) :: max_gradient = huge(1.0_dp)
      integer :: iterations = 0
      integer :: convergence = 1
   end type nspmix_result

   type, public :: hcnm_result
      real(dp), allocatable :: p(:)
      real(dp) :: ll = -huge(1.0_dp)
      real(dp) :: maxgrad = huge(1.0_dp)
      integer :: iterations = 0
      integer :: convergence = 1
   end type hcnm_result
end module nspmix_types
