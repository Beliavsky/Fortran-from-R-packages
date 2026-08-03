! SPDX-License-Identifier: GPL-3.0-or-later
module pinstimation_types
   use pinstimation_kinds, only : dp, i8
   implicit none
   private

   type, public :: trade_counts
      integer(i8), allocatable :: buys(:)
      integer(i8), allocatable :: sells(:)
   contains
      procedure :: size => trade_count_size
      procedure :: valid => trade_count_valid
   end type trade_counts

   type, public :: pin_parameters
      real(dp) :: alpha = 0.5_dp
      real(dp) :: delta = 0.5_dp
      real(dp) :: mu = 1.0_dp
      real(dp) :: eps_b = 1.0_dp
      real(dp) :: eps_s = 1.0_dp
   end type pin_parameters

   type, public :: pin_result
      type(pin_parameters) :: parameters
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: pin = 0.0_dp
      real(dp) :: pin_good = 0.0_dp
      real(dp) :: pin_bad = 0.0_dp
      real(dp), allocatable :: posteriors(:,:)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 1
      logical :: converged = .false.
   end type pin_result

   type, public :: mpin_parameters
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: mu(:)
      real(dp) :: eps_b = 1.0_dp
      real(dp) :: eps_s = 1.0_dp
   contains
      procedure :: layers => mpin_layer_count
   end type mpin_parameters

   type, public :: mpin_result
      type(mpin_parameters) :: parameters
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: mpin = 0.0_dp
      real(dp), allocatable :: mpin_layer(:)
      real(dp), allocatable :: good_layer(:)
      real(dp), allocatable :: bad_layer(:)
      real(dp), allocatable :: posteriors(:,:)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 1
      logical :: converged = .false.
   end type mpin_result

   type, public :: adjpin_parameters
      real(dp) :: alpha = 0.5_dp
      real(dp) :: delta = 0.5_dp
      real(dp) :: theta = 0.2_dp
      real(dp) :: theta_p = 0.2_dp
      real(dp) :: eps_b = 1.0_dp
      real(dp) :: eps_s = 1.0_dp
      real(dp) :: mu_b = 1.0_dp
      real(dp) :: mu_s = 1.0_dp
      real(dp) :: d_b = 1.0_dp
      real(dp) :: d_s = 1.0_dp
   end type adjpin_parameters

   type, public :: adjpin_restrictions
      logical :: equal_theta = .false.
      logical :: equal_eps = .false.
      logical :: equal_mu = .false.
      logical :: equal_d = .false.
   end type adjpin_restrictions

   type, public :: adjpin_result
      type(adjpin_parameters) :: parameters
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: adjpin = 0.0_dp
      real(dp) :: psos = 0.0_dp
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 1
      logical :: converged = .false.
   end type adjpin_result

   type, public :: vpin_result
      real(dp) :: volume_bucket_size = 0.0_dp
      real(dp), allocatable :: buy_volume(:)
      real(dp), allocatable :: sell_volume(:)
      real(dp), allocatable :: duration(:)
      real(dp), allocatable :: vpin(:)
      real(dp), allocatable :: ivpin(:)
   end type vpin_result

   type, public :: bayes_pin_result
      type(pin_parameters) :: posterior_mean
      real(dp), allocatable :: draws(:,:)
      real(dp) :: acceptance_rate = 0.0_dp
      integer :: status = 0
   end type bayes_pin_result

contains

   pure integer function trade_count_size(self) result(n)
      class(trade_counts), intent(in) :: self
      if (allocated(self%buys)) then
         n = size(self%buys)
      else
         n = 0
      end if
   end function trade_count_size

   pure logical function trade_count_valid(self) result(ok)
      class(trade_counts), intent(in) :: self
      ok = allocated(self%buys) .and. allocated(self%sells)
      if (ok) ok = size(self%buys) == size(self%sells) .and. all(self%buys >= 0_i8) .and. all(self%sells >= 0_i8)
   end function trade_count_valid

   pure integer function mpin_layer_count(self) result(n)
      class(mpin_parameters), intent(in) :: self
      if (allocated(self%alpha)) then
         n = size(self%alpha)
      else
         n = 0
      end if
   end function mpin_layer_count

end module pinstimation_types
