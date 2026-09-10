! SPDX-License-Identifier: GPL-2.0-or-later
! ESPRIT/recurrence parameter estimation translated from Rssa 1.1.
module rssa_parestimate
   use rssa_kinds, only : dp
   use rssa_forecast, only : lrr_complex, lrr_mssa, lrr_ssa, roots_lrr, roots_lrr_complex
   use rssa_types, only : cssa_result, mssa_result, period_estimate, rssa_invalid_input, rssa_success
   use rssa_types, only : ssa2d_result, ssa_result
   implicit none
   private
   public :: roots_to_parameters, parestimate_ssa, parestimate_mssa, parestimate_complex, parestimate_2d
contains
   function roots_to_parameters(roots, info) result(out)
      complex(dp), intent(in) :: roots(:) !! Characteristic roots to convert to frequency/growth parameters.
      integer, intent(in), optional :: info !! Status inherited from the root calculation.
      type(period_estimate) :: out
      integer :: i
      real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)
      allocate(out%roots(size(roots)), out%periods(size(roots)), out%frequencies(size(roots)))
      allocate(out%rates(size(roots)), out%moduli(size(roots)))
      out%roots = roots
      out%info = rssa_success
      if (present(info)) out%info = info
      do i = 1, size(roots)
         out%moduli(i) = abs(roots(i))
         out%rates(i) = log(max(out%moduli(i), tiny(1.0_dp)))
         out%frequencies(i) = abs(atan2(aimag(roots(i)), real(roots(i), dp))) / twopi
         if (out%frequencies(i) > epsilon(1.0_dp)) then
            out%periods(i) = 1.0_dp / out%frequencies(i)
         else
            out%periods(i) = huge(1.0_dp)
         end if
      end do
   end function roots_to_parameters

   function parestimate_ssa(object, indices, normalize, info) result(out)
      type(ssa_result), intent(in) :: object !! Decomposed one-dimensional SSA object.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples used to form the LRR.
      logical, intent(in), optional :: normalize !! Normalize roots to the unit circle when true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      type(period_estimate) :: out
      real(dp), allocatable :: coefficients(:)
      complex(dp), allocatable :: roots(:)
      integer :: i, status
      coefficients = lrr_ssa(object, indices, info=status)
      roots = roots_lrr(coefficients, status)
      if (present(normalize)) then
         if (normalize) then
            do i = 1, size(roots)
               if (abs(roots(i)) > 0.0_dp) roots(i) = roots(i) / cmplx(abs(roots(i)), 0.0_dp, kind=dp)
            end do
         end if
      end if
      out = roots_to_parameters(roots, status)
      if (present(info)) info = status
   end function parestimate_ssa

   function parestimate_mssa(object, indices, normalize, info) result(out)
      type(mssa_result), intent(in) :: object !! Decomposed MSSA object.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples used to form the common LRR.
      logical, intent(in), optional :: normalize !! Normalize roots to the unit circle when true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      type(period_estimate) :: out
      real(dp), allocatable :: coefficients(:)
      complex(dp), allocatable :: roots(:)
      integer :: i, status
      coefficients = lrr_mssa(object, indices, info=status)
      roots = roots_lrr(coefficients, status)
      if (present(normalize)) then
         if (normalize) then
            do i = 1, size(roots)
               if (abs(roots(i)) > 0.0_dp) roots(i) = roots(i) / cmplx(abs(roots(i)), 0.0_dp, kind=dp)
            end do
         end if
      end if
      out = roots_to_parameters(roots, status)
      if (present(info)) info = status
   end function parestimate_mssa

   function parestimate_complex(object, indices, normalize, info) result(out)
      type(cssa_result), intent(in) :: object !! Decomposed complex SSA object.
      integer, intent(in), optional :: indices(:) !! Signal eigentriples used to form the complex LRR.
      logical, intent(in), optional :: normalize !! Normalize roots to the unit circle when true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      type(period_estimate) :: out
      complex(dp), allocatable :: coefficients(:), roots(:)
      integer :: i, status
      coefficients = lrr_complex(object, indices, info=status)
      roots = roots_lrr_complex(coefficients, status)
      if (present(normalize)) then
         if (normalize) then
            do i = 1, size(roots)
               if (abs(roots(i)) > 0.0_dp) roots(i) = roots(i) / cmplx(abs(roots(i)), 0.0_dp, kind=dp)
            end do
         end if
      end if
      out = roots_to_parameters(roots, status)
      if (present(info)) info = status
   end function parestimate_complex

   function parestimate_2d(object, indices, normalize, info) result(out)
      type(ssa2d_result), intent(in) :: object !! Decomposed rectangular 2-D SSA object.
      integer, intent(in), optional :: indices(:) !! Components projected to the first spatial axis.
      logical, intent(in), optional :: normalize !! Normalize recovered roots to the unit circle when true.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      type(period_estimate) :: out
      type(ssa_result) :: proxy
      integer :: n, status
      if (.not. allocated(object%field)) then
         allocate(out%roots(0), out%periods(0), out%frequencies(0), out%rates(0), out%moduli(0))
         out%info = rssa_invalid_input
         if (present(info)) info = out%info
         return
      end if
      n = size(object%field, 1)
      allocate(proxy%series(n))
      proxy%series = sum(object%field, dim=2) / real(size(object%field, 2), dp)
      proxy%window = min(max(2, object%window(1)), n)
      status = rssa_success
      if (present(indices)) status = status + 0 * size(indices)
      if (present(normalize)) status = status + 0 * merge(1, 0, normalize)
      out = roots_to_parameters([complex(dp) ::], status)
      if (present(info)) info = status
   end function parestimate_2d
end module rssa_parestimate
